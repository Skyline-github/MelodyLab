"""降解剂评估模型训练：骨架 5 折交叉验证 + 三族模型（HGB / ExtraTrees / MLP 集成）。

- 校准：折外（OOF）残差的 q90 作为 conformal 区间下限；各族 R²/Spearman 报告真实技能。
- 最终模型在全部清洗数据上训练；区间宽度 = max(1.645×成员间 std, q90 折外残差)。
- 训练产物缓存到 <data_root>/cache/degrader/，键含配方版本 + 数据文件指纹。
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import joblib
import numpy as np
from rdkit import Chem

from . import config as cfg
from .features import (canonical_smiles, descriptor_matrix, ecfp_matrix,
                       featurize, load_training_data, murcko_scaffold)


# ---------------------------------------------------------------- MLP（纯 torch 复刻教程结构）
def _build_mlp(n_features: int):
    import torch.nn as nn

    h1, h2 = cfg.MLP_HIDDEN
    return nn.Sequential(
        nn.Linear(n_features, h1), nn.ReLU(), nn.Dropout(0.0),
        nn.Linear(h1, h2), nn.ReLU(), nn.Dropout(cfg.MLP_DROPOUT),
        nn.Linear(h2, 1),
    )


def _train_mlp(X: np.ndarray, y: np.ndarray, seed: int, epochs: int) -> dict[str, np.ndarray]:
    import torch

    torch.set_num_threads(cfg.TORCH_THREADS)
    torch.manual_seed(seed)
    model = _build_mlp(X.shape[1])
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    xt = torch.from_numpy(X.astype(np.float32))
    yt = torch.from_numpy(y.astype(np.float32))
    model.train()
    n = len(xt)
    for _ in range(epochs):
        perm = torch.randperm(n)
        for i in range(0, n, cfg.MLP_BATCH):
            idx = perm[i:i + cfg.MLP_BATCH]
            opt.zero_grad()
            loss = ((model(xt[idx]).ravel() - yt[idx]) ** 2).mean()
            loss.backward()
            opt.step()
    return {k: v.detach().cpu().numpy() for k, v in model.state_dict().items()}


def mlp_predict(weights_list: list[dict[str, np.ndarray]], X: np.ndarray) -> list[np.ndarray]:
    import torch

    xt = torch.from_numpy(X.astype(np.float32))
    out = []
    for weights in weights_list:
        model = _build_mlp(X.shape[1])
        model.load_state_dict({k: torch.from_numpy(v) for k, v in weights.items()})
        model.eval()
        with torch.no_grad():
            out.append(model(xt).ravel().numpy())
    return out


# ---------------------------------------------------------------- 树模型
def _train_hgb(X, y, seed):
    from sklearn.ensemble import HistGradientBoostingRegressor

    return HistGradientBoostingRegressor(random_state=seed)


def _train_et(X, y, seed):
    from sklearn.ensemble import ExtraTreesRegressor

    return ExtraTreesRegressor(n_estimators=500, min_samples_leaf=2,
                               random_state=seed, n_jobs=-1)


# ---------------------------------------------------------------- 骨架分折
def scaffold_folds(smiles: list[str], n_folds: int = cfg.N_FOLDS,
                   seed: int = cfg.SEED) -> np.ndarray:
    """Murcko 骨架分组轮转分配（同组分子必在同折，防骨架泄漏）。"""
    groups: dict[str, list[int]] = {}
    for i, smi in enumerate(smiles):
        groups.setdefault(murcko_scaffold(smi), []).append(i)
    rng = np.random.default_rng(seed)
    keys = sorted(groups.keys())
    rng.shuffle(keys)
    fold_of = np.zeros(len(smiles), dtype=int)
    for i, k in enumerate(keys):
        for idx in groups[k]:
            fold_of[idx] = i % n_folds
    return fold_of


# ---------------------------------------------------------------- 训练与校准
def _cv_metrics(y: np.ndarray, oof: np.ndarray) -> dict:
    from scipy.stats import spearmanr

    mask = np.isfinite(oof)
    yy, pp = y[mask], oof[mask]
    if len(yy) < 10:
        return {"r2": None, "spearman": None}
    ss_res = float(np.sum((yy - pp) ** 2))
    ss_tot = float(np.sum((yy - yy.mean()) ** 2))
    return {"r2": round(1 - ss_res / ss_tot, 4) if ss_tot > 0 else None,
            "spearman": round(float(spearmanr(yy, pp).statistic), 4)}


def train_all(control=None) -> dict:
    """全量训练 + 5 折 CV 校准，返回模型 bundle。"""

    def cp(progress, message):
        if control is not None:
            control.checkpoint(progress=progress, phase="train", message=message)

    cp(5, "读取并清洗 PROTAC-DB 快照…")
    df = load_training_data()
    smiles_all = df["Smiles"].tolist()
    valid = [canonical_smiles(s) is not None for s in smiles_all]
    n_dropped = int(len(df) - sum(valid))
    df = df[valid].reset_index(drop=True)
    if len(df) < cfg.MIN_TRAIN_ROWS:
        raise RuntimeError(f"训练数据不足：仅 {len(df)} 行可用（阈值 {cfg.MIN_TRAIN_ROWS}）")
    smiles = df["Smiles"].tolist()
    y = df["y_dc50_log"].to_numpy(dtype=np.float64)

    cp(12, f"计算分子特征（{len(smiles)} 行）…")
    X = featurize(smiles)

    # ---- DB 重叠映射（全库任何行，含无 DC50 的行）----
    cp(18, "建立 PROTAC-DB 实测值索引…")
    db_dc50_map = _build_db_map()

    # ---- 骨架 5 折 CV ----
    fold_of = scaffold_folds(smiles)
    families = {"hgb": _train_hgb, "et": _train_et}
    oof: dict[str, np.ndarray] = {f: np.full(len(y), np.nan) for f in families}
    mlp_oof_seeds = [np.full(len(y), np.nan) for _ in range(cfg.MLP_CV_SEEDS)]
    n_folds = cfg.N_FOLDS
    for k in range(n_folds):
        tr, va = np.where(fold_of != k)[0], np.where(fold_of == k)[0]
        cp(20 + int(40 * k / n_folds),
           f"骨架交叉验证 {k + 1}/{n_folds} 折（树模型 + 轻量 MLP）…")
        for fam, fn in families.items():
            model = fn(X[tr], y[tr], cfg.SEED + k)
            model.fit(X[tr], y[tr])
            oof[fam][va] = model.predict(X[va])
        for s in range(cfg.MLP_CV_SEEDS):
            w = _train_mlp(X[tr], y[tr], cfg.SEED + 100 + s, cfg.MLP_CV_EPOCHS)
            mlp_oof_seeds[s][va] = mlp_predict([w], X[va])[0]
    oof["mlp"] = np.nanmean(np.vstack(mlp_oof_seeds), axis=0)

    cv = {fam: _cv_metrics(y, oof[fam]) for fam in oof}
    # conformal 下限：共识（三族 OOF 均值）残差的 q90 —— 与展示的共识预测同口径
    consensus_oof = np.nanmean(np.vstack([oof[f] for f in oof]), axis=0)
    finite = np.isfinite(consensus_oof)
    q90 = float(np.quantile(np.abs(consensus_oof[finite] - y[finite]), cfg.CONFORMAL_Q))
    q90_families = {fam: float(np.quantile(np.abs(oof[fam][np.isfinite(oof[fam])]),
                                           cfg.CONFORMAL_Q)) for fam in oof}

    # ---- 最终模型（全部数据）----
    cp(65, "训练最终树模型（全量数据）…")
    final = {"hgb": _train_hgb(X, y, cfg.SEED), "et": _train_et(X, y, cfg.SEED)}
    for name, model in final.items():
        model.fit(X, y)
    cp(75, f"训练 MLP 集成（{cfg.MLP_FINAL_SEEDS} 成员 × {cfg.MLP_FINAL_EPOCHS} epochs）…")
    mlps = [_train_mlp(X, y, cfg.SEED + 1000 + s, cfg.MLP_FINAL_EPOCHS)
            for s in range(cfg.MLP_FINAL_SEEDS)]

    desc = descriptor_matrix([Chem.MolFromSmiles(s) for s in smiles])
    train_q95 = {"MW": float(np.quantile(desc[:, 0], 0.95)),
                 "TPSA": float(np.quantile(desc[:, 2], 0.95))}

    bundle = {
        "version": cfg.MODEL_VERSION,
        "signature": cfg.db_signature(),
        "trained_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "db_snapshot": cfg.DB_SNAPSHOT,
        "n_train": len(y),
        "n_dropped_invalid": n_dropped,
        "recipe": {"ecfp_radius": cfg.ECFP_RADIUS, "ecfp_bits": cfg.ECFP_BITS,
                   "descriptors": cfg.DESC_NAMES, "n_folds": n_folds,
                   "mlp_final": [cfg.MLP_FINAL_SEEDS, cfg.MLP_FINAL_EPOCHS],
                   "mlp_cv": [cfg.MLP_CV_SEEDS, cfg.MLP_CV_EPOCHS]},
        "calibration": {"q90_residual_log": round(q90, 4),
                        "q90_per_family": {k: round(v, 4) for k, v in q90_families.items()},
                        "cv": cv},
        "train_q95": train_q95,
        "hgb": final["hgb"],
        "et": final["et"],
        "mlps": mlps,
        "train_smiles": smiles,
        "db_dc50_map": db_dc50_map,
    }
    cp(92, "保存模型缓存…")
    _save_bundle(bundle)
    return bundle


def _build_db_map() -> dict[str, str]:
    import pandas as pd

    db = pd.read_csv(cfg.DB_FILE, low_memory=False)[["Smiles", "DC50 (nM)"]].dropna(subset=["Smiles"])
    out: dict[str, str] = {}
    for smi, dc50 in zip(db["Smiles"], db["DC50 (nM)"]):
        can = canonical_smiles(str(smi))
        if can and can not in out and pd.notna(dc50):
            out[can] = str(dc50)
    return out


def _save_bundle(bundle: dict) -> None:
    path = cfg.cache_path()
    joblib.dump(bundle, path, compress=3)
    meta = {k: bundle[k] for k in ("version", "signature", "trained_at", "db_snapshot",
                                   "n_train", "n_dropped_invalid", "recipe", "calibration",
                                   "train_q95")}
    Path(str(path) + ".meta.json").write_text(
        json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    # 清理旧签名的缓存文件，避免堆积
    for old in cfg.CACHE_DIR.glob("models_*.joblib"):
        if old != path:
            old.unlink(missing_ok=True)
            Path(str(old) + ".meta.json").unlink(missing_ok=True)


def load_bundle() -> dict | None:
    path = cfg.cache_path()
    if not path.is_file():
        return None
    try:
        bundle = joblib.load(path)
        if bundle.get("signature") != cfg.db_signature():
            return None
        return bundle
    except Exception:  # noqa: BLE001 —— 缓存损坏时重训
        return None


def get_bundle(control=None) -> tuple[dict, str]:
    """取模型：缓存优先；否则全量训练（首次约 1–2 分钟）。"""
    bundle = load_bundle()
    if bundle is not None:
        return bundle, "cached"
    return train_all(control), "trained"


def invalidate_cache() -> int:
    n = 0
    for p in cfg.CACHE_DIR.glob("models_*.joblib*"):
        p.unlink(missing_ok=True)
        n += 1
    return n
