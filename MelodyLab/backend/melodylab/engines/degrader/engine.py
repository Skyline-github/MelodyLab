"""降解剂评估排序引擎：三族模型共识排名 + conformal 90% 区间 + 适用域检查。

设计原则（继承自研究项目的结论）：
- 预测只有**排序**意义，绝对值不可信 → 区间与适用域是一等公民输出；
- 三族（MLP 集成 / ExtraTrees / HistGB）各自排名取均值作为共识排名；
- 在库分子（PROTAC-DB 已收录）显式标注，防止把记忆当预测。
"""
from __future__ import annotations

import time
import uuid

import numpy as np
from rdkit import Chem
from scipy.stats import rankdata

from . import config as cfg
from .features import canonical_smiles, featurize, tanimoto_to_train
from .models import Calibration, CvMetric, DegraderReport, RankedRow
from .training import get_bundle, mlp_predict

FAMILIES = ("mlp", "et", "hgb")


def _train_fingerprints(train_smiles: list[str]) -> list:
    """训练集适用域指纹，一次计算整个查询任务复用。"""
    from .features import _AD_GEN

    fps = []
    for smi in train_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            fps.append(_AD_GEN.GetFingerprint(mol))
    return fps


def _family_rank(pred: np.ndarray) -> np.ndarray:
    """DC50 越低（越强）排名越靠前；rankdata 从 1 开始。"""
    return rankdata(pred, method="average")


def rank_compounds(items: list[dict], bundle: dict, model_status: str = "cached") -> DegraderReport:
    """items: features.parse_* 产出的 [{compound_id, smiles, e3?, source_line}]。"""
    report = DegraderReport(
        id=f"deg_{time.strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:6]}",
        created_at=time.strftime("%Y-%m-%d %H:%M:%S"),
        model_status=model_status,
        n_query=len(items),
        calibration=Calibration(
            n_train=bundle["n_train"], n_folds=bundle["recipe"]["n_folds"],
            q90_residual_log=bundle["calibration"]["q90_residual_log"],
            cv={k: CvMetric(**v) for k, v in bundle["calibration"]["cv"].items()},
            trained_at=bundle["trained_at"], db_snapshot=bundle["db_snapshot"]),
    )

    # 规范化 + 无效/重复处理
    seen: dict[str, str] = {}
    clean_items: list[dict] = []
    for it in items:
        can = canonical_smiles(it["smiles"])
        if can is None:
            report.rows.append(RankedRow(
                compound_id=it.get("compound_id") or f"line_{it.get('source_line', 0)}",
                smiles=it["smiles"], e3=it.get("e3", ""),
                status="invalid_smiles",
                source_line=it.get("source_line", 0),
                warnings=["SMILES 无法解析"]))
            report.n_invalid += 1
            continue
        if can in seen:
            report.n_duplicates += 1
            continue
        seen[can] = it.get("compound_id") or f"mol_{len(clean_items) + 1:03d}"
        clean_items.append({**it, "smiles": can})
    report.n_query = len(items)

    if not clean_items:
        report.notes.append("没有可评估的有效分子。")
        return report

    X = featurize([it["smiles"] for it in clean_items])
    mlp_preds = mlp_predict(bundle["mlps"], X)
    preds = {
        "mlp": np.mean(mlp_preds, axis=0),
        "et": bundle["et"].predict(X),
        "hgb": bundle["hgb"].predict(X),
    }
    members = np.vstack([preds["hgb"], preds["et"], *mlp_preds])  # 10 个成员用于 std
    ranks = {f: _family_rank(preds[f]) for f in FAMILIES}
    consensus_rank = np.mean([ranks[f] for f in FAMILIES], axis=0)

    q90 = bundle["calibration"]["q90_residual_log"]
    train_fps = _train_fingerprints(bundle["train_smiles"])
    for i, it in enumerate(clean_items):
        pred_log = float(np.mean([preds[f][i] for f in FAMILIES]))
        std = float(np.std(members[:, i]))
        z = max(cfg.Z90 * std, q90)
        row = RankedRow(
            compound_id=it.get("compound_id") or seen[it["smiles"]],
            smiles=it["smiles"], e3=it.get("e3", ""),
            source_line=it.get("source_line", 0),
            pred_log_dc50=round(pred_log, 3),
            pred_dc50_nm=round(float(np.exp(pred_log)), 2),
            lo90_nm=round(float(np.exp(pred_log - z)), 2),
            hi90_nm=round(float(np.exp(pred_log + z)), 2),
            rank_consensus=round(float(consensus_rank[i]), 2),
            rank_mlp=round(float(ranks["mlp"][i]), 2),
            rank_et=round(float(ranks["et"][i]), 2),
            rank_hgb=round(float(ranks["hgb"][i]), 2),
        )
        # ---- 适用域与在库检查 ----
        row.max_tanimoto = round(tanimoto_to_train(it["smiles"], [], train_fps=train_fps), 3)
        db_dc50 = bundle["db_dc50_map"].get(it["smiles"], "")
        if db_dc50:
            row.in_db = True
            row.db_dc50 = db_dc50
            row.warnings.append("PROTAC-DB 已收录（有实测 DC50），预测含记忆效应")
        if row.max_tanimoto < 0.35:
            row.warnings.append("与训练集最近邻相似度低，属外推区")
        mw, tpsa = X[i, cfg.ECFP_BITS], X[i, cfg.ECFP_BITS + 2]
        if mw > bundle["train_q95"]["MW"] or tpsa > bundle["train_q95"]["TPSA"]:
            row.warnings.append("MW/TPSA 超出训练集 95% 分位，透膜性可能限制实际细胞活性")
        report.rows.append(row)

    report.rows.sort(key=lambda r: (r.status != "ok", r.rank_consensus if r.rank_consensus else 1e9))
    report.notes = _methodology_notes(report)
    return report


def _methodology_notes(report: DegraderReport) -> list[str]:
    q90 = report.calibration.q90_residual_log if report.calibration else 0
    return [
        f"训练集为 PROTAC-DB 快照中 DC50 可回归的 {report.calibration.n_train if report.calibration else '—'} 行"
        f"（log nM 目标），骨架 5 折交叉验证；90% 区间下限（conformal q90）达 {q90:.2f} log 单位。",
        "本工具输出的是**相对排序参考**，不是绝对 DC50 预测：共识排名（三族均值）比任何单点数值可信；"
        "区间普遍跨越 1–2 个数量级属正常现象。",
        "三族模型（MLP 集成 / ExtraTrees / HistGB）排名可能明显不一致——对训练集没有的骨架，"
        "排名是模型相关的；同时进入多族前列的分子才是稳健候选。",
        "MW/TPSA 超分布的分子（大极性 PROTAC 常见）实际细胞活性很可能受透膜性限制，"
        "生化效力与细胞降解效力要分开解读。",
        "模型不使用蛋白序列信息（ESM 嵌入属二期可选增强）；对训练集未覆盖的 POI 一律视为外推。",
    ]
