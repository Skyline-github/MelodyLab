"""降解剂评估特征工程：PROTAC-DB 清洗、ECFP/描述符、Murcko 骨架、查询解析。

清洗规则与研究项目一致：DC50 (nM) 去掉缺失与含 '<' '>' '/' '~' '-' 的截断值，
取自然对数作为回归目标（log nM）。
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

from . import config as cfg


# ---------------------------------------------------------------- 分子特征
_ECFP_GEN = rdFingerprintGenerator.GetMorganGenerator(
    radius=cfg.ECFP_RADIUS, fpSize=cfg.ECFP_BITS, includeChirality=True)
_AD_GEN = rdFingerprintGenerator.GetMorganGenerator(
    radius=cfg.TANIMOTO_RADIUS, fpSize=cfg.ECFP_BITS)


def canonical_smiles(smi: str) -> str | None:
    mol = Chem.MolFromSmiles(smi) if smi else None
    return Chem.MolToSmiles(mol) if mol else None


def ecfp_matrix(smiles: list[str]) -> np.ndarray:
    out = np.zeros((len(smiles), cfg.ECFP_BITS), dtype=np.float32)
    for i, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            out[i] = _ECFP_GEN.GetFingerprint(mol)
    return out


def descriptor_matrix(mols: list) -> np.ndarray:
    rows = []
    for mol in mols:
        rows.append([
            Descriptors.MolWt(mol),
            Crippen.MolLogP(mol),
            Descriptors.TPSA(mol),
            Lipinski.NumHAcceptors(mol),
            Lipinski.NumHDonors(mol),
            rdMolDescriptors.CalcNumRotatableBonds(mol),
            Descriptors.RingCount(mol),
            mol.GetNumHeavyAtoms(),
            rdMolDescriptors.CalcFractionCSP3(mol),
            rdMolDescriptors.CalcNumAmideBonds(mol),
        ])
    return np.nan_to_num(np.asarray(rows, dtype=np.float32),
                         nan=0.0, posinf=0.0, neginf=0.0)


def featurize(smiles: list[str]) -> np.ndarray:
    """ECFP4-chiral 与 10 个理化描述符横向拼接。"""
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    return np.hstack([ecfp_matrix(smiles), descriptor_matrix(mols)])


def murcko_scaffold(smi: str) -> str:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return ""
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    except Exception:  # noqa: BLE001 —— 个别价态异常分子退化为空骨架
        return ""


def tanimoto_to_train(query_smi: str, train_smiles: list[str],
                      train_fps: list | None = None) -> float:
    """查询分子与训练集最近邻的 Tanimoto 相似度（适用域）。"""
    mol = Chem.MolFromSmiles(query_smi)
    if mol is None:
        return 0.0
    from rdkit import DataStructs

    fp = _AD_GEN.GetFingerprint(mol)
    if train_fps is None:
        train_fps = [_AD_GEN.GetFingerprint(Chem.MolFromSmiles(s))
                     for s in train_smiles if Chem.MolFromSmiles(s)]
    sims = DataStructs.BulkTanimotoSimilarity(fp, train_fps)
    return float(max(sims)) if sims else 0.0


# ---------------------------------------------------------------- 训练数据
def load_training_data() -> pd.DataFrame:
    """读取打包的 PROTAC-DB 快照并清洗出可回归的 DC50 行（log nM）。"""
    db = pd.read_csv(cfg.DB_FILE, low_memory=False)
    raw = db["DC50 (nM)"].dropna().astype(str)
    raw = raw[~raw.str.contains("<|>|/|~|-", regex=True)]
    dc50_nm = raw.astype(float)
    dc50_nm = dc50_nm[dc50_nm > 0]
    clean = db.loc[dc50_nm.index].copy()
    clean["y_dc50_log"] = np.log(dc50_nm.values)
    clean = clean.dropna(subset=["Smiles"])
    return clean.reset_index(drop=True)


# ---------------------------------------------------------------- 查询解析
def parse_smiles_text(text: str) -> list[dict]:
    """纯文本解析：每行一个 SMILES，或 'ID<空格>SMILES'（ID 在前，如 J52 O=C(...)）。"""
    items: list[dict] = []
    for lineno, line in enumerate(text.splitlines(), start=1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2 and not _looks_like_smiles(parts[0]):
            items.append({"compound_id": parts[0], "smiles": " ".join(parts[1:]),
                          "source_line": lineno})
        else:
            items.append({"compound_id": f"mol_{lineno:03d}", "smiles": line,
                          "source_line": lineno})
    return items


def parse_query_csv(path) -> list[dict]:
    """CSV 解析：需 smiles 列（大小写不敏感），可选 id/compound_id/name 列与 E3 列。"""
    df = pd.read_csv(path)
    cols = {c.lower().strip(): c for c in df.columns}
    smi_col = next((cols[k] for k in ("smiles", "smi", "canonical_smiles") if k in cols), None)
    if smi_col is None:
        raise ValueError("CSV 未找到 smiles 列")
    id_col = next((cols[k] for k in ("compound_id", "compound id", "id", "name") if k in cols), None)
    e3_col = next((cols[k] for k in ("e3 ligase", "e3_ligase", "e3") if k in cols), None)
    items: list[dict] = []
    for i, row in df.iterrows():
        smi = str(row[smi_col]).strip()
        if not smi or smi.lower() == "nan":
            continue

        def _clean(val, default: str) -> str:
            s = str(val).strip()
            return s if s and s.lower() not in ("nan", "none") else default

        items.append({
            "compound_id": _clean(row[id_col], f"mol_{i + 1:03d}") if id_col else f"mol_{i + 1:03d}",
            "smiles": smi,
            "e3": _clean(row[e3_col], "") if e3_col else "",
            "source_line": i + 2,
        })
    return items


def _looks_like_smiles(token: str) -> bool:
    return bool(Chem.MolFromSmiles(token)) or any(ch in token for ch in "()=[]#@/\\")
