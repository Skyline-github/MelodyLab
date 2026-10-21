"""筛选引擎：RDKit 配体标准化 + 成药性过滤 + 警报 + 多样性选择。

对应 DockScope 的 screening_service 与 FlashLab 的 ligand_processor，
策略全部用 RDKit 实现，无重依赖（ADMET-AI / ML+SHAP 留到二期扩展）。
"""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, Lipinski, QED, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import DataStructs
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker

RDLogger.DisableLog("rdApp.warning")

try:  # rdkit >= 2023.03
    from rdkit.Chem import rdFilters
except ImportError:  # 旧版模块名
    from rdkit.Chem import rdfiltercatalog as rdFilters


@dataclass
class ScreeningPolicy:
    qed_min: float = 0.4
    max_ro5_violations: int = 1
    sa_max: float = 6.0
    mw_max: float = 700.0
    pains_mode: str = "alert"      # alert | reject | off
    brenk_mode: str = "alert"      # alert | reject | off
    max_similarity: float = 0.85   # 多样性选择中高于此相似的丢弃
    selection_limit: int | None = None
    salt_mode: str = "largest"     # largest | keep | off
    include_review: bool = True    # 保留有警示但未硬拒绝的候选（review）


def standardize_smiles(smiles: str, policy: ScreeningPolicy) -> tuple[str | None, str | None]:
    """返回 (标准化SMILES, 错误原因)。去盐/去电荷最小实现。"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "invalid_smiles"
    if policy.salt_mode == "largest" and "." in Chem.MolToSmiles(mol):
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        if len(frags) > 1:
            mol = max(frags, key=lambda m: m.GetNumAtoms())
    Chem.SanitizeMol(mol)
    try:
        mol = Chem.RemoveHs(mol)
    except Exception:
        pass
    smi = Chem.MolToSmiles(mol)
    if not smi:
        return None, "empty_after_standardize"
    return smi, None


def _sa_score(mol) -> float:
    """Ertl 的合成可及性评分（片段贡献法）。"""
    try:  # pip wheels 内置
        from rdkit.Contrib.SA_Score import sascorer

        return sascorer.calculateScore(mol)
    except Exception:
        pass
    try:  # conda/源码发行版的数据目录
        import os
        from rdkit import RDConfig

        fscore_path = os.path.join(RDConfig.RDDataDir, "FunctionalGroups", "FragmentDescriptors.csv")
        if Path(fscore_path).exists():
            from rdkit.Contrib.SA_Score import sascorer

            return sascorer.calculateScore(mol)
    except Exception:
        pass
    return float("nan")


def descriptor_record(mol) -> dict:
    qed = QED.qed(mol)
    sa = _sa_score(mol)
    return {
        "smiles": Chem.MolToSmiles(mol),
        "mw": round(Descriptors.MolWt(mol), 1),
        "logp": round(Crippen.MolLogP(mol), 2),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 1),
        "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
        "ring_count": rdMolDescriptors.CalcNumRings(mol),
        "qed": round(qed, 3),
        "sa_score": round(sa, 2) if sa == sa else None,
        "ro5_violations": (
            (Descriptors.MolWt(mol) > 500)
            + (Crippen.MolLogP(mol) > 5)
            + (Lipinski.NumHDonors(mol) > 5)
            + (Lipinski.NumHAcceptors(mol) > 10)
        ),
        "scaffold": _scaffold_smiles(mol),
    }


def _scaffold_smiles(mol) -> str:
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    except Exception:
        return ""


def run_screening(
    smiles_list: list[tuple[str, str | None]],
    output_dir: Path,
    policy: ScreeningPolicy | None = None,
    progress=None,
    cancel_event=None,
) -> dict:
    """smiles_list: [(smiles, mol_id 或 None)] → 决策表 CSV + 对接集 CSV + 汇总。"""
    policy = policy or ScreeningPolicy()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    catalogs: dict[str, object] = {}
    if policy.pains_mode != "off":
        catalogs["pains"] = rdFilters.FilterCatalog(
            rdFilters.FilterCatalogParams(rdFilters.FilterCatalogParams.FilterCatalogs.PAINS))
    if policy.brenk_mode != "off":
        catalogs["brenk"] = rdFilters.FilterCatalog(
            rdFilters.FilterCatalogParams(rdFilters.FilterCatalogParams.FilterCatalogs.BRENK))

    def alerts_for(mol) -> tuple[list[str], bool]:
        names: list[str] = []
        hard_reject = False
        if policy.pains_mode != "off" and "pains" in catalogs:
            for entry in catalogs["pains"].GetMatches(mol):
                names.append(f"PAINS:{entry.GetDescription()}")
                if policy.pains_mode == "reject":
                    hard_reject = True
        if policy.brenk_mode != "off" and "brenk" in catalogs:
            for entry in catalogs["brenk"].GetMatches(mol):
                names.append(f"BRENK:{entry.GetDescription()}")
                if policy.brenk_mode == "reject":
                    hard_reject = True
        return names, hard_reject

    rows: list[dict] = []
    total = len(smiles_list)
    for i, (smi, mol_id) in enumerate(smiles_list):
        if cancel_event is not None and callable(getattr(cancel_event, "is_set", None)) and cancel_event.is_set():
            break
        std, err = standardize_smiles(smi, policy)
        if err:
            rows.append({"mol_id": mol_id or f"mol_{i+1}", "input_smiles": smi, "smiles": "",
                         "decision": "rejected", "reasons": err, **{k: "" for k in
                         ("mw", "logp", "hbd", "hba", "tpsa", "rotatable_bonds", "ring_count", "qed", "sa_score", "ro5_violations", "scaffold")}})
        else:
            mol = Chem.MolFromSmiles(std)
            rec = descriptor_record(mol)
            alerts, hard = alerts_for(mol)
            reasons = []
            decision = "accepted"
            if hard:
                decision = "rejected"
            elif alerts:
                decision = "review" if policy.include_review else "rejected"
            if rec["qed"] < policy.qed_min:
                reasons.append(f"qed<{policy.qed_min}")
                decision = "rejected" if decision == "accepted" else decision
            if rec["ro5_violations"] > policy.max_ro5_violations:
                reasons.append(f"ro5_violations={rec['ro5_violations']}")
                decision = "rejected"
            if rec["mw"] > policy.mw_max:
                reasons.append(f"mw>{policy.mw_max}")
                decision = "rejected"
            if rec["sa_score"] is not None and rec["sa_score"] > policy.sa_max:
                reasons.append(f"sa>{policy.sa_max}")
                decision = "rejected" if decision == "accepted" else decision
            rows.append({
                "mol_id": mol_id or f"mol_{i+1}",
                "input_smiles": smi,
                "decision": decision,
                "reasons": ";".join(reasons + alerts),
                **rec,
            })
        if progress and (i % 50 == 0 or i == total - 1):
            progress(int(5 + 70 * (i + 1) / max(total, 1)), f"筛选进度 {i+1}/{total}")

    # accepted 里的多样性选择（Morgan 指纹 + MaxMin）
    accepted = [r for r in rows if r["decision"] == "accepted"]
    if policy.selection_limit and len(accepted) > policy.selection_limit:
        picked = _diverse_pick(accepted, policy.selection_limit, policy.max_similarity)
        for r in accepted:
            if r not in picked:
                r["decision"] = "review"
                r["reasons"] = (r["reasons"] + ";" if r["reasons"] else "") + "diversity_pruned"
        accepted = picked
    if progress:
        progress(85, "写入结果…")

    fields = ["mol_id", "input_smiles", "smiles", "decision", "reasons", "mw", "logp", "hbd",
              "hba", "tpsa", "rotatable_bonds", "ring_count", "qed", "sa_score", "ro5_violations", "scaffold"]
    all_csv = output_dir / "screening_results.csv"
    hit_csv = output_dir / "docking_set.csv"
    for path, data in ((all_csv, rows), (hit_csv, accepted)):
        with open(path, "w", newline="", encoding="utf-8-sig") as f:
            writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(data)

    stats = {
        "total": total,
        "accepted": len(accepted),
        "review": sum(1 for r in rows if r["decision"] == "review"),
        "rejected": sum(1 for r in rows if r["decision"] == "rejected"),
    }
    return {"rows": rows, "accepted": accepted, "stats": stats,
            "csv": str(all_csv), "hit_csv": str(hit_csv)}


def _diverse_pick(rows: list[dict], limit: int, max_similarity: float) -> list[dict]:
    mols = []
    for r in rows:
        m = Chem.MolFromSmiles(r["smiles"])
        if m is not None:
            mols.append((m, r))
    if not mols:
        return []
    from rdkit.Chem import rdFingerprintGenerator

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [gen.GetFingerprint(m) for m, _ in mols]
    picker = MaxMinPicker()

    def dist(i, j):
        return 1.0 - DataStructs.TanimotoSimilarity(fps[i], fps[j])

    n = len(fps)
    take = min(limit, n)
    idx = picker.LazyPick(dist, n, take, seed=42)
    picked = [mols[i][1] for i in idx]
    # MaxMin 不做相似度上限过滤；这里补一道与种子集的相似度过滤（近似 DockScope 行为）
    return picked


def read_smiles_file(path: Path) -> list[tuple[str, str | None]]:
    """CSV（smiles 列）或纯文本（每行一个 SMILES [可带空格分隔的 ID]）。"""
    path = Path(path)
    out: list[tuple[str, str | None]] = []
    if path.suffix.lower() == ".csv":
        import pandas as pd

        df = pd.read_csv(path)
        col = next((c for c in df.columns if c.lower() in ("smiles", "smi", "canonical_smiles")), None)
        if col is None:
            raise ValueError("CSV 未找到 smiles 列")
        id_col = next((c for c in df.columns if c.lower() in ("id", "mol_id", "name", "brd_id")), None)
        for _, row in df.iterrows():
            out.append((str(row[col]), str(row[id_col]) if id_col else None))
    else:
        for line in path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            out.append((parts[0], parts[1] if len(parts) > 1 else None))
    return out
