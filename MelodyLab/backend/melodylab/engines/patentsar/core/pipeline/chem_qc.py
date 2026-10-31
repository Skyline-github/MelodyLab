"""RDKit-based chemistry quality control for OCR'd structures.

Sanitizes / canonicalizes SMILES, computes formula & molecular weight, and --
when a patent LCMS ``[M+H]+`` value is available -- cross-checks the observed
m/z against the RDKit-predicted one. Structures that fail this check (or are
duplicates) are flagged so a human can review them, mirroring the blog's
"高亮标注疑似错误" QC step.
"""
from __future__ import annotations

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

_PROTON = 1.007276


def canonicalize(smiles: str) -> str | None:
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def validate(smiles: str, observed_mz: float | None = None, mz_tol: float = 1.5) -> dict:
    flags: list[str] = []
    out = {
        "input_smiles": smiles,
        "canonical_smiles": None,
        "formula": None,
        "mw_avg": None,
        "mw_exact": None,
        "mz_calc": None,
        "mz_observed": observed_mz,
        "delta_mz": None,
        "valid": False,
        "flags": flags,
    }
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        flags.append("INVALID_SMILES")
        return out

    mol = Chem.RemoveHs(mol)
    out["canonical_smiles"] = Chem.MolToSmiles(mol)
    out["formula"] = rdMolDescriptors.CalcMolFormula(mol)
    out["mw_avg"] = round(Descriptors.MolWt(mol), 2)
    exact = Descriptors.ExactMolWt(mol)
    out["mw_exact"] = round(exact, 3)
    out["valid"] = True

    if observed_mz is not None:
        mz_calc = exact + _PROTON
        out["mz_calc"] = round(mz_calc, 3)
        delta = observed_mz - mz_calc
        out["delta_mz"] = round(delta, 3)
        if abs(delta) > mz_tol:
            flags.append("LCMS_MISMATCH")
    return out


def annotate_duplicates(rows: list[dict]) -> None:
    """Flag rows sharing the same canonical SMILES (in place)."""
    seen: dict[str, str] = {}
    for r in rows:
        c = r["qc"].get("canonical_smiles")
        if not c:
            continue
        if c in seen:
            r["qc"]["flags"].append(f"DUPLICATE_OF:{seen[c]}")
        else:
            seen[c] = r["example_no"]
