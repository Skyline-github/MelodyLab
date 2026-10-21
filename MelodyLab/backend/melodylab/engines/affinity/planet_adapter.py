"""PLANET 亲和力预测适配器（复用 assets/vendor/PLANET，pKd）。"""
from __future__ import annotations

import os
import sys
from pathlib import Path

from ...config import settings


def availability() -> dict:
    vendor = settings.vendor_dir("PLANET")
    param = settings.checkpoint("planet")
    return {
        "engine": "planet",
        "vendor": str(vendor) if vendor else None,
        "checkpoint": str(param) if param else None,
        "available": bool(vendor and param),
        "hint": None if (vendor and param) else "运行 scripts/setup_assets.py 复制 PLANET 代码与 PLANET.param 权重",
    }


def predict_affinity(protein_pdb: Path, ligand_sdf: Path) -> dict:
    info = availability()
    if not info["available"]:
        raise RuntimeError(f"PLANET 资产未就绪：{info['hint']}")
    import pandas as pd

    vendor = settings.vendor_dir("PLANET")
    original_cwd = os.getcwd()
    try:
        if str(vendor) not in sys.path:
            sys.path.insert(0, str(vendor))
        for mod in [k for k in list(sys.modules) if k.split(".")[0] in ("chemutils", "PLANET_model", "pred")]:
            del sys.modules[mod]
        os.chdir(vendor)
        import pred as planet_pred

        import tempfile

        with tempfile.TemporaryDirectory() as td:
            out_csv = Path(td) / "affinity.csv"
            planet_pred.single_process(str(protein_pdb), str(ligand_sdf), str(ligand_sdf), str(out_csv))
            if not out_csv.is_file():
                raise RuntimeError("PLANET 未产生输出")
            df = pd.read_csv(out_csv)
    finally:
        os.chdir(original_cwd)

    for col in ("Binding_Affinity", "affinity", "pKd", "predicted_affinity"):
        if col in df.columns and len(df) > 0:
            return {"affinity": float(df[col].iloc[0]), "unit": "pKd", "engine": "planet"}
    raise RuntimeError("PLANET 输出中未找到亲和力列")
