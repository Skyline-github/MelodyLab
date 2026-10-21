"""FlashDiff 分子生成适配器：复用 assets/vendor/FlashDiff_app 的开源实现与权重。

vendor 缺失时 available=False，API 返回明确的引导信息而不是 500。
"""
from __future__ import annotations

import os
import sys
import tempfile
from datetime import datetime
from pathlib import Path

import numpy as np

from ...config import settings


def availability() -> dict:
    vendor = settings.vendor_dir("FlashDiff_app")
    ckpt = settings.checkpoint("flashdiff")
    return {
        "engine": "flashdiff",
        "vendor": str(vendor) if vendor else None,
        "checkpoint": str(ckpt) if ckpt else None,
        "available": bool(vendor and ckpt),
        "hint": None if (vendor and ckpt) else "运行 scripts/setup_assets.py 复制 FlashDiff 代码与权重到 assets/vendor/",
    }


def _import_flashdiff(vendor: Path):
    if str(vendor) not in sys.path:
        sys.path.insert(0, str(vendor))
    # FlashDiff 与其他项目存在 utils/sample 命名冲突：先清缓存再导入
    for mod in [k for k in list(sys.modules) if k.split(".")[0] in ("utils", "sample", "models", "train") and "torch" not in k]:
        del sys.modules[mod]
    os.chdir(vendor)  # FlashDiff 内部有相对路径依赖
    from train import DataModule  # noqa: F401  (确保依赖链完整)
    import sample as flashdiff_sample

    return flashdiff_sample


def generate(protein_path: Path, ligand_path: Path, out_dir: Path, *,
             num_molecules: int = 5, batch_size: int = 10, pocket_radius: float = 15.0,
             num_steps: int = 50, seed: int = 2024, mean_atoms: int = 28, std_atoms: int = 5,
             progress=None, control=None) -> dict:
    info = availability()
    if not info["available"]:
        raise RuntimeError(f"FlashDiff 资产未就绪：{info['hint']}")
    from rdkit import Chem

    out_dir.mkdir(parents=True, exist_ok=True)
    vendor = settings.vendor_dir("FlashDiff_app")

    # 口袋中心 = 参考配体质心
    mol = Chem.MolFromMolFile(str(ligand_path))
    if mol is None:
        raise ValueError("无法解析参考配体 SDF")
    conf = mol.GetConformer()
    coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    pocket_center = np.mean(coords, axis=0).tolist()

    import argparse

    config = {
        "data": {
            "protein_path": str(protein_path),
            "is_pep": False,
            "pocket_args": {"pocket_coord": pocket_center, "radius": pocket_radius},
            "pocmol_args": {"data_id": "melodylab", "pdbid": "job"},
        },
        "sample": {"seed": seed, "batch_size": batch_size, "num_mols": num_molecules, "save_traj_prob": 1.0},
        "transforms": {
            "featurizer_pocket": {"center": pocket_center},
            "variable_mol_size": {
                "name": "variable_mol_size",
                "num_atoms_distri": {
                    "strategy": "mol_atoms_based",
                    "mean": {"coef": 0, "bias": mean_atoms},
                    "std": {"coef": 0, "bias": std_atoms},
                    "min": 5,
                },
            },
        },
        "task": {"name": "sbdd", "transform": {"name": "sbdd"}},
        "noise": {
            "name": "sbdd", "num_steps": num_steps, "prior": "from_train",
            "level": {"name": "advance", "min": 0.0, "max": 1.0,
                      "step2level": {"scale_start": 0.99999, "scale_end": 0.00001, "width": 3}},
        },
        "model": {"checkpoint": str(settings.checkpoint("flashdiff"))},
    }
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = out_dir / f"flashdiff_{timestamp}"

    original_cwd = os.getcwd()
    try:
        flashdiff_sample = _import_flashdiff(vendor)
        args = argparse.Namespace(outdir=str(output_dir.resolve()), device="cpu",
                                  batch_size=batch_size, shuffle=False, num_workers=0)
        flashdiff_sample.main(args, config_dict=config)
    finally:
        os.chdir(original_cwd)

    # 收集结果
    from rdkit import Chem as C

    generated = []
    if output_dir.is_dir():
        for sub in sorted(output_dir.rglob("*_SDF")):
            for sdf in sorted(sub.glob("*.sdf")):
                if sdf.name.endswith(("-out.sdf", "-bad.sdf")):
                    continue
                m = C.MolFromMolFile(str(sdf))
                if m is None:
                    continue
                smi = C.MolToSmiles(m)
                if "." in smi:
                    continue
                dest = out_dir / f"mol_{len(generated)+1}.sdf"
                dest.write_text(sdf.read_text(encoding="utf-8"), encoding="utf-8")
                traj = sub.parent / "SDF" / f"{sdf.stem}-out.sdf"
                traj_dest = None
                if traj.is_file():
                    traj_dest = out_dir / f"mol_{len(generated)+1}_trajectory.sdf"
                    traj_dest.write_text(traj.read_text(encoding="utf-8"), encoding="utf-8")
                generated.append({"id": len(generated) + 1, "smiles": smi,
                                  "num_atoms": m.GetNumAtoms(),
                                  "sdf": str(dest), "trajectory": str(traj_dest) if traj_dest else None})
    if not generated:
        raise RuntimeError("FlashDiff 未生成有效分子，请调整参数重试")
    return {"molecules": generated, "pocket_center": pocket_center}
