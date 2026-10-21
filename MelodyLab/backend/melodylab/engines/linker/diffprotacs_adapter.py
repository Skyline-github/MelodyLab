"""DiffPROTACs Linker 生成适配器（复用 assets/vendor/DiffPROTACs）。"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np

from ...config import settings


def availability() -> dict:
    vendor = settings.vendor_dir("DiffPROTACs")
    ckpt = settings.checkpoint("diffprotacs")
    return {
        "engine": "diffprotacs",
        "vendor": str(vendor) if vendor else None,
        "checkpoint": str(ckpt) if ckpt else None,
        "available": bool(vendor and ckpt),
        "hint": None if (vendor and ckpt) else "运行 scripts/setup_assets.py 复制 DiffPROTACs 代码与权重（注意：仅限研究用途）",
    }


def _read_sdf_atoms(sdf_path: Path):
    from rdkit import Chem

    mol = Chem.MolFromMolFile(str(sdf_path), removeHs=True)
    if mol is None:
        raise ValueError(f"无法解析 SDF: {sdf_path}")
    conf = mol.GetConformer()
    elements, coords = [], []
    for i in range(mol.GetNumAtoms()):
        symbol = mol.GetAtomWithIdx(i).GetSymbol()
        if symbol == "*":
            continue
        pos = conf.GetAtomPosition(i)
        elements.append(symbol)
        coords.append([pos.x, pos.y, pos.z])
    return elements, np.array(coords)


def _xyz_to_sdf(xyz_path: Path, sdf_path: Path) -> None:
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds

    mol = Chem.MolFromXYZFile(str(xyz_path))
    if mol is None:
        raise ValueError(f"无法解析 XYZ: {xyz_path}")
    for attempt in ("full", "hueckel", "connectivity"):
        try:
            if attempt == "full":
                rdDetermineBonds.DetermineBonds(mol)
            elif attempt == "hueckel":
                mol = Chem.MolFromXYZFile(str(xyz_path))
                rdDetermineBonds.DetermineBonds(mol, useHueckel=True)
            else:
                mol = Chem.MolFromXYZFile(str(xyz_path))
                rdDetermineBonds.DetermineConnectivity(mol)
                try:
                    rdDetermineBonds.DetermineBondOrders(mol)
                except Exception:
                    pass
            break
        except Exception:
            continue
    writer = Chem.SDWriter(str(sdf_path))
    writer.write(mol)
    writer.close()


def generate_linkers(lig1: Path, lig2: Path, out_dir: Path, *,
                     linker_len: int = 12, n_samples: int = 5,
                     progress=None, control=None) -> dict:
    info = availability()
    if not info["available"]:
        raise RuntimeError(f"DiffPROTACs 资产未就绪：{info['hint']}")
    vendor = settings.vendor_dir("DiffPROTACs")
    out_dir.mkdir(parents=True, exist_ok=True)

    elems1, coords1 = _read_sdf_atoms(lig1)
    elems2, coords2 = _read_sdf_atoms(lig2)
    all_elements = elems1 + elems2
    all_coords = np.vstack([coords1, coords2])
    original_center = np.mean(all_coords, axis=0)

    tmp = Path(os.path.abspath(out_dir / "_tmp"))
    tmp.mkdir(parents=True, exist_ok=True)
    fragments_xyz = tmp / "fragments.xyz"
    with open(fragments_xyz, "w", encoding="utf-8") as f:
        for elem, (x, y, z) in zip(all_elements, all_coords):
            f.write(f"{elem} {x:.9f} {y:.9f} {z:.9f}\n")

    output_dir = tmp / "output"
    original_cwd = os.getcwd()
    try:
        # 清理可能的模块名冲突后导入
        if str(vendor) not in sys.path:
            sys.path.insert(0, str(vendor))
        for mod in [k for k in list(sys.modules) if k.split(".")[0] in ("utils", "sample", "scripts", "models", "datasets")]:
            del sys.modules[mod]
        os.chdir(vendor)
        import sample as diffprotacs_sample

        args = argparse.Namespace(
            model_path=str(settings.checkpoint("diffprotacs").resolve()),
            output_dir=str(output_dir), xyz_path=str(fragments_xyz),
            linker_len=str(linker_len), n_samples=n_samples, name="melodylab_linker",
            batch_size=256, num_workers=0, exp_name="", checkpoints="checkpoints",
            epochs=None, in_node_nf=9, diffusion_steps=500,
            diffusion_noise_schedule="polynomial_2", diffusion_loss_type="l2",
            n_layers=6, hidden_nf=128, ffn_embedding_dim=1024,
            attention_heads=32, tanh=False, coords_range=10.0,
            dropout=0, activation_dropout=0,
        )
        diffprotacs_sample.main(args)
    finally:
        os.chdir(original_cwd)

    gen_base = output_dir / "melodylab_linker"
    n_frag = len(all_elements)
    linkers = []
    for idx in range(n_samples):
        xyz = gen_base / f"{idx}_.xyz"
        if not xyz.is_file():
            continue
        elements, coords = [], []
        for line in xyz.read_text(encoding="utf-8").splitlines():
            parts = line.split()
            if len(parts) == 4:
                try:
                    elements.append(parts[0])
                    coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
                except ValueError:
                    continue
        if not elements:
            continue
        full_coords = np.array(coords) + original_center
        no = len(linkers) + 1
        for tag, elems_, coords_ in (("full", elements, full_coords),
                                     ("linker", elements[n_frag:], full_coords[n_frag:])):
            if len(elems_) == 0:
                continue
            xyz_out = tmp / f"{no}_{tag}.xyz"
            with open(xyz_out, "w", encoding="utf-8") as f:
                f.write(f"{len(elems_)}\n\n")
                for e, (x, y, z) in zip(elems_, coords_):
                    f.write(f"{e} {x:.9f} {y:.9f} {z:.9f}\n")
            _xyz_to_sdf(xyz_out, out_dir / f"{no}_{tag}.sdf")
        linkers.append({
            "id": no,
            "num_atoms_full": len(elements),
            "num_atoms_linker": len(elements) - n_frag,
            "full_sdf": str(out_dir / f"{no}_full.sdf"),
            "linker_sdf": str(out_dir / f"{no}_linker.sdf"),
        })
    if not linkers:
        raise RuntimeError("DiffPROTACs 未生成有效 linker，请调整 linker_len 或 n_samples")
    shutil_rmtree(tmp)
    return {"linkers": linkers}


def shutil_rmtree(p: Path) -> None:
    import shutil

    shutil.rmtree(p, ignore_errors=True)
