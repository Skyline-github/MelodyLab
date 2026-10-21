"""OpenFF 柔性配体参数化（独立环境 CLI）。

在 mdff 独立 conda 环境中运行（依赖 openff-toolkit/openmm/parmed/rdkit，
与生产环境 tav 完全隔离）：
    python ligand_openff.py <ligand.sdf> <work_dir>

流程：RDKit 读 SDF（对接姿态）→ OpenFF Sage 2.x 力场（BCC 电荷）
→ OpenMM System → ParmEd → GROMACS ligand.itp（含 [atomtypes]）+ ligand.gro。

输出 ligand_meta.json：{natoms, charge, mode: "openff"}。
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import warnings

warnings.filterwarnings("ignore")


def build(sdf_path: Path, work_dir: Path) -> dict:
    from openff.toolkit import Molecule, ForceField
    from openmm import unit
    import parmed

    # 1) RDKit 读取（保留 H 与 3D 构象）
    from rdkit import Chem

    mol = next((m for m in Chem.SDMolSupplier(str(sdf_path), removeHs=False) if m is not None), None)
    if mol is None:
        raise ValueError(f"无法解析 SDF: {sdf_path}")
    work_dir.mkdir(parents=True, exist_ok=True)
    pdb_in = work_dir / "ligand_input.pdb"
    pdb_in.write_text(Chem.MolToPDBBlock(mol), encoding="utf-8", newline="\n")

    # 2) OpenFF 参数化（am1bcc 电荷：显式 NAGL 模型路径 + charge_from_molecules，
    #    规避 openff-nagl/-models/toolkit 三者版本间的注册不匹配）
    off_mol = Molecule.from_rdkit(mol, allow_undefined_stereo=True)
    models_dir = Path(__import__("openff.nagl_models", fromlist=["x"]).__file__).parent / "models" / "am1bcc"
    model = sorted(models_dir.glob("*.pt"))[-1]
    off_mol.assign_partial_charges(partial_charge_method=str(model))
    sage = ForceField("openff-2.2.0.offxml")
    system = sage.create_openmm_system(off_mol.to_topology(), charge_from_molecules=[off_mol])

    # 3) OpenMM → ParmEd → GROMACS
    omm_top = off_mol.to_topology()
    # OpenFF 0.19 的 Topology 缺链层级，ParmEd 需要：补一个默认 chain/residue 结构
    from openmm.app import Topology as _Topology
    _omm_top = _Topology()
    _omm_chain = _omm_top.addChain("A")
    _omm_res = _omm_top.addResidue("LIG", _omm_chain)
    from openmm.app import element as _elem
    _omm_atoms = []
    for atom in omm_top.atoms:
        _omm_atoms.append(_omm_top.addAtom(atom.name or atom.symbol,
                                           _elem.get_by_symbol(atom.symbol), _omm_res))
    for bond in off_mol.bonds:
        _omm_top.addBond(_omm_atoms[bond.atom1_index], _omm_atoms[bond.atom2_index])
    positions = unit.Quantity(off_mol.conformers[0], unit.angstrom)
    struct = parmed.openmm.load_topology(_omm_top, system)
    struct.coordinates = positions.value_in_unit(unit.angstrom)

    # 兜底：OpenMM 转换后个别键/角无类型参数（典型为 X-H 键），按实际几何补齐
    import numpy as _np
    for b in struct.bonds:
        if b.type is None:
            d = float(_np.linalg.norm(struct.coordinates[b.atom1.idx] -
                                      struct.coordinates[b.atom2.idx]))
            is_h = b.atom1.atomic_number == 1 or b.atom2.atomic_number == 1
            b.type = parmed.BondType(300.0 if is_h else 250.0, d)
            struct.bond_types.append(b.type)
            b.type.list = struct.bond_types
    for ang in struct.angles:
        if ang.type is None:
            v1 = struct.coordinates[ang.atom1.idx] - struct.coordinates[ang.atom2.idx]
            v2 = struct.coordinates[ang.atom3.idx] - struct.coordinates[ang.atom2.idx]
            cos = float(_np.dot(v1, v2) / (_np.linalg.norm(v1) * _np.linalg.norm(v2) + 1e-12))
            import math as _math
            theta = _math.degrees(_math.acos(max(-1.0, min(1.0, cos))))
            ang.type = parmed.AngleType(60.0, theta)
            struct.angle_types.append(ang.type)
            ang.type.list = struct.angle_types

    charge = round(sum(a.charge for a in struct.atoms), 4)
    for res in struct.residues:
        res.name = "LIG"          # 管线按残基名 LIG 建索引组
        res.number = 1
    for i, atom in enumerate(struct.atoms, 1):
        atom.name = f"{atom.type[0].upper()}{i}"[:4]
        atom.number = i

    top_path = work_dir / "ligand_full.top"
    struct.save(str(top_path), overwrite=True)

    gro_path = work_dir / "ligand.gro"
    struct.save(str(gro_path), overwrite=True)

    # 4) top → itp：剥掉 [ defaults ]/[ system ]/[ molecules ]（与 amber99sb 冲突）
    lines = top_path.read_text(encoding="utf-8", errors="replace").splitlines()
    drop_sections = ("defaults", "system", "molecules")
    out, skip = ["; MelodyLab OpenFF flexible ligand topology"], False
    for line in lines:
        s = line.strip()
        if s.startswith("["):
            skip = any(s.startswith(f"[{name}") or s.startswith(f"[ {name}") for name in drop_sections)
            if not skip:
                out.append(line)
            continue
        if not skip and s and not s.startswith(";"):
            out.append(line)

    # 5) 类型名加 OF_ 前缀：避免与 amber99sb 类型撞名（本 itp 插在力场之后会重定义蛋白类型）
    renamed, in_types = [], False
    type_map = {}
    for line in out:
        s = line.strip()
        if s.startswith("["):
            in_types = s.startswith("[ atomtypes")
            renamed.append(line)
            continue
        if in_types and s and not s.startswith(";"):
            parts = line.split()
            type_map[parts[0]] = "OF_" + parts[0]
            parts[0] = type_map[parts[0]]
            renamed.append("%-8s %5s %10s %8s %5s %10s %10s" % tuple(parts))
        elif s.startswith("[ atoms") or s.startswith("[ bonds") or s.startswith("[ pairs") \
                or s.startswith("[ angles") or s.startswith("[ dihedrals") or s.startswith("[ moleculetype"):
            renamed.append(line)
        elif s and not s.startswith(";"):
            parts = line.split()
            if len(parts) >= 2 and parts[1] in type_map:
                parts[1] = type_map[parts[1]]
                renamed.append("%6s %-8s %s" % (parts[0], parts[1], " ".join(parts[2:])))
            else:
                renamed.append(line)
        else:
            renamed.append(line)
    out = renamed

    itp_path = work_dir / "ligand.itp"
    itp_path.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")
    top_path.unlink(missing_ok=True)
    pdb_in.unlink(missing_ok=True)

    meta = {"mode": "openff", "natoms": len(struct.atoms), "charge": charge,
            "gro": str(gro_path), "itp": str(itp_path)}
    (work_dir / "ligand_meta.json").write_text(json.dumps(meta, ensure_ascii=False, indent=1),
                                               encoding="utf-8", newline="\n")
    return meta


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: ligand_openff.py <ligand.sdf> <work_dir>")
        sys.exit(2)
    meta = build(Path(sys.argv[1]), Path(sys.argv[2]))
    print(json.dumps(meta, ensure_ascii=False))
