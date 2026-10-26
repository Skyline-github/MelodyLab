"""三元复合物 MD 模式：输入构建 + 三体接触分析。

与常规 MD 的差异（其余流水线完全复用）：
- 输入来自 ternary 引擎的报告：complex_{cid}.pdb（链 A=POI、链 B=E3、
  HETATM UNL=PROTAC）+ 报告行的 SMILES；
- 配体 SDF 不重新嵌构象：SMILES 连接性 + UNL 晶体/预测坐标逐原子赋值
  （DeepTernary 输出保持输入 lig.pdb 的原子顺序），H 用 addCoords 补位——
  配体因此恰好落在蛋白-蛋白界面上；
- 分析增加三体接触：LIG↔POI、LIG↔E3、POI↔E3 随时间的接触原子对数
  （gmx mindist -on），这是 PROTAC 三元复合物的核心稳定性指标。
"""
from __future__ import annotations

import numpy as np
import subprocess
from pathlib import Path

from .gmx import run_gmx, xvg_series


def _read_heavy_atoms(complex_pdb: Path):
    """读复合物重原子：[(line_index, xyz(Å), group)]，group ∈ lig/A/B。

    氢判定用元素列 + 名字首段（PDB 1HB/2HG 风格的氢元素列常缺）。
    """
    import re

    atoms = []
    for li, line in enumerate(
            Path(complex_pdb).read_text(encoding="utf-8", errors="replace").splitlines()):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        el = line[76:78].strip().upper()
        name = line[12:16].strip()
        if not el:
            el = re.sub(r"^\d+", "", name)[:1].upper()
        if el == "H" or re.match(r"^\d*H", name):
            continue
        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        group = "lig" if line.startswith("HETATM") else (
            line[21] if line[21] in ("A", "B") else None)
        if group is not None:
            atoms.append((li, (x, y, z), group))
    return atoms


def preflight_clash(complex_pdb: Path) -> dict:
    """MD 前置碰撞预检：配体-蛋白与蛋白-蛋白的最近重原子距离（nm）。

    DeepTernary 预测复合物普遍带亚埃级穿插（模型伪影），单精度 EM 在相应
    巨力下无法松弛——超阈值者先走 declash_complex 几何退火。
    """
    atoms = _read_heavy_atoms(complex_pdb)
    lig = [p for _, p, g in atoms if g == "lig"]
    cha = [p for _, p, g in atoms if g == "A"]
    chb = [p for _, p, g in atoms if g == "B"]

    from scipy.spatial import cKDTree

    def _min(a, b):
        if not a or not b:
            return 1e9
        d, _ = cKDTree(np.asarray(b)).query(np.asarray(a), k=1)
        return float(np.min(d)) * 0.1  # Å → nm

    lig_clash = 0
    if lig and (cha or chb):
        t_p = cKDTree(np.asarray(cha + chb))
        lig_clash = int((cKDTree(np.asarray(lig)).query(t_p.data, k=1)[0] < 2.0).sum())
    return {"min_lig_prot": _min(lig, cha + chb),
            "min_prot_prot": _min(cha, chb),
            "n_lig_clash": lig_clash}


def declash_complex(complex_pdb: Path, out_pdb: Path, target_nm: float = 0.20,
                    max_iter: int = 800, step: float = 0.03) -> dict:
    """几何解穿插：跨分子重原子对（lig↔蛋白、A↔B）迭代推到 target_nm 以上。

    目标取 0.20 nm：氢比重原子多伸出 ~1 Å，0.16 nm 的重原子距仍会让
    pdb2gmx 生成的氢相互重叠（LJ 1e19）。只动跨分子对；配体-蛋白对只推
    配体（蛋白几何保持，配体畸变由约束松弛恢复）。EM 接力恢复理想键长。
    """
    atoms = _read_heavy_atoms(complex_pdb)
    lines = Path(complex_pdb).read_text(encoding="utf-8", errors="replace").splitlines()
    idxs = [li for li, _, _ in atoms]
    coords = np.asarray([p for _, p, _ in atoms], dtype=float)
    groups = [g for _, _, g in atoms]

    cross = {frozenset(("lig", "A")), frozenset(("lig", "B")), frozenset(("A", "B"))}

    def is_cross(i, j):
        return frozenset((groups[i], groups[j])) in cross

    from scipy.spatial import cKDTree

    target = target_nm * 10  # Å
    n_moved = 0
    it = 0
    for it in range(max_iter):
        pairs = cKDTree(coords).query_pairs(target, output_type="ndarray")
        bad = [(i, j) for i, j in pairs if is_cross(i, j)]
        if not bad:
            break
        n_moved += len(bad)
        for i, j in bad:
            v = coords[i] - coords[j]
            d = float(np.linalg.norm(v))
            v = v / d if d > 1e-6 else np.asarray([1.0, 0.0, 0.0])
            gi, gj = groups[i], groups[j]
            if "lig" in (gi, gj):
                # 配体-蛋白：只推配体（蛋白几何保持；配体畸变由约束重嵌恢复）
                k = i if gi == "lig" else j
                coords[k] += v * step * (1 if k == i else -1)
            else:
                coords[i] += v * step * 0.5
                coords[j] -= v * step * 0.5
    for k, li in enumerate(idxs):
        x, y, z = coords[k]
        line = lines[li]
        lines[li] = line[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[54:]
    Path(out_pdb).write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    after = preflight_clash(out_pdb)
    return {"iterations": it + 1 if it else 0, "pairs_moved": n_moved, **after}


def _reembed_ligand(mol, coords, k: float = 50.0):
    """MMFF 距离约束松弛：所有重原子对约束到 declash 距离 ±0.2 Å 后最小化。

    declash 的推挤会压坏分子内几何（重碰撞分子可出现近零距的分子内对，
    对应 LJ 项能量 1e19 级）。约束松弛在保持界面布局的同时恢复合法构象。
    成功返回新坐标数组，失败（MMFF 不可用/结果偏离过大）返回 None。
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Geometry import Point3D

    n = mol.GetNumAtoms()
    try:
        work = Chem.Mol(mol)
        conf = Chem.Conformer(n)
        for i in range(n):
            conf.SetAtomPosition(i, Point3D(*coords[i]))
        work.AddConformer(conf, assignId=True)
        props = AllChem.MMFFGetMoleculeProperties(work)
        if props is None:
            return None
        ff = AllChem.MMFFGetMoleculeForceField(work, props)
        for i in range(n):
            for j in range(i + 1, n):
                d = float(np.linalg.norm(np.asarray(coords[i]) - np.asarray(coords[j])))
                tol = 0.2 if d < 3.0 else 0.35
                # 下限 1.2 Å：不把 declash 压坏的距离锚死（损伤对被强制拉开）
                lo = max(d - tol, 1.2)
                hi = max(d + tol, lo + 0.2)
                ff.AddDistanceConstraint(i, j, lo, hi, k)
        ff.Initialize()
        ff.Minimize(maxIts=1500)
        new = work.GetConformer(0).GetPositions()
        rmsd = float(np.sqrt(((new - np.asarray(coords)) ** 2).sum(1).mean()))
        if rmsd > 2.0:
            return None
        return new
    except Exception:  # noqa: BLE001 —— 任何几何病态都回退原始坐标
        return None


def make_ligand_sdf_from_complex(complex_pdb: Path, smiles: str, out_sdf: Path) -> dict:
    """从复合物 PDB 的 UNL 坐标 + SMILES 连接性构建配体 SDF（带氢、含 3D 坐标）。"""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    from rdkit.Geometry import Point3D

    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("SMILES 无法解析")
    coords = []
    for line in Path(complex_pdb).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("HETATM") and line[17:20].strip() == "UNL":
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if len(coords) != mol.GetNumAtoms():
        raise ValueError(
            f"复合物配体原子数 {len(coords)} ≠ SMILES 重原子数 {mol.GetNumAtoms()}"
            "（DeepTernary 输出与查询 SMILES 不对应，无法构建界面构象）")
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, (x, y, z) in enumerate(coords):
        conf.SetAtomPosition(i, Point3D(x, y, z))
    mol.AddConformer(conf, assignId=True)
    # declash 推挤可能压坏分子内几何（近零距 1-4 对 → LJ 爆炸），约束重嵌修复
    reembedded = _reembed_ligand(mol, np.asarray(coords, dtype=float))
    reembed_rmsd = None
    if reembedded is not None:
        reembed_rmsd = float(np.sqrt(((reembedded - np.asarray(coords)) ** 2).sum(1).mean()))
        for i in range(mol.GetNumAtoms()):
            conf.SetAtomPosition(i, Point3D(*reembedded[i]))
    molh = Chem.AddHs(mol, addCoords=True)
    # AddHs 的理想方向在退化角上会把 H 放到别的原子顶上，先无约束 MMFF
    # 快速松弛（H 归位、重原子几乎不动），再做近零距校验
    try:
        hprops = AllChem.MMFFGetMoleculeProperties(molh)
        if hprops is not None:
            hff = AllChem.MMFFGetMoleculeForceField(molh, hprops)
            hff.Initialize()
            hff.Minimize(maxIts=300)
    except Exception:  # noqa: BLE001 —— MMFF 失败则仅依赖下面的校验
        pass
    # 校验：任意原子对 < 0.05 nm 即拒绝（防拓扑 LJ 奇异）。
    # RDKit 坐标单位是 Å，阈值 (0.05 nm)² = 0.25 Å²
    posh = molh.GetConformer().GetPositions()
    d2 = ((posh[:, None, :] - posh[None, :, :]) ** 2).sum(-1)
    np.fill_diagonal(d2, np.inf)
    if float(d2.min()) < 0.25:
        raise ValueError("配体构象存在近零距原子对（加氢后），无法参数化——"
                         "该复合物姿态几何病态，请更换分子")
    w = Chem.SDWriter(str(out_sdf))
    w.write(molh)
    w.close()
    return {"natoms_heavy": mol.GetNumAtoms(), "natoms": molh.GetNumAtoms(),
            "reembed_rmsd": round(reembed_rmsd, 3) if reembed_rmsd is not None else None}


def _parse_moleculetypes(top_text: str) -> dict[str, int]:
    """解析拓扑文本中各 [moleculetype] 的原子数（数据行计数）。"""
    counts: dict[str, int] = {}
    cur_mtype, in_atoms, cur_count = None, False, 0

    def _flush():
        nonlocal cur_mtype, cur_count
        if cur_mtype is not None:
            counts.setdefault(cur_mtype, cur_count)
        cur_mtype, cur_count = None, 0

    for line in top_text.splitlines():
        s = line.split(";")[0].strip()
        if s.startswith("["):
            sec = s.strip("[] ").strip().lower()
            if sec == "moleculetype":
                _flush()
                in_atoms = False
            elif sec == "atoms":
                in_atoms = True
            else:
                _flush()
                in_atoms = False
            continue
        if not s or s.startswith("#"):
            continue
        if in_atoms:
            cur_count += 1
        elif cur_mtype is None and len(s.split()) >= 2:
            cur_mtype = s.split()[0]
    _flush()
    return counts


def build_ternary_index(work: Path, out_ndx: Path = None) -> Path:
    """构建三元复合物分组索引：Chain_A / Chain_B / LIGX（按拓扑原子顺序切区间）。

    原子顺序确定：pdb2gmx 按输入顺序输出蛋白，merge_coordinates 把配体拼在
    蛋白之后，溶剂化/加离子的水与离子在最后。链分组按 pdb2gmx 命名
    （Protein_chain_A* / Protein_chain_B*；断裂段为 _A2、_A3… 同属链 A）。
    """
    topol = work / "topol.top"
    text = Path(topol).read_text(encoding="utf-8", errors="replace")
    mtype_natoms = _parse_moleculetypes(text)
    # 链拓扑在 pdb2gmx 写出的独立 itp 里（#include "topol_Protein_chain_*.itp"）
    import re as _re

    for inc in _re.findall(r'#include\s+"([^"]+)"', text):
        if inc.endswith(".ff/") or "/" in inc or "\\" in inc:
            continue  # 系统力场文件
        inc_path = work / inc
        if inc_path.is_file():
            mtype_natoms.update(_parse_moleculetypes(
                inc_path.read_text(encoding="utf-8", errors="replace")))
    lig_itp = work / "ligand.itp"
    if lig_itp.is_file():
        mtype_natoms.update(_parse_moleculetypes(
            lig_itp.read_text(encoding="utf-8", errors="replace")))

    in_mols = False
    ranges: dict[str, list[tuple[int, int]]] = {}
    pos = 1
    for line in text.splitlines():
        s = line.split(";")[0].strip()
        if s.startswith("["):
            in_mols = s.strip("[] ").strip().lower() == "molecules"
            continue
        if not in_mols or not s or s.startswith("#"):
            continue
        parts = s.split()
        name, nmol = parts[0], int(parts[1]) if len(parts) > 1 else 1
        n = mtype_natoms.get(name)
        if n is None:
            continue
        for _ in range(nmol):
            ranges.setdefault(name, []).append((pos, pos + n - 1))
            pos += n

    import re

    a_ranges = sorted((r for name, rs in ranges.items()
                       if re.match(r"protein_chain_a\d*$", name.lower())
                       for r in rs), key=lambda r: r[0])
    b_ranges = sorted((r for name, rs in ranges.items()
                       if re.match(r"protein_chain_b\d*$", name.lower())
                       for r in rs), key=lambda r: r[0])
    lig_ranges = [(a, b) for name, rs in ranges.items()
                  if name.upper() == "LIG" for a, b in rs]
    if not a_ranges or not b_ranges:
        raise RuntimeError(
            "拓扑中未找到两条蛋白链（Protein_chain_A*/B*），三元复合物模式要求 "
            "POI(chain A) 与 E3(chain B)。实际: "
            f"{sorted(ranges)}")
    if not lig_ranges:
        raise RuntimeError("拓扑中未找到 LIG 分子")

    chain_a = [a for r in a_ranges for a in range(r[0], r[1] + 1)]
    chain_b = [a for r in b_ranges for a in range(r[0], r[1] + 1)]
    lig = [a for r in lig_ranges for a in range(r[0], r[1] + 1)]

    ndx = Path(out_ndx) if out_ndx else work / "ternary.ndx"
    with open(ndx, "w", encoding="utf-8", newline="\n") as f:
        for title, atoms in (("Chain_A", chain_a), ("Chain_B", chain_b), ("LIGX", lig)):
            f.write(f"[ {title} ]\n")
            for i in range(0, len(atoms), 15):
                f.write(" ".join(str(a) for a in atoms[i:i + 15]) + "\n")
    return ndx


def ternary_contacts(work: Path, tpr: str, xtc: str) -> dict | None:
    """三体接触曲线：各对基团每帧接触原子对数（< 0.6 nm）。返回 None 表示不可用。"""
    try:
        ndx = work / "ternary.ndx"
        if not ndx.is_file():
            build_ternary_index(work, ndx)
        pairs = [("lig_A", "LIGX", "Chain_A"), ("lig_B", "LIGX", "Chain_B"),
                 ("A_B", "Chain_A", "Chain_B")]
        out: dict = {}
        for key, g1, g2 in pairs:
            xvg = work / f"contact_{key}.xvg"
            run_gmx(["mindist", "-f", xtc, "-s", tpr, "-n", str(ndx),
                     "-on", str(xvg), "-d", "0.6"],
                    cwd=work, stdin_text=f"{g1}\n{g2}\n")
            s = xvg_series(xvg)
            if s and s.get("time"):
                out[key] = s
                from .analyze import _png

                _png(xvg, work / f"contact_{key}.png",
                     f"Contacts {g1} - {g2} (< 0.6 nm)", "Number of contacts")
        return out or None
    except Exception:  # noqa: BLE001 —— 接触分析失败不阻塞主报告
        return None
