"""配体参数化（兼容路线）：刚体配体 + amber99sb 类型 + Gasteiger 电荷。

设计取舍（见 docs/二期设计文档.md §1 的 RDKit 兼容回退路线）：
- 跳过 GAFF2/AM1-BCC（需 AmberTools，Windows 不可用）；
- 配体内部几何用「刚体」表达：选 3 个重原子作锚，constraints 三角固定，
  其余原子做成 vsite3/3fd 虚拟站点（纯几何，无任何力常数表依赖），
  因此不会因参数近似导致模拟发散；
- LJ/电荷走 amber99sb 原子类型映射 + Gasteiger 电荷（归一化到形式电荷）。
局限：配体无内部柔性，适用「对接姿态稳定性预检」；生产级柔性模拟
未来经 WSL AmberTools/ACPYPE 路线扩展（params.parameterization 保留钩子）。
"""
from __future__ import annotations

import os
from pathlib import Path

ELEMENT_MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
                "P": 30.974, "S": 32.06, "Cl": 35.45, "Br": 79.904, "I": 126.90}


def load_ligand(sdf_path: Path | None, smiles: str | None):
    """优先 SDF（对接姿态），回退 SMILES + ETKDG 3D。返回带 3D 构象的加氢分子。"""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = None
    if sdf_path is not None and Path(sdf_path).is_file():
        mol = next((m for m in Chem.SDMolSupplier(str(sdf_path), removeHs=False) if m is not None), None)
    if mol is None:
        mol = Chem.MolFromSmiles(smiles or "")
        if mol is None:
            raise ValueError("无法解析配体（SDF/SMILES 均无效）")
    if not any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
        mol = Chem.AddHs(mol)
    if mol.GetNumConformers() == 0:
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        if AllChem.EmbedMolecule(mol, params) == -1:
            raise ValueError("配体 3D 构象生成失败")
    AllChem.MMFFOptimizeMolecule(mol)  # 仅优化坐标，不影响力场选择
    return mol


def amber_type(atom) -> str:
    """RDKit 原子 → amber99sb 原子类型（LJ/质量来源；刚体方案无需键合参数）。"""
    sym = atom.GetSymbol()
    hyb = str(atom.GetHybridization())
    aromatic = atom.GetIsAromatic()
    if sym == "C":
        if aromatic:
            return "CA"
        # 任一双键指向 N/O/S 视为羰基/亚胺碳
        for b in atom.GetBonds():
            if b.GetBondTypeAsDouble() > 1.5:
                other = b.GetOtherAtom(atom)
                if other.GetSymbol() in ("O", "N", "S"):
                    return "C"
        if hyb == "SP2":
            return "CC"  # 烯烃/环外双键碳
        return "CT"
    if sym == "N":
        if aromatic:
            return "NA" if any(x.GetSymbol() == "H" for x in atom.GetNeighbors()) else "N*"
        if hyb == "SP3":
            return "N3"
        if hyb == "SP":
            return "N2"
        return "N"  # 酰胺/其他 sp2
    if sym == "O":
        if hyb == "SP2" and atom.GetDegree() == 1:
            return "O2"  # 羰基氧
        if aromatic or hyb == "SP2":
            return "OS"
        return "OH" if atom.GetDegree() == 1 else "OS"
    if sym == "S":
        return "SH" if (atom.GetDegree() == 1 and any(n.GetSymbol() == "H" for n in atom.GetNeighbors())) else "S"
    if sym == "H":
        parent = atom.GetNeighbors()[0] if atom.GetNeighbors() else None
        if parent is not None and parent.GetSymbol() in ("N", "O", "S"):
            return "H"
        return "HA" if (parent is not None and parent.GetIsAromatic()) else "HC"
    return {"F": "F", "CL": "Cl", "BR": "Br", "I": "I", "P": "P"}.get(sym.upper(), "CT")


def gasteiger_charges(mol) -> list[float]:
    from rdkit.Chem import AllChem

    AllChem.ComputeGasteigerCharges(mol)
    q = [float(a.GetDoubleProp("_GasteigerCharge")) for a in mol.GetAtoms()]
    # Gasteiger 求和与形式电荷可能有 ~1e-2 漂移，精确归一到净形式电荷
    formal = sum(a.GetFormalCharge() for a in mol.GetAtoms())
    drift = sum(q) - formal
    if q:
        q = [x - drift / len(q) for x in q]
    return q


def _anchor_atoms(mol) -> tuple[int, int, int]:
    """选刚体锚：最远重原子对 + 使三角形面积最大的第三个重原子。"""
    from rdkit import Chem
    from rdkit.Chem import rdMolTransforms as T

    conf = mol.GetConformer()
    pos = conf.GetPositions()
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1] or list(range(mol.GetNumAtoms()))
    if len(heavy) < 3:
        heavy = list(range(mol.GetNumAtoms()))
    best, dmax = (heavy[0], heavy[1]), -1.0
    for i in range(len(heavy)):
        for j in range(i + 1, len(heavy)):
            d = sum((pos[heavy[i]][k] - pos[heavy[j]][k]) ** 2 for k in range(3))
            if d > dmax:
                dmax, best = d, (heavy[i], heavy[j])
    if len(heavy) < 3:
        return (best + (heavy[0],))[:3]
    i, j = best
    kmax, amax = heavy[0], -1.0
    for k in heavy:
        if k in (i, j):
            continue
        area = T.GetAngleDeg(conf, i, k, j)  # 用夹角近似：越接近 90° 面积越大
        score = abs(90.0 - area)
        if score < 0 or amax < 0 or score < amax:
            amax, kmax = score, k
    return i, j, kmax


def write_ligand_files(mol, lig_name: str, out_dir: Path, charge: int | None = None) -> dict:
    """生成 ligand.gro / ligand.itp。返回 {natoms, charge, gro, itp}。"""
    out_dir.mkdir(parents=True, exist_ok=True)
    n = mol.GetNumAtoms()
    if n < 2:
        raise ValueError("配体过小（<2 原子），无法构建 MD 拓扑")
    conf = mol.GetConformer()
    pos = conf.GetPositions() / 10.0  # Å → nm
    types = [amber_type(a) for a in mol.GetAtoms()]
    q = gasteiger_charges(mol)
    if charge is not None:
        drift = sum(q) - float(charge)
        q = [x - drift / n for x in q]
    anchors = _anchor_atoms(mol)
    real = set(anchors[: n if n < 3 else 3])

    def atom_name(i: int) -> str:
        return f"{types[i][0]}{i+1}"[:5]

    # ---- ligand.gro ----
    gro_lines = [f"MelodyLab rigid ligand {lig_name}"]
    gro_lines.append(f"{n:5d}")
    for i in range(n):
        gro_lines.append(
            f"{1 % 100000:5d}{f'{lig_name[:5]}':<5s}{atom_name(i):>5s}{(i+1) % 100000:5d}"
            f"{pos[i][0]:8.3f}{pos[i][1]:8.3f}{pos[i][2]:8.3f}"
        )
    box_min = pos.min(axis=0)
    gro_lines.append(f"{2.4:10.5f}{2.4:10.5f}{2.4:10.5f}")
    gro_path = out_dir / "ligand.gro"
    gro_path.write_text("\n".join(gro_lines) + "\n", encoding="utf-8", newline="\n")

    # ---- ligand.itp ----
    L = [f"; MelodyLab 刚性配体拓扑: {lig_name}  (anchor={list(anchors)})"]
    L.append("[ moleculetype ]\n; name  nrexcl\nLIG     1\n")
    L.append("[ atoms ]\n; nr  type  resnr residue atom cgnr charge mass")
    elements = [a.GetSymbol().capitalize() for a in mol.GetAtoms()]
    for i in range(n):
        mass = ELEMENT_MASS.get(elements[i], 12.011) if i in real else 0.0
        L.append(f"{i+1:6d} {types[i]:<6s} {1:5d} {lig_name[:5]:<5s} {atom_name(i):>5s} "
                 f"{1:5d} {q[i]:.6f} {mass:.4f}")
    real_idx = sorted(real)
    # 2 原子：1 条 constraint；>=3 原子：锚三角 3 条 constraints（刚体）
    if n == 2:
        L.append("\n[ constraints ]\n; i j funct length(nm)")
        L.append(f"1 2 1 {_dist(pos[0], pos[1]):.5f}")
    elif len(real_idx) >= 3:
        L.append("\n[ constraints ]\n; i j funct length(nm)")
        for a in range(3):
            for b in range(a + 1, 3):
                i1, i2 = real_idx[a], real_idx[b]
                L.append(f"{i1+1:5d} {i2+1:5d} 1 {_dist(pos[i1], pos[i2]):.5f}")
    vsites = [i for i in range(n) if i not in real]
    if vsites:
        # vsite3 列格式: site ai aj ak funct params
        # funct 4 = 3out: rs = ri + a*(rj-ri) + b*(rk-ri) + c*((rj-ri)×(rk-ri))
        # （funct 1 的线性组合只有 2 个系数，张成平面——非平面配体会被投影导致原子重叠）
        L.append("\n[ virtual_sites3 ]\n; site ai aj ak funct a b c")
        i0, i1, i2 = real_idx[0], real_idx[1], real_idx[2]
        for v in vsites:
            a, b, c = _vsite3out_coeff(pos[i0], pos[i1], pos[i2], pos[v])
            L.append(f"{v+1:5d} {i0+1:5d} {i1+1:5d} {i2+1:5d} 4 {a:.6f} {b:.6f} {c:.6f}")
        L.append("\n; vsite 与构造原子间的 exclusions（GROMACS 要求）")
        L.append("[ exclusions ]")
        for v in vsites:
            L.append(f"{v+1:5d} {i0+1:5d} {i1+1:5d} {i2+1:5d}")
    itp_path = out_dir / "ligand.itp"
    itp_path.write_text("\n".join(L) + "\n", encoding="utf-8", newline="\n")
    return {"natoms": n, "n_vsites": len(vsites), "charge": round(sum(q), 4),
            "gro": str(gro_path), "itp": str(itp_path),
            "types": {f"{types[i]}": types.count(types[i]) for i in range(n)}}


def _dist(a, b) -> float:
    return float(sum((a[k] - b[k]) ** 2 for k in range(3)) ** 0.5)


def refresh_geometry_from_gro(itp_path: Path, lig_gro: Path) -> dict:
    """EM 后几何刷新：按实际坐标重写 itp 的 constraints 长度与 vsite 3out 系数。

    EM（steep/cg）不执行约束求解，锚原子几何会漂移；不刷新则 NVT 首步
    LINCS 因巨大约束违反而崩溃。
    """
    lines = itp_path.read_text(encoding="utf-8").splitlines()
    natoms = 0
    in_atoms = False
    meta: dict[int, tuple] = {}  # nr -> (type, resname, atomname, charge)
    for l in lines:
        s = l.strip()
        if s.startswith("["):
            in_atoms = s.startswith("[ atoms")
            continue
        if in_atoms and s and not s.startswith(";"):
            p = s.split()
            meta[int(p[0])] = (p[1], p[3], p[4], float(p[6]))
            natoms = max(natoms, int(p[0]))
    gro_lines = Path(lig_gro).read_text(encoding="utf-8").splitlines()
    n_gro = int(gro_lines[1].strip())
    if n_gro != natoms:
        raise RuntimeError(f"LIG gro 原子数 {n_gro} 与 itp {natoms} 不一致")
    pos = {}
    for l in gro_lines[2:2 + n_gro]:
        pos[int(l[15:20])] = np_array([float(l[20:28]), float(l[28:36]), float(l[36:44])])

    masses = {}
    in_atoms = False
    for l in lines:
        s = l.strip()
        if s.startswith("["):
            in_atoms = s.startswith("[ atoms")
            continue
        if in_atoms and s and not s.startswith(";"):
            p = s.split()
            masses[int(p[0])] = float(p[7])
    real_idx = sorted(i - 1 for i, m in masses.items() if m > 0)
    if len(real_idx) < 3:
        raise RuntimeError("itp 中未识别出 3 个锚原子（质量>0）")
    i0, i1, i2 = (k + 1 for k in real_idx[:3])
    out, in_cons, in_vs = [], False, False
    for l in lines:
        s = l.strip()
        if s.startswith("["):
            in_cons = s.startswith("[ constraints")
            in_vs = s.startswith("[ virtual_sites3")
            out.append(l)
            continue
        if in_cons and s and not s.startswith(";"):
            p = s.split()
            a, b = int(p[0]), int(p[1])
            out.append(f"{a:5d} {b:5d} 1 {_dist(pos[a], pos[b]):.5f}")
            continue
        if in_vs and s and not s.startswith(";"):
            p = s.split()
            v = int(p[0])
            a, b, c = _vsite3out_coeff(pos[i0], pos[i1], pos[i2], pos[v])
            out.append(f"{v:5d} {i0:5d} {i1:5d} {i2:5d} 4 {a:.6f} {b:.6f} {c:.6f}")
            continue
        out.append(l)
    itp_path.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")
    return {"anchors": [i0, i1, i2], "n_vsites": natoms - len(real_idx)}


def extract_ligand_gro(src_gro: Path, dst_gro: Path) -> None:
    """从全体系 gro 中抽出 LIG 原子（序号重排 1..N），供几何刷新用。"""
    lines = Path(src_gro).read_text(encoding="utf-8").splitlines()
    n = int(lines[1].strip())
    lig = [l for l in lines[2:2 + n] if l[5:10].strip() == "LIG"]
    if not lig:
        raise RuntimeError(f"{src_gro} 中没有 LIG 原子")
    out = ["MelodyLab LIG extract", f"{len(lig):5d}"]
    for i, l in enumerate(lig, 1):
        out.append(f"{l[:15]}{i:5d}{l[20:]}")
    out.append(lines[2 + n])
    dst_gro.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")


def np_array(x):
    import numpy as np

    return np.array(x)


def wrap_solute_atoms(gro_path: Path) -> None:
    """把蛋白原子逐个规范化到 [0, L)（LIG 保持 unwrap 状态不动，保护 vsite 构造原子）。

    grompp 会对出盒原子做确定性逐个卷绕——大分子（拖无序尾巴的蛋白）在
    -pbc mol 质心入盒后总有少数原子越界，卷绕随微小热运动随机翻转，
    导致轨迹与参考结构的"包装状态"永远对不齐，RMSD 恒定虚高。
    与其对抗，不如抢先做同样的规范化：所有 gro 在 grompp 前统一成
    蛋白逐原子入胞的规范形，grompp 成为空操作，参考与轨迹天然一致。
    要求立方盒（prep 已用 -bt cubic）。动力学不受影响（PBC 感知）。
    """
    lines = Path(gro_path).read_text(encoding="utf-8").splitlines()
    n = int(lines[1])
    v = [float(x) for x in lines[2 + n].split()]
    box = [v[0], v[1], v[2]]  # 立方盒对角项即盒长
    out = lines[:2]
    changed = False
    for l in lines[2:2 + n]:
        resn = l[5:10].strip()
        if resn in ("SOL", "NA", "CL", "LIG"):
            out.append(l)
            continue
        x, y, z = float(l[20:28]), float(l[28:36]), float(l[36:44])
        wx = x - box[0] * floor_div(x, box[0])
        wy = y - box[1] * floor_div(y, box[1])
        wz = z - box[2] * floor_div(z, box[2])
        if (wx, wy, wz) != (x, y, z):
            changed = True
            out.append(f"{l[:20]}{wx:8.3f}{wy:8.3f}{wz:8.3f}{l[44:]}")
        else:
            out.append(l)
    if changed:
        out.append(lines[2 + n])
        Path(gro_path).write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")


def floor_div(x: float, L: float) -> int:
    import math

    return math.floor(x / L + 1e-9)


def unwrap_ligand_in_gro(gro_path: Path) -> None:
    """把 LIG 分子整体放进主原胞（grompp 会对出界原子逐个卷绕，
    而 vsite 重构不是 PBC 不变的——锚原子被卷到不同镜像会撕碎刚体配体）。"""
    import numpy as np

    lines = Path(gro_path).read_text(encoding="utf-8").splitlines()
    n = int(lines[1].strip())
    v = [float(x) for x in lines[2 + n].split()]
    if len(v) == 3:
        # GROMACS gro 盒子顺序: v1x v2y v3z 或 9 数(含交叉项)
        box = np.array(np.diag(v))
    else:
        # v1x v2y v3z v1y v1z v2x v2z v3x v3y
        box = np.array([[v[0], v[3], v[4]],
                        [v[5], v[1], v[6]],
                        [v[7], v[8], v[2]]])
    inv = np.linalg.inv(box)
    serials, coords, is_lig = [], [], []
    for l in lines[2:2 + n]:
        serials.append(int(l[15:20]))
        coords.append([float(l[20:28]), float(l[28:36]), float(l[36:44])])
        is_lig.append(l[5:10].strip() == "LIG")
    if not any(is_lig):
        return
    C = np_array(coords)
    lig_idx = [i for i, f in enumerate(is_lig) if f]
    anchor = C[lig_idx[0]]
    frac = (C[lig_idx] - anchor) @ inv
    frac -= np.round(frac)          # 最小镜像展开到同一周期镜像
    C[lig_idx] = anchor + frac @ box
    # 整数晶格平移把锚移入主胞：周期等价操作，物理位置不变。
    # （非整数平移是物理搬运，会把配体挪进蛋白——严禁）
    f0 = anchor @ inv
    shift = -np.floor(f0 + 1e-9) @ box
    C[lig_idx] += shift
    out = lines[:2]
    for i, l in enumerate(lines[2:2 + n]):
        if is_lig[i]:
            out.append(f"{l[:20]}{C[i][0]:8.3f}{C[i][1]:8.3f}{C[i][2]:8.3f}{l[44:]}")
        else:
            out.append(l)
    out.append(lines[2 + n])
    Path(gro_path).write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")


def _vsite3out_coeff(ri, rj, rk, rp) -> tuple[float, float, float]:
    """3out (funct 4): rp = ri + a*(rj-ri) + b*(rk-ri) + c*((rj-ri)×(rk-ri))。

    3 个基向量一般线性无关 → 任意 3D 位置可精确表达。
    """
    import numpy as np

    r12 = np.asarray(rj) - np.asarray(ri)
    r13 = np.asarray(rk) - np.asarray(ri)
    M = np.stack([r12, r13, np.cross(r12, r13)], axis=1)
    coef, *_ = np.linalg.lstsq(M, np.asarray(rp) - np.asarray(ri), rcond=None)
    return float(coef[0]), float(coef[1]), float(coef[2])
