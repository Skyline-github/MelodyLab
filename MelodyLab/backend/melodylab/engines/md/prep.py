"""复合物制备：蛋白选链清洗 → pdb2gmx → 并入刚体配体 → 盒子/溶剂化/加离子。

产物（work/ 目录内）：
    protein_clean.pdb → protein.gro/topol.top → complex.gro/topol_fix.top
    → boxed.gro → solv.gro → complex_ionized.gro (topol.top 最终版)
"""
from __future__ import annotations

import os
import re
import shutil
from pathlib import Path

from .gmx import run_gmx
from . import ligand as LG


def filter_protein_chains(pdb_in: Path, pdb_out: Path, chains: str | None = None) -> dict:
    """保留 ATOM 记录（可选链），丢弃水/杂原子/氢；返回统计。"""
    keep_chains = {c.strip().upper() for c in chains.split(",") if c.strip()} if chains else None
    out_lines, n_atoms, seen_chains = [], 0, set()
    for line in pdb_in.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("TER", "END")):
            continue
        if not line.startswith("ATOM"):
            continue
        altloc = line[16]
        if altloc not in (" ", "A"):
            continue
        element = line[76:78].strip() if len(line) >= 78 else ""
        if element == "H":
            continue
        chain = line[21].strip().upper()
        if keep_chains and chain not in keep_chains:
            continue
        seen_chains.add(chain)
        out_lines.append(line)
        n_atoms += 1
    if n_atoms < 20:
        raise ValueError(f"蛋白有效重原子过少（{n_atoms}），请检查输入 PDB")
    content = "\n".join(out_lines) + "\nEND\n"
    pdb_out.write_text(content, encoding="utf-8", newline="\n")
    return {"atoms": n_atoms, "chains": sorted(seen_chains)}


BACKBONE = ("N", "CA", "C", "O")

# 20 种标准残基的完整重原子名单（OXT 等端基除外）
RESIDUE_ATOMS = {
    "ALA": "N CA C O CB", "ARG": "N CA C O CB CG CD NE CZ NH1 NH2",
    "ASN": "N CA C O CB CG OD1 ND2", "ASP": "N CA C O CB CG OD1 OD2",
    "CYS": "N CA C O CB SG", "GLN": "N CA C O CB CG CD OE1 NE2",
    "GLU": "N CA C O CB CG CD OE1 OE2", "GLY": "N CA C O",
    "HIS": "N CA C O CB CG ND1 CD2 CE1 NE2", "ILE": "N CA C O CB CG1 CG2 CD1",
    "LEU": "N CA C O CB CG CD1 CD2", "LYS": "N CA C O CB CG CD CE NZ",
    "MET": "N CA C O CB CG SD CE", "PHE": "N CA C O CB CG CD1 CD2 CE1 CE2 CZ",
    "PRO": "N CA C O CB CG CD", "SER": "N CA C O CB OG",
    "THR": "N CA C O CB OG1 CG2",
    "TRP": "N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
    "TYR": "N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH", "VAL": "N CA C O CB CG1 CG2",
}


def repair_protein(pdb_in: Path, pdb_out: Path) -> dict:
    """修蛋白：剔除骨架或侧链不完整的标准残基，链内断裂处插 TER（pdb2gmx 按段补端基）。

    另按空间连续性分段：相邻残基 CA 距离 > 1.5 nm（晶体伪影/离群片段）视为断裂，
    每条链只保留残基最多的主片段——离群片段会让门控 RMSD 永远对不上参考。
    """
    residues: dict[tuple, list[str]] = {}
    order: list[tuple] = []
    for line in pdb_in.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM"):
            continue
        key = (line[21].strip(), int(line[22:26]), line[26].strip())
        if key not in residues:
            residues[key] = []
            order.append(key)
        residues[key].append(line)

    def ca_coord(lines: list[str]) -> tuple[float, float, float] | None:
        for l in lines:
            if l[12:16].strip() == "CA":
                return (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        return None

    complete = []  # [(key, lines, ca)]
    dropped = []
    for key in order:
        lines = residues[key]
        resname = lines[0][17:20].strip().upper()
        names = {l[12:16].strip() for l in lines}
        expected = RESIDUE_ATOMS.get(resname)
        if expected is not None:
            missing = [n for n in expected.split() if n not in names]
            if missing:
                dropped.append(f"{key[0]}{key[1]}{resname}(-{','.join(missing)})")
                continue
        complete.append((key, lines, ca_coord(lines)))

    # 按链分段：相邻残基 CA 距离 > 15 Å（=1.5 nm，晶体伪影/离群片段）视为断裂。
    # 注意 PDB 坐标单位是埃（Å），不是 nm。
    frags: dict[str, list[list[tuple]]] = {}
    prev = None  # (key, ca)
    for key, lines, ca in complete:
        ch = key[0]
        if prev is None or prev[0][0] != ch:
            frags.setdefault(ch, []).append([])
        else:
            far = (prev[1] is not None and ca is not None
                   and sum((a - b) ** 2 for a, b in zip(prev[1], ca)) > 15.0 ** 2)
            if far:
                frags[ch].append([])
        frags[ch][-1].append(key)
        prev = (key, ca)

    kept_keys: set[tuple] = set()
    dropped_frags = []
    for ch, fl in frags.items():
        best = max(fl, key=len)
        kept_keys.update(best)
        for f in fl:
            if f is not best:
                dropped_frags.append(f"链{ch} {f[0][1]}-{f[-1][1]}（{len(f)} 残基）")

    out: list[str] = []
    prev_key = None
    for key, lines, _ in complete:
        if key not in kept_keys:
            continue
        if prev_key is not None and prev_key[0] == key[0] and key[1] != prev_key[1] + 1:
            out.append("TER")
        out.extend(lines)
        prev_key = key
    out.append("TER")
    pdb_out.write_text("\n".join(out) + "\nEND\n", encoding="utf-8", newline="\n")
    n_res = len(kept_keys)
    return {"dropped_residues": dropped, "dropped_fragments": dropped_frags,
            "n_residues": n_res}


def pdb2gmx(protein_pdb: Path, work: Path, forcefield: str, water: str) -> dict:
    """生成蛋白拓扑与坐标（非交互：-ff/-water 显式指定，stdin=DEVNULL）。"""
    run_gmx(["pdb2gmx", "-f", str(protein_pdb), "-o", "protein.gro", "-p", "topol.top",
             "-i", "posre.itp", "-ff", forcefield, "-water", water, "-ignh"], cwd=work)
    gro, top = work / "protein.gro", work / "topol.top"
    if not gro.is_file() or not top.is_file():
        raise RuntimeError("pdb2gmx 未生成拓扑，请检查蛋白结构（缺失残基/非标准残基）")
    natoms = int(gro.read_text(encoding="utf-8").splitlines()[1].strip())
    return {"natoms": natoms, "gro": str(gro), "top": str(top)}


def merge_topology(work: Path, lig_itp: Path, lig_natoms: int, mode: str = "rigid") -> Path:
    """把 ligand.itp 插入 topol.top，[molecules] 加 LIG。

    - rigid 模式：插在链拓扑之后、水拓扑之前；
    - openff 模式：itp 含 [atomtypes]，必须紧跟力场 include（早于任何 moleculetype）。
    """
    top = work / "topol.top"
    text = top.read_text(encoding="utf-8")
    lines = text.splitlines()
    out, inserted_itp, inserted_mol = [], False, False
    for line in lines:
        s = line.strip().lower()
        if not inserted_itp:
            if mode == "openff" and 'forcefield.itp"' in s:
                out.append(line)
                out.append('; Include ligand topology (MelodyLab OpenFF flexible route)')
                out.append(f'#include "{os.path.basename(lig_itp)}"')
                inserted_itp = True
                continue
            if mode != "openff" and (s.startswith('#include "amber99sb.ff/tip3p') or s == "[ system ]"):
                out.append('; Include ligand topology (MelodyLab rigid-ligand route)')
                out.append(f'#include "{os.path.basename(lig_itp)}"')
                inserted_itp = True
        out.append(line)
        if not inserted_mol and s.startswith("sol"):
            m = re.match(r"^(\s*)(\S+)", line)
            indent = m.group(1) if m else ""
            out.append(f"{indent}LIG                1")
            inserted_mol = True
    if not inserted_itp:
        raise RuntimeError("topol.top 未找到拓扑插入点，无法插入配体拓扑")
    if not inserted_mol:
        out.append("LIG                1")
    fixed = work / "topol_fix.top"
    fixed.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")
    return fixed


def merge_coordinates(protein_gro: Path, ligand_gro: Path, out_gro: Path) -> int:
    """拼接坐标文件（都无速度场；盒子取蛋白的）。"""
    plines = protein_gro.read_text(encoding="utf-8").splitlines()
    llines = ligand_gro.read_text(encoding="utf-8").splitlines()
    n_p, n_l = int(plines[1].strip()), int(llines[1].strip())
    body = plines[2:2 + n_p] + llines[2:2 + n_l]
    out_gro.write_text(
        "\n".join([plines[0], str(n_p + n_l), *body, plines[2 + n_p]]) + "\n",
        encoding="utf-8", newline="\n")
    return n_p + n_l


def build_system(work: Path, params: dict, control) -> dict:
    """完整体系构建：返回 {gro, top, natoms, n_water...}。"""
    ff = params.get("forcefield") or "amber99sb"
    water = params.get("water") or "tip3p"
    salt = float(params.get("salt_molar", 0.15))
    lig_gro = work / "ligand.gro"
    lig_itp = work / "ligand.itp"

    control.checkpoint(progress=8, phase="prep", message=f"pdb2gmx 蛋白拓扑（{ff}/{water}）…")
    prot = pdb2gmx(work / "protein_clean.pdb", work, ff, water)

    control.checkpoint(progress=14, phase="prep", message="并入刚体配体拓扑…")
    top_fix = merge_topology(work, lig_itp, 0, params.get("ligand_mode", "rigid"))
    complex_gro = work / "complex.gro"
    merge_coordinates(Path(prot["gro"]), lig_gro, complex_gro)

    control.checkpoint(progress=18, phase="prep", message="定义盒子与溶剂化…")
    # 立方盒：大蛋白在十二面体的三斜表示里必然有原子戳出盒外，每次 grompp
    # 会把它们逐个卷绕（整数盒向量平移），轨迹与参考的包装状态永远对不齐。
    # 立方盒无斜切边，整分子 + -ur compact 后蛋白可完整入胞，卷绕成为空操作。
    run_gmx(["editconf", "-f", "complex.gro", "-o", "boxed.gro", "-c", "-d", "1.0",
             "-bt", "cubic"], cwd=work)
    run_gmx(["solvate", "-cp", "boxed.gro", "-cs", "spc216.gro", "-o", "solv.gro",
             "-p", "topol_fix.top"], cwd=work)
    # 蛋白在三斜盒表示里必然跨越边界，grompp 会逐原子卷绕而撕开链——
    # 先整分子 + 装入紧凑原胞，之后所有 grompp 的卷绕成为空操作
    # 注意 -s 必须用溶剂化后的 solv.gro（含全部原子）；用 boxed.gro 会丢掉所有水
    run_gmx(["trjconv", "-s", "solv.gro", "-f", "solv.gro", "-o", "solv_c.gro",
             "-pbc", "mol", "-ur", "compact"], cwd=work, stdin_text="System\n")
    (work / "solv_c.gro").replace(work / "solv.gro")
    from .ligand import wrap_solute_atoms

    wrap_solute_atoms(work / "solv.gro")

    control.checkpoint(progress=24, phase="prep", message="加离子与中和…")
    from .mdps import write_mdp
    write_mdp(work / "ions.mdp", "ions")
    # grompp 用最终 topol：solvate 会更新 topol_fix.top（SOL 计数），这里直接引用
    run_gmx(["grompp", "-f", "ions.mdp", "-c", "solv.gro", "-p", "topol_fix.top",
             "-o", "ions.tpr", "-maxwarn", "5"], cwd=work)
    run_gmx(["genion", "-s", "ions.tpr", "-o", "ionized.gro", "-p", "topol_fix.top",
             "-pname", "NA", "-nname", "CL", "-neutral", "-conc", f"{salt}",
             "-n", str(_make_sol_index(work))], cwd=work, stdin_text="SOL\n")
    final_top = work / "topol_fix.top"
    text = final_top.read_text(encoding="utf-8")
    final_top = work / "topol.top"
    final_top.write_text(text, encoding="utf-8", newline="\n")
    natoms = int((work / "ionized.gro").read_text(encoding="utf-8").splitlines()[1].strip())
    return {"gro": "ionized.gro", "top": final_top.name, "natoms": natoms,
            "forcefield": ff, "water": water, "salt_molar": salt}


def _make_sol_index(work: Path) -> Path:
    """生成只含 SOL 组的 index 文件供 genion 用（-n 需要显式文件避免交互选组歧义）。"""
    ndx = work / "sol.ndx"
    if ndx.is_file():
        return ndx
    run_gmx(["make_ndx", "-f", "solv.gro", "-o", "sol.ndx"], cwd=work, stdin_text="q\n")
    return ndx
