"""对接引擎适配：配体/受体准备 + Vina / Smina / Vina-GPU 运行与打分解析。"""
from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

from ...config import settings

CREATE_NO_WINDOW = getattr(subprocess, "CREATE_NO_WINDOW", 0)


# ---------------- 配体准备（RDKit + Meeko） ----------------
def prepare_ligand_pdbqt(smiles: str | None, sdf_path: Path | None, out_pdbqt: Path) -> None:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    if sdf_path is not None:
        mol = next((m for m in Chem.SDMolSupplier(str(sdf_path), removeHs=False) if m is not None), None)
    else:
        mol = Chem.MolFromSmiles(smiles or "")
    if mol is None:
        raise ValueError("无法解析配体（SMILES/SDF 无效）")
    if not any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
        mol = Chem.AddHs(mol)
    if mol.GetNumConformers() == 0:
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        if AllChem.EmbedMolecule(mol, params) == -1:
            raise ValueError(f"配体 3D 构象生成失败: {smiles or sdf_path}")
        AllChem.MMFFOptimizeMolecule(mol)
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    setups = MoleculePreparation().prepare(mol)
    for setup in setups:
        text, ok, err = PDBQTWriterLegacy.write_string(setup)
        if ok:
            out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
            # newline="\n"：Windows 默认会把 \n 翻成 CRLF，Vina-GPU 等旧解析器不容忍 \r
            out_pdbqt.write_text(text, encoding="utf-8", newline="\n")
            return
    raise ValueError(f"Meeko 配体转换失败: {err}")


# ---------------- 受体准备 ----------------
_VINA_CHARGE = {"C": -0.15, "N": -0.39, "O": -0.38, "S": -0.43, "P": 0.0,
                "F": -0.35, "CL": -0.35, "BR": -0.35, "I": -0.35}


def write_fallback_receptor_pdbqt(protein_path: Path, out_pdbqt: Path) -> int:
    """无 Meeko Polymer 依赖的刚性受体 PDBQT（非极性氢合并近似）。

    直接从重原子生成：保留全部重原子 + AD4 原子类型 + 简化电荷。
    对 Vina 打分足够稳健（DockScope 同样保留此兜底路线）。
    """
    from .pockets import VDW, parse_pdb_atoms

    atoms = parse_pdb_atoms(protein_path, keep_hetatm=True, remove_waters=True)
    if len(atoms) < 20:
        raise ValueError(f"受体原子过少: {protein_path}")
    lines = []
    for i, a in enumerate(atoms, 1):
        el = a["element"]
        ad_type = el if el in VDW else "A" if el == "C" else el.capitalize()
        x, y, z = a["xyz"]
        q = _VINA_CHARGE.get(el, 0.0)
        # 标准 PDBQT 列位：occ 55-60, temp 61-66, 电荷 71-76, 原子类型 78-79
        lines.append(
            f"ATOM  {i % 100000:5d} {el:>2}   LIG A   1    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"  1.00  0.00    {q:6.3f} {ad_type:<2s}"
        )
    lines.append("TER")
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    out_pdbqt.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return len(atoms)


def prepare_receptor(protein_path: Path, out_pdbqt: Path) -> int:
    """优先 Meeko Polymer（正规路线），失败回退刚性重原子近似。"""
    try:
        from meeko import MoleculePreparation, Polymer, PDBQTWriterLegacy, ResidueChemTemplates

        pdb = Path(protein_path).read_text(encoding="utf-8", errors="replace")
        templates = ResidueChemTemplates.create_from_defaults()
        try:
            polymer = Polymer.from_pdb_string(pdb, chem_templates=templates,
                                              mk_prep=MoleculePreparation(), allow_bad_res=True)
        except TypeError:
            polymer = Polymer.from_pdb_string(pdb, chem_templates=templates,
                                              mk_prep=MoleculePreparation())
        # meeko 0.8+: write_string_from_polymer 返回 (text, error)；成功时 error 为空串
        text, err = PDBQTWriterLegacy.write_string_from_polymer(polymer)
        if text and not err:
            out_pdbqt.parent.mkdir(parents=True, exist_ok=True)
            out_pdbqt.write_text(text, encoding="utf-8", newline="\n")
            return sum(1 for l in text.splitlines() if l.startswith(("ATOM", "HETATM")))
    except Exception:
        pass
    return write_fallback_receptor_pdbqt(protein_path, out_pdbqt)


# ---------------- 引擎运行 ----------------
def _run(cmd: list[str], log_path: Path, timeout: int = 600, cwd: Path | None = None,
         cancel_event=None) -> subprocess.CompletedProcess:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w", encoding="utf-8", errors="replace") as log:
        proc = subprocess.Popen(
            cmd, stdout=log, stderr=subprocess.STDOUT, cwd=str(cwd) if cwd else None,
            creationflags=CREATE_NO_WINDOW if os.name == "nt" else 0,
        )
        try:
            rc = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            raise TimeoutError(f"引擎超时（>{timeout}s）: {cmd[0]}")
    if cancel_event is not None and getattr(cancel_event, "is_set", None) and cancel_event.is_set():
        raise InterruptedError("任务被取消")
    return subprocess.CompletedProcess(cmd, rc, "", "")


def select_engine(requested: str | None = None) -> str:
    """GPU 优先选择引擎；逐个 --version 试跑验证可执行性。"""
    def usable(name: str) -> bool:
        path = settings.tool(name)
        if path is None:
            return False
        try:
            proc = subprocess.run([str(path), "--version"], capture_output=True,
                                  timeout=8, creationflags=CREATE_NO_WINDOW)
            return proc.returncode in (0, 1)
        except Exception:
            return False

    order = ([requested] if requested else ["vina_gpu", "vina_cpu", "smina"])
    for name in order:
        if name and usable(name):
            return name
    raise RuntimeError("未找到可用对接引擎（vina_gpu/vina_cpu/smina）；"
                       "请确认引擎可执行文件与本机平台匹配并在设置中配置路径")


def _nvidia_ok() -> bool:
    from shutil import which

    smi = which("nvidia-smi")
    if not smi:
        return False
    try:
        return subprocess.run([smi], capture_output=True, timeout=6).returncode == 0
    except Exception:
        return False


def run_docking(engine: str, receptor_pdbqt: Path, ligand_pdbqt: Path,
                out_pdbqt: Path, log_path: Path, pocket: dict,
                *, exhaustiveness: int = 8, num_modes: int = 9,
                seed: int = 42, thread: int = 8192, search_depth: int = 1,
                timeout: int = 900, cancel_event=None) -> dict:
    """运行指定引擎，返回 {engine, score, pose_path, log}。"""
    center = {k: float(pocket[k]) for k in ("center_x", "center_y", "center_z")}
    size = float(pocket.get("size") or max(20.0, float(pocket.get("radius", 5.0)) * 2 + 8.0))
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)

    if engine == "vina_gpu":
        exe = settings.tool("vina_gpu")
        cmd = [str(exe),
               "--receptor", str(receptor_pdbqt), "--ligand", str(ligand_pdbqt),
               "--out", str(out_pdbqt),
               "--center_x", f"{center['center_x']:.3f}", "--center_y", f"{center['center_y']:.3f}",
               "--center_z", f"{center['center_z']:.3f}",
               "--size_x", f"{size:.1f}", "--size_y", f"{size:.1f}", "--size_z", f"{size:.1f}",
               "--thread", str(thread), "--search_depth", str(search_depth),
               "--exhaustiveness", "1", "--seed", str(seed), "--num_modes", str(num_modes)]
        _run(cmd, log_path, timeout=max(timeout, 1800), cancel_event=cancel_event)

    elif engine == "vina_cpu":
        exe = settings.tool("vina")
        cmd = [str(exe),
               "--receptor", str(receptor_pdbqt), "--ligand", str(ligand_pdbqt),
               "--out", str(out_pdbqt),
               "--center_x", f"{center['center_x']:.3f}", "--center_y", f"{center['center_y']:.3f}",
               "--center_z", f"{center['center_z']:.3f}",
               "--size_x", f"{size:.1f}", "--size_y", f"{size:.1f}", "--size_z", f"{size:.1f}",
               "--exhaustiveness", str(exhaustiveness), "--num_modes", str(num_modes),
               "--seed", str(seed), "--cpu", str(os.cpu_count() or 2)]
        _run(cmd, log_path, timeout=timeout, cancel_event=cancel_event)

    elif engine == "smina":
        exe = settings.tool("smina")
        cmd = [str(exe),
               "-r", str(receptor_pdbqt), "-l", str(ligand_pdbqt), "-o", str(out_pdbqt),
               "--center_x", f"{center['center_x']:.3f}", "--center_y", f"{center['center_y']:.3f}",
               "--center_z", f"{center['center_z']:.3f}",
               "--size_x", f"{size:.1f}", "--size_y", f"{size:.1f}", "--size_z", f"{size:.1f}",
               "--exhaustiveness", str(exhaustiveness), "--num_modes", str(num_modes),
               "--seed", str(seed), "--cpu", str(os.cpu_count() or 2)]
        _run(cmd, log_path, timeout=timeout, cancel_event=cancel_event)

    else:
        raise ValueError(f"未知引擎: {engine}")

    score = parse_score(out_pdbqt, log_path)
    if score is None:
        raise RuntimeError(f"{engine} 未产生有效打分（见 {log_path.name}）")
    return {"engine": engine, "score": score, "pose_path": str(out_pdbqt), "log": str(log_path)}


def parse_score(pose_pdbqt: Path, log_path: Path | None = None) -> float | None:
    m = re.search(r"REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)", Path(pose_pdbqt).read_text(errors="replace")) if Path(pose_pdbqt).is_file() else None
    if m:
        return float(m.group(1))
    if log_path and Path(log_path).is_file():
        text = Path(log_path).read_text(errors="replace")
        m = re.search(r"^\s+1\s+(-?\d+(?:\.\d+)?)", text, re.M)  # vina 表格
        if m:
            return float(m.group(1))
        m = re.search(r"Affinity:\s+(-?\d+(?:\.\d+)?)", text)  # smina
        if m:
            return float(m.group(1))
    return None


def pose_to_sdf(pose_pdbqt: Path, out_sdf: Path) -> bool:
    """Meeko 把对接姿态转回 SDF（供 3D 查看）；失败返回 False 不阻塞。"""
    try:
        from meeko import PDBQTMolecule, RDKitMolCreate

        pm = PDBQTMolecule.from_file(str(pose_pdbqt))
        sdf, _fail = RDKitMolCreate.write_sd_string(pm)
        if sdf:
            out_sdf.parent.mkdir(parents=True, exist_ok=True)
            out_sdf.write_text(sdf, encoding="utf-8")
            return True
    except Exception:
        pass
    return False
