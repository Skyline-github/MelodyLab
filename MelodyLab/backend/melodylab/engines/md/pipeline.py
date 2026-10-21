"""MD 流水线 worker：阶段状态机（prep→em→nvt→gate→npt→prod→analysis）。

每阶段完成即写 work/md_state.json；暂停/中断后 resume 跳过已完成阶段，
未完成阶段的 mdrun 若有 checkpoint（.cpt）则以 -cpi 续跑。
门控：NVT（预检）后蛋白骨架 RMSD 尾段均值 < rmsd_threshold 才进入 NPT/生产。
"""
from __future__ import annotations

import hashlib
import json
import os
import time
import zipfile
from pathlib import Path

from ...config import settings
from . import analyze, prep
from .gmx import gpu_available, run_gmx, xvg_last_average, xvg_series
from .ligand import load_ligand, write_ligand_files
from .mdps import write_mdp

PS_PER_STEP = 0.002  # dt


def _state_path(work: Path) -> Path:
    return work / "md_state.json"


def _load_state(work: Path) -> dict:
    f = _state_path(work)
    if f.is_file():
        try:
            return json.loads(f.read_text(encoding="utf-8"))
        except Exception:
            pass
    return {"done": {}}


def _save_state(work: Path, state: dict) -> None:
    _state_path(work).write_text(json.dumps(state, ensure_ascii=False, indent=1),
                                 encoding="utf-8", newline="\n")


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _locate_inputs(job_dir: Path, params: dict) -> tuple[Path, Path | None, str | None]:
    """返回 (蛋白 pdb, 配体 sdf|None, smiles|None)。优先对接任务引用，否则用上传。"""
    inputs = job_dir / "inputs"
    run_id, task_id = params.get("run_id"), params.get("task_id")
    if run_id and task_id:
        src_task = settings.jobs_dir / run_id / "work" / "tasks" / task_id
        if not src_task.is_dir():
            raise ValueError(f"对接任务不存在: {run_id}/{task_id}")
        sdf = src_task / "docked.sdf"
        proteins = [p for p in (settings.jobs_dir / run_id / "inputs").glob("*")
                    if p.suffix.lower() in (".pdb", ".ent")]
        if not proteins:
            raise ValueError("对接 run 目录中未找到受体 PDB")
        return proteins[0], (sdf if sdf.is_file() else None), None
    proteins = [p for p in inputs.glob("*") if p.suffix.lower() in (".pdb", ".ent")]
    if not proteins:
        raise ValueError("请上传蛋白 PDB（或引用一个对接 run_id+task_id）")
    lig = next((p for p in inputs.glob("*.sdf")), None)
    smi = None
    if lig is None and (inputs / "ligands.txt").is_file():
        line = (inputs / "ligands.txt").read_text(encoding="utf-8").strip().splitlines()
        smi = line[0].split()[0] if line else None
    return proteins[0], lig, smi


def _md_python() -> Path | None:
    """OpenFF 参数化子进程解释器：环境变量优先，缺省找同根 conda env mdff。"""
    import sys as _sys

    env = os.environ.get("MELODYLAB_MD_PYTHON")
    if env and Path(env).is_file():
        return Path(env)
    sibling = Path(_sys.executable).parents[1] / "mdff" / "python.exe"
    return sibling if sibling.is_file() else None


def _make_ligand_sdf(smiles: str, out_sdf: Path) -> Path:
    from rdkit import Chem
    from rdkit.Chem import AllChem

    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(m, randomSeed=42) == -1:
        raise ValueError(f"SMILES 3D 构象生成失败: {smiles}")
    AllChem.MMFFOptimizeMolecule(m)
    w = Chem.SDWriter(str(out_sdf))
    w.write(m)
    w.close()
    return out_sdf


def _ligand_openff(work: Path, lig_sdf: Path, lig_smiles: str | None) -> dict:
    """独立 mdff 环境（OpenFF Sage 柔性拓扑）子进程参数化。"""
    import subprocess

    md_py = _md_python()
    if lig_sdf is None:
        if not lig_smiles:
            raise ValueError("OpenFF 模式需要配体 SDF 或 SMILES")
        lig_sdf = _make_ligand_sdf(lig_smiles, work / "ligand_input.sdf")
    script = Path(__file__).parent / "ligand_openff.py"
    proc = subprocess.run(
        [str(md_py), str(script), str(lig_sdf), str(work)],
        capture_output=True, text=True, timeout=600,
        creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
    if proc.returncode != 0 or not (work / "ligand.itp").is_file():
        raise RuntimeError(f"OpenFF 配体参数化失败:\n{(proc.stderr or proc.stdout)[-600:]}")
    meta = json.loads((work / "ligand_meta.json").read_text(encoding="utf-8"))
    return meta


def _mdrun(work: Path, deffnm: str, use_gpu: bool, cancel_event, nt: int) -> None:
    args = ["mdrun", "-deffnm", deffnm, "-nt", str(nt), "-pin", "on"]
    state = _load_state(work)
    cpt = work / f"{deffnm}.cpt"
    resuming = bool(state.get("resumed_stage")) and cpt.is_file()
    if use_gpu:
        try:
            run_gmx([*args, "-nb", "gpu"], cwd=work, cancel_event=cancel_event, timeout=7200)
            return
        except RuntimeError as exc:
            if cancel_event is not None and getattr(cancel_event, "is_set", None) and cancel_event.is_set():
                raise
            _log(work, f"-nb gpu 失败，回退 CPU: {str(exc)[:160]}")
    # 该 Windows 构建存在瞬时段错误（0xC0000005），自动重试一次
    try:
        if resuming:
            run_gmx([*args, "-cpi", f"{deffnm}.cpt"], cwd=work, cancel_event=cancel_event, timeout=7200)
        else:
            run_gmx(args, cwd=work, cancel_event=cancel_event, timeout=7200)
    except RuntimeError as exc:
        if "3221225477" not in str(exc) or (cancel_event is not None and getattr(cancel_event, "is_set", None) and cancel_event.is_set()):
            raise
        _log(work, "mdrun 段错误，5 秒后重试一次")
        time.sleep(5)
        if resuming:
            run_gmx([*args, "-cpi", f"{deffnm}.cpt"], cwd=work, cancel_event=cancel_event, timeout=7200)
        else:
            run_gmx(args, cwd=work, cancel_event=cancel_event, timeout=7200)


def _log(work: Path, msg: str) -> None:
    with open(work / "melodylab_md.log", "a", encoding="utf-8") as f:
        f.write(f"[{time.strftime('%H:%M:%S')}] {msg}\n")


def _lig_group_number(work: Path) -> str:
    """从 index.ndx 找 LIG 组的序号（按组名选择在某些 trjconv 场景解析异常）。"""
    n = 0
    for line in (work / "index.ndx").read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if s.startswith("["):
            name = s.strip("[] ").strip()
            if name.upper() == "LIG":
                return str(n)
            n += 1
    return "LIG"


def _refresh_ligand(work: Path, source: str = "em.gro") -> None:
    """按 EM 后实际坐标刷新刚体配体的约束长度/3out 系数（防 NVT LINCS 崩溃）。

    关键：先做最小镜像展开再拟合——grompp 会把跨边界原子逐个卷绕，
    直接从撕裂几何拟合出的 vsite 系数会导致运行时内爆。
    """
    from .ligand import extract_ligand_gro, refresh_geometry_from_gro, unwrap_ligand_in_gro

    lig_gro = work / "lig_refresh.gro"
    extract_ligand_gro(work / source, lig_gro)
    unwrap_ligand_in_gro(lig_gro)
    refresh_geometry_from_gro(work / "ligand.itp", lig_gro)
    lig_gro.unlink(missing_ok=True)


def _make_whole(work: Path, stage: str) -> None:
    """mdrun 后把 PBC 拆开的分子拼回整分子，并装入紧凑原胞（-ur compact）。

    -ur compact 是关键：grompp 会对出盒原子逐个卷绕，多链蛋白在十二面体盒的
    三斜表示里必然跨越边界——不 compact 就会在下一次 grompp 被撕开，
    之后所有 tpr/xtc 与参考结构的坐标错位一个盒向量，RMSD 恒定虚高。
    """
    run_gmx(["trjconv", "-s", f"{stage}.tpr", "-f", f"{stage}.gro",
             "-o", f"{stage}_w.gro", "-pbc", "mol", "-ur", "compact"],
            cwd=work, stdin_text="System\n")
    (work / f"{stage}_w.gro").replace(work / f"{stage}.gro")
    # -pbc mol 按质心入盒，带无序尾巴的蛋白仍会有原子戳出盒外——
    # 统一做逐原子规范化，让 gro/tpr/xtc 共享同一包装状态
    from .ligand import wrap_solute_atoms

    wrap_solute_atoms(work / f"{stage}.gro")


def _compact(work: Path, tpr: str, src: str, out: str) -> str:
    """grompp 前把 LIG 平移入胞（vsites 要求构造原子同镜像）。

    纯 Python 实现，不走 trjconv：其 -pbc mol 对含虚拟站点分子会段错误
    （GROMACS 2020.6 已知问题）。蛋白无 vsites，逐原子卷绕无碍。
    """
    import shutil

    from .ligand import unwrap_ligand_in_gro, wrap_solute_atoms

    dst = work / out
    if src != out:
        shutil.copyfile(work / src, dst)
    unwrap_ligand_in_gro(dst)
    wrap_solute_atoms(dst)
    return out


def run_md_job(record: dict, control) -> dict:
    job_dir = settings.jobs_dir / record["id"]
    params = record["params"]
    work = job_dir / "work"
    outputs = job_dir / "outputs"
    outputs.mkdir(exist_ok=True)
    resuming = bool(params.get("resume", False))
    state = _load_state(work)
    if resuming:
        state["resumed_stage"] = True

    n_cpu = max(2, (os_cpu := __import__("os").cpu_count() or 4) - 2)
    use_gpu = bool(params.get("use_gpu", True)) and gpu_available()
    rmsd_threshold = float(params.get("rmsd_threshold", 0.3))
    seed = int(params.get("seed", 42))
    ref_t = float(params.get("ref_t", 300.0))

    def done(stage: str) -> bool:
        return state["done"].get(stage) is True

    # ---------- 阶段 1：体系构建 ----------
    if not done("prep"):
        control.checkpoint(progress=5, phase="prep", message="定位输入（蛋白/配体）…")
        protein, lig_sdf, lig_smiles = _locate_inputs(job_dir, params)
        _log(work, f"protein={protein.name} ligand={lig_sdf.name if lig_sdf else lig_smiles}")
        info = prep.filter_protein_chains(protein, work / "protein_clean.pdb",
                                          params.get("chains"))
        repair = prep.repair_protein(work / "protein_clean.pdb", work / "protein_clean.pdb")
        control.checkpoint(progress=6, phase="prep",
                           message=f"选链清洗完成（{info['atoms']} 重原子，链 {info['chains']}；"
                                   f"修复剔除 {repair['dropped_residues'] or '无'}）")
        prep.pdb2gmx(work / "protein_clean.pdb", work,
                     params.get("forcefield") or "amber99sb", params.get("water") or "tip3p")
        lig_mode = "openff" if _md_python() is not None else "rigid"
        state["ligand_mode"] = lig_mode
        if lig_mode == "openff":
            lig_meta = _ligand_openff(work, lig_sdf, lig_smiles)
            lig_info = {"natoms": lig_meta["natoms"], "charge": lig_meta["charge"],
                        "n_vsites": 0}
            control.checkpoint(progress=7, phase="prep",
                               message=f"配体 OpenFF 柔性拓扑完成（{lig_info['natoms']} 原子，"
                                       f"净电荷 {lig_info['charge']}）")
        else:
            mol = load_ligand(lig_sdf, lig_smiles)
            lig_info = write_ligand_files(mol, "LIG", work, params.get("ligand_charge"))
            control.checkpoint(progress=7, phase="prep",
                               message=f"配体刚体拓扑完成（{lig_info['natoms']} 原子，"
                                       f"{lig_info['n_vsites']} 虚拟站点，净电荷 {lig_info['charge']}）")
        params["ligand_mode"] = lig_mode
        sysinfo = prep.build_system(work, params, control)
        state["done"]["prep"] = True
        state["system"] = {**info, **{k: sysinfo[k] for k in ("natoms", "forcefield", "water")},
                           "ligand": {k: lig_info[k] for k in ("natoms", "charge")}}
        _save_state(work, state)
    else:
        sysinfo = {"natoms": state.get("system", {}).get("natoms"), "forcefield":
                   state.get("system", {}).get("forcefield")}

    natoms = state.get("system", {}).get("natoms")
    # 模拟时长（ns）→ 步数
    test_ns = float(params.get("test_ns", 1.0))
    npt_ns = float(params.get("npt_ns", 1.0))
    total_ns = float(params.get("total_ns", 5.0))

    def _grompp(stage: str, coord: str, extra: list[str] | None = None) -> None:
        write_mdp(work / f"{stage}.mdp", stage,
                  nsteps=int(nsteps_for(stage)), seed=seed, ref_t=ref_t)
        args = ["grompp", "-f", f"{stage}.mdp", "-c", coord, "-p", "topol.top",
                "-o", f"{stage}.tpr", "-maxwarn", "2"]
        if stage in ("nvt", "npt"):
            args += ["-r", coord]
        args += (extra or [])
        run_gmx(args, cwd=work, cancel_event=control.cancel_event)

    def nsteps_for(stage: str) -> float:
        if stage in ("em", "em2"):
            return 0  # write_mdp 内置步数
        ns = {"nvt": test_ns, "npt": npt_ns, "md": total_ns}[stage]
        return ns * 1000 / PS_PER_STEP

    def _run_stage(stage: str, coord: str, p_lo: int, p_hi: int, extra_grompp=None) -> None:
        if done(stage):
            return
        control.checkpoint(progress=p_lo, phase=stage,
                           message=f"{'续跑 ' if state.get('resumed_stage') and (work / f'{stage}.cpt').is_file() else ''}{stage} 模拟中…")
        if stage != "em":
            # em 输入 (ionized.gro) 由 solvate/genion 生成、无需整分子化；
            # 且 ions.tpr 是加离子前拓扑不能作 -s。em.tpr 出来后的阶段统一 compact。
            coord = _compact(work, "em.tpr", coord, f"{stage}_in.gro")
        _grompp(stage, coord, extra_grompp)
        _mdrun(work, stage, use_gpu, control.cancel_event, n_cpu)
        # mdrun 输出的分子会被 PBC 拆到不同镜像，拼回整分子供下一阶段 grompp
        _make_whole(work, stage)
        state["done"][stage] = True
        state["resumed_stage"] = False
        _save_state(work, state)
        control.checkpoint(progress=p_hi, phase=stage, message=f"{stage} 完成")

    _run_stage("em", "ionized.gro", 26, 27)
    rigid = state.get("ligand_mode", "rigid") == "rigid"
    if not done("em2"):
        if rigid:  # 刚体模式才需要几何刷新；OpenFF 柔性键合由力场参数自洽
            control.checkpoint(progress=27, phase="em2", message="刷新刚体几何后 CG 精修…")
            _refresh_ligand(work)
    _run_stage("em2", "em.gro", 28, 29)
    if not done("nvt") and rigid:
        control.checkpoint(progress=29, phase="nvt", message="刷新刚体几何（CG 后）…")
        _refresh_ligand(work, "em2.gro")

    # ---------- NVT 短模拟稳定性判定（预检门控，MD v0.2 的核心产品） ----------
    verdict = None
    if not done("gate"):
        stage_lo = "nvt_lo"
        if not done(stage_lo):
            control.checkpoint(progress=29, phase=stage_lo, message="NVT 稳定性判定（50K 短模拟）…")
            write_mdp(work / f"{stage_lo}.mdp", "nvt", nsteps=2000, seed=seed, ref_t=50.0)
            coord_lo = _compact(work, "em.tpr", "em2.gro", f"{stage_lo}_in.gro")
            run_gmx(["grompp", "-f", f"{stage_lo}.mdp", "-c", coord_lo, "-r", coord_lo,
                     "-p", "topol.top", "-o", f"{stage_lo}.tpr", "-maxwarn", "2"],
                    cwd=work, cancel_event=control.cancel_event)
            try:
                _mdrun(work, stage_lo, use_gpu, control.cancel_event, n_cpu)
                _make_whole(work, stage_lo)
                lincs_fail = _count_lincs_errors(work / f"{stage_lo}.log")
                if lincs_fail == 0 and (work / f"{stage_lo}.gro").is_file():
                    verdict = {"verdict": "stable", "n_steps": 2000, "ref_t": 50,
                               "note": "2000 步无约束失败，体系稳定"}
                else:
                    verdict = {"verdict": "unstable", "n_steps": 2000,
                               "lincs_errors": lincs_fail,
                               "note": "短模拟出现约束失败，体系不稳定（复合物/参数化问题）"}
            except (RuntimeError, TimeoutError) as exc:
                verdict = {"verdict": "engine_error", "note": f"引擎异常退出: {str(exc)[:200]}"}
            state["verdict"] = verdict
            state["done"][stage_lo] = True
            _save_state(work, state)
        else:
            verdict = state.get("verdict")
    else:
        verdict = {"verdict": "stable", "passed_full_gate": True}

    # verdict 为 stable 且用户要求长模拟 → 继续完整流水线；否则以判定结果收尾
    proceed_full = verdict is not None and verdict.get("verdict") == "stable" and total_ns > 0

    if proceed_full:
        _run_stage("nvt", "nvt_lo.gro", 30, 45)

        # ---------- 门控 ----------
        if not done("gate"):
            control.checkpoint(progress=45, phase="gate", message="NVT 预检门控：计算骨架 RMSD…")
            run_gmx(["rms", "-s", "em.gro", "-f", "nvt.xtc", "-n", str(_ensure_index(work)),
                     "-o", "rms_nvt.xvg", "-fit", "rot+trans"],
                    cwd=work, stdin_text="Backbone\nBackbone\n", cancel_event=control.cancel_event)
            rms = xvg_last_average(work / "rms_nvt.xvg", 0.2)
            if rms is None:
                raise RuntimeError("NVT 门控 RMSD 计算失败")
            if rms >= rmsd_threshold:
                state["done"]["gate"] = False
                state["gate_failed"] = {"rmsd": rms, "threshold": rmsd_threshold}
                _save_state(work, state)
                raise RuntimeError(
                    f"预检门控未通过：NVT 末段骨架 RMSD {rms:.3f} nm ≥ 阈值 {rmsd_threshold} nm。"
                    "体系不稳定（复合物/参数化问题），已停止，不出正式轨迹。")
            state["done"]["gate"] = True
            state["gate_rmsd"] = rms
            _save_state(work, state)

        _run_stage("npt", "nvt.gro", 60, 68)
        _run_stage("md", "npt.gro", 70, 82)

    # ---------- 分析 ----------
    if not done("analysis"):
        control.checkpoint(progress=85, phase="analysis", message="轨迹分析中…")
        if (work / "md.xtc").is_file():
            report = analyze.run_analysis(work, outputs, control, rmsd_threshold)
            report["gate"] = {"threshold_nm": rmsd_threshold,
                              "nvt_rmsd_tail_avg": state.get("gate_rmsd"),
                              "passed": True}
        else:
            # 烟雾判定模式：只用 nvt_lo 短轨迹做骨架/配体 RMSD
            report = _smoke_analysis(work)
            a_dir = outputs / "analysis"
            a_dir.mkdir(parents=True, exist_ok=True)
            for png in ("rms_smoke_bb.png", "rms_smoke_lig.png"):
                if (work / png).is_file():
                    (work / png).replace(a_dir / png)
        report["verdict"] = verdict
        report["system"] = state.get("system", {})
        (outputs / "md_report.json").write_text(
            json.dumps(report, ensure_ascii=False, indent=1), encoding="utf-8", newline="\n")
        if (work / "md.xtc").is_file():
            _trajectory_frames(work, outputs / "trajectory_preview.pdb", n_frames=PREVIEW_FRAMES)
        _traceability_zip(job_dir, outputs)
        state["done"]["analysis"] = True
        _save_state(work, state)
    else:
        report = json.loads((outputs / "md_report.json").read_text(encoding="utf-8"))

    # ---------- 汇总 ----------
    arts = [{"name": "md_report.json", "path": str(outputs / "md_report.json"), "kind": "json"}]
    if (outputs / "trajectory_preview.pdb").is_file():
        arts.append({"name": "trajectory_preview.pdb",
                     "path": str(outputs / "trajectory_preview.pdb"), "kind": "pdb"})
    if (work / "md.xtc").is_file():
        arts.append({"name": "production_trajectory.xtc", "path": str(work / "md.xtc"), "kind": "xtc"})
    for png in sorted((outputs / "analysis").glob("*.png")):
        arts.append({"name": png.name, "path": str(png), "kind": "png"})
    zip_path = outputs / "md_traceability.zip"
    if zip_path.is_file():
        arts.append({"name": "md_traceability.zip", "path": str(zip_path), "kind": "zip"})
    lig_rms = (report.get("ligand_rmsd") or {}).get("avg_tail")
    bb_rms = (report.get("backbone_rmsd") or {}).get("avg_tail")
    v = verdict.get("verdict") if verdict else "unknown"
    msg = {"stable": "稳定性判定：✅ 稳定",
           "unstable": "稳定性判定：⚠️ 不稳定（见报告诊断）",
           "engine_error": "稳定性判定：⚠️ 引擎异常（GROMACS 崩溃，见溯源包日志）"}.get(v, "完成")
    if bb_rms is not None:
        msg += f"；骨架 RMSD {bb_rms} nm"
    if lig_rms is not None:
        msg += f"，配体 RMSD {lig_rms} nm"
    return {"stats": {
                "natoms": natoms,
                "verdict": v,
                "gate_rmsd_nm": state.get("gate_rmsd"),
                "backbone_rmsd_tail_nm": bb_rms,
                "ligand_rmsd_tail_nm": lig_rms,
                "engine": "gpu" if use_gpu else "cpu",
            },
            "message": f"MD 完成（{'GPU' if use_gpu else 'CPU'}，{natoms} 原子）—— {msg}",
            "artifacts": arts}


def _count_lincs_errors(log_path: Path) -> int:
    if not log_path.is_file():
        return 0
    import re

    return len(re.findall(r"Constraint error|can not be settled",
                          Path(log_path).read_text(encoding="utf-8", errors="replace")))


def _smoke_analysis(work: Path) -> dict:
    """烟雾判定模式的轻量分析：nvt_lo 短轨迹的骨架/配体 RMSD。"""
    report: dict = {"plots": [], "curves": {}, "mode": "smoke"}
    if not (work / "nvt_lo.xtc").is_file():
        report["note"] = "短模拟未产出轨迹"
        return report
    ndx = _ensure_index(work)
    lig_no = _lig_group_number(work)
    try:
        run_gmx(["rms", "-s", "nvt_lo.tpr", "-f", "nvt_lo.xtc", "-n", ndx,
                 "-o", "rms_smoke_bb.xvg", "-fit", "rot+trans"],
                cwd=work, stdin_text="Backbone\nBackbone\n")
        run_gmx(["rms", "-s", "nvt_lo.tpr", "-f", "nvt_lo.xtc", "-n", ndx,
                 "-o", "rms_smoke_lig.xvg", "-fit", "rot+trans"],
                cwd=work, stdin_text=f"Backbone\n{lig_no}\n")
        analyze._png(work / "rms_smoke_bb.xvg", work / "rms_smoke_bb.png",
                     "Backbone RMSD (50K smoke)", "RMSD (nm)")
        analyze._png(work / "rms_smoke_lig.xvg", work / "rms_smoke_lig.png",
                     "Ligand RMSD (50K smoke)", "RMSD (nm)")
        (work / "rms_smoke_bb.png").replace(work / "rms_smoke_bb.png")
        report["backbone_rmsd"] = {"avg_tail": xvg_last_average(work / "rms_smoke_bb.xvg", 0.3)}
        report["ligand_rmsd"] = {"avg_tail": xvg_last_average(work / "rms_smoke_lig.xvg", 0.3)}
        report["curves"]["backbone_rmsd"] = xvg_series(work / "rms_smoke_bb.xvg")
    except Exception as exc:
        report["smoke_error"] = str(exc)[:200]
    return report


def _ensure_index(work: Path) -> str:
    ndx = work / "index.ndx"
    if not ndx.is_file():
        # em.gro 是第一个经过 -pbc mol -ur compact 的完整体系，
        # 原子数与后续所有 tpr/xtc 一致（ions.tpr 是加离子前的旧拓扑，勿用）
        run_gmx(["make_ndx", "-f", "em.gro", "-o", "index.ndx"], cwd=work, stdin_text="r LIG\nq\n")
    return "index.ndx"


PREVIEW_FRAMES = 30  # 轨迹预览导出帧数（分析阶段与预览路由共用）


def _ndx_first_group_numbers(work: Path) -> dict[str, int]:
    """index.ndx 组名 → 首次出现的组号（make_ndx 会追加重名组，按名选择不可靠）。"""
    groups: dict[str, int] = {}
    idx = 0
    for line in (work / "index.ndx").read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if s.startswith("["):
            name = s.strip("[] ").strip()
            groups.setdefault(name, idx)
            idx += 1
    return groups


def _unwrap_pdb_frame(lines: list[str]) -> str:
    """单帧 PDB 的可视化 PBC 整理：剥 MODEL/ENDMDL 与离子，再沿连接顺序
    做最小镜像展开（模拟 trjconv -pbc whole，但按文件顺序=拓扑连接顺序行走，
    不依赖 tpr 键表）。

    背景：坐标经逐原子规范化后，柔性尾巴的少数原子被卷绕到盒子对侧，
    动力学不受影响，但画出来会"断成几截"。行走展开把每个原子拉回与
    前一相邻原子最近镜像的位置，整条链在视觉上重连。
    """
    out: list[str] = []
    prev: tuple[float, float, float] | None = None
    box_hint = 0.0
    for l in lines:
        if l.startswith(("MODEL", "ENDMDL")):
            continue
        if l.startswith("CRYST1"):
            try:
                box_hint = float(l[6:15])
            except ValueError:
                pass
            out.append(l)
            continue
        if l.startswith(("ATOM", "HETATM")):
            if l[17:20].strip() in ("NA", "CL"):
                continue  # 离子不进可视化（散布全盒，会被误认为小分子）
            x, y, z = float(l[30:38]), float(l[38:46]), float(l[46:54])
            L = box_hint or 114.0
            if prev is not None:
                # 最小镜像：与前一个原子（键连）的距离超过 L/2 就平移一个晶格
                if x - prev[0] > L / 2:
                    x -= L
                elif prev[0] - x > L / 2:
                    x += L
                if y - prev[1] > L / 2:
                    y -= L
                elif prev[1] - y > L / 2:
                    y += L
                if z - prev[2] > L / 2:
                    z -= L
                elif prev[2] - z > L / 2:
                    z += L
            prev = (x, y, z)
            out.append(f"{l[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{l[54:]}")
            continue
        out.append(l)
    return "\n".join(out)


def _trajectory_frames(work: Path, out_pdb: Path, n_frames: int = 15) -> None:
    """XTC → 多 MODEL 的 PDB 快照（供前端 3Dmol 顺序播放）。

    只导出蛋白+配体：全体系 14.5 万原子 × 30 帧的 PDB 有 170MB+，
    浏览器端 3Dmol 根本消化不了。
    """
    xtc = work / "md.xtc"
    tpr = work / "md.tpr"
    if not xtc.is_file():
        return
    try:
        g = _ndx_first_group_numbers(work)
        # 本构建的 trjconv 多组选择只认第一个组；non-Water = 蛋白+配体+离子，一步到位
        sel = str(g["non-Water"]) if "non-Water" in g else str(g["Protein"])
    except Exception:
        sel = "Protein"
    s = xvg_series(work / "gyrate.xvg") if (work / "gyrate.xvg").is_file() else None
    total_ps = s["time"][-1] if s and s["time"] else None
    models = []
    for i in range(n_frames):
        dump = total_ps * i / max(1, n_frames - 1) if total_ps else None
        args = ["trjconv", "-s", "md.tpr", "-f", "md.xtc", "-o", f"_frame{i}.pdb", "-pbc", "no"]
        args += (["-dump", f"{dump:.1f}"] if dump is not None else ["-b", "-1", "-e", "-1"])
        try:
            run_gmx(args, cwd=work, stdin_text=f"{sel}\n")
            models.append((work / f"_frame{i}.pdb").read_text(encoding="utf-8"))
        except Exception:
            continue
    if models:
        merged = []
        for i, m in enumerate(models, 1):
            body = _unwrap_pdb_frame(m.strip().splitlines())
            merged.append(f"MODEL     {i:4d}\n{body}\nENDMDL")
        out_pdb.write_text("\n".join(merged) + "\n", encoding="utf-8", newline="\n")
    for i in range(n_frames):
        f = work / f"_frame{i}.pdb"
        if f.is_file():
            f.unlink()


def _traceability_zip(job_dir: Path, outputs: Path) -> None:
    """溯源包：输入/拓扑/mdp/分析图表 + SHA-256 清单。"""
    work = job_dir / "work"
    z = outputs / "md_traceability.zip"
    include = [job_dir / "inputs", work / "protein_clean.pdb", work / "topol.top",
               work / "ligand.itp", work / "ligand.gro", work / "md_state.json",
               work / "melodylab_md.log"]
    include += list(work.glob("*.mdp")) + list(work.glob("*.xvg"))
    include += list((outputs / "analysis").glob("*.png")) + [outputs / "md_report.json"]
    manifest = []
    with zipfile.ZipFile(z, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in include:
            if p.is_file():
                zf.write(p, f"md/{p.relative_to(job_dir)}")
                manifest.append({"file": str(p.relative_to(job_dir)), "sha256": _sha256(p),
                                 "bytes": p.stat().st_size})
            elif p.is_dir():
                for sub in sorted(p.rglob("*")):
                    if sub.is_file():
                        zf.write(sub, f"md/{sub.relative_to(job_dir)}")
                        manifest.append({"file": str(sub.relative_to(job_dir)),
                                         "sha256": _sha256(sub), "bytes": sub.stat().st_size})
        zf.writestr("md/SHA256SUMS.json", json.dumps(manifest, ensure_ascii=False, indent=1))
