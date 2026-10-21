"""MD 结果分析：gmx rms/rmsf/gyrate/energy → XVG → JSON + PNG + 报告。"""
from __future__ import annotations

import json
from pathlib import Path

from .gmx import run_gmx, xvg_series

WORK = "work"  # pipeline 以 gmx 工作目录为 work/，分析输出到 outputs/


def _png(xvg_path: Path, png_path: Path, title: str, ylabel: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    s = xvg_series(xvg_path)
    plt.figure(figsize=(7, 4), dpi=130)
    for name, vals in s["series"].items():
        t = s["time"]
        if len(t) == len(vals):
            plt.plot(t, vals, lw=1.0, label=name)
    plt.xlabel(s["x_label"] or "Time (ps)")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend(fontsize=7)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    png_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(png_path)
    plt.close()


def make_index(work: Path) -> Path:
    ndx = work / "index.ndx"
    if not ndx.is_file():
        run_gmx(["make_ndx", "-f", "md.tpr", "-o", "index.ndx"], cwd=work, stdin_text="r LIG\nq\n")
    return ndx


def run_analysis(work: Path, outputs: Path, control, gate_rmsd: float | None = None) -> dict:
    """执行全部分析，返回报告 dict（含图表路径与关键统计）。"""
    outputs.mkdir(parents=True, exist_ok=True)
    a_dir = outputs / "analysis"
    a_dir.mkdir(exist_ok=True)
    ndx = make_index(work)
    # index.ndx 可能含多个重名 LIG 组（make_ndx 每次运行都会追加），
    # 按名字选择会报 Multiple groups——按组号选择
    from .pipeline import _lig_group_number

    lig_sel = _lig_group_number(work)
    report: dict = {"plots": [], "curves": {}}

    def grp(*names: str) -> str:
        return "\n".join(names) + "\n"

    # 1) 蛋白骨架 RMSD（相对 em 结构）
    control.checkpoint(progress=88, phase="analysis", message="RMSD 分析…")
    run_gmx(["rms", "-s", "em.gro", "-f", "md.xtc", "-n", "index.ndx", "-o", "rms_backbone.xvg",
             "-fit", "rot+trans"], cwd=work, stdin_text=grp("Backbone", "Backbone"))
    _png(work / "rms_backbone.xvg", a_dir / "rmsd_backbone.png", "Backbone RMSD", "RMSD (nm)")
    bb = xvg_series(work / "rms_backbone.xvg")
    bb_final = _series_last(bb, 0.2)
    report["backbone_rmsd"] = {"avg_tail": bb_final, "max": _series_max(bb)}
    report["curves"]["backbone_rmsd"] = bb

    # 2) 配体 RMSD（以蛋白骨架拟合，衡量姿态稳定性）
    run_gmx(["rms", "-s", "em.gro", "-f", "md.xtc", "-n", "index.ndx", "-o", "rms_ligand.xvg",
             "-fit", "rot+trans"], cwd=work, stdin_text=grp("Backbone", lig_sel))
    _png(work / "rms_ligand.xvg", a_dir / "rmsd_ligand.png", "Ligand RMSD (fit on Backbone)",
         "RMSD (nm)")
    lig = xvg_series(work / "rms_ligand.xvg")
    lig_final = _series_last(lig, 0.2)
    report["ligand_rmsd"] = {"avg_tail": lig_final, "max": _series_max(lig)}
    report["curves"]["ligand_rmsd"] = lig

    # 3) RMSF（C-alpha）
    control.checkpoint(progress=91, phase="analysis", message="RMSF / 回转半径…")
    try:
        run_gmx(["rmsf", "-s", "md.tpr", "-f", "md.xtc", "-n", "index.ndx", "-o", "rmsf.xvg",
                 "-res"], cwd=work, stdin_text=grp("C-alpha"))
        _png(work / "rmsf.xvg", a_dir / "rmsf_ca.png", "C-alpha RMSF", "RMSF (nm)")
        report["curves"]["rmsf_ca"] = xvg_series(work / "rmsf.xvg")
    except Exception as exc:  # C-alpha 组缺失时跳过，不阻塞
        report["rmsf_error"] = str(exc)[:200]

    # 4) 回转半径（蛋白）
    try:
        run_gmx(["gyrate", "-s", "md.tpr", "-f", "md.xtc", "-n", "index.ndx", "-o", "gyrate.xvg"],
                cwd=work, stdin_text=grp("Protein"))
        _png(work / "gyrate.xvg", a_dir / "rgyration.png", "Radius of Gyration", "Rg (nm)")
        report["radius_gyration"] = {"avg_tail": _series_last(xvg_series(work / "gyrate.xvg"), 0.2)}
        report["curves"]["rgyration"] = xvg_series(work / "gyrate.xvg")
    except Exception as exc:
        report["gyrate_error"] = str(exc)[:200]

    # 5) 能量项（先 gmx check 拿项名顺序，再按名选号）
    control.checkpoint(progress=94, phase="analysis", message="能量分析…")
    try:
        terms = _energy_term_numbers(work)
        wanted = [n for n in ("Potential", "Kinetic-En.", "Total-Energy", "Temperature", "Pressure")
                  if n in terms]
        selection = " ".join(str(terms[n] + 1) for n in wanted)
        run_gmx(["energy", "-f", "md.edr", "-o", "energy.xvg"], cwd=work,
                stdin_text=f"{selection}\n0\n")
        _png(work / "energy.xvg", a_dir / "energy.png", "Energies", "")
        es = xvg_series(work / "energy.xvg")
        report["curves"]["energy"] = es
        report["energy_tail_avg"] = {k: _series_last({"time": es["time"], "series": {k: v}}, 0.2)
                                     for k in es["series"]}
    except Exception as exc:
        report["energy_error"] = str(exc)[:200]

    # 6) NVT 预检门控判定（rms_nvt.xvg 由 pipeline 在门控步生成）
    nvt_rms = work / "rms_nvt.xvg"
    if gate_rmsd is not None and nvt_rms.is_file():
        report["gate"] = {"threshold_nm": gate_rmsd,
                          "nvt_rmsd_tail_avg": _series_last(xvg_series(nvt_rms), 0.2)}
    return report


def _series_last(s: dict, frac: float) -> float | None:
    names = list(s["series"].keys())
    if not names:
        return None
    vals = s["series"][names[0]]
    if not vals:
        return None
    n = max(1, int(len(vals) * frac))
    return round(sum(vals[-n:]) / n, 4)


def _series_max(s: dict) -> float | None:
    names = list(s["series"].keys())
    if not names or not s["series"][names[0]]:
        return None
    return round(max(s["series"][names[0]]), 4)


def _energy_term_numbers(work: Path) -> dict[str, int]:
    """gmx check -f md.edr 输出里按名定位能量项序号。"""
    proc = run_gmx(["check", "-f", "md.edr"], cwd=work)
    names: dict[str, int] = {}
    started = False
    idx = 0
    for line in (proc.stdout or "").splitlines():
        if "Energy" in line and "Average" in line:
            started = True
            continue
        if started:
            if not line.strip():
                break
            parts = line.split()
            # 形如 "LJ (SR)   123.456 ..." —— 项名 = 行首数字列之前的部分
            first_num = next((i for i, p in enumerate(parts) if _is_float(p)), None)
            if first_num is None:
                continue
            idx += 1
            names[" ".join(parts[:first_num]).strip()] = idx - 1
    return names


def _is_float(s: str) -> bool:
    try:
        float(s.replace("-", "").replace(".", ""))
        return True
    except ValueError:
        return False
