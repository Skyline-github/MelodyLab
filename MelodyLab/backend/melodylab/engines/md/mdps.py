"""GROMACS mdp 参数文件生成（分阶段：ions/em/nvt/npt/md）。"""
from __future__ import annotations

BASE_PM = {
    "cutoff-scheme": "Verlet",
    "nstlist": 20,
    "coulombtype": "PME",
    "fourierspacing": "0.16",
    "pme_order": 4,
    "rcoulomb": 1.0,
    "rvdw": 1.0,
    "rlist": 1.0,
    "constraints": "h-bonds",
    "constraint-algorithm": "lincs",
    "lincs-iter": 1,
    "lincs-order": 4,
    "pbc": "xyz",
}


def _write(path, title: str, pairs: dict) -> None:
    lines = [f"; {title} — MelodyLab MD"]
    if "integrator" not in pairs:
        lines.append("integrator = md")
    lines += [f"{k} = {v}" for k, v in pairs.items()]
    path.write_text("\n".join(lines) + "\n\n", encoding="utf-8", newline="\n")


def write_mdp(path, stage: str, *, nsteps: int = 0, seed: int = 42, ref_t: float = 300.0) -> None:
    """stage ∈ {ions, em, nvt, npt, md}；nsteps 单位=步数。"""
    if stage == "ions":
        _write(path, "ions", {"integrator": "steep", "emtol": 1000.0, "emstep": 0.01,
                              "nstlist": 10, **{k: v for k, v in BASE_PM.items()
                                                if k not in ("constraints", "constraint-algorithm",
                                                             "lincs-iter", "lincs-order")}})
        return
    if stage in ("em", "em2"):
        # 最小化阶段不用约束（CG 遇到约束失败会直接报错；无约束让 H 自由松弛）
        pm = {k: v for k, v in BASE_PM.items()
              if k not in ("constraints", "constraint-algorithm", "lincs-iter", "lincs-order")}
        _write(path, stage, {"integrator": "steep" if stage == "em" else "cg",
                             "emtol": 2000.0 if stage == "em" else 500.0,
                             "emstep": 0.01, "nsteps": 50000,
                             "nstlist": 10, "nstenergy": 100, "nstlog": 100, **pm})
        return
    # 通用输出频率
    out = {"dt": 0.002, "nsteps": int(nsteps),
           "nstxout": 0, "nstvout": 0, "nstenergy": 500, "nstlog": 500,
           "nstxout-compressed": 1000, "compressed-x-precision": 1000,
           "tcoupl": "V-rescale", "tc-grps": "Protein Non-Protein",
           "tau-t": "0.1 0.1", "ref-t": f"{ref_t} {ref_t}"}
    if stage == "nvt":
        out.update({"pcoupl": "no", "gen-vel": "yes", "gen-temp": 50.0,
                    "gen-seed": int(seed), "define": "-DPOSRES", "continuation": "no"})
        if ref_t <= 60.0:  # 低温预平衡段：更小步长
            out["dt"] = 0.001
    elif stage == "npt":
        out.update({"pcoupl": "Berendsen", "pcoupltype": "isotropic", "tau-p": 2.0,
                    "ref-p": 1.0, "compressibility": "4.5e-5", "gen-vel": "no",
                    "define": "-DPOSRES", "continuation": "yes"})
    elif stage == "md":
        out.update({"pcoupl": "Parrinello-Rahman", "pcoupltype": "isotropic", "tau-p": 2.0,
                    "ref-p": 1.0, "compressibility": "4.5e-5", "gen-vel": "no",
                    "nstxout-compressed": 5000, "continuation": "yes"})
    else:
        raise ValueError(f"未知 mdp 阶段: {stage}")
    _write(path, stage, out)
