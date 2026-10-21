"""二期 · 分子动力学引擎（GROMACS）。

数据流（对齐 DockScope md_service）：
  对接最佳复合物 → 选链/修蛋白 → 配体参数化（刚体兼容路线，见 ligand.py 头注）
  → GROMACS 分阶段（minim → NVT 预检门控(RMSD<阈值) → NPT → 生产）
  → mdrun(GPU 可用 -nb gpu) → RMSD/RMSF/Rg/能量 → XTC 轨迹 → 溯源 ZIP。

**状态：实验性，默认停用。** Windows 首启验证中 GROMACS 2020.6 流水线结构已走通
（prep/参数化/分阶段/分析/溯源全部实现），但刚体配体在镜像搬运后的初始构象仍会
出现配体-蛋白近距接触，导致 NVT 发散（详见 docs/交接文档.md §7.14）。
设置环境变量 MELODYLAB_ENABLE_MD=1 可重新启用。
"""
import os
import sys as _sys
from pathlib import Path

from ...config import settings


def availability() -> dict:
    gmx = settings.tool("gromacs")
    gpu = False
    if gmx is not None:
        try:
            from .gmx import gpu_available

            gpu = gpu_available()
        except Exception:
            gpu = False
    md_py = os.environ.get("MELODYLAB_MD_PYTHON")
    mdff = (Path(md_py).is_file() if md_py else
            Path(_sys.executable).parents[1].joinpath("mdff", "python.exe").is_file())
    return {
        "engine": "gromacs",
        "available": gmx is not None and bool(mdff),
        "path": str(gmx) if gmx else None,
        "gpu": gpu,
        "parameterization": "OpenFF Sage 柔性配体（mdff 独立环境子进程）+ amber99sb 蛋白",
        "hint": None if (gmx is not None and mdff) else
        "需要内置 GROMACS 与 mdff 环境（openff-toolkit/openmm/parmed）",
    }
