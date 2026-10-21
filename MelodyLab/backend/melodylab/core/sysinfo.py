"""环境自检：依赖、引擎可执行文件、GPU、可选模型资产。"""
from __future__ import annotations

import importlib
import shutil
import subprocess

from ..config import settings


def _module_version(name: str) -> str | None:
    try:
        mod = importlib.import_module(name)
        return getattr(mod, "__version__", None) or "ok"
    except Exception:
        return None


def _probe_executable(path, args: list[str] | None = None, timeout: float = 8.0) -> dict:
    import os

    if path is None:
        return {"available": False, "path": None, "path_exists": False}
    path = str(path)
    if not os.path.isfile(path):
        return {"available": False, "path": path, "path_exists": False}
    try:
        proc = subprocess.run(
            [path, *(args or [])],
            capture_output=True,
            text=True,
            timeout=timeout,
            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
        )
        return {"available": proc.returncode in (0, 1), "path": path, "path_exists": True}
    except Exception:
        return {"available": False, "path": path, "path_exists": True}


def gpu_info() -> dict:
    nvidia_smi = shutil.which("nvidia-smi")
    if not nvidia_smi:
        return {"available": False, "driver": None}
    try:
        proc = subprocess.run(
            [nvidia_smi, "--query-gpu=name,driver_version,memory.total", "--format=csv,noheader"],
            capture_output=True, text=True, timeout=8,
        )
        first = proc.stdout.strip().splitlines()[0] if proc.stdout.strip() else ""
        name, driver, mem = (first.split(",") + ["", "", ""])[:3]
        return {"available": proc.returncode == 0, "driver": driver.strip(), "name": name.strip(), "memory": mem.strip()}
    except Exception:
        return {"available": False, "driver": None}


def java_home() -> str | None:
    import os
    from pathlib import Path

    home = os.environ.get("JAVA_HOME")
    if home and (Path(home) / "bin" / "java.exe").is_file():
        return home
    return None


def environment_report() -> dict:
    deps = {
        name: _module_version(mod)
        for name, mod in [
            ("fastapi", "fastapi"),
            ("rdkit", "rdkit"),
            ("pandas", "pandas"),
            ("numpy", "numpy"),
        ]
    }
    vina = settings.tool("vina")
    smina = settings.tool("smina")
    vina_gpu = settings.tool("vina_gpu")
    gnina = settings.tool("gnina")
    p2rank = settings.tool("p2rank")

    engines = {
        "vina_cpu": {**_probe_executable(vina, ["--version"]), "required": True, "label": "AutoDock Vina (CPU)"},
        "smina": {**_probe_executable(smina, ["--version"]), "required": False, "label": "Smina (CPU)"},
        "vina_gpu": {**_probe_executable(vina_gpu), "required": False, "label": "Vina-GPU (CUDA)"},
        "gnina": {**_probe_executable(gnina, ["--version"]), "required": False, "label": "GNINA (CNN 打分)"},
        "p2rank": {"available": p2rank is not None, "path": str(p2rank) if p2rank else None, "required": False, "label": "P2Rank 口袋预测"},
    }

    assets = {
        name: bool(settings.checkpoint(name))
        for name in ("flashdiff", "diffprotacs", "planet", "unimol_docking")
    }
    vendors = {
        name: bool(settings.vendor_dir(name))
        for name in ("FlashDiff_app", "DiffPROTACs", "PLANET", "Uni-Mol/unimol_docking_v2", "p2rank_2.5")
    }

    gpu = gpu_info()
    docking_modes = []
    if engines["vina_gpu"]["available"] and gpu["available"]:
        docking_modes.append("vina_gpu")
    if engines["vina_cpu"]["available"]:
        docking_modes.append("vina_cpu")
    if engines["smina"]["available"]:
        docking_modes.append("smina")
    if engines["gnina"]["available"]:
        docking_modes.append("gnina")

    # 二期引擎可用性预览（实现前即可见资产就绪状态）
    phase2 = {}
    try:
        from ..engines.md import availability as md_avail
        from ..engines.knowledge import availability as kb_avail
        from ..engines.ppi import availability as ppi_avail
        from ..engines.geo import availability as geo_avail
        from ..engines.retro import availability as retro_avail

        phase2 = {
            "md": md_avail(),
            "knowledge": kb_avail(),
            "ppi": ppi_avail(),
            "geo": geo_avail(),
            "retrosynthesis": retro_avail(),
        }
    except Exception:
        phase2 = {}

    core_ready = bool(docking_modes)
    return {
        "status": "ready" if core_ready else "incomplete",
        "python": {
            "version": _module_version("sys") or "",
            "deps": deps,
            "deps_ok": all(deps.values()),
        },
        "gpu": gpu,
        "engines": engines,
        "docking_modes": docking_modes,
        "selected_docking_mode": docking_modes[0] if docking_modes else None,
        "assets": assets,
        "vendors": vendors,
        "phase2": phase2,
        "java_home": java_home(),
        "data_root": str(settings.data_root),
        "assets_root": str(settings.assets),
    }
