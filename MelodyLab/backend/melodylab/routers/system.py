"""系统路由：/api/health、/api/environment、引擎路径配置。"""
from __future__ import annotations

import os
import platform

from fastapi import APIRouter, HTTPException

from ..config import settings
from ..core.sysinfo import environment_report

router = APIRouter(prefix="/api", tags=["system"])


@router.get("/health")
def health():
    return {"ok": True, "app": "MelodyLab", "pid": os.getpid()}


@router.get("/settings")
def get_settings():
    """当前设置：数据根目录及来源。"""
    from ..config import configured_data_root, data_root_marker

    env_override = os.environ.get("MELODYLAB_DATA_ROOT")
    return {
        "data_root": str(settings.data_root),
        "source": ("env" if env_override else
                   ("config" if configured_data_root() else "default")),
        "marker_file": str(data_root_marker()),
        "env_override": env_override,
        "note": "修改数据根目录后需重启后端生效",
    }


@router.post("/settings/data-root")
def set_data_root(payload: dict):
    """选择数据根目录（任务/报告/结果保存位置）。payload: {path: "..."}"""
    from ..config import data_root_marker

    raw = (payload or {}).get("path", "").strip()
    if not raw:
        raise HTTPException(400, "path 不能为空")
    p = Path(raw).expanduser()
    try:
        p.mkdir(parents=True, exist_ok=True)
        probe = p / ".melodylab-write-test"
        probe.write_text("ok", encoding="utf-8")
        probe.unlink()
    except OSError as exc:
        raise HTTPException(400, f"路径不可写: {exc}")
    marker = data_root_marker()
    marker.write_text(str(p.resolve()), encoding="utf-8")
    return {"ok": True, "data_root": str(p.resolve()), "marker": str(marker),
            "restart_required": True}


@router.get("/environment")
def environment():
    report = environment_report()
    report["platform"] = f"{platform.system()} {platform.release()}"
    return report


class ToolPathsPayload(dict):
    pass


@router.post("/tools")
def save_tool_paths(payload: dict):
    """持久化引擎路径覆盖（写入数据根 tools.env，进程内即时生效）。"""
    allowed = {"vina", "smina", "vina_gpu", "gnina", "p2rank"}
    lines = []
    for name, value in (payload or {}).items():
        if name not in allowed:
            continue
        os.environ[f"MELODYLAB_{name.upper()}_EXECUTABLE"] = value or ""
        if value:
            lines.append(f"{name}={value}")
    env_file = settings.data_root / "tools.env"
    env_file.write_text("\n".join(lines), encoding="utf-8")
    return {"ok": True, "saved": str(env_file)}


@router.get("/tools")
def get_tool_paths():
    out = {}
    for name in ("vina", "smina", "vina_gpu", "gnina", "p2rank", "gromacs"):
        p = settings.tool(name)
        out[name] = str(p) if p else None
    return out


@router.get("/environment/conda")
def conda_scan():
    """扫描本机 conda/系统 Python 环境与 MelodyLab 依赖就绪度（应用内版）。"""
    from ..core import envscan

    return envscan.scan_all()


@router.post("/environment/conda/install")
def conda_install(payload: dict):
    """向指定环境补装 MelodyLab 缺失依赖。payload: {python: "..."}"""
    from ..core import envscan

    python = (payload or {}).get("python")
    if not python:
        return {"ok": False, "error": "缺少 python 字段"}
    return envscan.install_missing(python)
