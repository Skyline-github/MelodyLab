"""MelodyLab 全局配置。

数据根目录的解析优先级：
  1. 环境变量 MELODYLAB_DATA_ROOT
  2. 项目根目录的 data_root.txt（设置页写入，跨重启持久）
  3. 桌面打包形态 / 用户主目录回退
引擎路径均支持环境变量覆盖，打包后（PyInstaller/Electron）与开发环境共用同一套解析逻辑。
"""
from __future__ import annotations

import os
import sys
from pathlib import Path


def data_root_marker() -> Path:
    """数据根目录的持久化标记文件（随应用存放，避免鸡生蛋问题）。"""
    return project_root() / "data_root.txt"


def configured_data_root() -> Path | None:
    """data_root.txt 中用户选择的路径；不存在或失效返回 None。"""
    marker = data_root_marker()
    if marker.is_file():
        try:
            p = Path(marker.read_text(encoding="utf-8").strip())
            if str(p) and p.is_dir():
                return p.resolve()
        except OSError:
            pass
    return None


def _default_data_root() -> Path:
    if env := os.environ.get("MELODYLAB_DATA_ROOT"):
        return Path(env).expanduser().resolve()
    if configured := configured_data_root():
        return configured
    if os.environ.get("MELODYLAB_DESKTOP") == "1":
        # 桌面打包形态：优先短路径，规避 Windows MAX_PATH 问题（同 DockScope 教训）
        system_drive = os.environ.get("SystemDrive", "C:")
        candidates = [
            Path(system_drive) / "MelodyLabData",
            Path.home() / "MelodyLabData",
        ]
    else:
        candidates = [Path.home() / "MelodyLabData"]
    for c in candidates:
        try:
            c.mkdir(parents=True, exist_ok=True)
            probe = c / f".write-test-{os.getpid()}"
            probe.write_text("ok")
            probe.unlink()
            return c
        except OSError:
            continue
    return candidates[-1]


def project_root() -> Path:
    """仓库根目录（backend/ 的上一级）；打包后为 assets 所在目录。"""
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent
    return Path(__file__).resolve().parents[2]


def assets_root() -> Path:
    env = os.environ.get("MELODYLAB_ASSETS")
    if env:
        return Path(env).expanduser().resolve()
    return project_root() / "assets"


def frontend_dist() -> Path | None:
    env = os.environ.get("MELODYLAB_FRONTEND_DIST")
    if env and Path(env).is_dir():
        return Path(env)
    dist = project_root() / "frontend" / "dist"
    return dist if dist.is_dir() else None


class Settings:
    def __init__(self) -> None:
        self.data_root = _default_data_root()
        self.assets = assets_root()
        self.jobs_dir = self.data_root / "jobs"
        self.logs_dir = self.data_root / "logs"
        self.cors_origins = [
            o.strip()
            for o in os.environ.get(
                "MELODYLAB_CORS_ORIGINS",
                "http://localhost:5173,http://127.0.0.1:5173",
            ).split(",")
            if o.strip()
        ]
        for d in (self.jobs_dir, self.logs_dir):
            d.mkdir(parents=True, exist_ok=True)

    # ---- 引擎可执行文件解析：env 覆盖 > assets 目录 > PATH ----
    def tool(self, name: str) -> Path | None:
        override = os.environ.get(f"MELODYLAB_{name.upper()}_EXECUTABLE")
        if override and Path(override).is_file():
            return Path(override)
        candidates = {
            "vina": ["tools/vina_1.2.7_win.exe", "tools/vina.exe", "tools/vina"],
            "smina": ["tools/smina-runtime/smina.exe", "tools/smina"],
            "vina_gpu": ["tools/Vina-GPU-CUDA.exe", "tools/Vina-GPU-v2.0-CUDA.exe"],
            "gnina": ["tools/gnina.exe"],
            "p2rank": ["vendor/p2rank_2.5/prank.bat", "vendor/p2rank_2.5/prank",
                       "p2rank_2.5/prank.bat", "p2rank_2.5/prank"],
            "gromacs": ["tools/gromacs/bin/gmx.exe", "tools/gromacs/bin/gmx"],
            "java": None,  # 只查 PATH/JAVA_HOME
        }.get(name, [])
        if candidates:
            for rel in candidates:
                p = self.assets / rel
                if p.is_file():
                    return p
        from shutil import which

        found = which(name)
        return Path(found) if found else None

    def vendor_dir(self, name: str) -> Path | None:
        """可选模型代码目录（FlashDiff / DiffPROTACs / PLANET / Uni-Mol）。"""
        p = self.assets / "vendor" / name
        return p if p.is_dir() else None

    def checkpoint(self, name: str) -> Path | None:
        table = {
            "flashdiff": "vendor/FlashDiff_app/data/trained_models/FlashDiff/checkpoints/FlashDiff.ckpt",
            "diffprotacs": "vendor/DiffPROTACs/checkpoints/geom_best.ckpt",
            "planet": "vendor/PLANET/PLANET.param",
            "unimol_docking": "vendor/Uni-Mol/unimol_docking_v2/unimol_docking_v2_240517.pt",
        }
        rel = table.get(name)
        if not rel:
            return None
        p = self.assets / rel
        return p if p.is_file() else None


settings = Settings()
