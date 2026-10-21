"""本地 Python 环境检索器（纯标准库，无任何第三方依赖）。

在目标机器上运行：枚举 conda 环境与系统 Python，逐个探测 MelodyLab
所需的依赖，给出最佳环境推荐与缺失包的补装命令。
被 scripts/doctor.py（命令行）与 /api/environment/conda（应用内）共用。
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

# import 名 → pip 包名
REQUIRED = {
    "fastapi": "fastapi",
    "uvicorn": "uvicorn",
    "multipart": "python-multipart",
    "pydantic": "pydantic",
    "rdkit": "rdkit",
    "numpy": "numpy",
    "scipy": "scipy",
    "pandas": "pandas",
    "yaml": "pyyaml",
    "meeko": "meeko",
    "gemmi": "gemmi",
}
MIN_PY = (3, 10)


def _run(python: str | Path, code: str, timeout: float = 30.0) -> tuple[int, str]:
    try:
        proc = subprocess.run(
            [str(python), "-c", code], capture_output=True, text=True, timeout=timeout,
            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
        )
        return proc.returncode, (proc.stdout or "").strip()
    except Exception as exc:
        return 1, str(exc)


def find_conda_envs() -> list[dict]:
    """返回 [{name, python, version}]；base 环境含在内。"""
    out: list[dict] = []
    candidates: list[tuple[str, Path]] = []
    conda = None
    from shutil import which

    conda = which("conda")
    env_paths: list[Path] = []
    if conda:
        try:
            proc = subprocess.run([conda, "env", "list", "--json"],
                                  capture_output=True, text=True, timeout=30)
            data = json.loads(proc.stdout or "{}")
            env_paths = [Path(p) for p in data.get("envs", [])]
        except Exception:
            env_paths = []
    if not env_paths:  # conda 不在 PATH：探测常见安装位置
        home = Path.home()
        for root in (home / "anaconda3", home / "miniconda3", home / "miniforge3",
                     Path("C:/ProgramData/anaconda3"), Path("C:/anaconda3"),
                     Path("D:/anaconda3"), Path("D:/miniconda3")):
            if root.is_dir():
                env_paths = [root] + sorted((root / "envs").glob("*")) if (root / "envs").is_dir() else [root]
                break

    for p in env_paths:
        if not p.is_dir():
            continue
        exe = p / ("python.exe" if os.name == "nt" else "bin/python")
        if not exe.is_file():
            continue
        name = p.name if p.parent.name.lower() == "envs" else ("base" if (p / ("python.exe" if os.name == "nt" else "bin/python")).is_file() else p.name)
        candidates.append((name, exe))

    # 当前解释器（非 conda 时也有值）
    try:
        cur = Path(sys.executable).resolve()
        if not any(ex.resolve() == cur for _, ex in candidates):
            candidates.append(("system-python", cur))
    except Exception:
        pass
    return [{"name": n, "python": str(ex)} for n, ex in candidates]


def probe_env(python: str | Path) -> dict:
    """探测单个环境：Python 版本 + 各依赖可用性。"""
    rc, ver = _run(python, "import sys; print('%d.%d.%d' % sys.version_info[:3])")
    version = ver if rc == 0 else "?"
    try:
        ok_py = rc == 0 and tuple(int(x) for x in version.split(".")[:2]) >= MIN_PY
    except Exception:
        ok_py = False
    packages = {}
    for mod, pip_name in REQUIRED.items():
        prc, pout = _run(python, f"import {mod}; print(getattr({mod}, '__version__', 'ok'))", timeout=45)
        packages[pip_name] = {"ok": prc == 0, "version": pout if prc == 0 else None}
    missing = [pip for pip, info in packages.items() if not info["ok"]]
    present = [pip for pip, info in packages.items() if info["ok"]]
    return {
        "python": str(python),
        "python_version": version,
        "python_ok": ok_py,
        "packages": packages,
        "missing": missing,
        "present_count": len(present),
        "ready": ok_py and not missing,
        "install_cmd": (f'"{python}" -m pip install {" ".join(missing)}' if missing else ""),
    }


def scan_all() -> dict:
    envs = []
    for e in find_conda_envs():
        info = probe_env(e["python"])
        info["name"] = e["name"]
        envs.append(info)
    # 排序：Python 版本达标优先 → 已装包多优先
    envs.sort(key=lambda x: (not x["python_ok"], -x["present_count"], x["name"]))
    best = envs[0] if envs and envs[0]["python_ok"] else None
    return {
        "scan_target": f"{os.name}/{'win' if os.name == 'nt' else 'unix'}",
        "envs": envs,
        "recommended": best,
        "any_ready": any(e["ready"] for e in envs),
    }


def install_missing(python: str | Path, extra_args: list[str] | None = None) -> dict:
    """向选定环境补装缺失包（pip）。"""
    info = probe_env(python)
    if info["ready"]:
        return {"ok": True, "message": "环境已就绪，无需安装", "installed": []}
    cmd = [str(python), "-m", "pip", "install", *(extra_args or []), *info["missing"]]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    tail = (proc.stdout or "").strip().splitlines()[-12:]
    if proc.returncode == 0:
        return {"ok": True, "installed": info["missing"], "output_tail": tail}
    return {"ok": False, "error": (proc.stderr or proc.stdout or "")[-1500:], "output_tail": tail}
