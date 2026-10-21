"""GROMACS 命令运行与 XVG 解析。

所有 gmx 子进程 stdin 置 DEVNULL：pdb2gmx 等命令若因参数不全进入交互
提问会立刻报错退出，而不是挂死任务线程。
"""
from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path

from ...config import settings


def gmx_exe() -> Path:
    p = settings.tool("gromacs")
    if p is None:
        raise RuntimeError("未找到 GROMACS（gmx.exe）；请运行 scripts/setup_assets.py --with-gromacs")
    return p


def run_gmx(args: list[str], cwd: Path, log_path: Path | None = None,
            stdin_text: str | None = None, timeout: int = 1800,
            cancel_event=None) -> subprocess.CompletedProcess:
    """运行 gmx <args>，stdout/stderr 落盘并返回结果。

    mdrun 等长命令：轮询等待，期间响应 cancel_event（暂停/停止即杀进程）。
    """
    import time as _time

    exe = gmx_exe()
    cmd = [str(exe), *args]
    cwd.mkdir(parents=True, exist_ok=True)
    log_path = log_path or (cwd / ("gmx_" + re.sub(r"[^a-z0-9]+", "_", args[0].lower()) + ".log"))
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w", encoding="utf-8", errors="replace") as log:
        proc = subprocess.Popen(
            cmd, cwd=str(cwd), stdout=log, stderr=subprocess.STDOUT,
            stdin=subprocess.PIPE if stdin_text is not None else subprocess.DEVNULL,
            text=(stdin_text is not None),
            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
        )
        try:
            if stdin_text is not None:
                proc.communicate(stdin_text, timeout=timeout)
            else:
                deadline = _time.monotonic() + timeout
                while True:
                    rc = proc.poll()
                    if rc is not None:
                        break
                    if cancel_event is not None and getattr(cancel_event, "is_set", None) and cancel_event.is_set():
                        proc.kill()
                        raise InterruptedError("MD 任务被取消（模拟进程已终止）")
                    if _time.monotonic() > deadline:
                        proc.kill()
                        raise TimeoutError(f"gmx 超时（>{timeout}s）: {' '.join(args[:3])}")
                    _time.sleep(1.0)
            rc = proc.returncode
        except subprocess.TimeoutExpired:
            proc.kill()
            raise TimeoutError(f"gmx 超时（>{timeout}s）: {' '.join(args[:3])}")
    if rc != 0:
        tail = Path(log_path).read_text(encoding="utf-8", errors="replace")[-1500:]
        raise RuntimeError(f"gmx {' '.join(args[:2])} 失败（rc={rc}）:\n{tail}")
    return subprocess.CompletedProcess(cmd, rc, "", "")


def parse_xvg(path: Path) -> dict:
    """解析 xvg → {legend, x_label, columns, data: [[...]]}。"""
    legends, data = [], []
    x_label, col_labels = "x", []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("@"):
            body = line[1:].strip()
            if body.startswith("title"):
                legends.append(body.split('"')[1] if '"' in body else body[5:].strip())
            elif body.startswith("xaxis label"):
                x_label = body.split('"')[1] if '"' in body else ""
            elif body.startswith(("s", "yaxis label")) and '"' in body:
                col_labels.append(body.split('"')[1])
            continue
        parts = line.split()
        try:
            data.append([float(p) for p in parts])
        except ValueError:
            continue
    return {"legend": legends[0] if legends else Path(path).stem,
            "x_label": x_label, "col_labels": col_labels, "data": data}


def xvg_series(path: Path) -> dict:
    """解析 xvg → {time: [...], series: {名称: [...]}}（用于 JSON 报告与画图）。"""
    d = parse_xvg(path)
    if not d["data"]:
        return {"time": [], "series": {}}
    cols = list(zip(*d["data"]))
    labels = d["col_labels"] or [f"s{i}" for i in range(len(cols) - 1)]
    return {"time": list(cols[0]), "x_label": d["x_label"],
            "series": {labels[i]: list(cols[i + 1]) for i in range(len(cols) - 1)}}


def xvg_last_average(path: Path, tail_frac: float = 0.2) -> float | None:
    """取 xvg 各列最后 tail_frac 区间的均值（第 0 列除外）。"""
    s = xvg_series(path)
    names = list(s["series"].keys())
    if not names:
        return None
    vals = s["series"][names[0]]
    n = max(1, int(len(vals) * tail_frac))
    return sum(vals[-n:]) / n


def gpu_available() -> bool:
    from shutil import which

    smi = which("nvidia-smi")
    if not smi:
        return False
    try:
        return subprocess.run([smi], capture_output=True, timeout=6).returncode == 0
    except Exception:
        return False
