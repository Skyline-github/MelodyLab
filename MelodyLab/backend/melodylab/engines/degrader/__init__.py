"""降解剂评估引擎（PROTAC DC50 共识排序）。

打包的 PROTAC-DB 快照 + 骨架 5 折交叉验证的三族模型（MLP 集成 / ExtraTrees / HistGB）；
模型缓存于 <data_root>/cache/degrader/，报告落 <data_root>/degrader_reports/。
定位是排序参考工具：输出共识排名 + 90% 区间 + 适用域，不是绝对值预测。
"""
from __future__ import annotations

import importlib.util

from ...config import settings
from . import config as cfg


def model_info() -> dict | None:
    """缓存模型元信息（sidecar JSON，避免加载完整 bundle）。"""
    import json

    try:
        sidecar = cfg.cache_path().with_name(cfg.cache_path().name + ".meta.json")
    except OSError:
        return None
    if sidecar.is_file():
        try:
            return json.loads(sidecar.read_text(encoding="utf-8"))
        except Exception:  # noqa: BLE001
            return None
    return None


def availability() -> dict:
    libs = {name: importlib.util.find_spec(name) is not None
            for name in ("rdkit", "sklearn", "torch", "joblib", "scipy", "pandas")}
    data_ok = cfg.DB_FILE.is_file()
    missing = [k for k, ok in libs.items() if not ok]
    info = model_info() or {}
    available = data_ok and not missing
    hint = None
    if missing:
        hint = f"缺少依赖：{', '.join(missing)}（pip install {' '.join(missing)}）"
    elif not data_ok:
        hint = f"训练数据缺失：{cfg.DB_FILE}"
    return {
        "engine": "degrader",
        "available": available,
        "hint": hint,
        "missing_libs": missing,
        "db_snapshot": cfg.DB_SNAPSHOT if data_ok else "",
        "db_file": str(cfg.DB_FILE) if data_ok else "",
        "db_rows": 5388,
        "model_ready": bool(info),
        "model": info,
        "cache_dir": str(cfg.CACHE_DIR),
        "reports_dir": str(cfg.REPORTS_DIR),
    }
