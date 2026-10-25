"""靶点立项评估引擎（立项五问，移植自 target-triage）。

可用性探测 + worker（service.run_triage_job）由 JobManager 调度；
外部数据源在 sources/，阈值判据在 rules.yaml。
"""
from __future__ import annotations

import importlib.util
import os

from ...config import settings


def availability() -> dict:
    has_httpx = importlib.util.find_spec("httpx") is not None
    has_yaml = importlib.util.find_spec("yaml") is not None
    hpa_ready = (settings.data_root / "cache" / "triage" / "normal_tissue_v22.done").is_file()
    token_file = (settings.data_root / "oncokb_token.txt").is_file()
    return {
        "engine": "triage",
        "available": has_httpx and has_yaml,
        "oncokb_configured": os.getenv("ONCOKB_TOKEN", "").strip() != "" or token_file,
        "hpa_bulk_ready": hpa_ready,
        "cache_dir": str(settings.data_root / "cache" / "triage"),
        "hint": None if (has_httpx and has_yaml) else "需要 httpx 与 pyyaml（pip install httpx pyyaml）",
    }
