"""简单 JSON 文件缓存：按 (service, key) 哈希落盘，带 TTL。

与原 target-triage 版的差异：cache_set 支持可选 ttl，
失败结果用短 TTL（CACHE_TTL_FAILURE），避免一次网络抖动把
"信息不足"钉死一周。
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path
from typing import Any

from .config import CACHE_DIR, CACHE_TTL, CACHE_TTL_FAILURE


def _cache_path(service: str, key: str) -> Path:
    digest = hashlib.sha256(key.encode()).hexdigest()[:16]
    return CACHE_DIR / f"{service}_{digest}.json"


def cache_get(service: str, key: str) -> Any | None:
    path = _cache_path(service, key)
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        return None
    if time.time() - payload.get("ts", 0) > payload.get("ttl", CACHE_TTL):
        return None
    return payload.get("data")


def cache_set(service: str, key: str, data: Any, ttl: int | None = None) -> None:
    path = _cache_path(service, key)
    tmp = path.with_suffix(".tmp")
    payload = {"ts": time.time(), "ttl": ttl or CACHE_TTL, "data": data}
    tmp.write_text(json.dumps(payload, ensure_ascii=False))
    tmp.replace(path)


# 失败场景统一用这个（语义显式，便于各 source 调用）
FAILURE_TTL = CACHE_TTL_FAILURE
