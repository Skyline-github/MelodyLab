"""共享异步 HTTP 客户端：超时、重试、UA。"""
from __future__ import annotations

import asyncio
from typing import Any

import httpx

from ..config import MAX_RETRIES, REQUEST_TIMEOUT, USER_AGENT


async def request_json(
    method: str,
    url: str,
    *,
    params: dict | None = None,
    json_body: dict | None = None,
    headers: dict | None = None,
) -> Any:
    """发起请求并解析 JSON；失败抛 httpx.HTTPError。"""
    merged = {"User-Agent": USER_AGENT, **(headers or {})}
    async with httpx.AsyncClient(timeout=REQUEST_TIMEOUT, follow_redirects=True) as client:
        last_err: Exception | None = None
        for attempt in range(MAX_RETRIES + 1):
            try:
                resp = await client.request(
                    method, url, params=params, json=json_body, headers=merged
                )
                resp.raise_for_status()
                return resp.json()
            except (httpx.HTTPError, ValueError) as exc:  # 网络错误 / 非 JSON 响应
                last_err = exc
                if attempt < MAX_RETRIES:
                    await asyncio.sleep(1.5 * (attempt + 1))
        raise last_err  # type: ignore[misc]
