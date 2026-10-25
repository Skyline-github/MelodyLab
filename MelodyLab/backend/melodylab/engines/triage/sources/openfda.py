"""openFDA drug/label —— 药物标签中提及该靶点/突变的产品线索。"""
from __future__ import annotations

from ..cache import FAILURE_TTL, cache_get, cache_set
from ..config import OPENFDA_DRUG_LABEL
from .http_client import request_json

_SERVICE = "openfda"


async def get_label_mentions(term: str, limit: int = 10) -> dict:
    """在药物标签（适应症/机制描述）全文中检索该基因的提及。"""
    key = f"{term}:{limit}"
    cached = cache_get(_SERVICE, key)
    if cached is not None:
        return cached
    query = f'indications_and_usage:"{term}"'
    out: dict = {"total": 0, "labels": [], "note": ""}
    try:
        data = await request_json(
            "GET", OPENFDA_DRUG_LABEL,
            params={"search": query, "limit": str(limit)},
        )
        out["total"] = data.get("meta", {}).get("results", {}).get("total", 0)
        for r in data.get("results", [])[:limit]:
            ofda = r.get("openfda", {})
            labels = r.get("indications_and_usage", [""])
            out["labels"].append({
                "brand": (ofda.get("brand_name") or [""])[0],
                "manufacturer": (ofda.get("manufacturer_name") or [""])[0],
                "indication_snippet": (labels[0][:160] + "…") if labels and len(labels[0]) > 160 else (labels[0] if labels else ""),
            })
    except Exception as exc:  # noqa: BLE001 —— openFDA 限流或无结果时降级
        out["note"] = f"openFDA 查询失败（限流 240 次/分钟）：{exc}"
    # 失败（限流/网络）只缓存 1 小时；"total=0 但查询成功"是真负结果，正常缓存
    failed = bool(out["note"])
    cache_set(_SERVICE, key, out, ttl=FAILURE_TTL if failed else None)
    return out
