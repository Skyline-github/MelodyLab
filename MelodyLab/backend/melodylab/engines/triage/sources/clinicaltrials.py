"""ClinicalTrials.gov API v2 —— 竞争格局：该靶点的临床试验统计。"""
from __future__ import annotations

import asyncio

from ..cache import cache_get, cache_set
from ..config import CLINICALTRIALS_API, REQUEST_TIMEOUT, USER_AGENT
from .http_client import request_json

_SERVICE = "clinicaltrials"


async def get_trials(gene_symbol: str, indication: str = "") -> dict:
    """查询该靶点（可选限定适应症）的临床试验：总数、状态分布、终止原因。"""
    key = f"{gene_symbol}|{indication}"
    cached = cache_get(_SERVICE, key)
    if cached is not None:
        return cached

    import httpx

    params = {
        "query.intr": gene_symbol,
        "pageSize": "100",
        "fields": "NCTId|BriefTitle|OverallStatus|WhyStopped|Condition|LeadSponsorName",
        "countTotal": "true",
    }
    if indication:
        params["query.cond"] = indication
    async with httpx.AsyncClient(timeout=REQUEST_TIMEOUT, follow_redirects=True) as client:
        resp = await client.get(CLINICALTRIALS_API, params=params,
                                headers={"User-Agent": USER_AGENT})
        resp.raise_for_status()
        data = resp.json()

    total = int(data.get("totalCount", 0))
    studies = data.get("studies", []) or []
    status_dist: dict[str, int] = {}
    terminated: list[dict] = []
    conditions: dict[str, int] = {}
    for s in studies:
        proto = s.get("protocolSection", {})
        status = (proto.get("statusModule", {}).get("overallStatus")
                  or "UNKNOWN").replace("_", " ")
        status_dist[status] = status_dist.get(status, 0) + 1
        for c in proto.get("conditionsModule", {}).get("conditions", []) or []:
            conditions[c] = conditions.get(c, 0) + 1
        why = proto.get("statusModule", {}).get("whyStopped", "")
        if why:
            terminated.append({
                "nct_id": proto.get("identificationModule", {}).get("nctId", ""),
                "title": proto.get("identificationModule", {}).get("briefTitle", ""),
                "why": why[:200],
            })
    top_conditions = sorted(conditions.items(), key=lambda kv: -kv[1])[:10]
    out = {
        "total": total,
        "returned": len(studies),
        "status_dist": status_dist,
        "terminated": terminated[:10],
        "terminated_total": sum(1 for k, v in status_dist.items()
                                if "TERMINATED" in k.upper() or "WITHDRAWN" in k.upper()),
        "top_conditions": [{"condition": c, "count": n} for c, n in top_conditions],
    }
    cache_set(_SERVICE, key, out)
    return out
