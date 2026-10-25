"""STRING API —— 蛋白互作网络（已知 + 预测）。"""
from __future__ import annotations

from ..cache import cache_get, cache_set
from ..config import STRING_API
from .http_client import request_json

_SERVICE = "string"
HUMAN_TAXID = 9606


async def get_interactions(gene_symbol: str, limit: int = 10) -> dict:
    """返回 top 互作伙伴及 combined_score（0-1000）。"""
    key = f"{gene_symbol}:{limit}"
    cached = cache_get(_SERVICE, key)
    if cached is not None:
        return cached
    params = {"identifiers": gene_symbol, "species": str(HUMAN_TAXID), "limit": str(limit)}
    text = (await _get_tsv("network", params)).strip().splitlines()
    partners = []
    for line in text:
        cols = line.split("\t")
        # network TSV: stringId_A, stringId_B, preferredName_A, preferredName_B,
        #              ncbiTaxonId, score(0~1)
        if len(cols) >= 6:
            score = float(cols[5])
            p1, p2 = cols[2], cols[3]
            partner = p2 if p1.upper() == gene_symbol.upper() else p1
            partners.append({"partner": partner, "score": round(score, 3)})
    partners.sort(key=lambda x: -x["score"])
    out = {"partners": partners[:limit]}
    cache_set(_SERVICE, key, out)
    return out


async def _get_tsv(endpoint: str, params: dict) -> str:
    import httpx
    from ..config import REQUEST_TIMEOUT, USER_AGENT

    url = f"{STRING_API}/tsv-no-header/{endpoint}"
    async with httpx.AsyncClient(timeout=REQUEST_TIMEOUT, follow_redirects=True) as client:
        last_err: Exception | None = None
        for attempt in range(3):
            try:
                resp = await client.get(url, params=params,
                                        headers={"User-Agent": USER_AGENT})
                resp.raise_for_status()
                return resp.text
            except httpx.HTTPError as exc:
                last_err = exc
                if attempt < 2:
                    import asyncio
                    await asyncio.sleep(1.5)
        raise last_err  # type: ignore[misc]
