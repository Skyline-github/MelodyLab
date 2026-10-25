"""GWAS Catalog V2 REST —— 人类遗传学信号（旧 REST 已弃用）。

V2 端点如不可用则优雅降级：Q1 的遗传学证据主要来自 Open Targets 的
genetic_association 分型，此模块是补充。
"""
from __future__ import annotations

from ..cache import FAILURE_TTL, cache_get, cache_set
from ..config import GWAS_CATALOG_V2
from .http_client import request_json

_SERVICE = "gwas"


async def get_associations_for_gene(gene_symbol: str, size: int = 10) -> dict:
    """按基因符号查 GWAS 关联（rsID、性状、p 值）。"""
    key = f"{gene_symbol}:{size}"
    cached = cache_get(_SERVICE, key)
    if cached is not None:
        return cached

    out: dict = {"available": False, "associations": [], "note": ""}
    try:
        params = {"q": f"{gene_symbol}", "size": str(size)}
        # V2 API：/associations 支持按基因搜索
        data = await request_json(
            "GET", f"{GWAS_CATALOG_V2}/associations", params=params
        )
        assocs = data.get("associations") or data.get("_embedded", {}).get("associations") or []
        rows = []
        for a in assocs[:size]:
            trait = ""
            loci = a.get("loci") or []
            if loci:
                genes = [g.get("geneName", "")
                         for au in loci for g in (au.get("authorReportedGenes") or [])]
                trait = ", ".join(x for x in genes if x)
            eff = a.get("effectorGenes") or []
            rows.append({
                "rsid": (a.get("variantId") or ""),
                "trait": (a.get("diseaseTrait") or {}).get("trait", ""),
                "pvalue": a.get("pvalue", {}).get("value") if isinstance(a.get("pvalue"), dict) else a.get("pvalue"),
                "genes": ", ".join(eff) or trait,
            })
        out = {"available": True, "associations": rows, "note": ""}
    except Exception as exc:  # noqa: BLE001 —— 优雅降级
        out["note"] = f"GWAS Catalog V2 查询失败，请通过链接手动核查：{exc}"
    # 失败结果只缓存 1 小时，网络恢复后可自动重查
    cache_set(_SERVICE, key, out, ttl=FAILURE_TTL if not out.get("available") else None)
    return out
