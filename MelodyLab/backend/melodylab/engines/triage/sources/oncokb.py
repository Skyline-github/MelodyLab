"""OncoKB —— 精准肿瘤学伴随诊断/治疗线索（可选，需免费学术 token）。

无 token 时调用方应降级为手动核查清单。
"""
from __future__ import annotations

from ..cache import FAILURE_TTL, cache_get, cache_set
from ..config import ONCOKB_API, ONCOKB_TOKEN
from .http_client import request_json

_SERVICE = "oncokb"


def is_configured() -> bool:
    return bool(ONCOKB_TOKEN)


async def get_gene_actionability(gene_symbol: str) -> dict:
    """查询基因在 OncoKB 中的治疗级别信息（伴随诊断/药物线索）。"""
    if not is_configured():
        return {"available": False,
                "note": "未配置 ONCOKB_TOKEN（oncoKB.org 免费申请学术 token 后启用）"}
    key = gene_symbol
    cached = cache_get(_SERVICE, key)
    if cached is not None:
        return cached
    try:
        genes = await request_json(
            "GET", f"{ONCOKB_API}/genes/lookup",
            params={"query": gene_symbol},
            headers={"Authorization": f"Bearer {ONCOKB_TOKEN}"},
        )
        if not genes:
            return {"available": False, "note": "OncoKB 未收录该基因"}
        gene_id = genes[0]["entrezGeneId"]
        # 治疗与伴随诊断相关变异（示例：取最高级别 therapeutic 变异）
        variants = await request_json(
            "GET", f"{ONCOKB_API}/genes/{gene_id}/variants",
            headers={"Authorization": f"Bearer {ONCOKB_TOKEN}"},
        )
        top = []
        for v in variants[:20]:
            top.append({"variant": v.get("name", ""),
                        "oncogenic": v.get("oncogenic", ""),
                        "highest_level": v.get("highestSensitiveLevel")
                        or v.get("highestResistanceLevel") or ""})
        out = {"available": True, "gene": genes[0].get("hugoSymbol", gene_symbol),
               "variants": top}
    except Exception as exc:  # noqa: BLE001
        out = {"available": False, "note": f"OncoKB 查询失败：{exc}",
               "_failed": True}
    cache_set(_SERVICE, key, out, ttl=FAILURE_TTL if out.get("_failed") else None)
    out.pop("_failed", None)
    return out
