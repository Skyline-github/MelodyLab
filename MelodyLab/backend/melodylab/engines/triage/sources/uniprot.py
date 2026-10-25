"""UniProt REST —— 靶点蛋白基础注释：功能、亚细胞定位、Ensembl 交叉引用。"""
from __future__ import annotations

from ..cache import cache_get, cache_set
from ..config import UNIPROT_API
from .http_client import request_json

_SERVICE = "uniprot"

_FIELDS = ("accession,id,protein_name,gene_names,cc_function,"
           "cc_subcellular_location,cc_disease,xref_ensembl")


async def get_protein(gene_symbol: str, organism_id: int = 9606) -> dict | None:
    """基因符号 → 人源 reviewed 蛋白条目（Swiss-Prot）。"""
    key = f"{gene_symbol}:{organism_id}"
    cached = cache_get(_SERVICE, key)
    if cached is not None:
        return cached
    params = {
        "query": f"(gene:{gene_symbol}) AND (organism_id:{organism_id}) AND reviewed:true",
        "format": "json",
        "size": "1",
        "fields": _FIELDS,
    }
    data = await request_json("GET", UNIPROT_API, params=params)
    results = data.get("results") or []
    if not results:
        return None
    r = results[0]
    comments = {c["commentType"]: c.get("texts", []) for c in r.get("comments", [])}

    def _texts(ctype: str) -> list[str]:
        out = []
        for t in comments.get(ctype, []):
            v = t.get("value", "")
            if v and "By similarity" not in v[:20]:
                out.append(v)
        return out[:3]

    ensembl = next(
        (x["id"] for x in r.get("uniProtKBCrossReferences", [])
         if x["database"] == "Ensembl" and x["id"].startswith("ENSG")),
        None,
    )
    info = {
        "uniprot_id": r["primaryAccession"],
        "protein_name": (r.get("proteinDescription", {}).get("recommendedName", {})
                         .get("fullName", {}).get("value", "")),
        "ensembl_id": ensembl,
        "function": _texts("FUNCTION"),
        "subcellular_location": _texts("SUBCELLULAR LOCATION"),
        "disease_notes": _texts("DISEASE"),
    }
    cache_set(_SERVICE, key, info)
    return info
