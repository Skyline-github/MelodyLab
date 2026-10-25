"""Open Targets Platform GraphQL API —— 评估主干数据源。

提供：靶点注释与 tractability、靶点-疾病关联（含遗传学证据分型）、
safety liabilities、遗传约束（LOF z 分数）、必需基因标记与
小鼠敲除表型（Genetics 独立 API 已退役，相关数据并入 Platform）。
"""
from __future__ import annotations

import json
from typing import Any

from ..cache import cache_get, cache_set
from ..config import OPEN_TARGETS_GRAPHQL
from .http_client import request_json

_SERVICE = "opentargets"


async def _graphql(query: str, variables: dict) -> Any:
    key = query + json.dumps(variables, sort_keys=True)
    cached = cache_get(_SERVICE, key)
    if cached is not None:
        return cached
    data = await request_json(
        "POST", OPEN_TARGETS_GRAPHQL, json_body={"query": query, "variables": variables}
    )
    if data.get("errors"):
        raise RuntimeError(f"Open Targets GraphQL error: {data['errors'][:1]}")
    result = data["data"]
    cache_set(_SERVICE, key, result)
    return result


async def search_target(symbol: str) -> dict | None:
    """基因符号 → {ensemblId, approvedName}；找不到返回 None。"""
    query = """
    query($q: String!) {
      search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 5}) {
        hits { id entity name }
      }
    }"""
    data = await _graphql(query, {"q": symbol})
    for hit in data["search"]["hits"]:
        if hit["entity"] == "target":
            return {"ensembl_id": hit["id"], "approved_name": hit["name"]}
    return None


async def get_target_overview(ensembl_id: str) -> dict:
    """靶点基础注释：符号、全名、生物类型、可药性、遗传约束、必需性。"""
    query = """
    query($id: String!) {
      target(ensemblId: $id) {
        id approvedSymbol approvedName biotype isEssential
        geneticConstraint { constraintType score }
        tractability { label modality value }
      }
    }"""
    data = await _graphql(query, {"id": ensembl_id})
    t = data["target"]
    constraint = {c["constraintType"]: c["score"]
                  for c in (t.get("geneticConstraint") or [])}
    return {
        "symbol": t["approvedSymbol"],
        "name": t["approvedName"],
        "biotype": t.get("biotype", ""),
        "is_essential": t.get("isEssential"),
        # LOF 约束 z 分数：越高越受约束（z >= 3.09 约等于经典高约束阈值）
        "lof_z": constraint.get("lof"),
        "tractability": [
            {"modality": x["modality"], "label": x["label"], "value": x["value"]}
            for x in (t.get("tractability") or [])
            if x.get("value")
        ],
    }


async def get_associations(ensembl_id: str, top_n: int = 15) -> list[dict]:
    """靶点-疾病关联分数（多源证据整合），按分排序。"""
    query = """
    query($id: String!, $n: Int!) {
      target(ensemblId: $id) {
        associatedDiseases(page: {index: 0, size: $n}) {
          rows {
            disease { id name }
            score
            datatypeScores { id score }
          }
        }
      }
    }"""
    data = await _graphql(query, {"id": ensembl_id, "n": top_n})
    rows = data["target"]["associatedDiseases"]["rows"]
    out = []
    for r in rows:
        dts = {d["id"]: d["score"] for d in (r.get("datatypeScores") or [])}
        out.append(
            {
                "disease": r["disease"]["name"],
                "efo_id": r["disease"]["id"],
                "overall": round(r["score"], 3),
                "genetic": round(dts.get("genetic_association", 0), 3),
                "somatic": round(dts.get("somatic_mutation", 0), 3),
                "drugs": round(dts.get("known_drug", 0), 3),
                "literature": round(dts.get("literature", 0), 3),
            }
        )
    return out


async def get_safety(ensembl_id: str) -> list[dict]:
    """已知安全事件（来自 PharmGKB、ToxCast、文献 curated 等）。"""
    query = """
    query($id: String!) {
      target(ensemblId: $id) {
        safetyLiabilities { event datasource url }
      }
    }"""
    data = await _graphql(query, {"id": ensembl_id})
    return [
        {"event": s["event"], "datasource": s.get("datasource", ""), "url": s.get("url", "")}
        for s in (data["target"].get("safetyLiabilities") or [])
    ]


_LETHAL_KEYWORDS = ("lethal", "lethality")


async def get_mouse_knockout(ensembl_id: str, max_rows: int = 200) -> dict:
    """小鼠敲除表型：是否记录致死类表型（必需基因的强证据）。"""
    query = """
    query($id: String!) {
      target(ensemblId: $id) {
        mousePhenotypes {
          modelPhenotypeId
          modelPhenotypeLabel
          modelPhenotypeClasses { id label }
          targetInModelEnsemblId
        }
      }
    }"""
    data = await _graphql(query, {"id": ensembl_id})
    rows = data["target"].get("mousePhenotypes") or []
    phenotypes: list[dict] = []
    lethal_labels: list[str] = []
    for row in rows[:max_rows]:
        classes = [c["label"] for c in (row.get("modelPhenotypeClasses") or [])]
        label = row.get("modelPhenotypeLabel", "")
        phenotypes.append({"label": label, "classes": classes})
        hay = (label + " " + " ".join(classes)).lower()
        if any(k in hay for k in _LETHAL_KEYWORDS):
            lethal_labels.append(label)
    return {"has_knockout_data": bool(rows), "lethal_phenotype": bool(lethal_labels),
            "lethal_examples": sorted(set(lethal_labels))[:5],
            "phenotype_count": len(rows),
            "phenotypes": phenotypes[:60]}
