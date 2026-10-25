"""ChEMBL REST —— 靶点对应药物、作用机制与研发阶段（同家族已上市药物线索）。"""
from __future__ import annotations

from ..cache import cache_get, cache_set
from ..config import CHEMBL_API
from .http_client import request_json

_SERVICE = "chembl"
_PHASE_LABEL = {4: "已上市", 3: "III 期", 2: "II 期", 1: "I 期", 0: "临床前/其他",
                -1: "撤市/终止"}


async def _resolve_target_chembl_id(gene_symbol: str) -> str | None:
    key = f"resolve:{gene_symbol}"
    cached = cache_get(_SERVICE, key)
    if cached is not None:
        return cached or None
    data = await request_json(
        "GET", f"{CHEMBL_API}/target/search",
        params={"q": gene_symbol, "format": "json", "limit": "20"},
    )
    sym = gene_symbol.upper()

    def _score(t: dict) -> int:
        """匹配优先级：人源 > pref_name 精确 > 组件描述精确 > 组件同义词精确；
        同分时 SINGLE PROTEIN 优先于 PPI/其他类型。"""
        best = 0
        if (t.get("organism") or "").lower() == "homo sapiens":
            best += 40
        if (t.get("pref_name") or "").upper() == sym:
            return best + 30
        for c in t.get("target_components") or []:
            if (c.get("component_description") or "").upper() == sym:
                best = max(best, 20)
            syns = {s.get("component_synonym", "").upper()
                    for s in (c.get("target_component_synonyms") or [])}
            if sym in syns:
                best = max(best, 10)
        if t.get("target_type") == "SINGLE PROTEIN":
            best += 2
        return best

    candidates = sorted(data.get("targets", []), key=lambda t: -_score(t))
    chosen = next((t for t in candidates if _score(t) > 0), None)
    chembl_id = chosen["target_chembl_id"] if chosen else ""
    cache_set(_SERVICE, key, chembl_id)
    return chembl_id or None


async def get_drugs(gene_symbol: str) -> dict:
    """作用于此靶点的药物及其最高研发阶段。"""
    chembl_id = await _resolve_target_chembl_id(gene_symbol)
    if not chembl_id:
        return {"resolved": False, "drugs": [], "note": "ChEMBL 未解析到该靶点"}
    data = await request_json(
        "GET", f"{CHEMBL_API}/mechanism",
        params={"target_chembl_id": chembl_id, "format": "json", "limit": "50"},
    )
    drugs: dict[str, dict] = {}
    for m in data.get("mechanisms", []):
        name = m.get("molecule_chembl_id", "")
        d = drugs.setdefault(name, {
            "drug": m.get("molecule_chembl_id", ""),
            "name": "",
            "action": m.get("action_type", ""),
            "moa": m.get("mechanism_of_action", ""),
            "max_phase": 0,
        })
        d["moa"] = d["moa"] or m.get("mechanism_of_action", "")
        d["action"] = d["action"] or m.get("action_type", "")
    # 拉取分子信息拿名称与研发阶段（单分子端点直接返回对象）
    for chembl_mol in list(drugs)[:30]:
        try:
            m = await request_json(
                "GET", f"{CHEMBL_API}/molecule/{chembl_mol}.json"
            )
            drugs[chembl_mol]["name"] = m.get("pref_name", chembl_mol)
            drugs[chembl_mol]["max_phase"] = int(float(m.get("max_phase") or 0))
            drugs[chembl_mol]["withdrawn"] = bool(m.get("withdrawn_flag"))
        except Exception:  # noqa: BLE001
            continue
    rows = sorted(drugs.values(), key=lambda d: -d["max_phase"])
    approved = [d for d in rows if d["max_phase"] == 4]
    return {"resolved": True, "chembl_id": chembl_id, "drugs": rows,
            "approved_count": len(approved), "total_count": len(rows)}
