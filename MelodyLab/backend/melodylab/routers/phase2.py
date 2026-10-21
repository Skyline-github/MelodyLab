"""二期路由：知识库 / PPI / GEO / 逆合成。"""
from __future__ import annotations

import json
from pathlib import Path

from fastapi import APIRouter, File, Form, HTTPException, UploadFile

from ..config import settings
from ..core.jobs import manager
from ..engines.geo import pipeline as geo_pipe
from ..engines.knowledge import service as knowledge
from ..engines.ppi import service as ppi
from ..engines.retro import service as retro

router = APIRouter(prefix="/api", tags=["phase2"])


# ================= 知识库检索 =================
@router.get("/knowledge/availability")
def knowledge_availability():
    import importlib.util

    has_http = importlib.util.find_spec("httpx") is not None
    return {"engine": "knowledge", "available": has_http,
            "sources": ["chembl", "pubchem", "uniprot", "rcsb", "alphafold", "zinc22", "geo"],
            "hint": None if has_http else "需要 httpx（pip install httpx）"}


@router.post("/knowledge/activity/search")
def knowledge_activity(
    query: str = Form(...),
    query_type: str = Form("gene"),
    activity_types: str = Form(""),
    max_records: int = Form(50),
):
    record = manager.create("knowledge", f"活性检索: {query}", {
        "task": "activity", "query": query, "query_type": query_type,
        "activity_types": [t.strip() for t in activity_types.split(",") if t.strip()],
        "max_records": max_records,
    })
    return manager.submit(record["id"], knowledge.run_knowledge_job)


@router.post("/knowledge/structures/search")
def knowledge_structures(
    query: str = Form(...),
    resolution: float | None = Form(None),
    require_ligand: bool = Form(False),
    method: str = Form(""),
    max_records: int = Form(20),
):
    record = manager.create("knowledge", f"结构检索: {query}", {
        "task": "structures", "query": query, "resolution": resolution,
        "require_ligand": require_ligand, "method": method or None,
        "max_records": max_records,
    })
    return manager.submit(record["id"], knowledge.run_knowledge_job)


@router.post("/knowledge/zinc/search")
def knowledge_zinc(smiles: str = Form(...), max_records: int = Form(50)):
    record = manager.create("knowledge", f"ZINC 相似检索", {
        "task": "zinc", "smiles": smiles, "max_records": max_records})
    return manager.submit(record["id"], knowledge.run_knowledge_job)


@router.post("/knowledge/geo/search")
def knowledge_geo(query: str = Form(...), max_records: int = Form(20)):
    record = manager.create("knowledge", f"GEO 检索: {query}", {
        "task": "geo", "query": query, "max_records": max_records})
    return manager.submit(record["id"], knowledge.run_knowledge_job)


@router.get("/knowledge/jobs/{job_id}")
def knowledge_job(job_id: str):
    rec = manager.get(job_id)
    if rec is None or rec["kind"] != "knowledge":
        raise HTTPException(404, "知识库任务不存在")
    return rec


# ================= PPI 网络 =================
@router.post("/ppi/networks")
def ppi_create(
    identifiers: str = Form(""),
    identifiers_file: UploadFile | None = File(None),
    species: int = Form(9606),
    required_score: int = Form(400),
    include_partners: bool = Form(False),
    max_partners: int = Form(20),
    network_type: str = Form("functional"),
    name: str = Form(""),
):
    ids = [x.strip() for x in identifiers.replace(",", "\n").splitlines() if x.strip()]
    if identifiers_file is not None and identifiers_file.filename:
        raw = identifiers_file.file.read().decode("utf-8", errors="replace")
        ids += [x.strip() for x in raw.replace(",", "\n").splitlines() if x.strip()]
    ids = list(dict.fromkeys(ids))
    if not ids:
        raise HTTPException(400, "请提供蛋白列表（文本或文件）")
    record = manager.create("ppi", name or f"PPI 网络（{len(ids)} 蛋白）", {
        "identifiers": ids, "species": species, "required_score": required_score,
        "include_partners": include_partners, "max_partners": max_partners,
        "network_type": network_type,
    })
    return manager.submit(record["id"], _run_ppi_job)


def _run_ppi_job(record: dict, control) -> dict:
    params = record["params"]
    job_dir = settings.jobs_dir / record["id"]
    outputs = job_dir / "outputs"
    outputs.mkdir(parents=True, exist_ok=True)
    control.checkpoint(progress=10, phase="string", message="查询 STRING-DB…")
    try:
        nodes, edges = ppi.fetch_network(
            params["identifiers"], params.get("species", 9606),
            params.get("required_score", 400),
            params.get("max_partners", 20) if params.get("include_partners") else 0,
            params.get("network_type", "functional"))
    except Exception as exc:
        # 离线降级：仅节点无边
        nodes = [{"id": i} for i in params["identifiers"]]
        edges = []
        degraded = str(exc)[:200]
    else:
        degraded = None
    control.checkpoint(progress=55, phase="analysis",
                       message=f"图分析（{len(nodes)} 节点 {len(edges)} 边）…")
    analysis = ppi.analyze_graph(nodes, edges)
    control.checkpoint(progress=80, phase="export", message="导出 GraphML / Cytoscape…")
    files = ppi.export_graph(nodes, edges, outputs / "network", analysis)
    control.checkpoint(progress=90, phase="enrichment", message="STRING 富集分析…")
    try:
        enrichment = ppi.fetch_enrichment([n["id"] for n in nodes],
                                          params.get("species", 9606)) if edges else []
    except Exception:
        enrichment = []
    result_path = outputs / "ppi_results.json"
    result_path.write_text(json.dumps({
        "nodes": nodes, "edges": edges, "analysis": analysis,
        "enrichment": enrichment, "degraded": degraded,
    }, ensure_ascii=False, indent=1), encoding="utf-8", newline="\n")
    return {"stats": {"n_nodes": len(nodes), "n_edges": len(edges),
                      "n_communities": analysis.get("stats", {}).get("n_communities", 0),
                      "degraded": degraded,
                      "hub_genes": _hubs(analysis)},
            "message": (f"PPI 完成：{len(nodes)} 节点 / {len(edges)} 边"
                        + ("（离线降级模式）" if degraded else "")),
            "artifacts": [{"name": "ppi_results.json", "path": str(result_path), "kind": "json"},
                          *[{"name": f.name, "path": str(f), "kind": "file"} for f in files]]}


def _hubs(analysis: dict) -> list[dict]:
    cent = analysis.get("centrality") or {}
    rows = [{"gene": g, **(v or {})} for g, v in cent.items()]
    rows.sort(key=lambda r: -(r.get("degree") or 0))
    return rows[:10]


@router.get("/ppi/networks")
def ppi_list():
    return {"networks": manager.list(kind="ppi")}


@router.get("/ppi/networks/{job_id}/results")
def ppi_results(job_id: str):
    rec = manager.get(job_id)
    if rec is None or rec["kind"] != "ppi":
        raise HTTPException(404, "PPI 任务不存在")
    f = manager.job_dir(job_id) / "outputs" / "ppi_results.json"
    if not f.is_file():
        raise HTTPException(404, "结果尚未生成")
    return json.loads(f.read_text(encoding="utf-8"))


# ================= GEO 差异分析 =================
@router.get("/geo/availability")
def geo_availability():
    import importlib.util

    ok = all(importlib.util.find_spec(m) for m in ("pandas", "numpy", "scipy", "openpyxl"))
    return {"engine": "geo_diff", "available": ok, "method": "welch t-test + BH-FDR"}


@router.post("/geo/analyze")
def geo_analyze(
    gse: str = Form(...),
    manual_groups: str = Form(""),
    name: str = Form(""),
):
    manual = None
    if manual_groups.strip():
        try:
            manual = json.loads(manual_groups)
        except json.JSONDecodeError:
            raise HTTPException(400, "manual_groups 需为 JSON：{\"control\": [...], \"treatment\": [...]}")
    record = manager.create("geo", name or f"GEO 差异分析: {gse}",
                            {"gse": gse.upper().strip(), "manual_groups": manual})
    return manager.submit(record["id"], geo_pipe.run_geo_job)


@router.get("/geo/jobs/{job_id}")
def geo_job(job_id: str):
    rec = manager.get(job_id)
    if rec is None or rec["kind"] != "geo":
        raise HTTPException(404, "GEO 任务不存在")
    return rec


# ================= 逆合成规划 =================
@router.get("/retro/availability")
def retro_availability():
    try:
        from rdkit.Chem import BRICS  # noqa: F401
        from rdkit.Chem import Recap  # noqa: F401

        return {"engine": "retrosynthesis", "available": True,
                "engines": ["brics", "recap", "hybrid"]}
    except Exception:
        return {"engine": "retrosynthesis", "available": False, "hint": "需要 rdkit"}


@router.post("/retro/plan")
def retro_plan(
    smiles: str = Form(...),
    engine: str = Form("hybrid"),
    max_routes: int = Form(5),
    max_depth: int = Form(2),
):
    try:
        result = retro.plan_routes(smiles.strip(), engine, max_routes, max_depth)
    except ValueError as exc:
        raise HTTPException(400, str(exc))
    return result
