"""靶点立项评估路由（立项五问，移植自 target-triage）。"""
from __future__ import annotations

import json

from fastapi import APIRouter, Form, HTTPException
from fastapi.responses import FileResponse

from ..core.jobs import manager
from ..engines.triage import availability
from ..engines.triage import service as triage_service
from ..engines.triage.config import DATA_SOURCES, REPORTS_DIR

router = APIRouter(prefix="/api/triage", tags=["triage"])

_DOWNLOAD_MEDIA = {
    "json": "application/json",
    "md": "text/markdown",
    "xlsx": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
}


@router.get("/availability")
def triage_availability():
    info = availability()
    uses = {
        "Open Targets Platform": "靶点-疾病关联、安全事件、遗传约束、小鼠敲除（主干）",
        "GWAS Catalog": "人类遗传学信号补充",
        "Human Protein Atlas": "44 种正常组织蛋白表达 + RNA 表达",
        "UniProt": "蛋白功能与亚细胞定位注释",
        "STRING": "蛋白-蛋白互作网络",
        "ChEMBL": "靶点药物、作用机制与研发阶段",
        "ClinicalTrials.gov": "竞争格局：试验数量、状态、终止原因",
        "openFDA": "药物标签中靶点/突变提及（上市产品线索）",
        "OncoKB（可选）": "伴随诊断与精准肿瘤学证据（需配置 ONCOKB_TOKEN）",
    }
    info["sources"] = [{**s, "use": uses.get(s["name"], "")} for s in DATA_SOURCES]
    return info


@router.post("/assess")
def triage_assess(gene: str = Form(...), indication: str = Form("")):
    gene = gene.strip().upper()
    if not gene:
        raise HTTPException(400, "请输入靶点基因符号")
    record = manager.create("triage", f"靶点评估: {gene}", {
        "gene": gene, "indication": (indication or "").strip()})
    return manager.submit(record["id"], triage_service.run_triage_job)


@router.get("/jobs")
def triage_jobs(limit: int = 50):
    return manager.list("triage", limit=limit)


@router.get("/jobs/{job_id}")
def triage_job(job_id: str):
    rec = manager.get(job_id)
    if rec is None or rec["kind"] != "triage":
        raise HTTPException(404, "评估任务不存在")
    return rec


@router.post("/jobs/{job_id}/resume")
def triage_job_resume(job_id: str):
    rec = manager.get(job_id)
    if rec is None or rec["kind"] != "triage":
        raise HTTPException(404, "评估任务不存在")
    return manager.resume(job_id, triage_service.run_triage_job)


@router.get("/reports")
def triage_reports(limit: int = 100):
    """评估归档列表（数据根 triage_reports/ 下的轻量索引）。"""
    items = []
    for path in sorted(REPORTS_DIR.glob("*.json"), key=lambda p: p.stat().st_mtime,
                       reverse=True):
        try:
            r = json.loads(path.read_text(encoding="utf-8"))
            items.append({"id": r.get("id"), "gene_symbol": r.get("gene_symbol"),
                          "gene_name": r.get("gene_name", ""),
                          "indication": r.get("indication", ""),
                          "overall_status": r.get("overall_status", "unknown"),
                          "risk_count": len(r.get("risk_flags") or []),
                          "watch_count": len(r.get("watch_flags") or []),
                          "created_at": r.get("created_at", "")})
        except Exception:  # noqa: BLE001
            continue
        if len(items) >= limit:
            break
    return items


@router.get("/reports/{report_id}")
def triage_report(report_id: str):
    path = REPORTS_DIR / f"{report_id}.json"
    if not path.is_file():
        raise HTTPException(404, "报告不存在")
    return json.loads(path.read_text(encoding="utf-8"))


@router.get("/reports/{report_id}/download/{ext}")
def triage_report_download(report_id: str, ext: str):
    if ext not in _DOWNLOAD_MEDIA:
        raise HTTPException(400, "仅支持 json / md / xlsx 下载")
    path = REPORTS_DIR / f"{report_id}.{ext}"
    if not path.is_file():
        raise HTTPException(404, "文件不存在")
    return FileResponse(path, media_type=_DOWNLOAD_MEDIA[ext], filename=path.name)
