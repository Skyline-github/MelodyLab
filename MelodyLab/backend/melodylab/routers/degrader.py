"""降解剂评估路由（PROTAC DC50 共识排序）。"""
from __future__ import annotations

import json

from fastapi import APIRouter, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse

from ..core.jobs import manager
from ..engines.degrader import availability, model_info
from ..engines.degrader import config as cfg
from ..engines.degrader import service as degrader_service

router = APIRouter(prefix="/api/degrader", tags=["degrader"])

_DOWNLOAD_MEDIA = {"json": "application/json", "csv": "text/csv", "md": "text/markdown"}


@router.get("/availability")
def degrader_availability():
    return availability()


@router.get("/example")
def degrader_example():
    """内置示例：PROTAC-DB 已发表 PROTAC（公开数据，含实测 DC50 可对照）。"""
    if not cfg.EXAMPLE_FILE.is_file():
        raise HTTPException(404, "示例文件缺失")
    from ..engines.degrader.features import parse_query_csv

    items = parse_query_csv(cfg.EXAMPLE_FILE)
    text = "\n".join(f"{it['compound_id']} {it['smiles']}" for it in items)
    return {"name": "PROTAC-DB 公开分子示例", "n": len(items), "text": text}


@router.post("/rank")
def degrader_rank(
    name: str = Form(""),
    smiles_text: str = Form(""),
    smiles_file: UploadFile | None = File(None),
):
    record = manager.create("degrader", name or "降解剂 DC50 排序", {"task": "rank"})
    inputs = manager.job_dir(record["id"]) / "inputs"
    if smiles_file and smiles_file.filename:
        dest = inputs / _sanitize(smiles_file.filename)
        with open(dest, "wb") as f:
            while chunk := smiles_file.file.read(1 << 20):
                f.write(chunk)
    if smiles_text.strip():
        (inputs / "queries.txt").write_text(smiles_text, encoding="utf-8")
    if not any(inputs.iterdir()):
        manager.delete(record["id"])
        raise HTTPException(400, "请粘贴 SMILES 或上传查询文件")
    return manager.submit(record["id"], degrader_service.run_degrader_job)


@router.post("/retrain")
def degrader_retrain():
    record = manager.create("degrader", "重训 DC50 排序模型", {"task": "retrain"})
    return manager.submit(record["id"], degrader_service.run_degrader_job)


@router.get("/jobs/{job_id}")
def degrader_job(job_id: str):
    rec = manager.get(job_id)
    if rec is None or rec["kind"] != "degrader":
        raise HTTPException(404, "任务不存在")
    return rec


@router.post("/jobs/{job_id}/resume")
def degrader_job_resume(job_id: str):
    rec = manager.get(job_id)
    if rec is None or rec["kind"] != "degrader":
        raise HTTPException(404, "任务不存在")
    return manager.resume(job_id, degrader_service.run_degrader_job)


@router.get("/reports")
def degrader_reports(limit: int = 100):
    """排序归档列表（数据根 degrader_reports/ 下的轻量索引）。"""
    items = []
    for path in sorted(cfg.REPORTS_DIR.glob("*.json"), key=lambda p: p.stat().st_mtime,
                       reverse=True):
        try:
            r = json.loads(path.read_text(encoding="utf-8"))
            items.append({
                "id": r.get("id"), "name": r.get("name", ""),
                "n_query": r.get("n_query", 0), "n_invalid": r.get("n_invalid", 0),
                "n_train": (r.get("calibration") or {}).get("n_train"),
                "model_status": r.get("model_status", ""),
                "top5": [row.get("compound_id") for row in (r.get("rows") or [])[:5]],
                "created_at": r.get("created_at", ""),
            })
        except Exception:  # noqa: BLE001
            continue
        if len(items) >= limit:
            break
    return items


@router.get("/reports/{report_id}")
def degrader_report(report_id: str):
    path = cfg.REPORTS_DIR / f"{report_id}.json"
    if not path.is_file():
        raise HTTPException(404, "报告不存在")
    return json.loads(path.read_text(encoding="utf-8"))


@router.get("/reports/{report_id}/download/{ext}")
def degrader_report_download(report_id: str, ext: str):
    if ext not in _DOWNLOAD_MEDIA:
        raise HTTPException(400, "仅支持 json / csv / md 下载")
    path = cfg.REPORTS_DIR / f"{report_id}.{ext}"
    if not path.is_file():
        raise HTTPException(404, "文件不存在")
    return FileResponse(path, media_type=_DOWNLOAD_MEDIA[ext], filename=path.name)


def _sanitize(name: str) -> str:
    keep = "".join(c for c in name.replace("\\", "/").split("/")[-1]
                   if c.isalnum() or c in "._-")
    return keep or "upload.bin"
