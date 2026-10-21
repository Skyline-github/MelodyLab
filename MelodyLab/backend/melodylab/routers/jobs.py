"""统一任务中心路由：/api/jobs（跨任务类型的列表/详情/控制/下载）。"""
from __future__ import annotations

import json
from pathlib import Path

from fastapi import APIRouter, HTTPException
from fastapi.responses import FileResponse

from ..core.jobs import manager

router = APIRouter(prefix="/api/jobs", tags=["jobs"])


@router.get("")
def list_jobs(kind: str | None = None, limit: int = 100):
    return {"jobs": manager.list(kind=kind, limit=limit)}


@router.get("/{job_id}")
def get_job(job_id: str):
    record = manager.get(job_id)
    if record is None:
        raise HTTPException(404, "任务不存在")
    return record


@router.delete("/{job_id}")
def delete_job(job_id: str):
    if not manager.delete(job_id):
        raise HTTPException(404, "任务不存在")
    return {"ok": True}


@router.post("/{job_id}/pause")
def pause_job(job_id: str):
    return manager.pause(job_id) or _404()


@router.post("/{job_id}/stop")
def stop_job(job_id: str):
    return manager.stop(job_id) or _404()


@router.get("/{job_id}/artifacts/{name}")
def download_artifact(job_id: str, name: str):
    record = manager.get(job_id)
    if record is None:
        raise HTTPException(404, "任务不存在")
    for a in record["artifacts"]:
        if a["name"] == name and Path(a["path"]).is_file():
            return FileResponse(a["path"], filename=Path(a["path"]).name)
    raise HTTPException(404, "产物不存在")


@router.get("/{job_id}/log")
def job_log(job_id: str, tail: int = 200):
    log = manager.job_dir(job_id) / "logs" / "run.log"
    if not log.is_file():
        return {"lines": []}
    lines = log.read_text(encoding="utf-8", errors="replace").splitlines()
    return {"lines": lines[-tail:]}


def _404():
    raise HTTPException(404, "任务不存在")
