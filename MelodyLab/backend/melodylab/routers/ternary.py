"""三元复合物预测路由（DeepTernary 三元复合物 + BSA 结构健康）。"""
from __future__ import annotations

import json

from fastapi import APIRouter, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, PlainTextResponse

from ..core.jobs import manager
from ..engines.ternary import availability
from ..engines.ternary import config as cfg
from ..engines.ternary import service as ternary_service

router = APIRouter(prefix="/api/ternary", tags=["ternary"])

_DOWNLOAD_MEDIA = {"json": "application/json", "csv": "text/csv", "md": "text/markdown"}


@router.get("/availability")
def ternary_availability():
    return availability()


@router.get("/example")
def ternary_example():
    """内置示例：合成演示分子（TD2 + PEG + 泊马度胺架构，公开来源，非任何实研系列）。"""
    if not cfg.EXAMPLE_FILE.is_file():
        raise HTTPException(404, "示例文件缺失")
    text = cfg.EXAMPLE_FILE.read_text(encoding="utf-8").strip()
    return {"name": "合成演示分子 SYN-01（Galectin-3 × CRBN 架构示例）",
            "n": 1, "text": text}


@router.post("/predict")
def ternary_predict(
    name: str = Form(""),
    smiles_text: str = Form(""),
    default_e3: str = Form(""),
    seeds: int = Form(2),
    smiles_file: UploadFile | None = File(None),
):
    if not availability()["available"]:
        raise HTTPException(400, "三元复合物引擎未就绪：见可用性接口的 hint")
    record = manager.create("ternary", name or "三元复合物预测", {
        "task": "predict", "default_e3": (default_e3 or "").strip().upper(),
        "seeds": max(1, min(int(seeds), 4))})
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
    return manager.submit(record["id"], ternary_service.run_ternary_job)


@router.get("/jobs/{job_id}")
def ternary_job(job_id: str):
    rec = manager.get(job_id)
    if rec is None or rec["kind"] != "ternary":
        raise HTTPException(404, "任务不存在")
    return rec


@router.post("/jobs/{job_id}/resume")
def ternary_job_resume(job_id: str):
    rec = manager.get(job_id)
    if rec is None or rec["kind"] != "ternary":
        raise HTTPException(404, "任务不存在")
    return manager.resume(job_id, ternary_service.run_ternary_job)


@router.get("/reports")
def ternary_reports(limit: int = 100):
    items = []
    for path in sorted(cfg.REPORTS_DIR.glob("*/report.json"),
                       key=lambda p: p.stat().st_mtime, reverse=True):
        try:
            r = json.loads(path.read_text(encoding="utf-8"))
            items.append({
                "id": r.get("id"), "name": r.get("name", ""),
                "n_query": r.get("n_query", 0), "n_ok": r.get("n_ok", 0),
                "seeds": r.get("seeds", 0),
                "top": [row.get("compound_id") for row in (r.get("rows") or [])
                        if row.get("status") == "ok"][:5],
                "created_at": r.get("created_at", ""),
            })
        except Exception:  # noqa: BLE001
            continue
        if len(items) >= limit:
            break
    return items


@router.get("/reports/{report_id}")
def ternary_report(report_id: str):
    path = cfg.REPORTS_DIR / report_id / "report.json"
    if not path.is_file():
        raise HTTPException(404, "报告不存在")
    return json.loads(path.read_text(encoding="utf-8"))


@router.get("/reports/{report_id}/complex/{cid}", response_class=PlainTextResponse)
def ternary_complex(report_id: str, cid: str):
    """最优 seed 的三元复合物 PDB（3Dmol 展示用）。"""
    if "/" in cid or "\\" in cid or ".." in cid:
        raise HTTPException(400, "非法分子 ID")
    path = cfg.REPORTS_DIR / report_id / f"complex_{cid}.pdb"
    if not path.is_file():
        raise HTTPException(404, "复合物文件不存在")
    return PlainTextResponse(path.read_text(encoding="utf-8"),
                             media_type="chemical/x-pdb")


@router.get("/reports/{report_id}/download/{ext}")
def ternary_report_download(report_id: str, ext: str):
    if ext not in _DOWNLOAD_MEDIA:
        raise HTTPException(400, "仅支持 json / csv / md 下载")
    path = cfg.REPORTS_DIR / report_id / f"report.{ext}"
    if not path.is_file():
        raise HTTPException(404, "文件不存在")
    return FileResponse(path, media_type=_DOWNLOAD_MEDIA[ext], filename=f"{report_id}.{ext}")


def _sanitize(name: str) -> str:
    keep = "".join(c for c in name.replace("\\", "/").split("/")[-1]
                   if c.isalnum() or c in "._-")
    return keep or "upload.bin"
