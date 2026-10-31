"""专利 SAR 挖掘路由：/api/patentsar（挖掘长任务 + 同步 SAR 分析/总结）。"""
from __future__ import annotations

import time
from pathlib import Path

from fastapi import APIRouter, File, Form, HTTPException, UploadFile

from ..core.jobs import manager
from ..engines.patentsar import service

router = APIRouter(prefix="/api/patentsar", tags=["patentsar"])


@router.get("/availability")
def patentsar_availability():
    return service.availability()


@router.get("/example")
def patentsar_example():
    return service.example_info()


@router.post("/mine")
def patentsar_mine(
    name: str = Form(""),
    pdf: UploadFile | None = File(None),
    example: bool = Form(False),
    path_type: str = Form("auto"),
    structure_pages: str = Form(""),
    assay_pages: str = Form(""),
    mz_pages: str = Form(""),
    markush_config: str = Form(""),
    markush_pages: str = Form(""),
    target: str = Form("TARGET"),
    kinds: str = Form("example"),
    activity_col: str = Form(""),
    core: str = Form(""),
    hand_drawn: bool = Form(False),
    zoom: float = Form(3.0),
    text_zoom_mult: float = Form(2.0),
    mz_tol: float = Form(1.5),
    vlm: str = Form("auto"),
    with_sar: bool = Form(True),
    with_summary: bool = Form(True),
):
    params = {
        "path_type": path_type, "structure_pages": structure_pages,
        "assay_pages": assay_pages, "mz_pages": mz_pages,
        "markush_config": markush_config, "markush_pages": markush_pages,
        "target": target, "kinds": kinds, "activity_col": activity_col,
        "core": core, "hand_drawn": hand_drawn, "zoom": zoom,
        "text_zoom_mult": text_zoom_mult, "mz_tol": mz_tol, "vlm": vlm,
        "with_sar": with_sar, "with_summary": with_summary, "example": example,
    }
    record = manager.create("patentsar", name or f"专利 SAR 挖掘: {target}", params)
    if not example:
        inputs = manager.job_dir(record["id"]) / "inputs"
        if not (pdf and pdf.filename):
            manager.delete(record["id"])
            raise HTTPException(400, "请上传专利 PDF（或使用内置示例）")
        dest = inputs / _sanitize(pdf.filename)
        with open(dest, "wb") as f:
            while chunk := pdf.file.read(1 << 20):
                f.write(chunk)
    return manager.submit(record["id"], service.run_mine_job)


@router.get("/jobs/{job_id}")
def patentsar_job(job_id: str):
    rec = manager.get(job_id)
    if rec is None or rec["kind"] != "patentsar":
        raise HTTPException(404, "任务不存在")
    return rec


@router.post("/jobs/{job_id}/resume")
def patentsar_job_resume(job_id: str):
    rec = manager.get(job_id)
    if rec is None or rec["kind"] != "patentsar":
        raise HTTPException(404, "任务不存在")
    return manager.resume(job_id, service.run_mine_job)


@router.post("/sar")
def patentsar_sar(
    dataset: UploadFile = File(...),
    activity_col: str = Form(""),
    target: str = Form("TARGET"),
    core: str = Form(""),
):
    """同步 R 基团分析 + Free-Wilson + 3D-QSAR；返回报告文本与产物路径。"""
    import csv

    outdir = service.OUTPUTS_DIR / f"sar_{time.strftime('%Y%m%d_%H%M%S')}"
    updir = outdir / "input"
    updir.mkdir(parents=True, exist_ok=True)
    ds_path = updir / _sanitize(dataset.filename or "dataset.csv")
    with open(ds_path, "wb") as f:
        while chunk := dataset.file.read(1 << 20):
            f.write(chunk)
    if not activity_col:
        with open(ds_path, newline="") as fh:
            header = next(csv.reader(fh), [])
        activity_col = next((c for c in header if c.lower().startswith("ec50")), "")
    if not activity_col:
        raise HTTPException(400, "数据集中未找到 ec50* 活性列，请在参数中指定")
    try:
        res = service.sar_analysis(ds_path, activity_col, target, outdir, core=core or None)
    except Exception as e:  # noqa: BLE001
        raise HTTPException(500, f"SAR 分析失败：{type(e).__name__}: {e}") from e
    report_md = (outdir / "sar_report.md").read_text(encoding="utf-8") if (outdir / "sar_report.md").is_file() else ""
    return {"outdir": str(outdir), "stats": res, "sar_report_md": report_md}


@router.post("/summarize")
def patentsar_summarize(
    dataset: UploadFile = File(...),
    rgroup_table: UploadFile | None = File(None),
    target: str = Form("TARGET"),
    lang: str = Form("zh"),
):
    """同步叙事总结（引用校验门）；返回 markdown + 结构化 JSON。"""
    outdir = service.OUTPUTS_DIR / f"summary_{time.strftime('%Y%m%d_%H%M%S')}"
    updir = outdir / "input"
    updir.mkdir(parents=True, exist_ok=True)
    ds_path = updir / _sanitize(dataset.filename or "dataset.csv")
    with open(ds_path, "wb") as f:
        while chunk := dataset.file.read(1 << 20):
            f.write(chunk)
    rg_path = None
    if rgroup_table and rgroup_table.filename:
        rg_path = updir / _sanitize(rgroup_table.filename)
        with open(rg_path, "wb") as f:
            while chunk := rgroup_table.file.read(1 << 20):
                f.write(chunk)
    try:
        s = service.summarize_dataset(ds_path, target, outdir,
                                      rgroup_table=rg_path, lang=lang)
    except Exception as e:  # noqa: BLE001
        raise HTTPException(500, f"总结失败：{type(e).__name__}: {e}") from e
    return {"outdir": str(outdir), "markdown": s["markdown"],
            "result": {"claims": s["result"].get("claims"),
                       "quarantined": s["result"].get("quarantined"),
                       "suggested_next": s["result"].get("suggested_next")}}


@router.post("/selftest")
def patentsar_selftest():
    return service.selftest()


def _sanitize(name: str) -> str:
    keep = "".join(c for c in name.replace("\\", "/").split("/")[-1]
                   if c.isalnum() or c in "._-")
    return keep or "upload.bin"
