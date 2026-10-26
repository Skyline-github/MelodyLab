"""MD 路由：任务创建 / 查询 / 结果 / 轨迹预览。"""
from __future__ import annotations

import threading
from pathlib import Path

from fastapi import APIRouter, File, Form, HTTPException, UploadFile

from ..config import settings
from ..core.jobs import manager
from ..engines.md import availability
from ..engines.md import pipeline as mdp

router = APIRouter(prefix="/api/md", tags=["md"])


def _save_upload(file: UploadFile, dest) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "wb") as f:
        while chunk := file.file.read(1 << 20):
            f.write(chunk)


@router.get("/availability")
def md_availability():
    return availability()


MD_PARAMS = {
    "chains": None, "total_ns": 5.0, "test_ns": 1.0, "npt_ns": 1.0,
    "forcefield": "amber99sb", "water": "tip3p", "ligand_charge": None,
    "salt_molar": 0.15, "rmsd_threshold": 0.3, "seed": 42, "use_gpu": True,
    "ref_t": 300.0, "run_id": None, "task_id": None,
}


@router.post("/jobs")
def create_md_job(
    run_id: str = Form(""),
    task_id: str = Form(""),
    chains: str = Form(""),
    protein_file: UploadFile | None = File(None),
    ligand_file: UploadFile | None = File(None),
    ligand_smiles: str = Form(""),
    total_ns: float = Form(5.0),
    test_ns: float = Form(1.0),
    npt_ns: float = Form(1.0),
    forcefield: str = Form("amber99sb"),
    water: str = Form("tip3p"),
    ligand_charge: int | None = Form(None),
    salt_molar: float = Form(0.15),
    rmsd_threshold: float = Form(0.3),
    seed: int = Form(42),
    use_gpu: bool = Form(True),
    ref_t: float = Form(300.0),
    ternary_report_id: str = Form(""),
    ternary_cid: str = Form(""),
):
    if not availability()["available"]:
        raise HTTPException(400, "GROMACS 不可用：请先运行 scripts/setup_assets.py --with-gromacs")
    is_ternary = bool(ternary_report_id.strip() and ternary_cid.strip())
    params = {
        "run_id": run_id or None, "task_id": task_id or None,
        "chains": chains or None,
        "total_ns": total_ns, "test_ns": test_ns, "npt_ns": npt_ns,
        "forcefield": forcefield, "water": water,
        "ligand_charge": ligand_charge, "salt_molar": salt_molar,
        "rmsd_threshold": rmsd_threshold, "seed": seed, "use_gpu": use_gpu,
        "ref_t": ref_t,
        "ternary_report_id": ternary_report_id or None,
        "ternary_cid": ternary_cid or None,
    }
    if not (run_id and task_id) and protein_file is None and not is_ternary:
        raise HTTPException(400, "请引用对接任务（run_id+task_id）、直传 protein_file，或引用三元复合物报告")
    record = manager.create(
        "md", f"三元复合物 MD: {ternary_cid.strip()}" if is_ternary else "分子动力学验证", params)
    job_dir = manager.job_dir(record["id"])
    if protein_file is not None and protein_file.filename:
        _save_upload(protein_file, job_dir / "inputs" / protein_file.filename)
    if ligand_file is not None and ligand_file.filename:
        _save_upload(ligand_file, job_dir / "inputs" / ligand_file.filename)
    if ligand_smiles.strip() and not is_ternary:
        (job_dir / "inputs" / "ligands.txt").write_text(ligand_smiles, encoding="utf-8")
    return manager.submit(record["id"], mdp.run_md_job)


@router.get("/jobs")
def list_md_jobs():
    return {"jobs": manager.list(kind="md")}


@router.get("/jobs/{job_id}")
def md_job(job_id: str):
    rec = manager.get(job_id)
    if rec is None or rec["kind"] != "md":
        raise HTTPException(404, "MD 任务不存在")
    return rec


@router.post("/jobs/{job_id}/resume")
def md_resume(job_id: str):
    rec = manager.get(job_id)
    if rec is None or rec["kind"] != "md":
        raise HTTPException(404, "MD 任务不存在")
    if rec["status"] != "paused":
        raise HTTPException(400, "只有已暂停的任务可恢复")
    rec["params"]["resume"] = True
    return manager.resume(job_id, mdp.run_md_job)


@router.get("/jobs/{job_id}/results")
def md_results(job_id: str):
    rec = manager.get(job_id)
    if rec is None or rec["kind"] != "md":
        raise HTTPException(404, "MD 任务不存在")
    out = manager.job_dir(job_id) / "outputs"
    report = None
    f = out / "md_report.json"
    if f.is_file():
        import json

        report = json.loads(f.read_text(encoding="utf-8"))
    return {"job": rec, "report": report}


@router.get("/jobs/{job_id}/trajectory-preview")
def md_trajectory(job_id: str, download: bool = False):
    from fastapi.responses import FileResponse

    p = _ensure_preview(job_id)
    return FileResponse(p, media_type="chemical/x-pdb",
                        filename="trajectory_preview.pdb" if download else None)


_preview_locks: dict[str, threading.Lock] = {}
_preview_locks_guard = threading.Lock()


def _ensure_preview(job_id: str) -> Path:
    """轨迹预览自愈：文件缺失或早于 md.xtc（旧版代码生成）时按当前逻辑重新生成。

    历史任务坏格式的预览在用户点开的那一刻自动修复，保证动画展示
    永远走新版导出逻辑，不依赖任务完成时用的是哪版代码。
    """
    rec = manager.get(job_id)
    if rec is None or rec["kind"] != "md":
        raise HTTPException(404, "MD 任务不存在")
    job_dir = manager.job_dir(job_id)
    xtc = job_dir / "work" / "md.xtc"
    p = job_dir / "outputs" / "trajectory_preview.pdb"
    if not xtc.is_file():
        raise HTTPException(404, "该任务没有生产轨迹（md.xtc），无预览可生成")

    def _stale() -> bool:
        return not p.is_file() or p.stat().st_mtime < xtc.stat().st_mtime

    if _stale():
        with _preview_locks_guard:
            lock = _preview_locks.setdefault(job_id, threading.Lock())
        with lock:
            if _stale():  # 双重检查：等锁期间可能已被并发请求修复
                p.parent.mkdir(parents=True, exist_ok=True)
                mdp._trajectory_frames(job_dir / "work", p, n_frames=mdp.PREVIEW_FRAMES)
    if not p.is_file():
        raise HTTPException(500, "轨迹预览生成失败，请查看后端日志")
    return p
