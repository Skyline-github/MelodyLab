"""业务任务路由：筛选 / 对接 / 分子生成 / Linker（创建 + 结果查询）。"""
from __future__ import annotations

import json
from pathlib import Path

from fastapi import APIRouter, File, Form, HTTPException, UploadFile

from ..config import settings
from ..core.jobs import JobControl, manager
from ..engines.docking import pipeline as dp
from ..engines.docking.pockets import detect_pockets
from ..engines.generation import flashdiff_adapter
from ..engines.linker import diffprotacs_adapter

router = APIRouter(prefix="/api", tags=["workloads"])


def _save_upload(file: UploadFile, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "wb") as f:
        while chunk := file.file.read(1 << 20):
            f.write(chunk)


def _sanitize(name: str) -> str:
    keep = "".join(c for c in name.replace("\\", "/").split("/")[-1] if c.isalnum() or c in "._-")
    return keep or "upload.bin"


# ================= 分子准入筛选 =================
@router.post("/screening")
def create_screening(
    name: str = Form(""),
    smiles_text: str = Form(""),
    smiles_file: UploadFile | None = File(None),
    qed_min: float = Form(0.4),
    max_ro5_violations: int = Form(1),
    sa_max: float = Form(6.0),
    pains_mode: str = Form("alert"),
    brenk_mode: str = Form("alert"),
    selection_limit: int | None = Form(None),
):
    record = manager.create("screening", name or "分子准入筛选", {
        "qed_min": qed_min, "max_ro5_violations": max_ro5_violations,
        "sa_max": sa_max, "pains_mode": pains_mode, "brenk_mode": brenk_mode,
        "selection_limit": selection_limit,
    })
    inputs = manager.job_dir(record["id"]) / "inputs"
    if smiles_file and smiles_file.filename:
        _save_upload(smiles_file, inputs / smiles_file.filename)
    if smiles_text.strip():
        (inputs / "ligands.txt").write_text(smiles_text, encoding="utf-8")
    if not any(inputs.iterdir()):
        manager.delete(record["id"])
        raise HTTPException(400, "请提供 SMILES 文本或上传文件")

    def worker(rec: dict, control: JobControl) -> dict:
        from ..engines.screening.rdkit_screening import ScreeningPolicy, run_screening

        job_dir = settings.jobs_dir / rec["id"]
        policy = ScreeningPolicy(
            qed_min=qed_min, max_ro5_violations=max_ro5_violations, sa_max=sa_max,
            pains_mode=pains_mode, brenk_mode=brenk_mode, selection_limit=selection_limit,
        )
        from ..engines.screening.rdkit_screening import read_smiles_file

        rows_in = []
        for p in sorted(inputs.iterdir()):
            rows_in.extend(read_smiles_file(p))
        if not rows_in:
            raise ValueError("未解析到任何 SMILES")
        control.checkpoint(progress=5, phase="screening", message=f"开始筛选 {len(rows_in)} 个分子")

        def progress_cb(pct, msg):
            control.heartbeat = control.heartbeat  # 保持引用
            manager.update(rec["id"], progress=pct, message=msg)

        result = run_screening(rows_in, job_dir / "outputs", policy,
                               progress=progress_cb, cancel_event=control.cancel_event)
        return {
            "stats": result["stats"],
            "message": f"入选 {result['stats']['accepted']}/{result['stats']['total']}",
            "artifacts": [
                {"name": "screening_results.csv", "path": result["csv"], "kind": "csv"},
                {"name": "docking_set.csv", "path": result["hit_csv"], "kind": "csv"},
            ],
        }

    return manager.submit(record["id"], worker)


# ================= 批量对接 =================
@router.post("/docking")
def create_docking(
    name: str = Form(""),
    ligand_smiles: str = Form(""),
    ligand_csv: UploadFile | None = File(None),
    ligand_files: list[UploadFile] = File(default=[]),
    protein_files: list[UploadFile] = File(default=[]),
    engine: str = Form("auto"),
    max_pockets: int = Form(3),
    exhaustiveness: int = Form(8),
    top_n: int = Form(10),
    seed: int = Form(42),
):
    if not protein_files or not any(f.filename for f in protein_files):
        raise HTTPException(400, "请上传至少一个蛋白结构（PDB/ENT）")
    record = manager.create("docking", name or "批量对接", {
        "engine": engine, "max_pockets": max_pockets,
        "exhaustiveness": exhaustiveness, "top_n": top_n, "seed": seed,
    })
    job_dir = manager.job_dir(record["id"])
    inputs = job_dir / "inputs"
    for f in protein_files:
        if f.filename:
            _save_upload(f, inputs / _sanitize(f.filename))
    if ligand_csv and ligand_csv.filename:
        dest = inputs / "ligands.csv" if ligand_csv.filename.lower().endswith(".csv") else inputs / _sanitize(ligand_csv.filename)
        _save_upload(ligand_csv, dest)
    if ligand_smiles.strip():
        (inputs / "ligands.txt").write_text(ligand_smiles, encoding="utf-8")
    for f in ligand_files or []:
        if f.filename and f.filename.lower().endswith(".sdf"):
            _save_upload(f, inputs / _sanitize(f.filename))
    if not (inputs / "ligands.csv").is_file() and not (inputs / "ligands.txt").is_file() \
            and not list(inputs.glob("*.sdf")):
        manager.delete(record["id"])
        raise HTTPException(400, "请提供配体：SMILES 文本 / CSV / SDF 文件")
    return manager.submit(record["id"], dp.run_docking_job)


@router.get("/docking/runs")
def list_runs():
    return {"runs": manager.list(kind="docking")}


@router.get("/docking/runs/{job_id}/results")
def run_results(job_id: str, limit: int = 100, sort: str = "asc"):
    record = manager.get(job_id)
    if record is None or record["kind"] != "docking":
        raise HTTPException(404, "运行不存在")
    work = manager.job_dir(job_id) / "work" / "tasks"
    rows = []
    for state_file in work.glob("*/state.json"):
        try:
            st = json.loads(state_file.read_text(encoding="utf-8"))
        except Exception:
            continue
        if st.get("status") == "completed":
            rows.append(st["result"])
        else:
            rows.append({"task_id": state_file.parent.name, "status": st.get("status"),
                         "error": st.get("error"), "score": None})
    rows.sort(key=lambda r: (r.get("score") is None, r.get("score") if r.get("score") is not None else 0),
              reverse=(sort == "desc"))
    return {"run": record, "results": rows[:limit], "total": len(rows)}


@router.post("/docking/runs/{job_id}/resume")
def resume_run(job_id: str):
    record = manager.get(job_id)
    if record is None:
        raise HTTPException(404, "运行不存在")
    record["params"]["resume"] = True
    return manager.resume(job_id, dp.run_docking_job) or HTTPException(400, "只有已暂停的任务可恢复")


# ================= 口袋预测（独立接口） =================
@router.post("/pockets")
async def predict_pockets(protein: UploadFile = File(...), max_pockets: int = Form(5)):
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / _sanitize(protein.filename or "protein.pdb")
        _save_upload(protein, p)
        try:
            pockets = detect_pockets(p, max_pockets=max_pockets)
        except Exception as exc:
            raise HTTPException(400, f"口袋检测失败: {exc}")
        return {"protein": p.name, "pockets": [pk.as_json() for pk in pockets]}


# ================= 3D 分子生成（FlashDiff） =================
@router.get("/generation/availability")
def generation_availability():
    return flashdiff_adapter.availability()


@router.post("/generation")
def create_generation(
    protein: UploadFile = File(...),
    ligand: UploadFile = File(...),
    num_molecules: int = Form(5),
    batch_size: int = Form(10),
    pocket_radius: float = Form(15.0),
    num_steps: int = Form(50),
    seed: int = Form(2024),
    mean_atoms: int = Form(28),
    std_atoms: int = Form(5),
):
    if not flashdiff_adapter.availability()["available"]:
        raise HTTPException(400, "FlashDiff 资产未就绪：请先运行 scripts/setup_assets.py")
    record = manager.create("generation", "3D 分子生成", {
        "num_molecules": num_molecules, "batch_size": batch_size,
        "pocket_radius": pocket_radius, "num_steps": num_steps, "seed": seed,
        "mean_atoms": mean_atoms, "std_atoms": std_atoms,
    })
    job_dir = manager.job_dir(record["id"])
    _save_upload(protein, job_dir / "inputs" / "protein.pdb")
    _save_upload(ligand, job_dir / "inputs" / "ligand.sdf")

    def worker(rec: dict, control: JobControl) -> dict:
        jd = settings.jobs_dir / rec["id"]
        p = rec["params"]
        control.checkpoint(progress=5, phase="generating", message="FlashDiff 生成中（CPU 较慢，请耐心等待）…")
        out = flashdiff_adapter.generate(
            jd / "inputs" / "protein.pdb", jd / "inputs" / "ligand.sdf",
            jd / "work", num_molecules=p["num_molecules"], batch_size=p["batch_size"],
            pocket_radius=p["pocket_radius"], num_steps=p["num_steps"], seed=p["seed"],
            mean_atoms=p["mean_atoms"], std_atoms=p["std_atoms"],
        )
        import shutil

        shutil.copy2(jd / "inputs" / "protein.pdb", jd / "outputs" / "protein.pdb")
        arts = [{"name": "protein.pdb", "path": str(jd / "outputs" / "protein.pdb"), "kind": "pdb"}]
        for m in out["molecules"]:
            arts.append({"name": f"mol_{m['id']}.sdf", "path": m["sdf"], "kind": "sdf"})
            if m["trajectory"]:
                arts.append({"name": f"mol_{m['id']}_trajectory.sdf", "path": m["trajectory"], "kind": "sdf"})
        return {"stats": {"num_generated": len(out["molecules"]),
                          "molecules": out["molecules"], "pocket_center": out["pocket_center"]},
                "message": f"生成 {len(out['molecules'])} 个有效分子", "artifacts": arts}

    return manager.submit(record["id"], worker)


# ================= Linker 生成（DiffPROTACs） =================
@router.get("/linker/availability")
def linker_availability():
    return diffprotacs_adapter.availability()


@router.post("/linker")
def create_linker(
    ligand1: UploadFile = File(...),
    ligand2: UploadFile = File(...),
    linker_len: int = Form(12),
    n_samples: int = Form(5),
):
    if not diffprotacs_adapter.availability()["available"]:
        raise HTTPException(400, "DiffPROTACs 资产未就绪：请先运行 scripts/setup_assets.py")
    record = manager.create("linker", "Linker 生成", {"linker_len": linker_len, "n_samples": n_samples})
    job_dir = manager.job_dir(record["id"])
    _save_upload(ligand1, job_dir / "inputs" / "e3_ligand.sdf")
    _save_upload(ligand2, job_dir / "inputs" / "warhead.sdf")

    def worker(rec: dict, control: JobControl) -> dict:
        jd = settings.jobs_dir / rec["id"]
        control.checkpoint(progress=5, phase="generating", message="DiffPROTACs 生成中…")
        out = diffprotacs_adapter.generate_linkers(
            jd / "inputs" / "e3_ligand.sdf", jd / "inputs" / "warhead.sdf",
            jd / "work", linker_len=rec["params"]["linker_len"],
            n_samples=rec["params"]["n_samples"])
        arts = []
        for lk in out["linkers"]:
            for key in ("full_sdf", "linker_sdf"):
                if Path(lk[key]).is_file():
                    arts.append({"name": Path(lk[key]).name, "path": lk[key], "kind": "sdf"})
        return {"stats": {"num_generated": len(out["linkers"]), "linkers": out["linkers"]},
                "message": f"生成 {len(out['linkers'])} 个 linker", "artifacts": arts}

    return manager.submit(record["id"], worker)
