"""批量对接流水线：输入标准化 → 口袋检测（缓存）→ 任务展开 → 引擎执行 → 汇总排名。

每个任务（蛋白 × 配体 × 口袋）写入 work/tasks/<task_id>/state.json，
暂停/中断后 resume 自动跳过已完成任务（job.json 断点恢复模式）。
"""
from __future__ import annotations

import csv
import json
import shutil
import zipfile
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from ...config import settings
from ..screening.rdkit_screening import read_smiles_file
from . import engines as E
from .pockets import detect_pockets

RDLogger.DisableLog("rdApp.error")


# ---------------- 输入标准化 ----------------
def collect_ligands(job_dir: Path, params: dict) -> list[dict]:
    """返回 [{id, smiles?, sdf_path?}]；来源优先级：csv > smiles 文本 > sdf 文件列表。"""
    inputs = job_dir / "inputs"
    ligands: list[dict] = []
    seen_ids: set[str] = set()

    def unique(base: str) -> str:
        i, name = 1, base
        while name in seen_ids:
            i += 1
            name = f"{base}_{i}"
        seen_ids.add(name)
        return name

    csv_path = inputs / "ligands.csv"
    txt_path = inputs / "ligands.txt"
    if csv_path.is_file() or txt_path.is_file():
        source = csv_path if csv_path.is_file() else txt_path
        for smi, mid in read_smiles_file(source):
            ligands.append({"id": unique(mid or f"lig_{len(ligands)+1}"), "smiles": smi})
    max_ligands = params.get("max_ligands")
    if max_ligands:
        ligands = ligands[: int(max_ligands)]
    for sdf in sorted(inputs.glob("*.sdf")):
        ligands.append({"id": unique(sdf.stem), "sdf_path": str(sdf)})
    if not ligands:
        raise ValueError("未提供配体：请上传 CSV/SMILES 文本或 SDF 文件")
    return ligands


def collect_proteins(job_dir: Path) -> list[Path]:
    inputs = job_dir / "inputs"
    exts = (".pdb", ".ent")
    proteins = [p for p in sorted(inputs.iterdir()) if p.suffix.lower() in exts and p.is_file()]
    if not proteins:
        raise ValueError("未提供蛋白结构（PDB/ENT）")
    return proteins


# ---------------- 口袋检测（按蛋白缓存） ----------------
def pockets_for_protein(protein: Path, work_dir: Path, params: dict) -> list[dict]:
    cache = work_dir / "pocket_cache"
    cache.mkdir(parents=True, exist_ok=True)
    out_json = cache / f"{protein.stem}.pockets.json"
    if out_json.is_file() and params.get("resume", False):
        try:
            return json.loads(out_json.read_text(encoding="utf-8"))["pockets"]
        except Exception:
            pass
    custom = (params.get("custom_pockets") or {}).get(protein.stem)
    if custom:
        pockets = [{**custom, "id": "custom", "score": 999.0, "residues": ""}]
    else:
        detected = detect_pockets(protein, max_pockets=int(params.get("max_pockets", 3)))
        pockets = [p.as_json() for p in detected]
    out_json.write_text(json.dumps({"protein": protein.name, "pockets": pockets}, ensure_ascii=False, indent=2),
                        encoding="utf-8")
    return pockets


# ---------------- 任务展开与执行 ----------------
def _task_id(protein_id: str, ligand_id: str, pocket_id: str) -> str:
    return f"{protein_id}__{ligand_id}__{pocket_id}"


def _task_state(task_dir: Path) -> dict | None:
    f = task_dir / "state.json"
    if f.is_file():
        try:
            return json.loads(f.read_text(encoding="utf-8"))
        except Exception:
            return None
    return None


def run_docking_job(record: dict, control) -> dict:
    """JobManager worker。control: JobControl（checkpoint/取消/暂停）。"""
    job_dir = settings.jobs_dir / record["id"]
    params = record["params"]
    work = job_dir / "work"
    outputs = job_dir / "outputs"
    work.mkdir(exist_ok=True)
    outputs.mkdir(exist_ok=True)

    engine_requested = params.get("engine") or "auto"
    engine = engine_requested if engine_requested != "auto" else None
    if engine is None:
        control.checkpoint(progress=2, phase="engine", message="选择对接引擎…")
        engine = E.select_engine(None)
    else:
        E.select_engine(engine)  # 校验可用性
    exhaustiveness = int(params.get("exhaustiveness", 8))

    # 1. 输入
    control.checkpoint(progress=5, phase="preparing", message="标准化输入…")
    ligands = collect_ligands(job_dir, params)
    proteins = collect_proteins(job_dir)

    # 2. 受体准备
    control.checkpoint(progress=8, phase="preparing_proteins", message=f"准备 {len(proteins)} 个受体…")
    receptor_dir = work / "receptors"
    receptor_dir.mkdir(exist_ok=True)
    receptors: dict[str, Path] = {}
    for protein in proteins:
        rec_pdbqt = receptor_dir / f"{protein.stem}.pdbqt"
        if not rec_pdbqt.is_file() or not params.get("resume"):
            E.prepare_receptor(protein, rec_pdbqt)
        receptors[protein.stem] = rec_pdbqt

    # 3. 口袋
    control.checkpoint(progress=15, phase="pockets", message="检测口袋…")
    pockets_map: dict[str, list[dict]] = {}
    for protein in proteins:
        pockets_map[protein.stem] = pockets_for_protein(protein, work, params)

    # 4. 任务展开
    tasks = []
    for protein in proteins:
        for lig in ligands:
            for pocket in pockets_map[protein.stem]:
                tasks.append({"task_id": _task_id(protein.stem, lig["id"], pocket["id"]),
                              "protein": protein.stem, "ligand": lig["id"],
                              "pocket": pocket, "ligand_ref": lig})
    total = len(tasks)
    if total == 0:
        raise ValueError("口袋检测无结果且未提供自定义口袋")

    # 5. 执行
    results, failures = [], []
    done = 0
    for t in tasks:
        task_dir = work / "tasks" / t["task_id"]
        state = _task_state(task_dir)
        if state and state.get("status") == "completed":
            results.append(state["result"])
            done += 1
            continue
        task_dir.mkdir(parents=True, exist_ok=True)
        control.checkpoint(
            progress=20 + int(70 * done / total),
            phase="docking",
            message=f"对接 {t['protein']} + {t['ligand']}（{t['pocket']['id']}） {done+1}/{total}",
            completed=done, total=total,
        )
        try:
            lig_pdbqt = task_dir / "ligand.pdbqt"
            if not lig_pdbqt.is_file():
                E.prepare_ligand_pdbqt(t["ligand_ref"].get("smiles"),
                                       Path(t["ligand_ref"]["sdf_path"]) if t["ligand_ref"].get("sdf_path") else None,
                                       lig_pdbqt)
            out_pdbqt = task_dir / "docked.pdbqt"
            result = E.run_docking(
                engine, receptors[t["protein"]], lig_pdbqt, out_pdbqt,
                task_dir / "engine.log", t["pocket"],
                exhaustiveness=exhaustiveness,
                cancel_event=control.cancel_event,
            )
            pose_sdf = task_dir / "docked.sdf"
            E.pose_to_sdf(out_pdbqt, pose_sdf)
            result.update({
                "protein": t["protein"], "ligand": t["ligand"],
                "pocket": t["pocket"]["id"],
                "pose_sdf": str(pose_sdf) if pose_sdf.is_file() else None,
                "pose_pdbqt": str(out_pdbqt),
                "task_id": t["task_id"],
            })
            (task_dir / "state.json").write_text(
                json.dumps({"status": "completed", "result": result}, ensure_ascii=False), encoding="utf-8")
            results.append(result)
        except Exception as exc:  # noqa: BLE001
            failures.append({"task": t["task_id"], "error": str(exc)[:300]})
            (task_dir / "state.json").write_text(
                json.dumps({"status": "failed", "error": str(exc)[:500]}, ensure_ascii=False), encoding="utf-8")
        finally:
            done += 1

    # 6. 汇总
    control.checkpoint(progress=95, phase="packaging", message="汇总排名与打包…")
    results.sort(key=lambda r: r["score"])
    csv_path = outputs / "docking_results.csv"
    fields = ["task_id", "protein", "ligand", "pocket", "engine", "score"]
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(results)

    top_n = int(params.get("top_n", 10))
    zip_path = outputs / f"{record['id']}_top{min(top_n, len(results))}.zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.write(csv_path, "docking_results.csv")
        for r in results[:top_n]:
            td = work / "tasks" / r["task_id"]
            for f in ("docked.pdbqt", "docked.sdf"):
                p = td / f
                if p.is_file():
                    zf.write(p, f"top{top_n}/{r['task_id']}_{f}")
            log = td / "engine.log"
            if log.is_file():
                zf.write(log, f"top{top_n}/{r['task_id']}_engine.log")

    scores = [r["score"] for r in results]
    stats = {
        "engine": engine,
        "total_tasks": total,
        "completed": len(results),
        "failed": len(failures),
        "best_score": scores[0] if scores else None,
        "avg_score": round(sum(scores) / len(scores), 2) if scores else None,
        "top_zip": str(zip_path),
    }
    return {
        "stats": stats,
        "message": f"完成 {len(results)}/{total}（失败 {len(failures)}），最佳 {stats['best_score']} kcal/mol",
        "artifacts": [
            {"name": "docking_results.csv", "path": str(csv_path), "kind": "csv"},
            {"name": f"TOP{top_n} 结果包", "path": str(zip_path), "kind": "zip"},
        ],
    }
