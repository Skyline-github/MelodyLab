"""三元复合物预测 worker：输入制备（主环境）→ DeepTernary 子进程推理 → BSA → 报告。

JobManager worker 约定：(record, control) -> dict。
制备在主环境（rdkit）；推理以子进程调用 deepternary 独立环境（dgl 与主环境
torch 不兼容，见侦察结论）；BSA 用 Bio.PDB 在主环境计算。
"""
from __future__ import annotations

import json
import shutil
import subprocess
import threading
import time
import uuid
from pathlib import Path

import pandas as pd

from ...config import settings as app_settings
from . import config as cfg
from .models import PrepInfo, TernaryReport, TernaryRow
from .references import build_inputs, get_pack


# ---------------------------------------------------------------- 查询解析
def parse_queries_text(text: str, default_e3: str = "") -> list[dict]:
    """每行 'ID SMILES [E3]'（E3 缺省用表单默认值），或纯 SMILES。"""
    items: list[dict] = []
    for lineno, line in enumerate(text.splitlines(), start=1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2 and len(parts[0]) <= 12 and not any(
                ch in parts[0] for ch in "()=[]#@/\\"):
            cid, smi = parts[0], " ".join(parts[1:])
            e3 = parts[-1].upper() if len(parts) >= 3 and parts[-1].upper() in cfg.E3_REF else ""
        else:
            cid, smi, e3 = f"mol_{lineno:03d}", line, ""
        items.append({"compound_id": cid, "smiles": smi,
                      "e3": e3 or default_e3.upper(), "source_line": lineno})
    return items


def parse_queries_csv(path, default_e3: str = "") -> list[dict]:
    df = pd.read_csv(path)
    cols = {c.lower().strip(): c for c in df.columns}
    smi_col = next((cols[k] for k in ("smiles", "smi", "canonical_smiles") if k in cols), None)
    if smi_col is None:
        raise ValueError("CSV 未找到 smiles 列")
    id_col = next((cols[k] for k in ("compound_id", "compound id", "id", "name") if k in cols), None)
    e3_col = next((cols[k] for k in ("e3 ligase", "e3_ligase", "e3") if k in cols), None)
    items: list[dict] = []
    for i, row in df.iterrows():
        smi = str(row[smi_col]).strip()
        if not smi or smi.lower() == "nan":
            continue

        def _clean(val, default):
            s = str(val).strip()
            return s if s and s.lower() not in ("nan", "none") else default

        e3 = _clean(row[e3_col], "").upper() if e3_col else ""
        if e3 not in cfg.E3_REF:
            e3 = ""
        items.append({"compound_id": _clean(row[id_col], f"mol_{i + 1:03d}") if id_col else f"mol_{i + 1:03d}",
                      "smiles": smi, "e3": e3 or default_e3.upper(), "source_line": i + 2})
    return items


# ---------------------------------------------------------------- BSA（原样移植 run_dt.py）
def compute_bsa(complex_path: str) -> float:
    from copy import deepcopy

    from Bio.PDB import PDBParser, ShrakeRupley

    parser = PDBParser(QUIET=True)
    sasa = ShrakeRupley()
    cx = parser.get_structure("c", complex_path)
    sasa.compute(cx, level="C")
    asa_cx = sum(cx[0][c].sasa for c in cfg.COMPLEX_CHAINS)
    tot = 0.0
    for keep in cfg.COMPLEX_CHAINS:
        st = deepcopy(cx)
        for c in list(st[0]):
            if c.id != keep:
                st[0].detach_child(c.id)
        sasa.compute(st, level="C")
        tot += st[0][keep].sasa
    return float(tot - asa_cx)


# ---------------------------------------------------------------- 子进程推理
def _run_inference(tasks: list[dict], work_dir: Path, control) -> list[dict]:
    interpreter = cfg.interpreter_path()
    if interpreter is None:
        raise RuntimeError(
            "未找到 deepternary 独立环境 python。解决：创建 conda env deepternary"
            "（torch==2.2.1 + dgl==2.2.1 + torchdata==0.7.1，见 docs），"
            "或把其 python.exe 路径写入 数据根/deepternary_python.txt，"
            "或设置环境变量 MELODYLAB_DEEPTERNARY_PYTHON")
    task_json = work_dir / "dt_task.json"
    result_json = work_dir / "dt_result.json"
    log_path = work_dir / "dt_run.log"
    task_json.write_text(json.dumps(
        {"repo": str(cfg.VENDOR_DIR),
         "tasks": [{"name": t["name"], "input_dir": str(t["input_dir"]),
                    "out_dir": str(t["out_dir"]), "seeds": t["seeds"]} for t in tasks],
         "result_path": str(result_json)}, ensure_ascii=False), encoding="utf-8")

    dt_script = Path(__file__).resolve().parent / "dt_task.py"
    cmd = [str(interpreter), str(dt_script), str(task_json)]
    box: dict = {}

    def _run():
        try:
            with open(log_path, "w", encoding="utf-8") as log:
                proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                                      timeout=3600, cwd=str(cfg.VENDOR_DIR))
            box["code"] = proc.returncode
        except Exception as e:  # noqa: BLE001
            box["error"] = f"{type(e).__name__}: {e}"

    th = threading.Thread(target=_run, name="ternary-dt", daemon=True)
    th.start()
    while th.is_alive():
        th.join(timeout=2.0)
        if th.is_alive():
            control.checkpoint(message=f"DeepTernary 推理中（{len(tasks)} 个分子，"
                                       f"模型已加载）…")
    if "error" in box:
        raise RuntimeError(f"DeepTernary 子进程异常：{box['error']}（日志 {log_path}）")
    if box.get("code") != 0:
        raise RuntimeError(f"DeepTernary 子进程退出码 {box['code']}（日志 {log_path}）")
    return json.loads(result_json.read_text(encoding="utf-8"))


# ---------------------------------------------------------------- 报告导出
def save_report(report: TernaryReport, report_dir: Path) -> list[dict]:
    report_dir.mkdir(parents=True, exist_ok=True)
    (report_dir / "report.json").write_text(report.model_dump_json(indent=2),
                                            encoding="utf-8")
    (report_dir / "report.csv").write_text(render_csv(report), encoding="utf-8")
    (report_dir / "report.md").write_text(render_markdown(report), encoding="utf-8")
    base = report.name or "三元复合物预测"
    return [
        {"name": f"{base}.json", "path": str(report_dir / "report.json"), "kind": "json"},
        {"name": f"{base}.csv", "path": str(report_dir / "report.csv"), "kind": "csv"},
        {"name": f"{base}.md", "path": str(report_dir / "report.md"), "kind": "markdown"},
    ]


def render_csv(report: TernaryReport) -> str:
    rows = []
    for r in report.rows:
        d = r.model_dump()
        d["prep"] = json.dumps(d["prep"], ensure_ascii=False)
        rows.append(d)
    return pd.DataFrame(rows).to_csv(index=False)


def render_markdown(report: TernaryReport) -> str:
    md = [
        f"# 三元复合物预测：{report.name or report.id}",
        "",
        f"- **时间**：{report.created_at} · seeds={report.seeds} · "
        f"成功 {report.n_ok}/{report.n_query}",
        f"- {cfg.DB_NOTE}",
        "",
        "| 化合物 | E3 | 状态 | pred_p2_rmsd | clash | BSA (Å²) | 装配 | 提示 |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in report.rows:
        if r.status == "ok":
            md.append(f"| {r.compound_id} | {r.e3} | ✅ | {r.pred_p2_rmsd:.2f} "
                      f"| {r.clash:.3f} | {r.bsa:.0f} | {r.prep.lig_kind} "
                      f"| {r.prep.prep_warning or '—'} |")
        else:
            md.append(f"| {r.compound_id} | {r.e3 or '—'} | ❌ {r.status} "
                      f"| — | — | — | — | {r.reason} |")
    md += ["", "## 解读边界", ""]
    md += [f"- {n}" for n in report.notes]
    md += ["", "---", "*仅供研究参考，不构成任何实验或立项决策。*"]
    return "\n".join(md)


def _methodology_notes() -> list[str]:
    return [
        "pred_p2_rmsd 是 DeepTernary 的姿态质量代理指标（对参照 POI 的预测 RMSD，"
        "越小越好）；clash 是重叠原子占比（>0.2 提示输入几何差，结果仅供参考）。",
        "BSA（蛋白-蛋白埋藏表面积）**只有相对比较意义**：装配协议偏紧使绝对值整体偏高"
        "（文献经验高降解带 1100–1500 Å²，本协议产出普遍高于此）。",
        "结果强烈依赖输入装配质量：POI/E3 warhead 对齐 RMSD 大、装配碰撞高或 "
        "lig_kind=relinked 的分子应人工检查复合物 3D 结构后再采信。",
        "POI 参照当前仅内置 galectin-3C（5E89）；E3 支持 CRBN（4CI2）/VHL（4W9H）。"
        "其它 POI/E3 需扩展参照包后可用。",
    ]


# ---------------------------------------------------------------- worker
def run_ternary_job(record: dict, control) -> dict:
    params = record["params"]
    job_dir = app_settings.jobs_dir / record["id"]
    work_dir = job_dir / "work"
    report_id = f"tern_{time.strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:6]}"
    report_dir = cfg.REPORTS_DIR / report_id

    # ---- 解析查询 ----
    control.checkpoint(progress=5, phase="parse", message="解析查询分子…")
    csv_files = sorted((job_dir / "inputs").glob("*.csv"))
    if csv_files:
        items = parse_queries_csv(csv_files[0], params.get("default_e3", ""))
    else:
        txt = sorted(list((job_dir / "inputs").glob("*.txt"))
                     + list((job_dir / "inputs").glob("*.smi")))
        items = parse_queries_text(txt[0].read_text(encoding="utf-8"),
                                   params.get("default_e3", "")) if txt else []
    if not items:
        raise RuntimeError("没有可评估的查询分子")
    if len(items) > cfg.MAX_QUERY:
        items = items[:cfg.MAX_QUERY]
    seeds = int(params.get("seeds") or cfg.DEFAULT_SEEDS)

    # ---- 输入制备（主环境，逐分子）----
    control.checkpoint(progress=8, phase="prepare", message="加载参照晶体…")
    pack = get_pack()
    rows: list[TernaryRow] = []
    infer_tasks: list[dict] = []
    n = len(items)
    for i, it in enumerate(items):
        if control.stop_requested:
            control.checkpoint(message="已请求停止")
        control.checkpoint(progress=8 + int(37 * i / max(n, 1)),
                           message=f"输入制备 {i + 1}/{n}：{it['compound_id']}…")
        cid = it["compound_id"]
        input_dir = work_dir / "inputs" / cid
        input_dir.mkdir(parents=True, exist_ok=True)
        prep = build_inputs(it["smiles"], it.get("e3", ""), input_dir, pack, cid=cid)
        if prep["status"] != "ok":
            rows.append(TernaryRow(compound_id=cid, smiles=it["smiles"],
                                   e3=it.get("e3", ""), status="prep_failed",
                                   reason=prep["reason"]))
            continue
        rows.append(TernaryRow(
            compound_id=cid, smiles=it["smiles"], e3=prep["e3"], status="ok",
            prep=PrepInfo(mcs_poi=prep["mcs_poi"], mcs_e3=prep["mcs_e3"],
                          poi_rmsd=prep["poi_rmsd"], e3_rmsd=prep["e3_rmsd"],
                          lig_kind=prep["lig_kind"], best_clear=prep["best_clear"],
                          prep_warning=prep.get("warning") or "")))
        infer_tasks.append({"name": cid, "input_dir": input_dir,
                            "out_dir": work_dir / "predictions" / cid,
                            "seeds": seeds})
        rows[-1].complex_file = ""

    if not infer_tasks:
        report = TernaryReport(id=report_id, created_at=time.strftime("%Y-%m-%d %H:%M:%S"),
                               name=str(record.get("name") or ""), seeds=seeds,
                               n_query=len(items), n_ok=0, rows=rows,
                               notes=_methodology_notes())
        artifacts = save_report(report, report_dir)
        return {"stats": {"report_id": report_id, "n_query": len(items), "n_ok": 0},
                "message": f"全部 {len(items)} 个分子制备失败，无推理可执行",
                "artifacts": artifacts}

    # ---- DeepTernary 子进程批量推理 ----
    control.checkpoint(progress=48, phase="infer",
                       message=f"DeepTernary 推理 {len(infer_tasks)} 个分子"
                               f"（首次加载模型约 20–60 秒）…")
    for t in infer_tasks:
        t["out_dir"].mkdir(parents=True, exist_ok=True)
    results = _run_inference(infer_tasks, work_dir, control)
    by_name = {r["name"]: r for r in results}

    # ---- BSA + 汇总 ----
    report_dir.mkdir(parents=True, exist_ok=True)
    n_ok = 0
    for i, row in enumerate(rows):
        if row.status != "ok":
            continue
        control.checkpoint(progress=85 + int(10 * i / max(len(rows), 1)),
                           message=f"BSA 计算与汇总 {i + 1}/{len(rows)}…")
        res = by_name.get(row.compound_id)
        if res is None or "error" in res:
            row.status = "infer_failed"
            row.reason = (res or {}).get("error", "推理结果缺失")
            continue
        summary = pd.read_csv(res["summary_csv"]).sort_values("pred_p2_rmsd")
        best = summary.iloc[0]
        row.pred_p2_rmsd = round(float(best["pred_p2_rmsd"]), 3)
        row.clash = round(float(best["clash_ratio"]), 3)
        row.bsa = round(compute_bsa(best["complex_pred"]), 1)
        dest = report_dir / f"complex_{row.compound_id}.pdb"
        shutil.copyfile(best["complex_pred"], dest)
        row.complex_file = dest.name
        if row.prep.best_clear is not None and row.prep.best_clear <= 2.8:
            row.prep.prep_warning = (row.prep.prep_warning or
                                     f"装配最小间距 {row.prep.best_clear} Å，建议人工检查")
        if row.clash > 0.2:
            row.prep.prep_warning = (row.prep.prep_warning + "；" if row.prep.prep_warning else "")
            row.prep.prep_warning += f"clash={row.clash} 偏高，输入几何可能欠佳"
        n_ok += 1

    report = TernaryReport(id=report_id, created_at=time.strftime("%Y-%m-%d %H:%M:%S"),
                           name=str(record.get("name") or ""), seeds=seeds,
                           n_query=len(items), n_ok=n_ok, rows=rows,
                           notes=_methodology_notes())
    control.checkpoint(progress=97, phase="export", message="导出报告（JSON / CSV / Markdown）…")
    artifacts = save_report(report, report_dir)
    top_ids = [r.compound_id for r in report.rows if r.status == "ok"][:5]
    top = "、".join(top_ids)
    return {
        "stats": {"report_id": report_id, "n_query": report.n_query,
                  "n_ok": n_ok, "seeds": seeds,
                  "reports_dir": str(cfg.REPORTS_DIR)},
        "message": (f"三元复合物预测完成：{n_ok}/{len(items)} 成功"
                    + (f"，含 {top}" if top else "")),
        "artifacts": artifacts,
    }
