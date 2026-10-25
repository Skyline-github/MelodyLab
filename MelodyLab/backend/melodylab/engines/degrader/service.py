"""降解剂评估 worker：模型就绪 → 查询解析 → 共识排序 → 报告落盘。

JobManager worker 约定：(record, control) -> dict。
训练/predict 均为同步 CPU 计算，用 checkpoint 心跳让停止请求生效。
"""
from __future__ import annotations

import json
import time
from pathlib import Path

from ...config import settings as app_settings
from . import config as cfg
from .engine import rank_compounds
from .features import parse_query_csv, parse_smiles_text
from .models import DegraderReport
from .training import get_bundle, invalidate_cache


def _load_queries(inputs_dir: Path) -> list[dict]:
    csv_files = sorted(inputs_dir.glob("*.csv"))
    if csv_files:
        return parse_query_csv(csv_files[0])
    txt = sorted(list(inputs_dir.glob("*.txt")) + list(inputs_dir.glob("*.smi")))
    if txt:
        return parse_smiles_text(txt[0].read_text(encoding="utf-8"))
    raise RuntimeError("任务输入目录为空：没有可评估的查询分子")


def save_report(report: DegraderReport) -> list[dict]:
    json_path = cfg.REPORTS_DIR / f"{report.id}.json"
    csv_path = cfg.REPORTS_DIR / f"{report.id}.csv"
    md_path = cfg.REPORTS_DIR / f"{report.id}.md"
    json_path.write_text(report.model_dump_json(indent=2), encoding="utf-8")
    csv_path.write_text(render_csv(report), encoding="utf-8")
    md_path.write_text(render_markdown(report), encoding="utf-8")
    base = report.name or "降解剂评估"
    return [
        {"name": f"{base}_排序.json", "path": str(json_path), "kind": "json"},
        {"name": f"{base}_排序.csv", "path": str(csv_path), "kind": "csv"},
        {"name": f"{base}_报告.md", "path": str(md_path), "kind": "markdown"},
    ]


_CSV_COLS = ["compound_id", "e3", "pred_dc50_nm", "lo90_nm", "hi90_nm",
             "rank_consensus", "rank_mlp", "rank_et", "rank_hgb",
             "max_tanimoto", "in_db", "db_dc50", "status", "smiles", "warnings"]


def render_csv(report: DegraderReport) -> str:
    import pandas as pd

    rows = []
    for r in report.rows:
        d = r.model_dump()
        d["warnings"] = "；".join(d["warnings"])
        rows.append(d)
    return pd.DataFrame(rows, columns=_CSV_COLS).to_csv(index=False)


def render_markdown(report: DegraderReport) -> str:
    cal = report.calibration
    md = [
        f"# 降解剂 DC50 共识排序：{report.name or report.id}",
        "",
        f"- **评估时间**：{report.created_at}",
        f"- **模型状态**：{('首次训练' if report.model_status == 'trained' else '缓存模型')}"
        f"（训练集 {cal.n_train} 行，{cal.db_snapshot}）",
        f"- **查询分子**：{report.n_query}（无效 {report.n_invalid}，重复 {report.n_duplicates}）",
        f"- **骨架 {cal.n_folds} 折交叉验证**：" +
        " · ".join(f"{f} R²={m.r2}, ρ={m.spearman}" for f, m in cal.cv.items()) +
        f"，conformal q90 = {cal.q90_residual_log:.2f} log 单位",
        "",
        "> 本工具输出**相对排序参考**，不是绝对 DC50 预测。共识排名（三族均值）"
        "比单点数值可信；区间跨越数量级属正常。MD/TPSA 超分布的分子，实际细胞活性"
        "很可能受透膜性限制。",
        "",
        "## 共识排序（前 10）",
        "",
        "| 排名 | 化合物 | E3 | 预测 DC50 (nM) | 90% 区间 | MLP/ET/HGB 排名 | Tanimoto | 提示 |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in report.rows[:10]:
        md.append(
            f"| {r.rank_consensus} | {r.compound_id} | {r.e3 or '—'} "
            f"| {r.pred_dc50_nm} | {r.lo90_nm}–{r.hi90_nm} "
            f"| {r.rank_mlp}/{r.rank_et}/{r.rank_hgb} "
            f"| {r.max_tanimoto} | {'；'.join(r.warnings) or '—'} |")
    md += ["", "## 完整排序", "",
           "| 排名 | 化合物 | 预测 DC50 (nM) | 90% 区间 | 在库 | 提示 |", "|---|---|---|---|---|---|"]
    for r in report.rows:
        md.append(f"| {r.rank_consensus} | {r.compound_id} | {r.pred_dc50_nm} "
                  f"| {r.lo90_nm}–{r.hi90_nm} | {'✅ ' + r.db_dc50 if r.in_db else '—'} "
                  f"| {'；'.join(r.warnings) or '—'} |")
    md += ["", "## 方法说明", ""]
    md += [f"- {n}" for n in report.notes]
    md += ["", "---", f"*数据来源：{cfg.DB_SNAPSHOT}；仅供研究参考，不构成任何实验或立项决策。*"]
    return "\n".join(md)


def run_degrader_job(record: dict, control) -> dict:
    """JobManager worker（kind=degrader）。params: task=rank|retrain, name。"""
    params = record["params"]
    task = params.get("task", "rank")
    job_dir = app_settings.jobs_dir / record["id"]

    if task == "retrain":
        control.checkpoint(progress=5, phase="retrain", message="清除旧模型缓存…")
        invalidate_cache()
        control.checkpoint(progress=10, phase="train", message="重新训练（首次约 3–5 分钟）…")
        bundle, status = get_bundle(control)
        cal = bundle["calibration"]
        return {
            "stats": {"model_status": status, "n_train": bundle["n_train"],
                      "q90_residual_log": cal["q90_residual_log"],
                      "cache_dir": str(cfg.CACHE_DIR)},
            "message": (f"模型重训完成（{status}）：训练集 {bundle['n_train']} 行，"
                        f"q90 = {cal['q90_residual_log']:.2f} log 单位"),
            "artifacts": [],
        }

    control.checkpoint(progress=5, phase="load", message="解析查询分子…")
    items = _load_queries(job_dir / "inputs")
    if not items:
        raise RuntimeError("没有可评估的查询分子")
    if len(items) > cfg.MAX_QUERY:
        items = items[:cfg.MAX_QUERY]
    control.checkpoint(progress=12, phase="model",
                       message=f"准备模型（{len(items)} 个分子；若为首次运行将先训练，约 3–5 分钟）…")
    bundle, status = get_bundle(control)

    control.checkpoint(progress=60, phase="rank", message="三族模型推理与共识排名…")
    report = rank_compounds(items, bundle, model_status=status)
    report.name = str(record.get("name") or "")

    control.checkpoint(progress=85, phase="export", message="导出报告（JSON / CSV / Markdown）…")
    artifacts = save_report(report)
    top = "、".join(r.compound_id for r in report.rows[:5])
    return {
        "stats": {
            "report_id": report.id,
            "n_query": report.n_query, "n_invalid": report.n_invalid,
            "n_duplicates": report.n_duplicates,
            "n_train": report.calibration.n_train if report.calibration else None,
            "q90_residual_log": report.calibration.q90_residual_log if report.calibration else None,
            "model_status": status,
            "top5": top,
            "reports_dir": str(cfg.REPORTS_DIR),
        },
        "message": (f"排序完成：{report.n_query - report.n_invalid} 个有效分子，"
                    f"共识前 5：{top}"),
        "artifacts": artifacts,
    }
