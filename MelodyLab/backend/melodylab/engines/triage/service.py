"""靶点立项评估 worker：执行五问评估并导出 JSON / Markdown / Excel 报告。

JobManager worker 约定：(record, control) -> dict。
评估本体 engine.assess 是 asyncio 协程，在线程内用独立事件循环跑；
跑的时候另起线程 + 心跳轮询，让停止/暂停请求能在取数间隙生效。
"""
from __future__ import annotations

import asyncio
import json
import threading

from ...config import settings as app_settings
from .engine import assess
from .models import AssessmentReport, STATUS_LABEL

STATUS_BADGE = {"green": "🟢 证据充分", "yellow": "🟡 需注意",
                "red": "🔴 风险信号", "unknown": "⚪ 信息不足"}


# ---------------------------------------------------------------- 导出
def _ev_md(ev) -> list[str]:
    icon = {"green": "🟢", "yellow": "🟡", "red": "🔴", "unknown": "⚪"}[ev.status]
    lines = [f"### {icon} {ev.source} · {ev.check}", ""]
    lines.append(ev.summary)
    if ev.detail:
        lines.append("")
        lines.append("```json")
        lines.append(json.dumps(ev.detail, ensure_ascii=False, indent=2)[:2000])
        lines.append("```")
    if ev.links:
        lines.append("")
        lines.append(" · ".join(f"[{l.label}]({l.url})" for l in ev.links))
    lines.append("")
    return lines


def render_markdown(report: AssessmentReport) -> str:
    md = [
        f"# 靶点立项五问评估：{report.gene_symbol}",
        "",
        f"- **靶点全名**：{report.gene_name or '—'}（Ensembl: {report.ensembl_id}, "
        f"UniProt: {report.uniprot_id or '—'}）",
        f"- **限定适应症**：{report.indication or '不限'}",
        f"- **评估时间**：{report.created_at}",
        f"- **总体状态**：{STATUS_BADGE[report.overall_status]}",
        "",
    ]
    if report.risk_flags:
        md += ["## 🚩 风险信号（数据筛查提示，需人工判断，不等于结论）", ""]
        md += [f"- {v}" for v in report.risk_flags] + [""]
    if report.watch_flags:
        md += ["## 👀 建议进一步核查", ""]
        md += [f"- {v}" for v in report.watch_flags] + [""]
    md += ["> 五个维度是乘法关系：立项阶段真正要找的不是哪个维度得分最高，"
           "而是哪个维度会成为短板。红/黄/绿为阈值规则自动生成的筛查信号，不是结论。", ""]
    for q in report.questions:
        md += [f"## {q.title}：{q.question}", ""]
        for ev in q.evidences:
            md += _ev_md(ev)
        if q.manual_checklist:
            md += ["## 手动核查清单", ""]
            for ev in q.manual_checklist:
                md += [f"- **{ev.source} · {ev.check}**：{ev.summary}  "
                       + " · ".join(f"[{l.label}]({l.url})" for l in ev.links)]
            md.append("")
    md += ["---", "*本报告由开源数据自动生成，仅供研究参考，不构成任何投资或立项建议。*"]
    return "\n".join(md)


def export_excel(report: AssessmentReport, path) -> bool:
    """Excel 证据表（pandas + openpyxl）；openpyxl 缺失时返回 False 不阻塞。"""
    try:
        import pandas as pd
    except ImportError:
        return False
    auto_rows: list[dict] = []
    manual_rows: list[dict] = []
    for q in report.questions:
        for ev in q.evidences:
            auto_rows.append({
                "五问": q.title, "来源": ev.source, "检查点": ev.check,
                "状态": STATUS_BADGE[ev.status], "结论": ev.summary,
                "详情(JSON)": json.dumps(ev.detail, ensure_ascii=False)[:1000] if ev.detail else "",
                "链接": " ; ".join(l.url for l in ev.links),
            })
        for ev in q.manual_checklist:
            manual_rows.append({
                "五问": q.title, "来源": ev.source, "检查点": ev.check,
                "说明": ev.summary,
                "链接": " ; ".join(l.url for l in ev.links),
            })
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        pd.DataFrame(auto_rows).to_excel(writer, sheet_name="自动证据", index=False)
        if manual_rows:
            pd.DataFrame(manual_rows).to_excel(writer, sheet_name="手动核查清单", index=False)
        if report.risk_flags:
            pd.DataFrame({"风险信号（需人工判断）": report.risk_flags}).to_excel(
                writer, sheet_name="风险信号", index=False)
        if report.watch_flags:
            pd.DataFrame({"建议进一步核查": report.watch_flags}).to_excel(
                writer, sheet_name="建议核查", index=False)
    return True


def save_report(report: AssessmentReport) -> list[dict]:
    """报告三格式落盘到数据根 triage_reports/，返回 artifacts 描述。"""
    from .config import REPORTS_DIR

    rid = report.id
    json_path = REPORTS_DIR / f"{rid}.json"
    md_path = REPORTS_DIR / f"{rid}.md"
    xlsx_path = REPORTS_DIR / f"{rid}.xlsx"
    json_path.write_text(report.model_dump_json(indent=2), encoding="utf-8")
    md_path.write_text(render_markdown(report), encoding="utf-8")
    artifacts = [
        {"name": f"{report.gene_symbol}_评估报告.json", "path": str(json_path), "kind": "json"},
        {"name": f"{report.gene_symbol}_评估报告.md", "path": str(md_path), "kind": "markdown"},
    ]
    if export_excel(report, xlsx_path):
        artifacts.append({"name": f"{report.gene_symbol}_评估报告.xlsx",
                          "path": str(xlsx_path), "kind": "xlsx"})
    return artifacts


# ---------------------------------------------------------------- worker
def run_triage_job(record: dict, control) -> dict:
    """JobManager worker（kind=triage）。params: gene, indication。"""
    params = record["params"]
    gene = str(params["gene"]).strip().upper()
    indication = str(params.get("indication") or "").strip()

    control.checkpoint(progress=5, phase="assess",
                       message=f"开始评估 {gene}（并发查询 9 个公开数据库，约 10–40 秒）…")

    box: dict = {}
    error: list[str] = []

    def _run():
        try:
            box["report"] = asyncio.run(assess(gene, indication))
        except BaseException as exc:  # noqa: BLE001 —— 转交主线程抛出
            error.append(f"{exc}")

    t = threading.Thread(target=_run, name="triage-assess", daemon=True)
    t.start()
    while t.is_alive():
        t.join(timeout=2.0)
        if t.is_alive():
            # 心跳：保持任务记录活跃，同时让停止/暂停请求及时生效
            control.checkpoint(message=f"正在并发取数与评估 {gene} …")
    if error:
        raise RuntimeError(f"评估失败：{error[0]}")

    report: AssessmentReport = box["report"]
    control.checkpoint(progress=80, phase="export",
                       message="评估完成，正在导出报告（JSON / Markdown / Excel）…")

    artifacts = save_report(report)
    n_red, n_yellow = len(report.risk_flags), len(report.watch_flags)
    return {
        "stats": {
            "report_id": report.id,
            "gene_symbol": report.gene_symbol,
            "gene_name": report.gene_name,
            "ensembl_id": report.ensembl_id,
            "overall_status": report.overall_status,
            "overall_label": STATUS_LABEL[report.overall_status],
            "risk_count": n_red,
            "watch_count": n_yellow,
            "app_dir": str(app_settings.data_root / "triage_reports"),
        },
        "message": (f"评估完成：总体筛查信号 {STATUS_LABEL[report.overall_status]}；"
                    f"风险信号 {n_red} 项 / 建议核查 {n_yellow} 项"),
        "artifacts": artifacts,
    }
