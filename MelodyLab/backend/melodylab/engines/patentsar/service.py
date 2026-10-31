"""专利 SAR 挖掘引擎服务层（patent-sar-mining 并入）。

流程：专利 PDF（born-digital 文本层 / 扫描件）→ 结构化分子/活性数据集 →
R 基团分析 + 3D-QSAR → LLM 叙事总结（引用校验门）。

纯本地计算；LLM/VLM 为可选增强——未配置 SAR_LLM_API_KEY 时叙事总结自动降级为
确定性 digest 模式，VLM 标签消歧自动跳过，功能完整、批处理不中断。
"""
from __future__ import annotations

import json
import os
from pathlib import Path

from ...config import settings

DATA_DIR = Path(__file__).parent / "data"
OUTPUTS_DIR = settings.data_root / "patentsar"

_DEPENDENCIES = ("rdkit", "pymupdf", "cv2", "easyocr", "DECIMER",
                 "sklearn", "matplotlib", "numpy")


def _has(module: str) -> bool:
    import importlib.util
    try:
        return importlib.util.find_spec(module) is not None
    except (ImportError, ValueError):
        return False


def llm_config() -> dict:
    """LLM/VLM 配置：环境变量 SAR_LLM_*（base_url/model/api_key）+ SAR_VLM_MODEL。"""
    return {
        "base_url": os.environ.get("SAR_LLM_BASE_URL",
                                   "https://open.bigmodel.cn/api/paas/v4"),
        "model": os.environ.get("SAR_LLM_MODEL", "glm-4.6"),
        "vlm_model": os.environ.get("SAR_VLM_MODEL", "glm-4v-flash"),
        "api_key": os.environ.get("SAR_LLM_API_KEY", ""),
    }


def availability() -> dict:
    deps = {m: _has(m) for m in _DEPENDENCIES}
    born_digital = deps["rdkit"] and deps["pymupdf"] and deps["DECIMER"]
    scanned = born_digital and deps["cv2"] and deps["easyocr"]
    sar_ok = deps["rdkit"] and deps["sklearn"] and deps["matplotlib"] and deps["numpy"]
    llm = llm_config()
    missing = [m for m, ok in deps.items() if not ok]
    return {
        "engine": "patentsar",
        "available": bool(born_digital or scanned),
        "paths": {"born_digital": bool(born_digital), "scanned": bool(scanned)},
        "sar_analysis": bool(sar_ok),
        "llm": {"configured": bool(llm["api_key"]), "model": llm["model"],
                "vlm_model": llm["vlm_model"], "base_url": llm["base_url"]},
        "missing_dependencies": missing,
        "hint": ("" if born_digital else
                 "挖掘依赖缺失（pip install pymupdf decimer opencv-python easyocr）"
                 if missing else ""),
    }


def has_text_layer(pdf: Path) -> bool:
    """born-digital / 扫描件自动判别：前几页有足够文本层即 born-digital。"""
    import pymupdf
    doc = pymupdf.open(str(pdf))
    try:
        n = min(3, doc.page_count)
        chars = sum(len(doc[i].get_text()) for i in range(n))
    finally:
        doc.close()
    return chars > 100


def example_info() -> dict:
    """内置示例（合成硬基准 fixture：成对异构体 + 立体化学 + Markush R 表）。"""
    pdf = DATA_DIR / "fixture_patent_hard.pdf"
    return {
        "pdf": str(pdf) if pdf.is_file() else None,
        "name": "合成硬基准专利（born-digital：9a/9b 对映体同页、双手性中心、Markush R 表）",
        "params": {
            "path_type": "born_digital",
            "structure_pages": "1:5",
            "assay_pages": "6",
            "mz_pages": "7",
            "markush_pages": "5",
            "target": "AMY3R",
        },
    }


# ---------------------------------------------------------------------------
# 同步轻接口
# ---------------------------------------------------------------------------

def sar_analysis(dataset: Path, activity_col: str, target: str, outdir: Path,
                 core: str | None = None) -> dict:
    """R 基团分析 + Free-Wilson + 3D-QSAR（秒级，同步返回）。"""
    from .core.sar import run as sar_run

    res = sar_run.run(str(dataset), smiles_col="canonical_smiles",
                      id_col="example_no", activity_col=activity_col,
                      target=target, outdir=Path(outdir),
                      core=core)
    return res


def summarize_dataset(dataset: Path, target: str, outdir: Path,
                      rgroup_table: Path | None = None,
                      activity_col: str | None = None,
                      lang: str = "zh") -> dict:
    """叙事总结（LLM 或 digest），返回 markdown 文本供前端直接渲染。"""
    from .core.sar import summarize as su

    res = su.summarize(str(dataset), target, str(outdir),
                       rgroup_table=str(rgroup_table) if rgroup_table else None,
                       activity_col=activity_col, lang=lang)
    md_path = Path(res["outdir"]) / "sar_summary.md"
    return {"result": res, "markdown": md_path.read_text(encoding="utf-8")}


def selftest() -> dict:
    """离线自检：VLM 响应解析 + 引用校验门（不访问网络）。"""
    from .core.pipeline import vlm
    from .core.sar import summarize as su

    vlm_fails = vlm._selftest()
    rows = su.load_activity_rows([DATA_DIR / "sar_series.csv"], None, None)
    pack = su.build_pack(rows, "AMY3R")
    gate_fails = su._selftest(pack)
    return {"vlm_parse_fails": vlm_fails, "citation_gate_fails": gate_fails,
            "ok": vlm_fails == 0 and gate_fails == 0}


# ---------------------------------------------------------------------------
# 长任务 worker（/mine）
# ---------------------------------------------------------------------------

def run_mine_job(record: dict, control) -> dict:
    """任务中心 worker：挖掘 → SAR 分析 → 叙事总结，产物全部注册到 artifacts。"""
    from ...core.jobs import manager  # 局部导入避免循环

    p = dict(record.get("params") or {})
    job_dir = manager.job_dir(record["id"])
    inputs, outdir = job_dir / "inputs", job_dir / "outputs"
    outdir.mkdir(parents=True, exist_ok=True)

    pdf = None
    for cand in sorted(inputs.glob("*.pdf")):
        pdf = cand
        break
    if p.get("example"):
        bundled = DATA_DIR / "fixture_patent_hard.pdf"
        if bundled.is_file():
            pdf = bundled
    if pdf is None or not Path(pdf).is_file():
        raise ValueError("未找到输入 PDF（请上传或使用内置示例）")

    path_type = p.get("path_type") or ("born_digital" if has_text_layer(pdf)
                                       else "scanned")
    base = Path(pdf).stem
    control.checkpoint(progress=5, phase="mine",
                       message=f"path_type={path_type}, pdf={Path(pdf).name}")
    stats: dict = {"path_type": path_type, "dataset": None}

    # ---- 1) 挖掘 ----
    if path_type == "scanned":
        from .core.pipeline import run_scanned
        assay = [tuple(x.rsplit(":", 1)) for x in
                 (p.get("assay_pages") or "").split(",") if ":" in x]
        assay = [(int(pg), [t for t in tg.split("+") if t]) for pg, tg in assay]
        mk_cfg = _write_markush_config(p, job_dir)

        def _locate_done(n_cand, n_total):
            control.checkpoint(progress=15, phase="locate",
                               message=f"页级定位完成：{n_total} 页中 {n_cand} 页含结构图，开始精确 OCR")
            stats["located_structure_pages"] = n_cand

        res = run_scanned.run_scanned(
            str(pdf),
            run_scanned._parse_range(p.get("structure_pages") or "0:10000"),
            [k.strip() for k in (p.get("kinds") or "example").split(",") if k.strip()],
            assay, outdir, job_dir / "work",
            zoom=float(p.get("zoom") or 3.0),
            text_zoom_mult=float(p.get("text_zoom_mult") or 2.0),
            mz_tol=float(p.get("mz_tol") or 1.5),
            hand_drawn=bool(p.get("hand_drawn")),
            vlm_cfg=_vlm_cfg(p),
            mz_pages=[int(x) for x in (p.get("mz_pages") or "").split(",") if x.strip()],
            markush_config=str(mk_cfg) if mk_cfg else None,
            markush_pages=[int(x) for x in (p.get("markush_pages") or "").split(",") if x.strip()],
            locate_progress=_locate_done)
    else:
        from .core.pipeline import run
        mk_cfg = _write_markush_config(p, job_dir)
        res = run.run(
            str(pdf),
            run._parse_range(p.get("structure_pages") or "0:10000"),
            run._page_list(p.get("assay_pages") or ""),
            run._page_list(p.get("mz_pages") or ""),
            outdir, job_dir / "work",
            hand_drawn=bool(p.get("hand_drawn")),
            mz_tol=float(p.get("mz_tol") or 1.5),
            markush_config=str(mk_cfg) if mk_cfg else None,
            markush_pages=[int(x) for x in (p.get("markush_pages") or "").split(",") if x.strip()])
    dataset = outdir / f"{base}_dataset.csv"
    if not dataset.is_file():
        raise RuntimeError(f"挖掘未产出数据集：{dataset}")
    control.checkpoint(progress=70, phase="sar", message="挖掘完成，进行 R 基团分析")

    # ---- 2) SAR 分析（无活性数据时自动跳过：R 基团/3D-QSAR 无活性无意义）----
    artifacts = [{"name": dataset.name, "path": str(dataset), "kind": "file"}]
    stats["dataset"] = dataset.name
    for ext in ("structures.sdf", "check.html", "report.md"):
        f = outdir / f"{base}_{ext}"
        if f.is_file():
            artifacts.append({"name": f.name, "path": str(f), "kind": "file"})
    sar_dir = outdir / "sar"
    sar_stats = {}
    has_act, act_msg = _has_activity_data(dataset)
    if not has_act:
        sar_stats["skipped_sar"] = act_msg
        stats.update(sar_stats)
        control.checkpoint(progress=95, phase="done",
                           message=f"挖掘完成（无活性数据，已跳过 SAR/总结）：{dataset.name}")
        return {"message": f"完成：{path_type} 路径，结构数据集 {dataset.name}（无活性数据，跳过 SAR）",
                "stats": stats, "artifacts": artifacts}
    if p.get("with_sar", True):
        try:
            activity_col = p.get("activity_col") or _first_ec50_col(dataset)
            sres = sar_analysis(dataset, activity_col,
                                p.get("target") or "TARGET", sar_dir,
                                core=p.get("core") or None)
            sar_stats = {"rga_positions": sres.get("positions"),
                         "n_compounds": sres.get("n")}
            for f in sorted(sar_dir.glob("*")):
                if f.is_file():
                    artifacts.append({"name": f"sar/{f.name}", "path": str(f), "kind": "file"})
        except Exception as e:  # noqa: BLE001 - SAR 失败不拖垮挖掘结果
            sar_stats = {"error": f"{type(e).__name__}: {e}"}
        control.checkpoint(progress=85, phase="summarize", message="SAR 分析完成，生成叙事总结")

    # ---- 3) 叙事总结 ----
    summary_md = None
    if p.get("with_summary", True):
        try:
            s = summarize_dataset(dataset, p.get("target") or "TARGET", sar_dir,
                                  rgroup_table=sar_dir / "rgroup_table.csv"
                                  if (sar_dir / "rgroup_table.csv").is_file() else None,
                                  activity_col=p.get("activity_col") or None)
            summary_md = s["markdown"]
            for name in ("sar_summary.md", "sar_summary.json"):
                f = sar_dir / name
                if f.is_file():
                    artifacts.append({"name": f"sar/{name}", "path": str(f), "kind": "file"})
        except Exception as e:  # noqa: BLE001
            sar_stats["summary_error"] = f"{type(e).__name__}: {e}"

    stats.update(sar_stats)
    if summary_md:
        (outdir / "sar_summary_preview.md").write_text(summary_md, encoding="utf-8")
        stats["summary_headline"] = _headline(summary_md)
    return {"message": f"完成：{path_type} 路径，数据集 {dataset.name}",
            "stats": stats, "artifacts": artifacts}


def _headline(md: str) -> str:
    for line in md.splitlines():
        s = line.strip()
        if s and not s.startswith(("#", ">", "-", "|", "*")):
            return s[:200]
    return ""


def _first_ec50_col(dataset: Path) -> str:
    import csv
    with open(dataset, newline="") as fh:
        header = next(csv.reader(fh))
    return next((c for c in header if c.lower().startswith("ec50")), "ec50_amy3r")


def _has_activity_data(dataset: Path) -> tuple[bool, str]:
    """数据集是否含可用的活性数据：存在 ec50/ic50/ki/emax 列且至少一个真实数值。

    注意 write_csv 把缺失活性写成字面 "ND"——那是"无数据"标记，不算有活性。
    """
    import csv
    blanks = ("", "nd", "n/a", "na", "none", "-")
    with open(dataset, newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh)
        header = next(reader, [])
        act_cols = [i for i, c in enumerate(header)
                    if any(k in c.lower() for k in ("ec50", "ic50", "ki_", "emax", "activity"))]
        if not act_cols:
            return False, "数据集无活性列（ec50/ic50/ki/emax）"
        for row in reader:
            if any(row[i].strip().lower() not in blanks for i in act_cols if i < len(row)):
                return True, ""
    return False, "活性列全部为空（专利未附活性数据）"


def _write_markush_config(p: dict, job_dir: Path) -> Path | None:
    """markush 配置三来源：内嵌 JSON 文本 / bundled 示例 / None。"""
    text = (p.get("markush_config") or "").strip()
    if text:
        cfg_path = job_dir / "work" / "markush_config.json"
        cfg_path.parent.mkdir(parents=True, exist_ok=True)
        json.loads(text)  # 尽早暴露 JSON 错误
        cfg_path.write_text(text, encoding="utf-8")
        return cfg_path
    if p.get("example"):
        bundled = DATA_DIR / "markush_fixture.json"
        if bundled.is_file():
            return bundled
    return None


def _vlm_cfg(p: dict):
    from .core.pipeline import vlm
    if (p.get("vlm") or "auto") == "off":
        return None
    return vlm.default_config()
