"""GEO 表达差异分析流水线（对齐 DockScope geo_array_pipeline v6.2 思路）。

GSE id → NCBI 下载 series matrix → 解析表达矩阵 + 样本元数据 → 分组（自动/手动）
→ QC（PCA/样本相关性）→ 差异分析（welch t-test + BH-FDR）→ 火山图/热图 PNG
→ HTML 报告 + CSV/Excel。
"""
from __future__ import annotations

import gzip
import html
import io
import json
import re
import urllib.request
from pathlib import Path

import numpy as np

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"


def _http(url: str, timeout: float = 60.0, retries: int = 5) -> bytes:
    """带重试的下载（代理/NCBI 常见间歇性 SSL 重置）。"""
    import ssl
    import time

    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "MelodyLab/0.2"})
            ctx = ssl.create_default_context()
            with urllib.request.urlopen(req, timeout=timeout, context=ctx) as r:
                return r.read()
        except Exception as exc:
            last = exc
            time.sleep(2.0 * (attempt + 1))
    raise last


def _esummary(gse: str) -> dict:
    q = _http(f"{EUTILS}/esearch.fcgi?db=gds&term={gse}%5BACCN%5D&retmode=json&tool=MelodyLab")
    ids = json.loads(q).get("esearchresult", {}).get("idlist", [])
    if not ids:
        raise ValueError(f"GEO 未找到 {gse}")
    s = _http(f"{EUTILS}/esummary.fcgi?db=gds&id={ids[0]}&retmode=json&tool=MelodyLab")
    doc = json.loads(s).get("result", {})
    d = doc.get(ids[0], {})
    return {"accession": d.get("accession") or gse, "title": d.get("title"),
            "summary": d.get("summary"), "organism": d.get("taxon"),
            "n_samples": d.get("n_samples"), "platform": d.get("gpl")}


def download_series_matrix(gse: str, dest: Path) -> Path:
    """HTTPS 优先下载 GSExxx_series_matrix.txt(.gz)。"""
    gse = gse.upper().strip()
    prefix = f"{FTP}/{gse[:-3]}nnn/{gse}/matrix/"
    for name in (f"{gse}_series_matrix.txt.gz", f"{gse}_series_matrix.txt"):
        try:
            data = _http(prefix + name, timeout=120.0)
            p = dest / name
            p.write_bytes(data)
            return p
        except Exception:
            continue
    raise RuntimeError(f"series matrix 下载失败：{gse}（检查网络/GSE 是否为表达谱芯片）")


def parse_series_matrix(path: Path) -> dict:
    """解析 !Series/!Sample 元数据与表达矩阵。"""
    op = gzip.open if path.suffix == ".gz" else open
    series_meta: dict[str, list] = {}
    sample_meta: dict[str, dict] = {}
    matrix_rows: list[list[str]] = []
    in_matrix = False
    with op(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            # 注意顺序：table 标记以 "!series_" 开头，必须先于元数据分支判断
            if line.startswith("!series_matrix_table_begin"):
                in_matrix = True
            elif line.startswith("!series_matrix_table_end"):
                in_matrix = False
            elif in_matrix and line.strip():
                matrix_rows.append(line.split("\t"))
            elif line.startswith("!series_"):
                key = line.split("\t", 1)[0]
                series_meta.setdefault(key, []).append([c.strip('"') for c in line.split("\t")[1:]])
            elif line.startswith("!Sample_"):
                parts = line.split("\t")
                sample_meta.setdefault(parts[0], []).append(
                    [c.strip('"') for c in parts[1:]])
    if not matrix_rows:
        raise ValueError("series matrix 中没有表达矩阵")
    header = matrix_rows[0]
    probes = [r[0].strip('"') for r in matrix_rows[1:]]
    samples = [c.strip('"') for c in header[1:]]
    values = np.array([[float(x) if x not in ("", "null", "NaN") else np.nan
                        for x in r[1:len(samples) + 1]] for r in matrix_rows[1:]],
                      dtype=np.float64)
    sdata = {}
    for k, rows in sample_meta.items():
        # !Sample_* 同名键可出现多行（如多个 characteristics）；保留为行列表
        sdata[k] = [r[:len(samples)] + [""] * max(0, len(samples) - len(r)) for r in rows]
    return {"samples": samples, "probes": probes, "values": values,
            "sample_meta": sdata, "series_meta": series_meta}


def platform_annotation(gpl: str, cache: Path) -> dict[str, str]:
    """GPL → 探针-基因映射（优先 soft 文件中的 gene assignment 列）。缓存。"""
    cache.mkdir(parents=True, exist_ok=True)
    cfile = cache / f"{gpl}.json"
    if cfile.is_file():
        return json.loads(cfile.read_text(encoding="utf-8"))
    gpl_upper = gpl.upper()
    url = f"{FTP.replace('/series', '/platforms')}/{gpl_upper[:-3]}nnn/{gpl_upper}/annot/{gpl_upper}_annot.gz"
    mapping: dict[str, str] = {}
    try:
        raw = _http(url, timeout=180.0)
        text = gzip.decompress(raw).decode("utf-8", errors="replace")
        header = None
        gene_col = None
        id_col = None
        for line in text.splitlines():
            if line.startswith("!platform_table_begin"):
                continue
            if line.startswith("!platform_table_end"):
                break
            if header is None:
                cols = [c.strip('"') for c in line.split("\t")]
                header = cols
                id_col = 0
                for i, c in enumerate(cols):
                    if "gene assignment" in c.lower() or "gene symbol" in c.lower():
                        gene_col = i
                        break
                if gene_col is None:
                    break
                continue
            parts = line.split("\t")
            if len(parts) <= max(gene_col, id_col):
                continue
            sym = parts[gene_col]
            if "//" in sym:  # "symbol // desc // location" 取第一段
                sym = sym.split("//")[0].strip()
            m = re.findall(r"([A-Z][A-Z0-9]{1,9})\b", sym)
            if m:
                mapping[parts[id_col].strip('"')] = m[0]
    except Exception:
        mapping = {}
    cfile.write_text(json.dumps(mapping, ensure_ascii=False), encoding="utf-8", newline="\n")
    return mapping


def assign_groups(samples: list[str], sample_meta: dict, manual: dict | None = None) -> dict[str, list[str]]:
    """分组：手动覆盖（{"control": [...], "treatment": [...]}）或按来源关键词自动。"""
    if manual and manual.get("control") and manual.get("treatment"):
        return {"control": manual["control"], "treatment": manual["treatment"]}
    # 只用短分类字段做分组——protocol 长段落往往同时含两组关键词，会污染判定
    GROUP_KEYS = ("!Sample_characteristics_ch1", "!Sample_title",
                  "!Sample_source_name_ch1")
    texts = [""] * len(samples)
    for key in GROUP_KEYS:
        for row in sample_meta.get(key, []):
            for i in range(min(len(row), len(samples))):
                texts[i] += " " + str(row[i])
    texts = [x.lower() for x in texts]
    kw_ctl = ("control", "normal", "healthy", "non-tumor", "benign", "sham", "wild type", "wt", "nc", "n.c.")
    kw_trt = ("knockdown", "kd", "sirna", "overexpression", "oe", "ko", "infected", "stimulated",
              "tumor", "cancer", "disease", "treated")
    import re as _re

    def _match(text: str, kws: tuple) -> bool:
        for k in kws:
            if len(k) <= 3:
                if _re.search(r"\b" + _re.escape(k) + r"\d*\b", text):
                    return True
            elif k in text:
                return True
        return False

    # control 优先：命中 control 关键词的样本（如 "Control: NC siRNA"）不再参与 treatment 判定
    control = [samples[i] for i, t in enumerate(texts) if _match(t, kw_ctl)]
    treatment = [samples[i] for i, t in enumerate(texts) if _match(t, kw_trt) and samples[i] not in control]
    control = [s for s in dict.fromkeys(control) if s not in treatment]
    if not control or not treatment:
        raise ValueError(
            f"自动分组失败（control={len(control)}, treatment={len(treatment)}）；"
            "请用 manual_groups 显式指定 control/treatment 样本列表")
    return {"control": control, "treatment": treatment}


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR 校正。"""
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def differential_expression(values: np.ndarray, probes: list[str],
                            samples: list[str], groups: dict[str, list[str]],
                            gene_map: dict[str, str]) -> list[dict]:
    """Welch t-test + BH-FDR；返回按 p 排序的基因表。"""
    from scipy import stats

    ci = [samples.index(s) for s in groups["control"] if s in samples]
    ti = [samples.index(s) for s in groups["treatment"] if s in samples]
    # 尺度检测：RMA/limma 数据是 log2 尺寸（值 < 100），log2FC 用差值；线性数据用比值
    finite = values[np.isfinite(values)]
    is_log = bool(finite.size and finite.max() < 100 and finite.min() >= -20)
    rows = []
    for r, probe in enumerate(probes):
        x_c = values[r, ci]
        x_t = values[r, ti]
        x_c = x_c[~np.isnan(x_c)]
        x_t = x_t[~np.isnan(x_t)]
        if len(x_c) < 2 or len(x_t) < 2:
            continue
        t, p = stats.ttest_ind(x_t, x_c, equal_var=False)
        mc, mt = float(x_c.mean()), float(x_t.mean())
        fc = (mt - mc) if is_log else float(np.log2((mt + 1e-9) / (mc + 1e-9)))
        rows.append({
            "probe": probe,
            "gene": gene_map.get(probe, probe),
            "log2fc": float(fc),
            "mean_control": mc, "mean_treatment": mt,
            "t": float(t), "p_value": float(p),
        })
    if rows:
        ps = np.array([r["p_value"] for r in rows])
        for r, q in zip(rows, bh_fdr(ps)):
            r["p_adj"] = float(q)
        rows.sort(key=lambda r: r["p_value"])
    return rows


def qc_plots(values: np.ndarray, samples: list[str], groups: dict, out_dir: Path) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy import stats as st

    files = []
    X = np.nan_to_num(values, nan=0.0)
    # PCA
    Xc = X - X.mean(axis=1, keepdims=True)
    U, S, _ = np.linalg.svd(Xc.T, full_matrices=False)
    P = U[:, :2] * S[:2]
    fig, ax = plt.subplots(figsize=(6, 5), dpi=130)
    for gname, color in (("control", "#4477aa"), ("treatment", "#ee6677")):
        idx = [samples.index(s) for s in groups.get(gname, []) if s in samples]
        if idx:
            ax.scatter(P[idx, 0], P[idx, 1], label=gname, c=color, alpha=0.85)
    for i, s in enumerate(samples):
        ax.annotate(s[:10], (P[i, 0], P[i, 1]), fontsize=5, alpha=0.7)
    ax.set_title("PCA (samples)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    p = out_dir / "qc_pca.png"
    fig.savefig(p); plt.close(fig)
    files.append(p)
    # 样本相关性热图
    corr = np.corrcoef(X.T)
    fig, ax = plt.subplots(figsize=(6 + 0.12 * len(samples), 5 + 0.12 * len(samples)), dpi=110)
    im = ax.imshow(corr, cmap="RdBu_r", vmin=0.4, vmax=1.0)
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels([s[:8] for s in samples], rotation=90, fontsize=5)
    ax.set_yticks(range(len(samples)))
    ax.set_yticklabels([s[:8] for s in samples], fontsize=5)
    fig.colorbar(im, shrink=0.8)
    ax.set_title("Sample correlation")
    fig.tight_layout()
    p = out_dir / "qc_corr.png"
    fig.savefig(p); plt.close(fig)
    files.append(p)
    return files


def volcano_ma_plots(de_table: list[dict], out_dir: Path) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    files = []
    x = np.array([r["log2fc"] for r in de_table])
    p = np.array([r["p_value"] for r in de_table])
    p = np.clip(p, 1e-300, 1.0)
    y = -np.log10(p)
    sig = (np.array([r["p_adj"] for r in de_table]) < 0.05) & (np.abs(x) > 1)
    fig, ax = plt.subplots(figsize=(6.5, 5), dpi=130)
    ax.scatter(x[~sig], y[~sig], s=4, c="#bbbbbb", label="NS")
    ax.scatter(x[sig], y[sig], s=6, c="#cc3311", label="FDR<0.05 & |log2FC|>1")
    ax.axhline(-np.log10(0.05), ls="--", lw=0.8, c="#888888")
    ax.set_xlabel("log2 FC"); ax.set_ylabel("-log10 p")
    ax.set_title("Volcano"); ax.legend(fontsize=8)
    fig.tight_layout()
    fp = out_dir / "volcano.png"
    fig.savefig(fp); plt.close(fig)
    files.append(fp)
    # MA
    a = 0.5 * (np.array([r["mean_control"] for r in de_table]) +
               np.array([r["mean_treatment"] for r in de_table]))
    fig, ax = plt.subplots(figsize=(6.5, 5), dpi=130)
    ax.scatter(a[~sig], x[~sig], s=4, c="#bbbbbb")
    ax.scatter(a[sig], x[sig], s=6, c="#cc3311")
    ax.axhline(0, ls="--", lw=0.8, c="#888888")
    ax.set_xlabel("mean expression"); ax.set_ylabel("log2 FC")
    ax.set_title("MA plot")
    fig.tight_layout()
    fp = out_dir / "ma.png"
    fig.savefig(fp); plt.close(fig)
    files.append(fp)
    return files


def top_heatmap(values: np.ndarray, probes: list[str], samples: list[str],
                groups: dict, gene_map: dict, de_table: list[dict], out_dir: Path,
                top_n: int = 50) -> Path | None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    top = de_table[:top_n]
    if not top:
        return None
    idx = [probes.index(r["probe"]) for r in top if r["probe"] in probes]
    cols = groups["control"] + groups["treatment"]
    ci = [samples.index(s) for s in cols if s in samples]
    M = values[np.ix_(idx, ci)]
    M = np.nan_to_num(M, nan=0.0)
    M = M - M.mean(axis=1, keepdims=True)
    sd = M.std(axis=1, keepdims=True)
    M = M / np.where(sd == 0, 1, sd)
    fig, ax = plt.subplots(figsize=(1 + 0.16 * len(ci), 0.18 * len(idx) + 1.5), dpi=110)
    im = ax.imshow(M, aspect="auto", cmap="bwr", vmin=-2.5, vmax=2.5)
    labels = [gene_map.get(probes[i], probes[i]) for i in idx]
    ax.set_yticks(range(len(idx)))
    ax.set_yticklabels(labels, fontsize=5)
    ax.set_xticks(range(len(ci)))
    ax.set_xticklabels([samples[i][:8] for i in ci], rotation=90, fontsize=5)
    for j, s in enumerate(cols):
        if j < len(cols):
            gname = "control" if s in groups["control"] else "treatment"
            ax.text(j, -1.5, gname[0].upper(), ha="center", fontsize=6,
                    color="#4477aa" if gname == "control" else "#ee6677")
    fig.colorbar(im, shrink=0.6)
    ax.set_title(f"Top {len(idx)} DE genes (z-score)")
    fig.tight_layout()
    fp = out_dir / "heatmap_top.png"
    fig.savefig(fp); plt.close(fig)
    return fp


def write_html_report(out_dir: Path, meta: dict, groups: dict, de_table: list[dict],
                      stats: dict) -> Path:
    sig = [r for r in de_table if r.get("p_adj", 1) < 0.05]
    rows = "\n".join(
        f"<tr><td>{html.escape(str(r['gene']))}</td><td>{r['probe']}</td>"
        f"<td>{r['log2fc']:.3f}</td><td>{r['p_value']:.2e}</td><td>{r['p_adj']:.2e}</td></tr>"
        for r in sig[:100])
    doc = f"""<!doctype html><html lang="zh"><head><meta charset="utf-8">
<title>GEO 差异分析报告 — {html.escape(meta.get('accession', ''))}</title>
<style>body{{font-family:system-ui,sans-serif;margin:24px;max-width:960px}}
table{{border-collapse:collapse;width:100%;font-size:13px}}
td,th{{border:1px solid #ddd;padding:4px 8px;text-align:left}}
img{{max-width:100%;margin:12px 0}}h2{{margin-top:28px}}</style></head><body>
<h1>GEO 差异分析报告</h1>
<p><b>数据集：</b>{html.escape(str(meta.get('accession')))} — {html.escape(str(meta.get('title')))}<br>
<b>物种：</b>{html.escape(str(meta.get('organism')))} |
<b>平台：</b>{html.escape(str(meta.get('platform')))} |
<b>样本数：</b>{meta.get('n_samples')}<br>
<b>分组：</b>control = {len(groups.get('control', []))} 样本，treatment = {len(groups.get('treatment', []))} 样本</p>
<h2>统计摘要</h2>
<pre>{html.escape(json.dumps(stats, ensure_ascii=False, indent=1))}</pre>
<h2>QC</h2><img src="qc_pca.png"><img src="qc_corr.png">
<h2>差异结果</h2><img src="volcano.png"><img src="ma.png"><img src="heatmap_top.png">
<h2>显著差异基因（FDR&lt;0.05，前 100）</h2>
<table><tr><th>gene</th><th>probe</th><th>log2FC</th><th>p</th><th>p_adj</th></tr>
{rows}</table>
<p style="color:#888">方法：Welch t-test + Benjamini-Hochberg FDR（MelodyLab 0.2 GEO 引擎）</p>
</body></html>"""
    fp = out_dir / "geo_report.html"
    fp.write_text(doc, encoding="utf-8", newline="\n")
    return fp


def run_geo_job(record: dict, control) -> dict:
    """JobManager worker（kind=geo）。params: gse / manual_groups / top_n。"""
    import csv

    import pandas as pd

    from ...config import settings

    params = record["params"]
    job_dir = settings.jobs_dir / record["id"]
    work, outputs = job_dir / "work", job_dir / "outputs"
    work.mkdir(parents=True, exist_ok=True)
    outputs.mkdir(exist_ok=True)

    gse = params["gse"]
    control.checkpoint(progress=5, phase="download", message=f"获取 {gse} 元数据…")
    meta = _esummary(gse)
    control.checkpoint(progress=10, phase="download", message="下载 series matrix…")
    mat_path = download_series_matrix(gse, work)
    control.checkpoint(progress=25, phase="parse", message="解析矩阵与样本元数据…")
    data = parse_series_matrix(mat_path)
    samples, probes, values = data["samples"], data["probes"], data["values"]
    control.checkpoint(progress=35, phase="annotation",
                       message=f"平台注释映射（{meta.get('platform')}）…")
    gene_map = platform_annotation(meta.get("platform") or "", work / "gpl_cache")
    control.checkpoint(progress=45, phase="groups", message="样本分组…")
    groups = assign_groups(samples, data["sample_meta"], params.get("manual_groups"))
    control.checkpoint(progress=55, phase="qc", message="QC（PCA/相关性）…")
    qc_plots(values, samples, groups, outputs)
    control.checkpoint(progress=70, phase="de", message="差异分析（Welch t-test + BH-FDR）…")
    table = differential_expression(values, probes, samples, groups, gene_map)
    volcano_ma_plots(table, outputs)
    top_heatmap(values, probes, samples, groups, gene_map, table, outputs)

    csv_path = outputs / "de_table.csv"
    with open(csv_path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(table[0].keys()) if table else
                           ["probe", "gene", "log2fc", "mean_control", "mean_treatment", "t", "p_value", "p_adj"])
        w.writeheader()
        w.writerows(table)
    excel_path = outputs / "de_table.xlsx"
    pd.DataFrame(table).to_excel(excel_path, index=False)
    stats = {"n_probes_tested": len(table),
             "n_sig_fdr05": sum(1 for r in table if r.get("p_adj", 1) < 0.05),
             "n_sig_fdr05_fc1": sum(1 for r in table if r.get("p_adj", 1) < 0.05 and abs(r["log2fc"]) > 1),
             "groups": {k: len(v) for k, v in groups.items()}}
    report = write_html_report(outputs, meta, groups, table, stats)
    control.checkpoint(progress=95, phase="report", message="报告完成")
    top_genes = [{"gene": r["gene"], "log2fc": round(r["log2fc"], 3),
                  "p_adj": r["p_adj"]} for r in table[:50]]
    return {"stats": {**stats, "top_genes": top_genes},
            "message": f"{gse} 差异分析完成：{stats['n_probes_tested']} 探针，"
                       f"显著 {stats['n_sig_fdr05']}（FDR<0.05）",
            "artifacts": [
                {"name": "geo_report.html", "path": str(report), "kind": "html"},
                {"name": "de_table.csv", "path": str(csv_path), "kind": "csv"},
                {"name": "de_table.xlsx", "path": str(excel_path), "kind": "excel"},
                *[{"name": p.name, "path": str(p), "kind": "png"}
                  for p in sorted(outputs.glob("*.png"))],
            ]}
