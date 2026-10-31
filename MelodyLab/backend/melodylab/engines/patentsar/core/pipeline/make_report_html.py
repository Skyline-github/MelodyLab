"""Build a self-contained, shareable CARD-style HTML report for one mined patent.

Embeds SAR figures + Free-Wilson contributions + every compound as a card
(structure crop as base64, SMILES, activity badges, QC flag ribbon). One
portable .html file.

Usage:
  python -m pipeline.make_report_html --patent-dir output_scan/WO_2025240893_A1
"""
from __future__ import annotations

import argparse
import base64
import csv
import html
import math
import re
from pathlib import Path


def _b64(path: Path) -> str:
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode() if path and path.exists() else ""


def _pact(v):
    try:
        x = float(v)
        return -math.log10(x / 1e9) if x > 0 else None  # nM -> M -> pEC50
    except (TypeError, ValueError):
        return None


FW_BULLET = re.compile(r"-\s+(R\d+)\s+`([^`]+)`\s*:\s*([+\-0-9.]+)")


def parse_sar(sar_md: Path):
    """Return (fw_stats dict, contributions list[(pos,r,coef)], fw_r2, fw_q2)."""
    stats, contribs = {}, []
    if not sar_md.exists():
        return stats, contribs
    txt = sar_md.read_text()
    m = re.search(r"\|\s*Terms\s*\|.*?\n\|[-| ]+\n\|\s*(\d+)\s*\|\s*([\d.]+)\s*\|\s*(-?[\d.]+)\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)", txt)
    if m:
        stats = {"terms": m.group(1), "R2": m.group(2), "Q2": m.group(3),
                 "RMSE_fit": m.group(4), "RMSE_loo": m.group(5)}
    for mm in FW_BULLET.finditer(txt):
        contribs.append((mm.group(1), mm.group(2), float(mm.group(3))))
    contribs.sort(key=lambda t: -t[2])
    return stats, contribs


def build(patent_dir: Path) -> Path:
    base = patent_dir.name
    rows = list(csv.DictReader(open(patent_dir / f"{base}_dataset.csv")))
    sar_dir = patent_dir / "sar"
    struct_dir = patent_dir / "structures"
    img_by_label: dict[str, Path] = {}
    if struct_dir.exists():
        for p in sorted(struct_dir.glob("*.png")):
            parts = p.name.split("_")
            if len(parts) >= 2:
                img_by_label.setdefault(parts[1], p)

    n = len(rows)
    valid = [r for r in rows if r["canonical_smiles"]]
    amy = [r for r in rows if _pact(r["ec50_amy3r"]) is not None]
    ctr = [r for r in rows if _pact(r["ec50_ctr"]) is not None]
    flagged = [r for r in rows if r["flags"]]
    stats, contribs = parse_sar(sar_dir / "sar_report.md")

    # ---- SAR figures (core + boxplots + structural heatmap + 3D-QSAR) ----
    figs = ""
    core_p = sar_dir / "core_scaffold.png"
    if core_p.exists():
        figs += (f"<figure class='core'><img src='{_b64(core_p)}' width='560'>"
                 f"<figcaption>RGA 母核（FMCS 公共骨架，R1–Rn 为可变取代位点）</figcaption></figure>")
    for p in sorted(sar_dir.glob("rga_box_*.png"))[:3]:
        figs += (f"<figure><img src='{_b64(p)}' width='560'>"
                 f"<figcaption>{html.escape(p.stem.replace('rga_box_',''))} 位点箱线图（按 pEC50 降序）</figcaption></figure>")
    for p in sorted(sar_dir.glob("rga_heat_*.png"))[:2]:
        figs += (f"<figure><img src='{_b64(p)}' width='620'>"
                 f"<figcaption>{html.escape(p.stem.replace('rga_heat_','').replace('_','×'))} 组合热图（结构轴，格值=平均 pEC50）</figcaption></figure>")
    coeff = sar_dir / "qsar3d_coeff.png"
    if coeff.exists():
        figs += (f"<figure><img src='{_b64(coeff)}' width='480'>"
                 f"<figcaption>Field-based 3D-QSAR 系数场</figcaption></figure>")

    # ---- Free-Wilson top contributions with R-group structure thumbnails ----
    from ..sar import plots as _p
    tmp = sar_dir / "_contrib_imgs"; tmp.mkdir(exist_ok=True)

    def rgroup_thumb(rsmi):
        p = tmp / f"rg{abs(hash(rsmi)) % 10**9}.png"
        if not p.exists():
            if not _p._render_rgroup(rsmi, p, size=80):
                return ""
        return _b64(p)

    contrib_rows = ""
    for pos, rsmi, c in contribs[:16]:
        thumb = rgroup_thumb(rsmi)
        img = f"<img src='{thumb}' width='72'>" if thumb else html.escape(rsmi)
        contrib_rows += (f"<tr><td>{html.escape(pos)}</td>"
                         f"<td class='rgcell'>{img}</td>"
                         f"<td class='{'pos' if c>0 else 'neg'}'>{c:+.2f}</td></tr>")

    # ---- compound cards (valid, sorted by potency) ----
    def card(r):
        img = img_by_label.get(str(r["example_no"]))
        flag = r["flags"]
        cls = "flag" if flag else ""
        pec = _pact(r["ec50_amy3r"])
        b64 = _b64(img)
        img_html = f"<img src='{b64}'>" if b64 else "<span class='muted'>no image</span>"
        if r["ec50_amy3r"] not in ("", "ND"):
            badge = f"<span class='act'>{html.escape(r['ec50_amy3r'])} nM</span>"
        else:
            badge = "<span class='act nd'>ND</span>"
        ctrb = f"<span class='act ctr'>CTR {html.escape(r['ec50_ctr'])}</span>" if r["ec50_ctr"] not in ("", "ND") else ""
        pec_html = f"<span class='pec'>pEC50 {pec:.2f}</span>" if pec is not None else ""
        ribbon = f"<span class='ribbon'>{html.escape(flag)}</span>" if flag else ""
        return (f"<div class='cmp {cls}'>"
                f"<div class='cmpimg'>{img_html}</div>"
                f"<div class='cmpbody'>"
                f"<div class='cmpid'>Example {html.escape(str(r['example_no']))}{ribbon}</div>"
                f"<div class='smi'><code>{html.escape(r['canonical_smiles'] or '')}</code></div>"
                f"<div class='badges'>{badge}{ctrb}{pec_html}</div>"
                f"<div class='meta'>{html.escape(r['formula'] or '')} · MW {r['mw_avg'] or ''}"
                f" · m/z {r['mz_observed'] or '—'}/{r['mz_calc'] or '—'} · OCR {r['ocr_confidence'] or ''}</div>"
                f"</div></div>")

    cards_sorted = sorted(valid, key=lambda r: (-( _pact(r["ec50_amy3r"]) or -99)))
    cards_html = "".join(card(r) for r in cards_sorted)

    doc = f"""<!doctype html><html lang="zh"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(base)} — 专利挖掘 & SAR</title>
<style>
:root{{--brand:#0b3d5c;--bg:#f4f6f8}}
body{{font-family:-apple-system,Segoe UI,Roboto,"Microsoft YaHei",sans-serif;margin:0;background:var(--bg);color:#1b1b1b}}
header{{background:linear-gradient(120deg,#0b3d5c,#12628f);color:#fff;padding:28px 40px}}
header h1{{margin:0;font-size:22px}} header p{{margin:4px 0 0;opacity:.9;font-size:13px}}
main{{max-width:1180px;margin:0 auto;padding:22px 28px}}
section{{background:#fff;border:1px solid #e3e6e9;border-radius:12px;padding:20px 24px;margin:18px 0;box-shadow:0 1px 3px rgba(0,0,0,.04)}}
h2{{font-size:16px;margin:0 0 14px;color:var(--brand);border-left:4px solid var(--brand);padding-left:10px}}
.cards{{display:flex;flex-wrap:wrap;gap:12px}}
.metric{{flex:1 1 150px;background:#eef5f9;border-radius:10px;padding:16px;text-align:center}}
.metric b{{display:block;font-size:28px;color:var(--brand)}} .metric span{{font-size:12px;color:#555}}
.figrow{{display:flex;flex-wrap:wrap;gap:18px;align-items:flex-start}}
figure{{margin:0}} figcaption{{font-size:12px;color:#555;margin-top:4px;text-align:center}}
table.fw{{border-collapse:collapse;font-size:12px;width:100%}}
table.fw th,table.fw td{{border-bottom:1px solid #eee;padding:4px 8px;text-align:left;vertical-align:middle}}
.rgcell img{{display:block}}
.pos{{color:#0a7d3f;font-weight:600}} .neg{{color:#b00}}
.grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(300px,1fr));gap:14px}}
.cmp{{display:flex;gap:10px;background:#fff;border:1px solid #e3e6e9;border-radius:10px;padding:10px;align-items:stretch}}
.cmp.flag{{border-color:#e0a;background:#fff5f8}}
.cmpimg{{flex:0 0 110px;display:flex;align-items:center;justify-content:center;background:#fafafa;border-radius:6px}}
.cmpimg img{{max-width:108px;max-height:130px}}
.cmpbody{{flex:1;min-width:0}}
.cmpid{{font-weight:700;font-size:13px;color:var(--brand)}}
.ribbon{{float:right;font-size:10px;background:#e0a;color:#fff;border-radius:4px;padding:1px 5px}}
.smi{{margin:4px 0}} .smi code{{font-size:10px;word-break:break-all;color:#333}}
.badges{{margin:3px 0}} .act{{display:inline-block;font-size:11px;background:#0a7d3f;color:#fff;border-radius:4px;padding:1px 6px;margin-right:5px}}
.act.ctr{{background:#7a4fb0}} .act.nd{{background:#999}} .pec{{font-size:11px;color:#0a7d3f;font-weight:600}}
.meta{{font-size:10px;color:#777}} .muted{{color:#bbb;font-size:11px}}
.note{{font-size:13px;color:#444}} ul{{margin:6px 0}}
</style></head><body>
<header><h1>{html.escape(base)} — AI-Agent 专利挖掘 → SAR 汇总</h1>
<p>Lilly DACRA (AMY3R / CTR) · WIPO 扫描版 A1 · 开源管线 (easyocr + OpenCV + DECIMER + RDKit + sklearn) · 自包含单文件</p></header>
<main>
<section><h2>1 · 关键指标</h2><div class="cards">
<div class="metric"><b>{n}</b><span>提取结构</span></div>
<div class="metric"><b>{100*len(valid)//n}%</b><span>有效 SMILES</span></div>
<div class="metric"><b>{len(amy)}</b><span>AMY3R 活性</span></div>
<div class="metric"><b>{len(ctr)}</b><span>CTR 活性</span></div>
<div class="metric"><b>{len(flagged)}</b><span>QC 标记待复核</span></div>
</div></section>

<section><h2>2 · SAR 结果（母核=马库什核心，RGA 7 位点）</h2>
<p class="note"><b>Free-Wilson 加和 QSAR</b>：Terms {stats.get('terms','—')} · <b>R²={stats.get('R2','—')}</b> · <b>Q²(LOO)={stats.get('Q2','—')}</b> · RMSE_loo={stats.get('RMSE_loo','—')}（n={len(amy)} 真实化合物）</p>
<div class="figrow">{figs}</div>
<h3 style="font-size:14px;margin:16px 0 6px">Top R-基团贡献（pEC50，正=提升活性）</h3>
<table class="fw"><tr><th>位点</th><th>R 基团</th><th>贡献</th></tr>{contrib_rows}</table>
</section>

<section><h2>3 · 化合物卡片（{len(valid)} 个，按 AMY3R 活性降序；粉框=QC 待复核）</h2>
<div class="grid">{cards_html}</div></section>

<section><h2>4 · 方法与边界</h2><ul class="note">
<li>整页 easyocr(zoom 3.0) → OpenCV 分子切图(含表格行按网格线拆分) → DECIMER 图→SMILES → RDKit 质控(MW/LCMS [M+H]+ 交叉校验+去重) → 按 Example 号 join 活性表 → RGA/Free-Wilson/3D-QSAR。</li>
<li>QC 标记 {len(flagged)} 行（INVALID / LCMS_MISMATCH / DUPLICATE）以粉框+红标呈现，供人工在 Flare 内核对修正。</li>
<li>3D-QSAR 为开源近似（单构象+粗网格），精度低于 Flare 力场/FEP；异构系列需指定母核才能得干净 RGA。</li>
<li>结果未经人工核验。</li>
</ul></section>
</main></body></html>"""

    out = patent_dir / f"{base}_REPORT.html"
    out.write_text(doc)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patent-dir", required=True)
    a = ap.parse_args()
    out = build(Path(a.patent_dir))
    print(f"wrote {out} ({out.stat().st_size//1024} KB)")


if __name__ == "__main__":
    main()
