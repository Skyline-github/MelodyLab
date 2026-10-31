"""Scanned-patent orchestration: page -> OCR text + structure crops -> dataset.

This is the real-WO-A1 path (no text layer). For each page in the structure
range we: OCR the page (cached), segment & crop molecule drawings, run DECIMER
on each crop, and pull the Example/Intermediate label and LCMS m/z from the OCR
text. Separately, assay-table pages are parsed from OCR word boxes into an
activity map. Everything is joined by Example number, RDKit-QC'd, and written
with the same deliverables as the born-digital path (pipeline.report).

Usage:
  python -m pipeline.run_scanned --pdf ../Patent_test/WO_2025240893_A1.pdf \
      --structure-pages 247:260 --kinds example --outdir output_scan
"""
from __future__ import annotations

import argparse
import re
import gc
import shutil
from pathlib import Path

from . import assay_text, chem_qc, markush, ocr, ocr_text, report, seg, vlm

# "Preparation of Example 214" / "Example 214:" / "Example 225a"; "ref" covers
# reference-compound captions (fixture convention, kept out of numeric ids)
EX_NUM = r"[0-9]{1,4}[a-z]?|ref"
EXAMPLE_RE = re.compile(rf"\bExamples?\s*[:\-]?\s*({EX_NUM})\b")
INTERMEDIATE_RE = re.compile(rf"\bIntermediates?\s*[:\-]?\s*({EX_NUM})\b")
MZ_RE = re.compile(r"m/?z\s*=?\s*([0-9]{2,4}\.?[0-9]{0,2})", re.IGNORECASE)


def parse_label(text: str) -> tuple[str | None, str | None]:
    """Return (kind, number). Prefer an Example label over Intermediate."""
    m = EXAMPLE_RE.search(text)
    if m:
        return "example", m.group(1)
    m = INTERMEDIATE_RE.search(text)
    if m:
        return "intermediate", m.group(1)
    return None, None


def parse_mz(text: str) -> float | None:
    m = MZ_RE.search(text)
    return float(m.group(1)) if m else None


_LABEL_LINE = re.compile(rf"[Ee]xamples?\s*[:\-]?\s*({EX_NUM})\b")
_INTERM_LINE = re.compile(rf"[Ii]ntermediates?\s*[:\-]?\s*({EX_NUM})\b")
_MZ_LINE = re.compile(r"m/?z\s*=?\s*([0-9]{2,4}\.?[0-9]{0,2})", re.IGNORECASE)


def _label_positions(rec: dict) -> list[dict]:
    """{y, cx, kind, number} for each line carrying an Example/Intermediate label.

    easyocr returns line-level text boxes, so we regex inside each line and use
    the line box top as the label's vertical position; cx (line center x) is
    kept for below-caption matching.
    """
    out = []
    for w in rec["words"]:
        text = w["text"]
        m = _INTERM_LINE.search(text)
        # prefer whichever label type appears first in the line
        e = _LABEL_LINE.search(text)
        if e and (not m or e.start() <= m.start()):
            kind, num, mm = "example", e.group(1), e
        elif m:
            kind, num, mm = "intermediate", m.group(1), m
        else:
            continue
        x0, _, x1, _ = w["bbox"]
        out.append({"y": w["bbox"][1], "cx": (x0 + x1) / 2, "kind": kind,
                    "number": num})
    return out


def _mz_positions(rec: dict) -> list[tuple[int, float]]:
    """(y, m/z value) for each line carrying an LCMS 'm/z NNN'."""
    out = []
    for w in rec["words"]:
        for m in _MZ_LINE.finditer(w["text"]):
            try:
                out.append((w["bbox"][1], float(m.group(1))))
            except ValueError:
                pass
    return out


def _assign_label(top_y: int, bot_y: int, cx: float, labels) -> tuple[str | None, str | None]:
    """Label for a structure region, mirroring pdf_images' heuristic:
    1. caption directly BELOW the region within a horizontal window (patent
       convention: 'Example N' under each drawing) - nearest below wins;
    2. else nearest heading ABOVE the region top (WO convention:
       'Preparation of Example N' above the structures);
    3. else the globally nearest label."""
    if not labels:
        return None, None
    below = [l for l in labels if l["y"] >= bot_y - 20 and abs(l["cx"] - cx) < 350]
    if below:
        best = min(below, key=lambda l: l["y"])
        return best["kind"], best["number"]
    above = [l for l in labels if l["y"] <= top_y + 30]
    if above:
        best = max(above, key=lambda l: l["y"])
        return best["kind"], best["number"]
    best = min(labels, key=lambda l: abs(l["y"] - top_y))
    return best["kind"], best["number"]


def _assign_mz(top_y: int, bot_y: int, mzs):
    """Nearest m/z below the region (LCMS text follows the structure)."""
    below = [m for m in mzs if m[0] >= bot_y - 20]
    if below:
        return min(below, key=lambda m: m[0])[1]
    if mzs:
        return min(mzs, key=lambda m: abs(m[0] - bot_y))[1]
    return None



def parse_scanned_table(rec: dict, targets: list[str]) -> dict[str, dict]:
    """Cluster OCR words into rows; read id + 2 values per target.

    Column layout per target: EC50 then Emax (numeric, graded letter, or ND),
    taken positionally after the Example id — so a two-target table
    (AMY3R + CTR) with 4 value columns is parsed fully instead of dropping
    columns 3-4."""
    words = rec["words"]
    rows: list[list] = []
    for w in sorted(words, key=lambda w: (w["bbox"][1], w["bbox"][0])):
        cy = (w["bbox"][1] + w["bbox"][3]) / 2
        if rows and abs(cy - rows[-1][0]) < 22:
            rows[-1][1].append(w)
        else:
            rows.append([cy, [w]])
    n_vals = 2 * len(targets)
    out: dict[str, dict] = {}
    for _, ws in rows:
        ws.sort(key=lambda w: w["bbox"][0])
        toks = [w["text"].replace(",", "").strip() for w in ws]
        toks = [t for t in toks if t]
        if not toks:
            continue
        idm = re.fullmatch(rf"({EX_NUM})", toks[0])
        if not idm or len(toks) < n_vals + 1:
            continue
        vals = toks[1:1 + n_vals]
        if not all(assay_text.NUM_RE.match(v) or assay_text.GRADE_RE.match(v)
                   or v.lower() == "nd" for v in vals):
            continue
        rec_vals = {}
        for i, t in enumerate(targets):
            e, m = vals[2 * i], vals[2 * i + 1]
            rec_vals[f"ec50_{t}"] = None if e.lower() == "nd" else e
            rec_vals[f"emax_{t}"] = None if m.lower() == "nd" else m
        out[idm.group(1)] = rec_vals
    return out


def process_structure_page(pdf, idx, workdir, zoom=3.0, hand_drawn=False,
                           vlm_cfg: dict | None = None):
    rec = ocr_text.get_page(pdf, idx, workdir, zoom=zoom)
    labels = _label_positions(rec)
    mzs = _mz_positions(rec)
    regions = seg.find_structure_regions(rec["image"], Path(workdir) / "crops")
    ocr_res = ocr.ocr_images([r.crop_path for r in regions],
                             cache_json=Path(workdir) / "decimer_cache.json",
                             confidence=True, hand_drawn=hand_drawn)
    out = []
    for r in regions:
        cx = (r.bbox[0] + r.bbox[2]) / 2
        kind, num = _assign_label(r.bbox[1], r.bbox[3], cx, labels)
        mz = _assign_mz(r.bbox[1], r.bbox[3], mzs)
        res = ocr_res[str(r.crop_path)]
        conf = _mean_conf(res.get("confidence"))
        out.append({"page": idx, "kind": kind, "number": num, "region_bbox": r.bbox,
                    "crop": r.crop_path, "smiles": res.get("smiles"),
                    "mz_observed": mz, "ocr_confidence": conf})
    # VLM disambiguation: fix OCR-garbled labels / one-heading-two-structures
    # / multi-molecule crops. Silent no-op without an API key.
    if vlm_cfg is not None and out:
        assigns = vlm.disambiguate_page(
            rec["image"],
            [{"crop_path": r.crop_path, "label": s["number"]}
             for r, s in zip(regions, out)],
            vlm_cfg)
        if assigns:
            for s, a in zip(out, assigns):
                if a["kind"] == "none":
                    continue
                if a["example_no"] and a["example_no"] != s["number"]:
                    s["number"] = a["example_no"]
                    s["vlm_override"] = True
                if a["kind"] != s["kind"]:
                    s["kind"] = a["kind"]
                    s["vlm_override"] = True
                if a["multi_molecule"]:
                    s["vlm_multi"] = True
    return out, rec


def _mean_conf(conf):
    if not conf:
        return None
    try:
        vals = [float(p) for (_, p) in conf] if isinstance(conf[0], (list, tuple)) else [float(conf)]
    except Exception:
        return None
    return round(sum(vals) / len(vals), 4) if vals else None


def run_scanned(pdf, structure_pages, kinds, assay_pages, outdir, workdir,
                zoom=3.0, text_zoom_mult: float = 2.0, mz_tol=1.5,
                hand_drawn=False, vlm_cfg: dict | None = None,
                mz_pages: list[int] | None = None,
                markush_config: str | None = None,
                markush_pages: list[int] | None = None,
                auto_locate: bool = True, locate_progress=None) -> dict:
    outdir = Path(outdir); workdir = Path(workdir)
    outdir.mkdir(parents=True, exist_ok=True); workdir.mkdir(parents=True, exist_ok=True)
    import pymupdf
    doc = pymupdf.open(str(pdf))
    total = doc.page_count
    doc.close()
    # text-heavy pages (tables, R-definitions) need higher dpi than structure
    # pages: at zoom 3.0 a lone "1" id fuses with the next column's "12"
    tzoom = zoom * text_zoom_mult

    # ---- 两遍法：范围是全文档默认值时，先用廉价启发式定位结构页 ----
    # （locate 渲染低分辨率图 + OpenCV 形态学，不跑 OCR/DECIMER，百页秒级；
    #   全文档盲目精确 OCR 每页要数秒~数十秒，绝大多数页纯浪费）
    # 用户显式指定了小范围（<60 页）时不定位、照旧全扫——那是人工兜底通道。
    start, end = structure_pages
    candidate_pages = None
    if auto_locate and (end - start) >= 60:
        from .locate import scan_pdf
        located = scan_pdf(str(pdf))
        candidate_pages = [p for p in located["structure_pages"] if start <= p < end]
        if locate_progress is not None:
            locate_progress(len(candidate_pages), total)
        if not candidate_pages:
            # 假阴性防护：启发式一页没定位到（小图/非常规排版都可能）就回退
            # 全范围精确扫——宁可慢，不可静默交出空数据集
            candidate_pages = None

    # activity map from assay pages
    activity: dict[str, dict] = {}
    for spec in assay_pages:  # (page_idx, [targets]) - e.g. 6:amy3r+ctr
        p, targets = spec
        rec = ocr_text.get_page(pdf, p, workdir, zoom=tzoom)
        for k, v in parse_scanned_table(rec, targets).items():
            activity.setdefault(k, {}).update(v)

    # optional dedicated LCMS pages: "Example N: LCMS m/z X" -> {N: X}
    mz_map: dict[str, float] = {}
    for p in (mz_pages or []):
        if p >= total:
            continue
        rec = ocr_text.get_page(pdf, p, workdir, zoom=tzoom)
        for w in rec["words"]:
            me = re.search(rf"[Ee]xamples?\s*[:\-]?\s*({EX_NUM})\b", w["text"])
            mm = MZ_RE.search(w["text"])
            if me and mm:
                try:
                    mz_map[me.group(1)] = float(mm.group(1))
                except ValueError:
                    pass

    # structures
    img_dir = outdir / "structures"; img_dir.mkdir(exist_ok=True)
    rows: list[dict] = []
    if candidate_pages is not None:
        page_iter = iter(candidate_pages)  # 定位后的候选页（可能为空=整本无结构图）
    else:
        page_iter = iter(range(max(0, start), min(end, total)))
    for idx in page_iter:
        structs, _ = process_structure_page(pdf, idx, workdir, zoom=zoom,
                                            hand_drawn=hand_drawn, vlm_cfg=vlm_cfg)
        for ri, s in enumerate(structs):
            if s["kind"] is None:
                # 无标签结构（中文专利/Markush 通式页无 "Example N" 标题）：
                # 标记为 unlabeled，example_no 用页号+区域内序号
                s["kind"] = "unlabeled"
                if not s["number"]:
                    s["number"] = f"p{idx}_r{ri}"
            if kinds and s["kind"] not in kinds and s["kind"] != "unlabeled":
                continue
            mzv = s["mz_observed"]
            if mzv is None and s["number"]:
                mzv = mz_map.get(s["number"])
            qc = chem_qc.validate(s["smiles"], observed_mz=mzv, mz_tol=mz_tol)
            if s["kind"] == "unlabeled":
                qc.setdefault("flags", []).append("UNLABELED_STRUCTURE")
            if s.get("vlm_multi"):
                qc.setdefault("flags", []).append("VLM_MULTI_MOLECULE")
            dst = img_dir / f"p{s['page']}_{s['number']}_{Path(s['crop']).name}"
            shutil.copy(s["crop"], dst)
            rows.append({
                "example_no": s["number"],
                "kind": s["kind"],
                "image_rel": f"structures/{dst.name}",
                "structure_page": s["page"],
                "raw_smiles": s["smiles"],
                "qc": qc,
                "assay": activity.get(s["number"]) if s["number"] else None,
                "ocr_confidence": s["ocr_confidence"],
            })
        gc.collect()  # bound memory on large scanned docs
    chem_qc.annotate_duplicates(rows)
    # kinds 过滤语义保持：unlabeled 行只有在"过滤后一行不剩"时才保留
    # （正常专利不受影响；Markush/无标签专利保底出全量结构）
    if kinds and "unlabeled" not in kinds:
        kept = [r for r in rows if r["kind"] != "unlabeled"]
        if kept:
            rows = kept

    # Markush R-table-defined examples: enumerate concrete SMILES from the
    # explicitly configured core template + R definitions, full QC + join
    n_markush = 0
    if markush_config and markush_pages:
        cfg_mk = markush.load_config(markush_config)
        lines = markush.harvest_scanned_lines(pdf, [p for p in markush_pages
                                                    if p < total],
                                              workdir, zoom=tzoom)
        for ex, combo in markush.parse_combinations(lines).items():
            if ex in {str(r["example_no"]) for r in rows}:
                continue
            smi, err = markush.enumerate_smiles(cfg_mk, combo)
            if smi is None:
                print(f"  [markush] Example {ex} skipped: {err}")
                continue
            qc = chem_qc.validate(smi, observed_mz=mz_map.get(ex), mz_tol=mz_tol)
            qc.setdefault("flags", []).append("MARKUSH_ENUMERATED")
            rows.append({"example_no": ex, "kind": "markush_example",
                         "image_rel": None, "structure_page": None, "qc": qc,
                         "assay": activity.get(ex), "ocr_confidence": None})
            n_markush += 1

    base = Path(pdf).stem
    report.write_csv(rows, outdir / f"{base}_dataset.csv")
    report.write_sdf(rows, outdir / f"{base}_structures.sdf")
    report.write_html(rows, outdir / f"{base}_check.html", title=f"{base} (scanned) - check table")
    meta = {"pdf": str(pdf), "structure_pages": structure_pages, "assay_pages": assay_pages,
            "lcms_pages": "in-text", "n_structures": len(rows), "hand_drawn": hand_drawn,
            "auto_located": candidate_pages is not None,
            "located_structure_pages": candidate_pages if candidate_pages is not None else None}
    report.write_markdown_report(rows, meta, outdir / f"{base}_report.md")
    return {"rows": rows, "activity": activity, "outdir": outdir, "base": base,
            "auto_located": candidate_pages is not None,
            "n_candidates": len(candidate_pages) if candidate_pages is not None else None}


def _parse_range(spec):
    a, b = spec.split(":")
    return int(a), int(b)


def _parse_assay(spec):
    # "378:amy3r,385:ctr" or one page with both targets: "6:amy3r+ctr"
    out = []
    for part in filter(None, spec.split(",")):
        p, t = part.split(":")
        out.append((int(p), [x for x in t.split("+") if x]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdf", required=True)
    ap.add_argument("--structure-pages", default="247:255")
    ap.add_argument("--kinds", default="example", help="comma list: example,intermediate")
    ap.add_argument("--assay-pages", default="", help="e.g. 378:amy3r,385:ctr")
    ap.add_argument("--outdir", default="output_scan")
    ap.add_argument("--workdir", default="output_scan/work")
    ap.add_argument("--zoom", type=float, default=3.0)
    ap.add_argument("--text-zoom-mult", type=float, default=2.0,
                    help="zoom multiplier for text-heavy pages (assay / LC-MS / "
                         "Markush R-table); tables need higher dpi to keep "
                         "columns separated")
    ap.add_argument("--mz-tol", type=float, default=1.5)
    ap.add_argument("--hand-drawn", action="store_true")
    ap.add_argument("--vlm", choices=["auto", "off"], default="auto",
                    help="VLM label disambiguation: 'auto' activates when "
                         "SAR_LLM_API_KEY is set (model: SAR_VLM_MODEL, "
                         "default glm-4v-flash); 'off' keeps heuristic labels")
    ap.add_argument("--mz-pages", default="",
                    help="page indices holding LCMS lines, e.g. '7' (m/z read "
                         "from OCR text and joined by Example number)")
    ap.add_argument("--markush-config", default=None,
                    help="JSON config: core_template with {rN} placeholders + "
                         "per-position label->SMILES table (pipeline/markush.py)")
    ap.add_argument("--markush-pages", default="",
                    help="page indices holding the generic structure + R-table "
                         "text, e.g. '5'")
    a = ap.parse_args()
    kinds = [k.strip() for k in a.kinds.split(",") if k.strip()] or None
    vlm_cfg = vlm.default_config() if a.vlm == "auto" else None
    mz_pages = [int(x) for x in a.mz_pages.split(",") if x.strip()]
    mk_pages = [int(x) for x in a.markush_pages.split(",") if x.strip()]
    res = run_scanned(a.pdf, _parse_range(a.structure_pages), kinds,
                      _parse_assay(a.assay_pages), Path(a.outdir), Path(a.workdir),
                      zoom=a.zoom, text_zoom_mult=a.text_zoom_mult,
                      mz_tol=a.mz_tol, hand_drawn=a.hand_drawn,
                      vlm_cfg=vlm_cfg, mz_pages=mz_pages,
                      markush_config=a.markush_config, markush_pages=mk_pages)
    rows = res["rows"]
    print(f"structures: {len(rows)} | valid SMILES: {sum(1 for r in rows if r['qc']['valid'])} "
          f"| flagged: {sum(1 for r in rows if r['qc']['flags'])} | activity rows: {len(res['activity'])}")
    for r in rows:
        print(f"  p{r['structure_page']} {r['kind']} {r['example_no']}: "
              f"{(r['qc']['canonical_smiles'] or 'INVALID')[:48]} "
              f"flags={r['qc']['flags']}")


if __name__ == "__main__":
    main()
