"""End-to-end CLI for stage III: patent PDF -> structure/activity dataset + QC.

Pipeline (mirrors the blog's AI-Agent workflow, built on open-source tools):
  1. extract structure images from a page range        (PyMuPDF)
  2. chemical OCR image -> SMILES                        (DECIMER, cached)
  3. parse AMY3R/CTR assay table + LCMS m/z from text    (heuristics)
  4. RDKit validation, MW / LCMS cross-check, dedup QC   (RDKit)
  5. join everything by Example number
  6. deliverables: SDF + CSV + HTML check table + MD report

Usage:
  python -m pipeline.run --pdf data/fixture_patent.pdf \
      --structure-pages 1:6 --assay-pages 6 --lcms-pages 7 --outdir output
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pymupdf

from . import assay_text, chem_qc, markush, ocr, pdf_images, report


def _parse_range(spec: str) -> tuple[int, int]:
    a, b = spec.split(":")
    return int(a), int(b)


def _page_list(spec: str) -> list[int]:
    return [int(x) for x in spec.replace(" ", "").split(",") if x != ""]


def _mean_conf(conf) -> float | None:
    if not conf:
        return None
    try:
        vals = [float(p) for (_, p) in conf] if isinstance(conf[0], (list, tuple)) else [float(conf)]
    except Exception:
        return None
    return round(sum(vals) / len(vals), 4) if vals else None


def run(pdf: str, structure_pages: tuple[int, int], assay_pages: list[int],
        lcms_pages: list[int], outdir: Path, workdir: Path,
        hand_drawn: bool = False, mz_tol: float = 1.5,
        markush_config: str | None = None,
        markush_pages: list[int] | None = None) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    workdir.mkdir(parents=True, exist_ok=True)

    # 1) structures -> images
    recs = pdf_images.extract_structure_images(pdf, workdir, page_range=structure_pages)
    # 2) OCR
    ocr_res = ocr.ocr_images([r.image_path for r in recs],
                             cache_json=workdir / "ocr_cache.json",
                             confidence=True, hand_drawn=hand_drawn)
    # 3) text tables
    doc = pymupdf.open(pdf)
    assay_txt = "\n".join(doc[p].get_text() for p in assay_pages if p < doc.page_count)
    lcms_txt = "\n".join(doc[p].get_text() for p in lcms_pages if p < doc.page_count)
    doc.close()
    assay_map = assay_text.parse_assay_table(assay_txt)
    lcms_map = assay_text.parse_lcms(lcms_txt)

    # 4) QC + join; copy images next to HTML for relative embedding
    img_dir = outdir / "structures"
    img_dir.mkdir(exist_ok=True)
    rows: list[dict] = []
    for r in recs:
        res = ocr_res[str(r.image_path)]
        smiles = res.get("smiles")
        observed = lcms_map.get(r.label) if r.label else None
        qc = chem_qc.validate(smiles, observed_mz=observed, mz_tol=mz_tol)
        dst = img_dir / r.image_path.name
        shutil.copy(r.image_path, dst)
        rows.append({
            "example_no": r.label or Path(r.image_path).stem,
            "image_rel": f"structures/{dst.name}",
            "structure_page": r.page_index,
            "qc": qc,
            "assay": assay_map.get(r.label) if r.label else None,
            "ocr_confidence": _mean_conf(res.get("confidence")),
        })
    chem_qc.annotate_duplicates(rows)

    # Markush R-table-defined examples (e.g. "Example 12: R1 = Et, R2 = H"):
    # enumerate concrete SMILES from the explicit core template + R table,
    # full QC incl. LC-MS cross-check, joined to the assay map like any row
    n_markush = 0
    if markush_config and markush_pages:
        mk_cfg = markush.load_config(markush_config)
        d2 = pymupdf.open(pdf)
        lines: list[str] = []
        for p in markush_pages:
            if p < d2.page_count:
                lines += d2[p].get_text("text").splitlines()
        d2.close()
        existing = {str(r["example_no"]) for r in rows}
        for ex, combo in markush.parse_combinations(lines).items():
            if ex in existing:
                continue
            smi, err = markush.enumerate_smiles(mk_cfg, combo)
            if smi is None:
                print(f"  [markush] Example {ex} skipped: {err}")
                continue
            qc = chem_qc.validate(smi, observed_mz=lcms_map.get(ex), mz_tol=mz_tol)
            qc.setdefault("flags", []).append("MARKUSH_ENUMERATED")
            rows.append({"example_no": ex, "kind": "markush_example",
                         "image_rel": None, "structure_page": None, "qc": qc,
                         "assay": assay_map.get(ex), "ocr_confidence": None})
            n_markush += 1

    # 6) deliverables
    base = Path(pdf).stem
    report.write_csv(rows, outdir / f"{base}_dataset.csv")
    report.write_sdf(rows, outdir / f"{base}_structures.sdf")
    report.write_html(rows, outdir / f"{base}_check.html",
                      title=f"{base} - structure/activity check table")
    meta = {"pdf": str(pdf), "structure_pages": structure_pages,
            "assay_pages": assay_pages, "lcms_pages": lcms_pages,
            "n_structures": len(recs), "hand_drawn": hand_drawn}
    report.write_markdown_report(rows, meta, outdir / f"{base}_report.md")
    return {"rows": rows, "meta": meta, "outdir": outdir, "base": base}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdf", default="data/fixture_patent.pdf")
    ap.add_argument("--structure-pages", default="1:6")
    ap.add_argument("--assay-pages", default="6")
    ap.add_argument("--lcms-pages", default="7")
    ap.add_argument("--outdir", default="output")
    ap.add_argument("--workdir", default="output/work")
    ap.add_argument("--hand-drawn", action="store_true")
    ap.add_argument("--mz-tol", type=float, default=1.5)
    ap.add_argument("--markush-config", default=None,
                    help="JSON config: core_template with {rN} placeholders + "
                         "per-position label->SMILES table (pipeline/markush.py)")
    ap.add_argument("--markush-pages", default="",
                    help="page indices holding the generic structure + R-table "
                         "text, e.g. '5'")
    a = ap.parse_args()

    summary = run(a.pdf, _parse_range(a.structure_pages), _page_list(a.assay_pages),
                  _page_list(a.lcms_pages), Path(a.outdir), Path(a.workdir),
                  hand_drawn=a.hand_drawn, mz_tol=a.mz_tol,
                  markush_config=a.markush_config,
                  markush_pages=_page_list(a.markush_pages))
    rows = summary["rows"]
    flagged = [r for r in rows if r["qc"]["flags"]]
    print(f"structures: {len(rows)} | valid SMILES: {sum(1 for r in rows if r['qc']['valid'])} "
          f"| flagged: {len(flagged)}")
    for r in flagged:
        print(f"  - Example {r['example_no']}: {r['qc']['flags']}")
    print(f"outputs -> {summary['outdir']}/{summary['base']}_dataset.csv"
          f", _structures.sdf, _check.html, _report.md")


if __name__ == "__main__":
    main()
