"""Detect assay-table pages and their target (AMY3R vs CTR) in a scanned patent.

For each candidate page it OCRs (cached), counts parseable [number, EC50, Emax]
rows, and looks for a target keyword in the page text. Run AFTER the structure
batch (it spins up easyocr, so avoid concurrent heavy jobs). The reported pages
feed --assay-pages into pipeline.run_scanned to regenerate joined deliverables.

Usage:
  python -m pipeline.detect_assay --pdf ../Patent_test/WO_2025240893_A1.pdf \
      --range 378:395
"""
from __future__ import annotations

import argparse
import re

from . import ocr_text, run_scanned

KEY = re.compile(r"(AMY3R|AMYLIN|CTR|CALCITONIN|CGRP|GIP)", re.IGNORECASE)


def detect(pdf, rng, workdir, zoom=3.0):
    doc_pages = range(rng[0], rng[1])
    for p in doc_pages:
        rec = ocr_text.get_page(pdf, p, workdir, zoom=zoom)
        rows = run_scanned.parse_scanned_table(rec, "x")
        kws = sorted({m.group(1).upper() for m in KEY.finditer(rec["text"])})
        if rows or kws:
            print(f"p{p}: rows={len(rows):>3} keywords={kws} "
                  f"| first={list(rows)[:4]}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdf", required=True)
    ap.add_argument("--range", required=True, help="start:end")
    ap.add_argument("--workdir", default="output_scan/work")
    a = ap.parse_args()
    s, e = (int(x) for x in a.range.split(":"))
    detect(a.pdf, (s, e), a.workdir)


if __name__ == "__main__":
    main()
