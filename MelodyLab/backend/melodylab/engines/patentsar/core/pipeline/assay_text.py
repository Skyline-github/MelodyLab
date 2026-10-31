"""Parse AMY3R / CTR activity tables and LCMS data from patent text.

Heuristic but tolerant parser for the layout the blog describes: an Example
number followed by EC50 / Emax columns (numeric, graded letter, or ``ND``),
and ``Example n: LCMS m/z X`` lines. Text is taken from the PDF so no external
OCR of the tables is needed.
"""
from __future__ import annotations

import re

NUM_RE = re.compile(r"^-?\d+(\.\d+)?$")
GRADE_RE = re.compile(r"^[A-Za-z]$")
LCMS_RE = re.compile(r"[Ee]xample[\s.\-]*([A-Za-z0-9\-]+).*?m/?z\s*=?\s*([0-9]+\.?[0-9]*)",
                     re.DOTALL)
COL_RE = re.compile(
    r"^\s*(?:[Ee]xample[\s.\-]*)?([A-Za-z0-9\-]+)[ \t]+"
    r"(\S+)[ \t]+(\S+)[ \t]+(\S+)[ \t]+(\S+)\s*$")


def _is_value(tok: str) -> bool:
    return bool(NUM_RE.match(tok) or GRADE_RE.match(tok) or tok.lower() == "nd")


def _norm(tok: str):
    """Return None for ND, else the raw string (kept so grades survive)."""
    return None if tok.lower() == "nd" else tok


def parse_assay_table(text: str) -> dict[str, dict]:
    """Map Example -> {ec50_amy3r, emax_amy3r, ec50_ctr, emax_ctr}.

    Column order assumed (per prompt): AMY3R EC50, AMY3R Emax, CTR EC50, CTR Emax.
    """
    rows: dict[str, dict] = {}
    for line in text.splitlines():
        m = COL_RE.search(line)
        if not m:
            continue
        ex, c1, c2, c3, c4 = m.groups()
        # skip header/garbage rows where the four 'values' aren't plausible
        if not all(_is_value(t) for t in (c1, c2, c3, c4)):
            continue
        rows[ex] = {
            "ec50_amy3r": _norm(c1), "emax_amy3r": _norm(c2),
            "ec50_ctr": _norm(c3), "emax_ctr": _norm(c4),
            "raw": line.strip(),
        }
    return rows


def parse_lcms(text: str) -> dict[str, float]:
    out: dict[str, float] = {}
    for m in LCMS_RE.finditer(text):
        ex, mz = m.group(1), m.group(2)
        try:
            out[ex] = float(mz)
        except ValueError:
            pass
    return out


if __name__ == "__main__":
    import pymupdf, sys
    doc = pymupdf.open(sys.argv[1] if len(sys.argv) > 1 else "data/fixture_patent.pdf")
    assay_txt = doc[6].get_text()
    lcms_txt = doc[7].get_text()
    print("assay:", parse_assay_table(assay_txt))
    print("lcms:", parse_lcms(lcms_txt))
