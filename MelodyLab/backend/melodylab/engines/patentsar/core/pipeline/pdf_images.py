"""Locate and crop chemical-structure images from a patent PDF (PyMuPDF).

For every embedded image inside the requested page range we render just its
placed rectangle to a PNG (so vector/raster structures are captured at print
resolution) and attach the nearest ``Example <n>`` text label found on the same
page. Labels drive the later join with the activity table.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import pymupdf

EXAMPLE_RE = re.compile(r"[Ee]xample[\s.\-]*([A-Za-z0-9\-]+)")


@dataclass
class StructureImage:
    image_path: Path
    page_index: int
    label: str | None
    bbox: tuple[float, float, float, float]


def _iter_pages(doc: pymupdf.Document, page_range: tuple[int, int] | None):
    start, end = page_range if page_range else (0, doc.page_count)
    for idx in range(max(0, start), min(end, doc.page_count)):
        yield doc[idx]


def _page_labels(page: pymupdf.Page) -> list[tuple[str, pymupdf.Rect]]:
    """Return (label, rect) for every 'Example n' string on the page."""
    out = []
    for block in page.get_text("dict").get("blocks", []):
        for line in block.get("lines", []):
            text = "".join(s["text"] for s in line["spans"])
            m = EXAMPLE_RE.search(text)
            if m:
                out.append((m.group(1), pymupdf.Rect(line["bbox"])))
    return out


def _nearest_label(cx: float, cy: float, labels) -> str | None:
    if not labels:
        return None
    # prefer labels below-left of the image center (patent caption convention),
    # tie-broken by euclidean distance.
    def score(item):
        lab, rect = item
        lx = (rect.x0 + rect.x1) / 2
        ly = (rect.y0 + rect.y1) / 2
        dy = ly - cy
        below_bonus = -1000 if (dy >= 0 and abs(lx - cx) < 160) else 0
        return below_bonus + ((lx - cx) ** 2 + (ly - cy) ** 2) ** 0.5
    return min(labels, key=score)[0]


def extract_structure_images(
    pdf_path: str | Path,
    workdir: str | Path,
    page_range: tuple[int, int] | None = None,
    zoom: float = 2.0,
    min_area: float = 900.0,
) -> list[StructureImage]:
    pdf_path = Path(pdf_path)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    out_dir = workdir / "structures"
    out_dir.mkdir(parents=True, exist_ok=True)

    doc = pymupdf.open(pdf_path)
    records: list[StructureImage] = []
    for page in _iter_pages(doc, page_range):
        labels = _page_labels(page)
        for info in page.get_image_info(xrefs=False):
            rect = pymupdf.Rect(info["bbox"])
            if rect.get_area() < min_area:
                continue
            cx, cy = (rect.x0 + rect.x1) / 2, (rect.y0 + rect.y1) / 2
            label = _nearest_label(cx, cy, labels)
            pix = page.get_pixmap(matrix=pymupdf.Matrix(zoom, zoom), clip=rect)
            key = f"p{page.number}_{'x' if label is None else label}.png"
            path = out_dir / key
            pix.save(path)
            records.append(StructureImage(path, page.number, label, tuple(rect)))
    doc.close()
    return records


if __name__ == "__main__":
    import sys
    pdf = sys.argv[1] if len(sys.argv) > 1 else "data/fixture_patent.pdf"
    recs = extract_structure_images(pdf, "output/work", page_range=(1, 6))
    for r in recs:
        print(f"page={r.page_index} label={r.label} -> {r.image_path.name}")
    print(f"extracted {len(recs)} structure images")
