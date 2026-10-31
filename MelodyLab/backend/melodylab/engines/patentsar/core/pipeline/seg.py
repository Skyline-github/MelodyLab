"""Locate chemical-structure drawings inside a full scanned patent page (OpenCV).

A scanned page is one big image mixing text lines and molecule drawings. We
binarise to ink, dilate so that a molecule's atoms/bonds merge into a single
blob (while a text line merges only into a thin, wide strip), then keep blobs
that look like structures: tall enough and not text-line-shaped. Each kept blob
is cropped (with padding) to a PNG for DECIMER.

Prep pages ("Preparation of Example/Intermediate N:") typically yield one blob;
table pages can yield several. Header/footer bands are ignored.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import cv2
import numpy as np


@dataclass
class Region:
    bbox: tuple[int, int, int, int]
    crop_path: Path


def _binarise(img: np.ndarray) -> np.ndarray:
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY) if img.ndim == 3 else img
    # adaptive threshold -> ink = 255 on black = 0
    inv = 255 - gray
    _, binary = cv2.threshold(inv, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    return binary


def _split_components(binary, x0, y0, x1, y1, min_area, dil):
    """Split a wide blob (compound-table row: number|name|structure|m/z) into
    cells at vertical grid lines, then keep the molecule cell. A molecule is one
    large connected drawing, so its biggest connected component fills most of the
    cell; a text cell shatters into many small letter components."""
    sub = binary[y0:y1, x0:x1]
    if sub.size == 0:
        return []
    H, W = sub.shape[:2]
    col_ink = (sub > 0).mean(axis=0)
    line_cols = np.flatnonzero(col_ink > 0.5)
    bounds = [0]
    i = 0
    while i < len(line_cols):
        j = i
        while j + 1 < len(line_cols) and line_cols[j + 1] - line_cols[j] <= 2:
            j += 1
        mid = int((line_cols[i] + line_cols[j]) / 2)
        if mid - bounds[-1] > 40:
            bounds.append(mid)
        i = j + 1
    bounds.append(W)

    k = cv2.getStructuringElement(cv2.MORPH_RECT, (3, 3))
    boxes = []
    for a, b in zip(bounds, bounds[1:]):
        if b - a < 60:
            continue
        block = sub[:, a:b]
        rows = np.flatnonzero((block > 0).sum(axis=1) > 0)
        if rows.size == 0:
            continue
        ry0, ry1 = rows[0], rows[-1]
        cell = block[ry0:ry1 + 1]
        bh, bw = cell.shape[:2]
        if bw * bh < min_area or bh < 120 or not (0.4 <= bw / max(bh, 1) <= 2.2):
            continue
        m = cv2.dilate(cell, k, iterations=1)
        n, labels, stats, _ = cv2.connectedComponentsWithStats(m, connectivity=8)
        if n <= 1:
            continue
        largest = stats[1:, cv2.CC_STAT_AREA].max()
        ink = (cell > 0).sum()
        if ink and largest / max(ink, 1) > 0.5:  # one dominant connected drawing
            boxes.append((x0 + a, y0 + ry0, x0 + b, y0 + ry1 + 1))
    return boxes


def find_structure_regions(
    image_path: str | Path,
    out_dir: str | Path,
    min_height: int = 150,
    min_width: int = 120,
    min_area: int = 20000,
    max_aspect: float = 6.0,
    dilate: int = 9,
    pad: int = 12,
    band: float = 0.10,
    max_regions: int = 12,
) -> list[Region]:
    image_path = Path(image_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    img = cv2.imread(str(image_path))
    if img is None:
        return []
    h, w = img.shape[:2]
    binary = _binarise(img)
    # blank the header/footer bands so running numbers / page no. aren't picked
    binary[: int(h * band), :] = 0
    binary[int(h * (1 - band)):, :] = 0

    kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (dilate, dilate))
    merged = cv2.dilate(binary, kernel, iterations=1)
    contours, _ = cv2.findContours(merged, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    boxes: list[tuple[int, int, int, int]] = []
    for c in contours:
        x, y, bw, bh = cv2.boundingRect(c)
        area = bw * bh
        aspect = bw / max(bh, 1)
        if bh < min_height or bw < min_width or area < min_area:
            continue
        if aspect > max_aspect or aspect < 1 / max_aspect:
            continue
        # a very wide blob is often a compound-table row (text+molecule+m/z) or
        # two molecules -> strip grid lines and split; use the parts if found
        if aspect > 1.5:
            subs = _split_components(binary, x, y, x + bw, y + bh, min_area // 3, dil=5)
            if subs:
                boxes.extend(subs)
                continue
        boxes.append((x, y, x + bw, y + bh))

    regions: list[Region] = []
    for (x, y, bw_end, bh_end) in boxes:
        x0 = max(0, x - pad); y0 = max(0, y - pad)
        x1 = min(w, bw_end + pad); y1 = min(h, bh_end + pad)
        crop = img[y0:y1, x0:x1]
        out = out_dir / f"{image_path.stem}_r{len(regions)}.png"
        cv2.imwrite(str(out), crop)
        regions.append(Region((x0, y0, x1, y1), out))
        if len(regions) >= max_regions:
            break
    # top-to-bottom, left-to-right reading order
    regions.sort(key=lambda r: (r.bbox[1] // 100, r.bbox[0]))
    return regions


if __name__ == "__main__":
    import sys
    img = sys.argv[1] if len(sys.argv) > 1 else "output/work_scan/pages/page_247.png"
    regs = find_structure_regions(img, "output/probe/crops")
    for r in regs:
        print(r.bbox, "->", r.crop_path.name)
    print(f"{len(regs)} region(s)")
