"""Cheap page-type locator for scanned patents (no OCR / no DECIMER).

Renders each page at low resolution and uses OpenCV heuristics to flag:
  * structure pages  -> contain a large ink blob (a chemical drawing) that is
    tall relative to the page (text lines are short)
  * assay-table pages -> contain long horizontal grid lines (a bordered table)
Used to pick --structure-pages / --assay-pages ranges for run_scanned without
OCR'ing the whole document first.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import cv2
import numpy as np
import pymupdf


def _page_arrays(doc, idx, zoom):
    pix = doc[idx].get_pixmap(matrix=pymupdf.Matrix(zoom, zoom))
    arr = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
    img = cv2.cvtColor(arr[..., :3], cv2.COLOR_RGB2BGR)
    return img


def _has_structure(img) -> bool:
    h, w = img.shape[:2]
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    inv = 255 - gray
    _, b = cv2.threshold(inv, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    b[: int(h * 0.10), :] = 0
    b[int(h * 0.90):, :] = 0
    k = cv2.getStructuringElement(cv2.MORPH_RECT, (7, 7))
    merged = cv2.dilate(b, k, iterations=1)
    contours, _ = cv2.findContours(merged, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    for c in contours:
        x, y, bw, bh = cv2.boundingRect(c)
        if bh >= 0.15 * h and bw >= 0.12 * w and 0.25 <= bw / bh <= 4.0:
            # 排除文本段落块：整栏文字经膨胀也会连成高大 blob（WO 正文页
            # 185/248 全误报的教训）。化学结构图无文本行纹理——blob 内水平
            # 投影的行峰间距（pitch）远大于正文行距；双栏正文 pitch≈0.024h，
            # 结构图 ≥0.1h。取 0.05h 为界。
            proj = b[y:y + bh, x:x + bw].sum(axis=1) / 255
            th = max(1.0, float(proj.max()) * 0.3)
            rows = proj > th
            peaks, prev = 0, False
            for v in rows:
                if v and not prev:
                    peaks += 1
                prev = v
            pitch = bh / peaks if peaks else bh
            if pitch >= 0.05 * h:
                return True
    return False


def _has_table(img) -> bool:
    h, w = img.shape[:2]
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    binary = cv2.adaptiveThreshold(~gray, 255, cv2.ADAPTIVE_THRESH_MEAN_C,
                                   cv2.THRESH_BINARY, 15, -2)
    hk = cv2.getStructuringElement(cv2.MORPH_RECT, (max(20, w // 8), 1))
    vk = cv2.getStructuringElement(cv2.MORPH_RECT, (1, max(20, h // 8)))
    horiz = cv2.erode(cv2.dilate(binary, hk), hk)
    vert = cv2.erode(cv2.dilate(binary, vk), vk)
    # a bordered table has >=3 long horizontal rules and >=2 long vertical rules
    n_h = _count_lines(horiz, axis=0, min_frac=0.4, h=h)
    n_v = _count_lines(vert, axis=1, min_frac=0.4, h=w)
    return n_h >= 3 and n_v >= 2


def _count_lines(mask, axis, min_frac, h):
    proj = mask.sum(axis=axis)
    thresh = min_frac * h * 255
    lines = 0
    prev = False
    for v in proj:
        cur = v > thresh
        if cur and not prev:
            lines += 1
        prev = cur
    return lines


def scan_pdf(pdf: str, zoom: float = 0.8) -> dict:
    doc = pymupdf.open(pdf)
    n = doc.page_count
    struct_pages, table_pages = [], []
    for idx in range(n):
        img = _page_arrays(doc, idx, zoom)
        if _has_table(img):
            table_pages.append(idx)
        elif _has_structure(img):
            struct_pages.append(idx)
    doc.close()
    return {"pdf": pdf, "pages": n, "structure_pages": struct_pages,
            "table_pages": table_pages}


def _ranges(pages: list[int], gap: int = 3) -> list[tuple[int, int]]:
    if not pages:
        return []
    out, start, prev = [], pages[0], pages[0]
    for p in pages[1:]:
        if p - prev > gap:
            out.append((start, prev + 1))
            start = p
        prev = p
    out.append((start, prev + 1))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdfs", nargs="+")
    ap.add_argument("--zoom", type=float, default=0.8)
    a = ap.parse_args()
    for pdf in a.pdfs:
        r = scan_pdf(pdf, a.zoom)
        doc = pymupdf.open(pdf); n = doc.page_count; doc.close()
        sr = _ranges(r["structure_pages"])
        print(f"\n=== {Path(pdf).name} ({n} pages) ===")
        print(f"structure pages: {len(r['structure_pages'])} | ranges(start:end): {sr}")
        print(f"table pages: {len(r['table_pages'])} -> {r['table_pages'][:40]}"
              f"{' ...' if len(r['table_pages']) > 40 else ''}")


if __name__ == "__main__":
    main()
