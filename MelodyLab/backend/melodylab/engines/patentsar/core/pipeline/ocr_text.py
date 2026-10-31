"""Whole-page text OCR for scanned patents (easyocr), with on-disk cache.

Scanned WO-A1 PDFs have no text layer, so we render each page to a PNG at high
zoom (default 3.0 ~ 172 dpi equivalent -> chemical text becomes reliably
readable) and run easyocr. Per-page results (joined text + word boxes) are
cached to JSON so long batches are resumable and re-runs are cheap.

OCR quality is very zoom-sensitive on these dense chemical pages: at 1.4x the
numbers garble ("416"->"4[6"); at 3.0x "Preparation of Intermediate 416" and
"ES-MS m/z 683 (M+H)" come out clean.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pymupdf

_READER = None


def _easyocr_model_dir() -> str:
    """easyocr 权重落工作区（默认 ~/.EasyOCR 会写 C 盘）。

    仓库根 = ocr_text.py 向上 7 级（pipeline→core→patentsar→engines→
    melodylab→backend→MelodyLab），与 run.py 的 PYSTOW_HOME 同一个 .model_cache。
    """
    d = Path(__file__).resolve().parents[6] / ".model_cache" / "EasyOCR"
    d.mkdir(parents=True, exist_ok=True)
    return str(d)


def get_reader():
    global _READER
    if _READER is None:
        import easyocr
        _READER = easyocr.Reader(["en"], gpu=False, verbose=False,
                                 model_storage_directory=_easyocr_model_dir())
    return _READER


@dataclass
class OcrWord:
    text: str
    conf: float
    bbox: tuple[int, int, int, int]  # x0,y0,x1,y1 in image pixels


def _boxes_to_words(result) -> list[OcrWord]:
    words = []
    for box, text, conf in result:
        xs = [p[0] for p in box]
        ys = [p[1] for p in box]
        words.append(OcrWord(text, round(float(conf), 3),
                             (int(min(xs)), int(min(ys)), int(max(xs)), int(max(ys)))))
    return words


def render_page(pdf: str | Path, idx: int, out_png: Path, zoom: float = 3.0) -> tuple[int, int]:
    out_png.parent.mkdir(parents=True, exist_ok=True)
    doc = pymupdf.open(str(pdf))
    pix = doc[idx].get_pixmap(matrix=pymupdf.Matrix(zoom, zoom))
    out_png.write_bytes(pix.tobytes("png"))
    w, h = pix.width, pix.height
    doc.close()
    return w, h


def get_page(pdf: str | Path, idx: int, workdir: str | Path, zoom: float = 3.0) -> dict:
    """Return {'text', 'words':[{text,conf,bbox}], 'image', 'size'} for a page, cached."""
    workdir = Path(workdir)
    cache = workdir / "ocr" / f"page_{idx}_z{int(zoom*10)}.json"
    if cache.exists():
        return json.loads(cache.read_text())

    img = workdir / "pages" / f"page_{idx}.png"
    w, h = render_page(pdf, idx, img, zoom=zoom)
    reader = get_reader()
    raw = reader.readtext(str(img), detail=1)
    words = _boxes_to_words(raw)
    text = " ".join(wd.text for wd in words)
    rec = {"page": idx, "zoom": zoom, "image": str(img), "size": [w, h],
           "text": text, "words": [{"text": wd.text, "conf": wd.conf,
                                    "bbox": wd.bbox} for wd in words]}
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(rec))
    return rec


def words_from_rec(rec: dict) -> list[OcrWord]:
    return [OcrWord(w["text"], w["conf"], tuple(w["bbox"])) for w in rec["words"]]


if __name__ == "__main__":
    import sys
    pdf = sys.argv[1] if len(sys.argv) > 1 else "../Patent_test/WO_2025240893_A1.pdf"
    idx = int(sys.argv[2]) if len(sys.argv) > 2 else 247
    rec = get_page(pdf, idx, "output/work_scan", zoom=3.0)
    print("page", rec["page"], "size", rec["size"], "words", len(rec["words"]))
    print(rec["text"][:400])
