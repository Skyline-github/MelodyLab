"""DECIMER optical chemical structure recognition (OCSR) wrapper.

Turns structure PNGs into SMILES. Results are cached on disk by image-content
hash so re-running the pipeline (or resuming an interrupted batch) does not
re-infer already-processed structures. The heavy TF/DECIMER import is done
lazily on first use.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")

_PREDICT = None


def _jsonable(o):
    import numpy as _np
    if isinstance(o, (_np.floating, _np.integer)):
        return o.item()
    if isinstance(o, _np.ndarray):
        return o.tolist()
    return str(o)


def _dump(obj) -> str:
    return json.dumps(obj, indent=2, default=_jsonable)


def _load_predictor():
    global _PREDICT
    if _PREDICT is None:
        import warnings
        warnings.filterwarnings("ignore")
        from DECIMER import predict_SMILES
        _PREDICT = predict_SMILES
    return _PREDICT


def _sha1(path: Path) -> str:
    return hashlib.sha1(path.read_bytes()).hexdigest()


def _load_cache(cache_json: Path) -> dict:
    if cache_json.exists():
        return json.loads(cache_json.read_text())
    return {}


def ocr_images(
    image_paths: list[Path],
    cache_json: str | Path = "output/work/ocr_cache.json",
    hand_drawn: bool = False,
    confidence: bool = False,
) -> dict[str, dict]:
    """Return {image_path_str: {smiles, confidence?, from_cache}}."""
    cache_json = Path(cache_json)
    cache_json.parent.mkdir(parents=True, exist_ok=True)
    cache = _load_cache(cache_json)
    predict = _load_predictor()
    results: dict[str, dict] = {}

    for img in map(Path, image_paths):
        key = _sha1(img)
        if key in cache:
            results[str(img)] = {**cache[key], "from_cache": True}
            continue
        payload = predict(str(img), confidence=confidence, hand_drawn=hand_drawn)
        if confidence:
            smiles, conf = payload
            rec = {"smiles": str(smiles), "confidence": conf}
        else:
            rec = {"smiles": str(payload)}
        cache[key] = rec
        results[str(img)] = {**rec, "from_cache": False}
        # flush after every inference so a crash does not lose completed work
        cache_json.write_text(_dump(cache))

    cache_json.write_text(_dump(cache))
    return results


if __name__ == "__main__":
    import sys
    from .pdf_images import extract_structure_images
    recs = extract_structure_images(sys.argv[1] if len(sys.argv) > 1 else "data/fixture_patent.pdf",
                                    "output/work", page_range=(1, 6))
    res = ocr_images([r.image_path for r in recs], confidence=True)
    for r in recs:
        rr = res[str(r.image_path)]
        print(f"Example {r.label}: {rr['smiles']}  conf={rr.get('confidence')}  cached={rr['from_cache']}")
