"""VLM label disambiguation for scanned structure pages.

The scanned path (run_scanned) assigns Example labels to structure crops
purely geometrically (nearest heading above the crop). That breaks on the
known weak spots: an enantiomer pair drawn under one heading, OCR-garbled
label text ("9b" -> "96"), or a crop that swallowed two molecules. This module
sends the full page render plus every crop to a vision-language model and asks
it to (re-)assign each crop to the Example/Intermediate label printed nearest
to it in the image, and to report crops that contain more than one molecule.

Any OpenAI-compatible VLM endpoint works (Zhipu GLM-4V by default). Configure
via SAR_LLM_BASE_URL / SAR_VLM_MODEL / SAR_LLM_API_KEY. With no API key (or on
any VLM error) callers keep the heuristic labels — degradation is silent by
design so batch runs never crash on the VLM.

Usage (debug one page):
  python -m pipeline.vlm --page output_scan/work/pages/page_3.png \
      --crop output_scan/work/crops/a.png --crop .../b.png
"""
from __future__ import annotations

import argparse
import base64
import re
import sys

from ..sar.summarize import _chat, _extract_json  # shared HTTP + JSON extraction

EX_NUM = r"[0-9]{1,4}[a-z]?"
VALID_KINDS = {"example", "intermediate", "none"}

SYSTEM_PROMPT = """你是专利文档版面分析助手。输入是一页扫描专利（第 1 张图）和从这页
切出的候选结构图（第 2 张起，按顺序编号 0..n-1）。任务：
1. 对每个候选图，找出页面上离它最近的 "Example <编号>" 或 "Intermediate <编号>"
   印刷标签（通常在结构图正下方或左下方），把标签编号原样读出来——注意 OCR 可能把
   编号读错（如 9b 读成 96），你要以图像为准。
2. 如果某张候选图里画了不止一个分子（成对异构体常见），把 multi_molecule 设为 true。
3. 如果某张候选图不是分子结构（表头、图表、碎片），kind 填 "none"。
只输出一个 JSON 对象，不要解释：
{"assignments": [{"crop_index": 0, "kind": "example|intermediate|none",
  "example_no": "<编号原样，kind 为 none 时填空字符串>", "multi_molecule": false}]}
assignments 必须恰好覆盖 0..n-1 每个候选图各一条。"""


def default_config(base_url: str | None = None, model: str | None = None,
                   api_key: str | None = None) -> dict:
    import os
    return {"base_url": base_url or os.environ.get("SAR_LLM_BASE_URL",
            "https://open.bigmodel.cn/api/paas/v4"),
            "model": model or os.environ.get("SAR_VLM_MODEL", "glm-4v-flash"),
            "api_key": api_key or os.environ.get("SAR_LLM_API_KEY", "")}


def enabled(cfg: dict) -> bool:
    return bool(cfg.get("api_key"))


def _b64_image(path: str) -> dict:
    data = base64.b64encode(open(path, "rb").read()).decode()
    return {"type": "image_url",
            "image_url": {"url": f"data:image/png;base64,{data}"}}


def build_messages(page_png: str, crops: list[dict]) -> list[dict]:
    crops_desc = "\n".join(
        f"- crop_index {i}: heuristic_label={c.get('label')!r}" for i, c in enumerate(crops))
    user = (f"页面：{page_png}\n候选结构图共 {len(crops)} 张，启发式标签（OCR+几何，可能错）：\n"
            f"{crops_desc}\n请按系统指令输出 assignments JSON。")
    content: list[dict] = [{"type": "text", "text": user},
                           _b64_image(page_png)]
    for c in crops:
        content.append(_b64_image(c["crop_path"]))
    return [{"role": "system", "content": SYSTEM_PROMPT},
            {"role": "user", "content": content}]


def parse_assignments(text: str, n_crops: int) -> list[dict] | None:
    """Validate the VLM's assignments; None if unusable."""
    parsed = _extract_json(text)
    if not isinstance(parsed, dict):
        return None
    items = parsed.get("assignments")
    if not isinstance(items, list):
        return None
    out: dict[int, dict] = {}
    for it in items:
        try:
            idx = int(it.get("crop_index"))
        except (TypeError, ValueError):
            return None
        if not 0 <= idx < n_crops or idx in out:
            return None
        kind = str(it.get("kind", "none")).lower()
        if kind not in VALID_KINDS:
            return None
        num = str(it.get("example_no", "") or "").strip()
        if kind == "none":
            num = ""
        elif not re.fullmatch(EX_NUM, num):
            return None
        out[idx] = {"crop_index": idx, "kind": kind, "example_no": num,
                    "multi_molecule": bool(it.get("multi_molecule"))}
    if len(out) != n_crops:
        return None
    return [out[i] for i in range(n_crops)]


def disambiguate_page(page_png: str, crops: list[dict], cfg: dict,
                      timeout: int = 300) -> list[dict] | None:
    """Return validated assignments (one per crop) or None on any failure."""
    if not enabled(cfg):
        return None
    try:
        text = _chat(build_messages(page_png, crops), cfg["base_url"],
                     cfg["api_key"], cfg["model"], temperature=0.1, timeout=timeout)
        return parse_assignments(text, len(crops))
    except Exception as e:  # noqa: BLE001 - VLM must never crash the batch
        print(f"[vlm] disambiguation failed ({type(e).__name__}: {e}); "
              f"keeping heuristic labels", file=sys.stderr)
        return None


def _selftest() -> int:
    good = ('```json\n{"assignments": ['
            '{"crop_index": 0, "kind": "example", "example_no": "9a", "multi_molecule": false},'
            '{"crop_index": 1, "kind": "example", "example_no": "9b", "multi_molecule": true}]}'
            '\n```')
    cases = [
        (good, 2, [{"kind": "example", "example_no": "9a", "multi_molecule": False},
                   {"kind": "example", "example_no": "9b", "multi_molecule": True}]),
        ('{"assignments": [{"crop_index": 0, "kind": "example", "example_no": "9a"}]}', 2,
         None),  # incomplete coverage
        ('{"assignments": [{"crop_index": 0, "kind": "example", "example_no": "9ab"},'
         '{"crop_index": 1, "kind": "example", "example_no": "9b"}]}', 2, None),  # bad number
        ('{"assignments": [{"crop_index": 0, "kind": "none", "example_no": ""},'
         '{"crop_index": 1, "kind": "example", "example_no": "214"}]}', 2,
         [{"kind": "none", "example_no": "", "multi_molecule": False},
          {"kind": "example", "example_no": "214", "multi_molecule": False}]),
        ('没有 JSON', 2, None),
    ]
    fails = 0
    for text, n, expect in cases:
        got = parse_assignments(text, n)
        if expect is None:
            ok = got is None
        else:
            ok = (got is not None and len(got) == len(expect)
                  and all(all(g[k] == e[k] for k in e) for g, e in zip(got, expect)))
        status = "PASS" if ok else "FAIL"
        fails += 0 if ok else 1
        print(f"  [{status}] vlm selftest parse case")
    return fails


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--page", help="rendered page PNG")
    ap.add_argument("--crop", action="append", default=[],
                    help="crop PNG (repeatable, in segmentation order)")
    ap.add_argument("--label", action="append", default=[],
                    help="heuristic label per crop (repeatable, same order)")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(1 if _selftest() else 0)
    if not a.page or not a.crop:
        ap.error("--page and at least one --crop required (or --selftest)")
    crops = [{"crop_path": c, "label": a.label[i] if i < len(a.label) else None}
             for i, c in enumerate(a.crop)]
    cfg = default_config()
    res = disambiguate_page(a.page, crops, cfg)
    print(res if res is not None else "no VLM available / unusable response")


if __name__ == "__main__":
    main()
