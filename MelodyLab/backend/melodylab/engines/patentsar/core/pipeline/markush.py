"""Markush R-table example enumeration (P2 MVP).

Real medicinal-chemistry patents define many series members only as R-group
combinations against a generic structure ("Example 12: R1 = Et, R2 = H")
without drawing each concrete molecule - the structure-extraction path cannot
see them, which is why the hard fixture scored 0/3 on Markush-defined
examples. Given an explicit core template with {rN} placeholders plus a
label -> SMILES-snippet table per position (an honest, user-specified Markush
core, as the README's known-limitations note recommends), this module
enumerates the concrete SMILES for those examples, RDKit-validates them, and
the orchestrators append them as fully joined dataset rows
(kind=markush_example, same QC incl. LC-MS cross-check as drawn structures).

Template semantics: each {rN} slot is replaced verbatim by the position's
SMILES snippet, so prefix attachment ("{r1}NC(=O)..."), ring substitution
("...c2{r2}ncc2c1" with value "[n](C)") and anything else SMILES can express
are all supported. Every enumerated molecule is validated with RDKit; a
combination whose label is missing from the table is skipped (and reported),
never guessed.
"""
from __future__ import annotations

import re

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

EXAMPLE_HEAD = re.compile(r"\bExamples?\s*[:\-]?\s*(\d{1,4})\b")
# [=:]* tolerant of easyocr dropping the '=' sign across word boxes
R_PAIR = re.compile(r"\bR\s*(\d)\s*[=:]*\s*([A-Za-z0-9]+)")
_PLACEHOLDER = re.compile(r"\{(\w+)\}")


def load_config(path: str) -> dict:
    import json
    cfg = json.load(open(path))
    if not cfg.get("core_template") or not cfg.get("positions"):
        raise ValueError("markush config needs 'core_template' and 'positions'")
    return cfg


def parse_combinations(lines: list[str]) -> dict[str, dict[str, str]]:
    """'Example 12: R1 = Et, R2 = H' -> {'12': {'R1': 'Et', 'R2': 'H'}}.

    Lines without an Example header (e.g. the R-group definition table itself)
    are ignored. Later lines for the same Example number win.
    """
    out: dict[str, dict[str, str]] = {}
    current = None
    for line in lines:
        m = EXAMPLE_HEAD.search(line)
        if m:
            current = m.group(1)
            pairs = R_PAIR.findall(line[m.end():])
            if pairs:
                out[current] = {f"R{p}": v for p, v in pairs}
            else:
                out.pop(current, None)
                current = None
            continue
        if current:
            pairs = R_PAIR.findall(line)
            if pairs:
                out[current].update({f"R{p}": v for p, v in pairs})
    return out


def _lookup(defs: dict, label: str) -> str | None:
    for k, v in defs.items():
        if k == label or k.lower() == label.lower():
            return v
    return None


def enumerate_smiles(config: dict, combo: dict[str, str],
                     ) -> tuple[str | None, str | None]:
    """Return (canonical_smiles, error). error is None on success; otherwise a
    short reason ('unknown R label X' / 'invalid SMILES ...') for the report."""
    smi = config["core_template"]
    pos_map = {k.lower(): v for k, v in config["positions"].items()}
    for rkey, label in combo.items():
        defs = pos_map.get(rkey.lower())
        if defs is None:
            return None, f"position {rkey} not in config"
        frag = _lookup(defs, label)
        if frag is None:
            return None, f"unknown R label {rkey}={label!r}"
        smi = smi.replace("{" + rkey.lower() + "}", frag)
    leftover = _PLACEHOLDER.search(smi)
    if leftover:
        return None, f"unfilled placeholder {leftover.group(0)}"
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None, f"invalid SMILES: {smi}"
    return Chem.MolToSmiles(mol), None


def harvest_scanned_lines(pdf, pages: list[int], workdir, zoom: float = 3.0,
                          y_tol: float = 22.0) -> list[str]:
    """Rebuild visual text lines from easyocr word boxes on the given pages
    (easyocr boxes are line-level already, but JPEG noise can split one visual
    line into several boxes at the same height)."""
    from . import ocr_text

    def line_texts(rec) -> list[str]:
        rows: list[list] = []
        for w in sorted(rec["words"], key=lambda w: (w["bbox"][1], w["bbox"][0])):
            cy = (w["bbox"][1] + w["bbox"][3]) / 2
            if rows and abs(cy - rows[-1][0]) < y_tol:
                rows[-1][1].append(w)
            else:
                rows.append([cy, [w]])
        out = []
        for _, ws in rows:
            ws.sort(key=lambda w: w["bbox"][0])
            t = " ".join(w["text"] for w in ws).strip()
            if t:
                out.append(t)
        return out

    lines: list[str] = []
    for p in pages:
        lines += line_texts(ocr_text.get_page(pdf, p, workdir, zoom=zoom))
    return lines


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print("usage: python -m pipeline.markush <config.json> <text-or-lines...>")
        sys.exit(2)
    cfg = load_config(sys.argv[1])
    combos = parse_combinations(sys.argv[2:])
    for ex, combo in combos.items():
        smi, err = enumerate_smiles(cfg, combo)
        print(f"Example {ex} {combo} -> {smi or f'SKIP ({err})'}")
