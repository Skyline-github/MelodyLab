"""Write the stage-III deliverables: SDF, CSV, HTML check table, Markdown report.

Each 'row' is a joined record:
    {example_no, image_rel, qc(chem_qc.validate dict), assay(dict|None),
     ocr_confidence(float|None), structure_page(int|None)}
Flagged rows (INVALID_SMILES / LCMS_MISMATCH / DUPLICATE) are highlighted in
the HTML for human review, as the blog's workflow does.
"""
from __future__ import annotations

import csv
import html
from pathlib import Path

from rdkit import Chem


def _sort_key(example_no: str):
    s = str(example_no)
    return (0, int(s)) if s.isdigit() else (1, s)


def _flags(row) -> list[str]:
    return row["qc"]["flags"]


def _v(val):
    """Render a missing assay value as 'ND' (never a blank cell)."""
    return "ND" if val is None or val == "" else val


def write_csv(rows: list[dict], path: Path) -> None:
    cols = ["example_no", "canonical_smiles", "smiles_raw", "formula", "mw_avg",
            "mz_calc", "mz_observed", "delta_mz", "ocr_confidence", "kind",
            "ec50_amy3r", "emax_amy3r", "ec50_ctr", "emax_ctr", "flags"]
    rows = sorted(rows, key=lambda r: _sort_key(r["example_no"]))
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            q, a = r["qc"], (r.get("assay") or {})
            w.writerow([
                r["example_no"], q.get("canonical_smiles"), r.get("raw_smiles"),
                q.get("formula"),
                q.get("mw_avg"), q.get("mz_calc"), q.get("mz_observed"),
                q.get("delta_mz"), r.get("ocr_confidence"), r.get("kind"),
                _v(a.get("ec50_amy3r")), _v(a.get("emax_amy3r")),
                _v(a.get("ec50_ctr")), _v(a.get("emax_ctr")),
                ";".join(_flags(r)),
            ])


def write_sdf(rows: list[dict], path: Path) -> None:
    rows = sorted(rows, key=lambda r: _sort_key(r["example_no"]))
    w = Chem.SDWriter(str(path))
    for r in rows:
        smi = r["qc"].get("canonical_smiles")
        if not smi:
            continue
        mol = Chem.MolFromSmiles(smi)
        a = r.get("assay") or {}
        mol.SetProp("_Name", f"Example {r['example_no']}")
        mol.SetProp("example_no", str(r["example_no"]))
        for k, v in (("formula", r["qc"].get("formula")),
                     ("mw_avg", r["qc"].get("mw_avg")),
                     ("ocr_confidence", r.get("ocr_confidence")),
                     ("ec50_amy3r", _v(a.get("ec50_amy3r"))),
                     ("emax_amy3r", _v(a.get("emax_amy3r"))),
                     ("ec50_ctr", _v(a.get("ec50_ctr"))),
                     ("emax_ctr", _v(a.get("emax_ctr"))),
                     ("flags", ";".join(_flags(r)))):
            if v is not None:
                mol.SetProp(k, str(v))
        w.write(mol)
    w.close()


def write_html(rows: list[dict], path: Path, title: str) -> None:
    rows = sorted(rows, key=lambda r: _sort_key(r["example_no"]))
    head = ("<tr>" + "".join(f"<th>{c}</th>" for c in
            ["Example", "Structure", "Canonical SMILES", "Formula", "MW",
             "m/z(calc/obs)", "conf", "AMY3R EC50/Emax", "CTR EC50/Emax",
             "Flags"]) + "</tr>")
    body = []
    for r in rows:
        q, a = r["qc"], (r.get("assay") or {})
        flagged = bool(_flags(r))
        cls = ' style="background:#ffe0e0"' if flagged else ""
        img = (f'<img src="{html.escape(str(r["image_rel"]))}" width="150">'
               if r.get("image_rel") else "(no structure)")
        mz = (f"{q.get('mz_calc')} / {q.get('mz_observed')}"
              if q.get("mz_observed") is not None else "-")
        body.append(
            f"<tr{cls}>"
            f"<td>{html.escape(str(r['example_no']))}</td>"
            f"<td>{img}</td>"
            f"<td><code>{html.escape(str(q.get('canonical_smiles') or ''))}</code></td>"
            f"<td>{q.get('formula')}</td><td>{q.get('mw_avg')}</td><td>{mz}</td>"
            f"<td>{r.get('ocr_confidence')}</td>"
            f"<td>{_v(a.get('ec50_amy3r'))} / {_v(a.get('emax_amy3r'))}</td>"
            f"<td>{_v(a.get('ec50_ctr'))} / {_v(a.get('emax_ctr'))}</td>"
            f"<td>{html.escape(';'.join(_flags(r)))}</td></tr>")
    doc = (f"<!doctype html><meta charset='utf-8'><title>{html.escape(title)}</title>"
           f"<style>body{{font-family:sans-serif}}table{{border-collapse:collapse}}"
           f"td,th{{border:1px solid #ccc;padding:4px 8px}}</style>"
           f"<h2>{html.escape(title)}</h2><table>{head}{''.join(body)}</table>")
    path.write_text(doc)


def write_markdown_report(rows: list[dict], meta: dict, path: Path) -> None:
    rows = sorted(rows, key=lambda r: _sort_key(r["example_no"]))
    flagged = [r for r in rows if _flags(r)]
    lines = [f"# Extraction report - {meta.get('pdf', '')}", ""]
    lines += [f"- PDF: `{meta.get('pdf')}`",
              f"- Structure page range: {meta.get('structure_pages')}",
              f"- Assay page(s): {meta.get('assay_pages')}",
              f"- LCMS page(s): {meta.get('lcms_pages')}",
              f"- Structures extracted: {meta.get('n_structures')}",
              f"- OCR model: DECIMER ({meta.get('hand_drawn','no')} hand_drawn)",
              f"- Valid SMILES: {sum(1 for r in rows if r['qc']['valid'])}",
              f"- Flagged for review: {len(flagged)}", ""]
    lines.append("| Example | canonical SMILES | formula | m/z obs/calc | flags |")
    lines.append("|---|---|---|---|---|")
    for r in rows:
        q = r["qc"]
        lines.append(f"| {r['example_no']} | `{q.get('canonical_smiles')}` | "
                     f"{q.get('formula')} | {q.get('mz_observed')}/{q.get('mz_calc')} | "
                     f"{';'.join(_flags(r)) or '-'} |")
    lines.append("")
    lines.append("## Notes")
    lines.append("- `ND` = not determined; graded EC50 letters (A/B/...) kept as text.")
    lines.append("- `LCMS_MISMATCH`: RDKit [M+H]+ vs patent m/z differs beyond tolerance.")
    lines.append("- Structures are OCR predictions and must be human-verified.")
    path.write_text("\n".join(lines))
