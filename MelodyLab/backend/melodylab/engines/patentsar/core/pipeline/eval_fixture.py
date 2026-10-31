"""Score a mined dataset CSV against fixture ground truth.

Turns the fixtures into an actual benchmark: per-example structure match
(canonical), stereochemistry match (isomeric), assay and LC-MS agreement,
plus explicit accounting for example kinds the pipeline does not support yet
(Markush R-table-defined examples are expected to be missing - that is the
honest baseline the P2 Markush module has to beat).

Usage:
  python -m pipeline.run --pdf data/fixture_patent_hard.pdf --structure-pages 1:5 \
      --assay-pages 6 --lcms-pages 7 --outdir output_hard
  python -m pipeline.eval_fixture --dataset output_hard/fixture_patent_hard_dataset.csv \
      --truth data/fixture_ground_truth_hard.json
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

MZ_TOL = 0.02  # [M+H]+ agreement window, Da


def _canon(smi: str | None, isomeric: bool = True) -> str | None:
    if not smi:
        return None
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(m, isomericSmiles=isomeric) if m else None


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def score(dataset_csv: str | Path, truth_json: str | Path) -> dict:
    with open(dataset_csv, newline="") as fh:
        rows = {r["example_no"]: r for r in csv.DictReader(fh)}
    truth = json.load(open(truth_json))

    counts = {"pass": 0, "partial_stereo": 0, "wrong_structure": 0, "missing": 0,
              "unsupported_expected": 0, "unexpected_present": 0}
    assay_ok = assay_n = mz_ok = mz_n = 0
    per: list[dict] = []
    for t in truth:
        tid, kind = str(t["example_no"]), t["kind"]
        row = rows.get(tid)
        e: dict = {"example_no": tid, "kind": kind}
        if kind == "intermediate":
            e["status"] = "pass" if row is None else "unexpected_present"
        elif kind == "markush_example":
            # defined by R-combination only: no drawn structure, so the
            # current pipeline cannot extract it (P2 Markush module target)
            if row is None:
                e["status"] = "unsupported_expected"
            else:
                e["status"] = "pass"
                e["note"] = "Markush-defined example extracted"
        else:
            dsmi = (row or {}).get("canonical_smiles") or ""
            tsmi = t.get("canonical_smiles") or ""
            if row is None or not dsmi:
                e["status"] = "missing"
            elif _canon(dsmi) == _canon(tsmi):
                e["status"] = "pass"
            elif _canon(dsmi, isomeric=False) == _canon(tsmi, isomeric=False):
                e["status"] = "partial_stereo"
                e["detail"] = f"got {dsmi} vs truth {tsmi}"
            else:
                e["status"] = "wrong_structure"
                e["detail"] = f"got {dsmi} vs truth {tsmi}"
            if row is not None:
                for col in ("ec50_amy3r", "emax_amy3r", "ec50_ctr", "emax_ctr"):
                    tv, dv = t.get(col), (row.get(col) or "").strip()
                    if tv is None:
                        continue
                    assay_n += 1
                    tn, dn = _num(tv), _num(dv)
                    if tv == dv or (tn is not None and dn is not None
                                    and abs(tn - dn) <= 1e-9 * max(1.0, abs(tn))):
                        assay_ok += 1
                    else:
                        e.setdefault("assay_mismatch", []).append(
                            f"{col}: {dv!r} vs {tv!r}")
                tmz, dmz = t.get("lcms_mz"), _num(row.get("mz_observed"))
                if tmz is not None and dmz is not None:
                    mz_n += 1
                    if abs(dmz - tmz) <= MZ_TOL:
                        mz_ok += 1
                    else:
                        e["mz_mismatch"] = f"{dmz} vs {tmz}"
        counts[e["status"]] = counts.get(e["status"], 0) + 1
        per.append(e)

    truth_ids = {str(t["example_no"]) for t in truth}
    extras = sorted(i for i in rows if i not in truth_ids)
    return {"dataset": str(dataset_csv), "truth": str(truth_json),
            "counts": counts, "assay": {"ok": assay_ok, "n": assay_n},
            "lcms": {"ok": mz_ok, "n": mz_n}, "extra_rows": extras,
            "per_example": per}


def _print_report(res: dict) -> None:
    def n_of(kind, status):
        return sum(1 for e in res["per_example"]
                   if e["kind"] == kind and e["status"] == status)
    n_drawn = sum(1 for e in res["per_example"] if e["kind"] in ("example", "reference"))
    print(f"drawn examples (example/reference, n={n_drawn}): "
          f"pass={n_of('example', 'pass') + n_of('reference', 'pass')} "
          f"| stereo-only={res['counts'].get('partial_stereo', 0)} "
          f"| wrong={res['counts'].get('wrong_structure', 0)} "
          f"| missing={res['counts'].get('missing', 0)}")
    print(f"markush-defined examples: extracted={n_of('markush_example', 'pass')} "
          f"unsupported_expected={res['counts'].get('unsupported_expected', 0)}")
    print(f"intermediate handling: unexpected_present={res['counts'].get('unexpected_present', 0)}")
    print(f"assay agreement: {res['assay']['ok']}/{res['assay']['n']}")
    print(f"LCMS m/z agreement: {res['lcms']['ok']}/{res['lcms']['n']}")
    if res["extra_rows"]:
        print(f"extra rows (no truth match): {', '.join(res['extra_rows'])}")
    for e in res["per_example"]:
        if e["status"] != "pass":
            line = f"  - Ex {e['example_no']} [{e['status']}]"
            if e.get("detail"):
                line += f" {e['detail']}"
            if e.get("assay_mismatch"):
                line += f" assay: {e['assay_mismatch']}"
            if e.get("mz_mismatch"):
                line += f" mz: {e['mz_mismatch']}"
            print(line)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", required=True, help="mined dataset CSV to score")
    ap.add_argument("--truth", required=True, help="fixture ground truth JSON")
    ap.add_argument("--report", default=None,
                    help="where to write the JSON score file "
                         "(default: <dataset stem>_eval.json next to the dataset)")
    a = ap.parse_args()
    res = score(a.dataset, a.truth)
    _print_report(res)
    out = a.report or str(Path(a.dataset).with_name(
        Path(a.dataset).stem + "_eval.json"))
    Path(out).write_text(json.dumps(res, indent=1))
    print(f"score file -> {out}")


if __name__ == "__main__":
    main()
