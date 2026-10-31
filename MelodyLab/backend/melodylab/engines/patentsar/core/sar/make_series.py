"""Generate a synthetic congeneric SAR series to demonstrate the sar module.

4-hydroxybenzamide scaffold with two DISTINCT attachment points (R1 on the
amide N, R2 on the phenol O) so R-group positions are unambiguous. Activity is
a designed additive (Free-Wilson-like) pEC50, so the RGA boxplots, the R1xR2
combination heatmap, and the Free-Wilson model all show clean, recoverable
signal. Clearly synthetic - for capability demo only.
"""
from __future__ import annotations

import csv
import random

CORE = "NC(=O)c1ccc(O)cc1"  # R1 on amide N, R2 on phenol O (asymmetric)

R1 = {"Me": "C", "Et": "CC", "nPr": "CCC", "iPr": "C(C)C", "cPr": "C1CC1",
      "Bn": "Cc1ccccc1", "Ph": "c1ccccc1", "HEt": "CCO"}
R2 = {"Me": "C", "Et": "CC", "iPr": "C(C)C", "All": "CC=C", "Bn": "Cc1ccccc1",
      "cPrM": "CC1CC1"}

# designed pEC50 contributions
A1 = {"Me": 0.0, "Et": 0.4, "nPr": 0.9, "iPr": 0.7, "cPr": 1.1, "Bn": 1.6,
      "Ph": 1.3, "HEt": 0.2}
A2 = {"Me": 0.0, "Et": 0.5, "iPr": 1.0, "All": 0.6, "Bn": 1.4, "cPrM": 1.2}
BASE = 6.0


def _ring(smi: str, digit: str) -> str:
    """Remap single-digit ring-closure numbers so concatenated fragments don't
    collide (our fragments use only digit '1')."""
    return smi.replace("1", digit)


def build(r1: str, r2: str) -> str:
    # core ring=1, R1 ring=2, R2 ring=3
    return f"{_ring(R1[r1], '2')}NC(=O)c1ccc(O{_ring(R2[r2], '3')})cc1"


def main(out: str = "data/sar_series.csv"):
    random.seed(7)
    rows, n = [], 0
    for r1 in R1:
        for r2 in R2:
            if random.random() < 0.28:
                continue
            n += 1
            pec50 = BASE + A1[r1] + A2[r2] + random.gauss(0, 0.15)
            ec50_nm = 10 ** (9 - pec50)
            rows.append({"example_no": str(n), "canonical_smiles": build(r1, r2),
                         "ec50_amy3r": f"{ec50_nm:.4g}", "R1": r1, "R2": r2})
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["example_no", "canonical_smiles",
                                            "ec50_amy3r", "R1", "R2"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} compounds -> {out}")


if __name__ == "__main__":
    main()
