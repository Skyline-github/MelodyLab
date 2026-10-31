"""Batch structure extraction across the three Lilly DACRA patents.

Robust to OOM on big scanned docs (091 is 697 pages): each patent is run in a
SUBPROCESS, and if it gets OOM-killed mid-range we simply retry. Per-page OCR
and DECIMER caches live on disk, so a restart skips already-done pages cheaply
(few new pages processed before the memory ceiling), and the loop converges
until one attempt finishes the whole range and writes the deliverables.
"""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

JOBS = [
    ("WO_2025240893_A1", "../Patent_test/WO_2025240893_A1.pdf", (247, 385)),
    ("WO_2026193091_A1", "../Patent_test/WO_2026193091_A1.pdf", (221, 600)),
    ("WO_2026193092_A1", "../Patent_test/WO_2026193092_A1.pdf", (66, 154)),
]

ROOT = Path("output_scan")
MAX_ATTEMPTS = 15


def run_patent(name, pdf, rng):
    outdir = ROOT / name
    workdir = outdir / "work"
    outdir.mkdir(parents=True, exist_ok=True)
    csv = outdir / f"{name}_dataset.csv"
    if csv.exists():
        print(f"=== {name} already done -> skip ===", flush=True)
        return
    for attempt in range(1, MAX_ATTEMPTS + 1):
        cached = len(list((workdir / "ocr").glob("*.json"))) if (workdir / "ocr").exists() else 0
        print(f"=== {name} attempt {attempt} (cached pages={cached}) ===", flush=True)
        t0 = time.time()
        rc = subprocess.run([
            sys.executable, "-m", "pipeline.run_scanned",
            "--pdf", pdf, "--structure-pages", f"{rng[0]}:{rng[1]}",
            "--kinds", "example", "--outdir", str(outdir), "--workdir", str(workdir),
        ]).returncode
        if csv.exists():
            print(f"{name}: DONE on attempt {attempt} ({time.time()-t0:.0f}s) -> {csv}",
                  flush=True)
            return
        print(f"{name}: attempt {attempt} ended rc={rc} (likely OOM), retrying...",
              flush=True)
        time.sleep(3)
    print(f"{name}: FAILED after {MAX_ATTEMPTS} attempts", flush=True)


def main():
    for name, pdf, rng in JOBS:
        run_patent(name, pdf, rng)
    print("\nALL DONE", flush=True)


if __name__ == "__main__":
    main()
