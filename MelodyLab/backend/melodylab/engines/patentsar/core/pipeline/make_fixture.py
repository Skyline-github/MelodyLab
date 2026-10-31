"""Deterministic, synthetic 'patent-like' PDF fixture for the mining pipeline.

Mirrors the structure of the real WO2025240893 case that the pipeline targets:
  * a Markush-like core scaffold with per-Example R-group variation
  * Example compound structure drawings on their own page region
  * an AMY3R / CTR activity table (Example no, EC50 nM, Emax %) incl. a reference
  * LCMS m/z text lines usable for a molecular-weight cross-check
It is generated with RDKit + PyMuPDF so the whole pipeline can be validated
end-to-end without downloading a 400-page real patent.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pymupdf
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import rdMolDraw2D

# Shared 1H-indazole-5-carboxamide core. The leading R string is the amide
# N-substituent, standing in for the patent's varied substituents.
CORE = "NC(=O)c1ccc2[nH]ncc2c1"


@dataclass
class Compound:
    example_no: str          # e.g. "1", "12", "ref" (reference cmpd)
    smiles: str              # full SMILES (empty string => intermediate w/o drawn structure)
    kind: str                # "example" | "intermediate" | "reference"
    ec50_amy3r: str | None   # numeric nM as text, graded grade (A/B/...), or None => ND
    emax_amy3r: str | None
    ec50_ctr: str | None
    emax_ctr: str | None
    lcms_mz: float | None    # observed [M+H]+ , m/z

    @property
    def mol(self):
        return Chem.MolFromSmiles(self.smiles) if self.smiles else None


def build_dataset() -> list[Compound]:
    # N-substituent fragments -> Example numbers (odd = examples, plus an
    # intermediate that has no drawn structure, plus a reference compound).
    r_frags = [
        ("1", "CC"),
        ("2", "CCO"),
        ("3", "CCCO"),
        ("4", "C(C)C"),
        ("5", "C1CC1"),
        ("6", "CCOCC"),
        ("7", "c1ccccc1"),
        ("8", "FC(F)(F)C"),
    ]
    compounds: list[Compound] = []
    # EC50 (nM) and Emax (%) chosen to look like a real SAR series.
    amy = {"1": (12.0, 88), "2": (35.0, 81), "3": (8.0, 92), "4": (140.0, 63),
           "5": (22.0, 79), "6": (None, None), "7": (76.0, 70), "8": (310.0, 41)}
    ctr = {"1": (1800.0, 55), "2": (950.0, 61), "3": (2400.0, 48), "4": (4100.0, 33),
           "5": (1200.0, 66), "6": (None, None), "7": (150.0, 72), "8": (8800.0, 22)}
    for ex, frag in r_frags:
        smiles = f"{frag}{CORE}"
        a = amy[ex]
        c = ctr[ex]
        compounds.append(Compound(ex, smiles, "example", _fmt(a[0]), _fmt(a[1]),
                                  _fmt(c[0]), _fmt(c[1]), None))
    # an intermediate: named in text but NOT an Example compound (no structure to draw)
    compounds.append(Compound("int-A", "", "intermediate", None, None, None, None, None))
    # a reference compound with graded activity
    compounds.append(Compound("ref", "CCNC(=O)c1ccc2cn[nH]c2c1", "reference",
                              "B", "72", "D", "40", None))
    # compute LCMS [M+H]+ from exact monoisotopic mass for every drawn structure
    for cpd in compounds:
        if cpd.mol is not None:
            mono = Descriptors.ExactMolWt(cpd.mol)
            cpd.lcms_mz = round(mono + 1.007276, 2)
    return compounds


def _fmt(v):
    if v is None:
        return None
    if isinstance(v, float):
        return ("%g" % v)
    return str(v)


def render_structure_png(cpd: Compound, out_png: Path) -> None:
    mol = cpd.mol
    d = rdMolDraw2D.MolDraw2DCairo(360, 260)
    d.drawOptions().addStereoAnnotation = False
    d.DrawMolecule(mol)
    d.FinishDrawing()
    out_png.write_bytes(d.GetDrawingText())


def make_pdf(compounds: list[Compound], png_dir: Path, pdf_path: Path) -> dict:
    """Build the fixture PDF and return the ground-truth layout metadata."""
    png_dir.mkdir(parents=True, exist_ok=True)
    structure_examples = [c for c in compounds if c.kind in ("example", "reference")]

    doc = pymupdf.open()
    layout = {"structure_pages": [], "assay_page": None, "intermediate_pages": []}

    # --- Page 1: title + Markush core (no per-compound structure) ---
    page = doc.new_page()
    page.insert_text((72, 90), "WO2025XXXXXX A1  (synthetic fixture)",
                     fontsize=14)
    page.insert_text((72, 120),
                     "4-aryl-1H-pyrazole-fused carboxamide derivatives as AMY3R agonists",
                     fontsize=10)

    # --- Structure pages: two Example compounds per page, image + caption ---
    per_page = 2
    x0 = 60
    for i in range(0, len(structure_examples), per_page):
        page = doc.new_page()
        pngs = []
        for j, cpd in enumerate(structure_examples[i:i + per_page]):
            png = png_dir / f"example_{cpd.example_no}.png"
            render_structure_png(cpd, png)
            rect = pymupdf.Rect(x0 + j * 250, 120, x0 + j * 250 + 240, 300)
            page.insert_image(rect, filename=str(png))
            page.insert_text((x0 + j * 250, 320), f"Example {cpd.example_no}", fontsize=10)
            pngs.append(png.name)
        layout["structure_pages"].append({"page_index": page.number,
                                          "images": pngs})

    # --- Assay page: text table (AMY3R + CTR), incl. ND and graded values ---
    page = doc.new_page()
    layout["assay_page"] = page.number
    y = 90
    page.insert_text((60, y), "Biological Example - AMY3R / CTR Assay (EC50 nM, Emax %)",
                     fontsize=11)
    y += 26
    header = "Example    AMY3R_EC50(nM)   AMY3R_Emax(%)    CTR_EC50(nM)   CTR_Emax(%)"
    page.insert_text((60, y), header, fontsize=9)
    for cpd in compounds:
        if cpd.kind == "intermediate":
            continue
        y += 18
        row = (f"{cpd.example_no:<9} {_fill(cpd.ec50_amy3r):<16} {_fill(cpd.emax_amy3r):<15} "
               f"{_fill(cpd.ec50_ctr):<13} {_fill(cpd.emax_ctr)}")
        page.insert_text((60, y), row, fontsize=9)

    # --- LCMS page: m/z per example for MW cross-check ---
    page = doc.new_page()
    page.insert_text((60, 90), "Analytical data (LCMS, [M+H]+)", fontsize=11)
    y = 116
    for cpd in structure_examples:
        if cpd.lcms_mz is not None:
            page.insert_text((60, y), f"Example {cpd.example_no}: LCMS m/z {cpd.lcms_mz}",
                             fontsize=9)
            y += 16

    doc.save(pdf_path)
    doc.close()
    return layout


def _fill(v):
    return "ND" if v is None else str(v)


# ---------------------------------------------------------------------------
# Hard fixture: a drawn enantiomer pair sharing one page (9a/9b), a
# two-centre diastereomer (10), and Markush R-table-defined examples (12-14,
# generic structure + R definitions, never drawn). These probe the pipeline's
# known weak spots; score with pipeline/eval_fixture.py.
# ---------------------------------------------------------------------------

MARKUSH_CORE_T = "NC(=O)c1ccc2{r2}ncc2c1"


def markush_smiles(r1: str, r2: str) -> str:
    """Concrete SMILES for an R1 (amide N) / R2 (indazole N) combination."""
    r2_frag = "[nH]" if not r2 else f"[n]({r2})"
    return f"{r1}{MARKUSH_CORE_T.format(r2=r2_frag)}"


def build_dataset_hard() -> list[Compound]:
    """Benchmark dataset: plain examples (1-4), drawn enantiomer pair 9a/9b,
    two-centre diastereomer 10, intermediate, reference, and Markush-defined
    examples 12-14 (concrete SMILES known for ground truth, but only given as
    R1/R2 combinations in the patent - no drawn structure to extract)."""
    soft = {"1": "CC", "2": "CCO", "3": "CCCO", "4": "C(C)C"}
    hard_smiles = {
        "9a": "CC[C@H](C)NC(=O)c1ccc2[nH]ncc2c1",
        "9b": "CC[C@@H](C)NC(=O)c1ccc2[nH]ncc2c1",
        "10": "C[C@H](O)[C@@H](C)NC(=O)c1ccc2[nH]ncc2c1",
        "12": markush_smiles("CC", ""),
        "13": markush_smiles("c1ccccc1", "C"),
        "14": markush_smiles("C1CC1", "C"),
    }
    # (example_no, smiles, kind) - drawn order = page layout order
    items = [(ex, f"{frag}{CORE}", "example") for ex, frag in soft.items()]
    items += [(ex, hard_smiles[ex], "example") for ex in ("9a", "9b", "10")]
    items += [(ex, hard_smiles[ex], "markush_example") for ex in ("12", "13", "14")]
    # EC50/Emax (AMY3R, CTR) - enantiomers differ, as they do in real patents
    acts = {"1": (12.0, 88, 1800.0, 55), "2": (35.0, 81, 950.0, 61),
            "3": (8.0, 92, 2400.0, 48), "4": (140.0, 63, 4100.0, 33),
            "9a": (45.0, 84, 1500.0, 58), "9b": (130.0, 66, 5200.0, 31),
            "10": (18.0, 90, 900.0, 64),
            "12": (55.0, 78, 2100.0, 52), "13": (6.5, 95, 3000.0, 45),
            "14": (15.0, 88, 1100.0, 61)}
    compounds: list[Compound] = []
    for ex, smiles, kind in items:
        a, c = acts[ex][:2], acts[ex][2:]
        compounds.append(Compound(ex, smiles, kind, _fmt(a[0]), _fmt(a[1]),
                                  _fmt(c[0]), _fmt(c[1]), None))
    compounds.append(Compound("int-A", "", "intermediate", None, None, None, None, None))
    compounds.append(Compound("ref", "CCNC(=O)c1ccc2cn[nH]c2c1", "reference",
                              "B", "72", "D", "40", None))
    for cpd in compounds:
        if cpd.mol is not None:
            assert Chem.MolFromSmiles(cpd.smiles), cpd.smiles
            mono = Descriptors.ExactMolWt(cpd.mol)
            cpd.lcms_mz = round(mono + 1.007276, 2)
    return compounds


def render_markush_png(out_png: Path) -> None:
    """The generic structure: shared core with R1 (amide N) and R2 (indazole
    N) as labelled dummy attachment points."""
    from rdkit.Chem import rdDepictor
    rw = Chem.RWMol(Chem.MolFromSmiles(CORE))
    amide_n = next(a.GetIdx() for a in rw.GetAtoms()
                   if a.GetSymbol() == "N" and not a.GetIsAromatic())
    ring_nh = next(a.GetIdx() for a in rw.GetAtoms()
                   if a.GetSymbol() == "N" and a.GetIsAromatic()
                   and a.GetTotalNumHs() == 1)
    nh = rw.GetAtomWithIdx(ring_nh)
    nh.SetNoImplicit(True)
    nh.SetNumExplicitHs(0)
    for idx, label in ((amide_n, "R1"), (ring_nh, "R2")):
        d = Chem.Atom(0)
        d.SetProp("dummyLabel", label)
        d.SetNoImplicit(True)
        i = rw.AddAtom(d)
        rw.AddBond(idx, i, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass
    rdDepictor.Compute2DCoords(mol)
    d = rdMolDraw2D.MolDraw2DCairo(360, 260)
    d.drawOptions().addStereoAnnotation = False
    d.DrawMolecule(mol)
    d.FinishDrawing()
    out_png.write_bytes(d.GetDrawingText())


def make_hard_pdf(compounds: list[Compound], png_dir: Path, pdf_path: Path) -> dict:
    """Hard-fixture PDF: drawn structures two per page (the 9a/9b enantiomer
    pair shares a page), a Markush page (generic structure + R definitions +
    combination-only Examples), then assay and LCMS pages."""
    png_dir.mkdir(parents=True, exist_ok=True)
    drawn = [c for c in compounds if c.kind in ("example", "reference")]
    doc = pymupdf.open()
    layout = {"structure_pages": [], "assay_page": None,
              "intermediate_pages": [], "markush_page": None}

    page = doc.new_page()
    page.insert_text((72, 90), "WO2025XXXXXX A1  (synthetic fixture, hard)",
                     fontsize=14)
    page.insert_text((72, 120),
                     "4-aryl-1H-pyrazole-fused carboxamide derivatives as AMY3R agonists",
                     fontsize=10)

    per_page, x0 = 2, 60
    for i in range(0, len(drawn), per_page):
        page = doc.new_page()
        pngs = []
        for j, cpd in enumerate(drawn[i:i + per_page]):
            png = png_dir / f"example_{cpd.example_no}.png"
            render_structure_png(cpd, png)
            rect = pymupdf.Rect(x0 + j * 250, 120, x0 + j * 250 + 240, 300)
            page.insert_image(rect, filename=str(png))
            page.insert_text((x0 + j * 250, 320), f"Example {cpd.example_no}", fontsize=10)
            pngs.append(png.name)
        layout["structure_pages"].append({"page_index": page.number,
                                          "images": pngs})

    # Markush page: generic structure + R-group definitions + examples given
    # as R combinations only (no drawn structures to extract)
    page = doc.new_page()
    layout["markush_page"] = page.number
    mk = png_dir / "markush_generic.png"
    render_markush_png(mk)
    page.insert_image(pymupdf.Rect(60, 90, 300, 270), filename=str(mk))
    y = 300
    page.insert_text((60, y), "R1: Me (methyl), Et (ethyl), Ph (phenyl), cPr (cyclopropyl)",
                     fontsize=9)
    page.insert_text((60, y + 14), "R2: H, Me (methyl)", fontsize=9)
    combos = [("12", "Et", "H"), ("13", "Ph", "Me"), ("14", "cPr", "Me")]
    for k, (ex, r1, r2) in enumerate(combos):
        page.insert_text((60, y + 32 + k * 14),
                         f"Example {ex}: R1 = {r1}, R2 = {r2}", fontsize=9)

    page = doc.new_page()
    layout["assay_page"] = page.number
    y = 90
    page.insert_text((60, y), "Biological Example - AMY3R / CTR Assay (EC50 nM, Emax %)",
                     fontsize=11)
    y += 26
    header = "Example    AMY3R_EC50(nM)   AMY3R_Emax(%)    CTR_EC50(nM)   CTR_Emax(%)"
    page.insert_text((60, y), header, fontsize=9)
    for cpd in compounds:
        if cpd.kind == "intermediate":
            continue
        y += 18
        row = (f"{cpd.example_no:<9} {_fill(cpd.ec50_amy3r):<16} {_fill(cpd.emax_amy3r):<15} "
               f"{_fill(cpd.ec50_ctr):<13} {_fill(cpd.emax_ctr)}")
        page.insert_text((60, y), row, fontsize=9)

    page = doc.new_page()
    page.insert_text((60, 90), "Analytical data (LCMS, [M+H]+)", fontsize=11)
    y = 116
    for cpd in compounds:
        if cpd.lcms_mz is not None:
            page.insert_text((60, y), f"Example {cpd.example_no}: LCMS m/z {cpd.lcms_mz}",
                             fontsize=9)
            y += 16

    doc.save(pdf_path)
    doc.close()
    return layout


def make_scan_pdf(src_pdf: Path, out_pdf: Path, zoom: float = 2.0,
                  seed: int = 7, jpeg_quality: int = 70) -> None:
    """Rasterize a fixture PDF into a scan-like PDF (no text layer): render at
    zoom, add mild sensor noise, JPEG-compress. Deterministic via seed."""
    import cv2
    import numpy as np
    rng = np.random.default_rng(seed)
    doc = pymupdf.open(str(src_pdf))
    out = pymupdf.open()
    for page in doc:
        pix = page.get_pixmap(matrix=pymupdf.Matrix(zoom, zoom), alpha=False)
        arr = np.frombuffer(pix.samples, dtype=np.uint8).reshape(
            pix.height, pix.width, pix.n)[..., :3]
        noise = rng.normal(0.0, 6.0, arr.shape).astype(np.int16)
        arr = np.clip(arr.astype(np.int16) + noise, 0, 255).astype(np.uint8)
        ok, buf = cv2.imencode(".jpg", arr,
                               [int(cv2.IMWRITE_JPEG_QUALITY), jpeg_quality])
        if not ok:
            raise RuntimeError("JPEG encode failed")
        new_page = out.new_page(width=page.rect.width, height=page.rect.height)
        new_page.insert_image(page.rect, stream=buf.tobytes())
    out.save(str(out_pdf))
    out.close()
    doc.close()


if __name__ == "__main__":
    import json
    import sys
    base = Path(__file__).resolve().parents[1]
    data = base / "data"
    if "--scan" in sys.argv:
        # scanned variant of the hard fixture (no text layer, noise + JPEG)
        make_scan_pdf(data / "fixture_patent_hard.pdf",
                      data / "fixture_patent_hard_scan.pdf")
        print(f"wrote {data / 'fixture_patent_hard_scan.pdf'}")
        sys.exit(0)
    if "--hard" in sys.argv:
        compounds = build_dataset_hard()
        layout = make_hard_pdf(compounds, data / "fixtures_png",
                               data / "fixture_patent_hard.pdf")
        out = data / "fixture_ground_truth_hard.json"
    else:
        compounds = build_dataset()
        layout = make_pdf(compounds, data / "fixtures_png", data / "fixture_patent.pdf")
        out = data / "fixture_ground_truth.json"
    truth = [
        {
            "example_no": c.example_no, "kind": c.kind,
            "canonical_smiles": (Chem.MolToSmiles(c.mol) if c.mol else None),
            "ec50_amy3r": c.ec50_amy3r, "emax_amy3r": c.emax_amy3r,
            "ec50_ctr": c.ec50_ctr, "emax_ctr": c.emax_ctr,
            "lcms_mz": c.lcms_mz,
        }
        for c in compounds
    ]
    out.write_text(json.dumps(truth, indent=2))
    print(f"wrote {out} with {len(compounds)} compounds")
