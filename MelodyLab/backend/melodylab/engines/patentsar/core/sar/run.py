"""SAR CLI: R-group analysis + field-based 3D-QSAR on a mined dataset.

Input is any CSV produced by the mining pipeline (needs a SMILES column and an
activity column, e.g. ec50_amy3r in nM). Outputs:
  - rgroup_table.csv   : compound x R-group-position matrix
  - rga_box_*.png      : per-position activity boxplots
  - rga_heat_*.png     : two-position combination heatmap (if applicable)
  - qsar3d_coeff.png   : top steric/electrostatic coefficient field points
  - sar_report.md      : R-group inventory + 3D-QSAR stats (Flare Table-1 style)

Usage:
  python -m sar.run --dataset output/fixture_patent_dataset.csv \
      --smiles-col canonical_smiles --id-col example_no --activity-col ec50_amy3r \
      --target AMY3R --outdir output/sar
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from . import plots, q3d, rgroup


def render_core(core_smiles, rows, out_png, positions=None, names=None):
    """Draw the RGA core scaffold with each variable attachment point labelled
    R1, R2, ... (a dummy stub carrying the position name), so the母核 is intuitive."""
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Chem import rdDepictor
    if core_smiles:
        q = Chem.MolFromSmiles(core_smiles) or Chem.MolFromSmarts(core_smiles)
    else:
        mols = [Chem.MolFromSmiles(r["smiles"]) for r in rows if r["smiles"]]
        mols = [m for m in mols if m]
        from rdkit.Chem import rdFMCS
        q = Chem.MolFromSmarts(rdFMCS.FindMCS(mols, maximizeBonds=True,
                          ringMatchesRingOnly=True, completeRingsOnly=True).smartsString)
    if q is None:
        raise ValueError("cannot build core mol")
    rw = Chem.RWMol()
    mapping = {}
    for a in q.GetAtoms():
        na = Chem.Atom(a.GetAtomicNum())
        na.SetIsAromatic(a.GetIsAromatic())
        mapping[a.GetIdx()] = rw.AddAtom(na)
    for b in q.GetBonds():
        rw.AddBond(mapping[b.GetBeginAtomIdx()], mapping[b.GetEndAtomIdx()],
                   b.GetBondType())
    # add a labelled dummy stub at each variable attachment position
    if positions and names:
        for p in positions:
            if p not in mapping:
                continue
            dummy = Chem.Atom(0)
            dummy.SetProp("dummyLabel", names[p])
            didx = rw.AddAtom(dummy)
            rw.AddBond(mapping[p], didx, Chem.BondType.SINGLE)
    m = rw.GetMol()
    try:
        Chem.SanitizeMol(m)
    except Exception:
        pass
    rdDepictor.Compute2DCoords(m)
    d = rdMolDraw2D.MolDraw2DCairo(760, 520)
    d.drawOptions().addStereoAnnotation = False
    d.DrawMolecule(m)
    d.FinishDrawing()
    Path(out_png).write_bytes(d.GetDrawingText())


def write_rgroup_csv(rows, positions, names, path: Path):
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "activity_ec50_nm"] + [names[p] for p in positions] + ["smiles"])
        for r in rows:
            w.writerow([r["cid"], r["activity"]] +
                       [r["groups"].get(p, "") for p in positions] + [r["smiles"]])


def plot_coefficients(model: q3d.Model3D, out: Path, k=60):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cs, ce = model.stats["_coef_steric"], model.stats["_coef_elec"]
    mag = np.abs(cs) + np.abs(ce)
    idx = np.argsort(-mag)[:k]
    pts = model.grid[idx]
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")
    sc = cs[idx]
    colors = np.where(sc >= 0, "g", "r")
    ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], c=colors, marker="o", label="steric")
    ec = ce[idx]
    ecolors = np.where(ec >= 0, "b", "m")
    ax.scatter(pts[:, 0] + 0.3, pts[:, 1], pts[:, 2], c=ecolors, marker="^", label="electrostatic")
    ax.set_title("3D-QSAR top coefficient fields\n(g/r steric +/-, b/m electrostatic +/-)")
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    fig.tight_layout(); fig.savefig(out, dpi=120); plt.close(fig)


def run(dataset, smiles_col, id_col, activity_col, target, outdir,
        spacing=2.0, n_comp=2, core=None):
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    comps = rgroup.load_dataset(dataset, smiles_col, id_col, activity_col)
    report = [f"# SAR report - {target}", "", f"- dataset: `{dataset}`",
              f"- core: {core or 'FMCS auto'}",
              f"- compounds with SMILES: {len(comps)}", ""]

    # R-group analysis
    rows, positions = rgroup.rgroup_table(comps, core_smiles=core)
    names = rgroup.position_names(positions)
    write_rgroup_csv(rows, positions, names, outdir / "rgroup_table.csv")
    inv = rgroup.inventory(rows, positions)
    report += ["## R-group analysis", "",
               f"- core positions detected (FMCS): {[names[p] for p in positions]}", ""]
    for pos_name, groups in inv.items():
        report.append(f"### {pos_name} ({len(groups)} distinct R-groups)")
        for g, cids in sorted(groups.items(), key=lambda kv: -len(kv[1]))[:15]:
            report.append(f"- `{g}` : n={len(cids)}  {', '.join(cids[:8])}")
        report.append("")
    made = plots.boxplot_by_position(rows, positions, names, outdir)
    report += [f"- RGA boxplots written: {len(made)}", ""]
    # clear stale heatmaps so the report shows only current, informative pairs
    for old in outdir.glob("rga_heat_*.png"):
        old.unlink()
    if len(positions) >= 2:
        # rank positions by number of distinct R-groups with activity data
        pop = sorted(positions, key=lambda p: -len(inv.get(names[p], {})))
        pairs = [(pop[0], pop[1]), (pop[0], pop[2])] if len(pop) >= 3 else [(pop[0], pop[1])]
        for pa, pb in pairs:
            hm = plots.combination_heatmap_structural(rows, pa, pb, names, outdir) \
                or plots.combination_heatmap(rows, pa, pb, names, outdir)
            if hm:
                report.append(f"- RGA heatmap: {hm.name}")
    # render the core scaffold used for RGA, with R1..Rn attachment labels
    try:
        render_core(core, rows, outdir / "core_scaffold.png", positions=positions, names=names)
        report.append("- core scaffold: core_scaffold.png")
    except Exception as e:
        report.append(f"- core render skipped: {e}")

    # Free-Wilson additive model (from the R-group table)
    report += ["", "## Free-Wilson additive QSAR (from R-group table)", ""]
    if positions:
        try:
            st, coefs = q3d.free_wilson(rows, positions)
            if "note" in st:
                report.append(f"- skipped: {st['note']} (n={st['n']})")
            else:
                report += ["| Terms | R2 | Q2(LOO) | RMSE_fit | RMSE_loo |",
                           "|---|---|---|---|---|",
                           f"| {st['n_terms']} | {st['R2']} | {st['Q2_loo']} "
                           f"| {st['RMSE_fit']} | {st['RMSE_loo']} |", "",
                           f"- n = {st['n']} compounds; R-group contributions "
                           "(pEC50, vs each position's reference R):"]
                for (p, g), c in sorted(coefs, key=lambda t: -t[1]):
                    report.append(f"  - {names[p]} `{g}` : {c:+.3f}")
        except Exception as e:
            report.append(f"- Free-Wilson skipped: {type(e).__name__}: {e}")

    # 3D-QSAR
    labeled = [c for c in comps if plots.to_pactivity(c.activity) is not None]
    report += ["", "## Field-based 3D-QSAR (open-source analog)", ""]
    if len(labeled) >= 4:
        try:
            model = q3d.build_model([c.cid for c in labeled], [c.smiles for c in labeled],
                                    [plots.to_pactivity(c.activity) for c in labeled],
                                    spacing=spacing, n_comp=n_comp)
            st = model.stats
            report += [
                "| Components | R2 | Q2(LOO) | RMSE_fit | RMSE_loo | tau_loo |",
                "|---|---|---|---|---|---|",
                f"| {st['n_comp']} | {st['R2']} | {st['Q2_loo']} | {st['RMSE_fit']} "
                f"| {st['RMSE_loo']} | {st['tau_loo']} |", "",
                f"- n = {st['n']} compounds, {model.X.shape[1]} field features "
                f"({len(model.grid)} grid pts x steric+electrostatic)",
                "- NOTE: single aligned conformer, coarse grid; a functional analog of "
                "Flare Field-based QSAR, not a validated predictive model.",
            ]
            plot_coefficients(model, outdir / "qsar3d_coeff.png")
            report.append(f"- coefficient field plot: qsar3d_coeff.png")
        except Exception as e:
            report.append(f"- 3D-QSAR skipped: {type(e).__name__}: {e}")
    else:
        report.append(f"- skipped: need >=4 compounds with numeric {activity_col} "
                      f"(have {len(labeled)})")

    (outdir / "sar_report.md").write_text("\n".join(report))
    return {"n": len(comps), "positions": [names[p] for p in positions],
            "boxplots": len(made), "outdir": outdir}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--smiles-col", default="canonical_smiles")
    ap.add_argument("--id-col", default="example_no")
    ap.add_argument("--activity-col", default="ec50_amy3r")
    ap.add_argument("--target", default="AMY3R")
    ap.add_argument("--outdir", default="output/sar")
    ap.add_argument("--spacing", type=float, default=2.5)
    ap.add_argument("--n-comp", type=int, default=3)
    ap.add_argument("--core", default=None,
                    help="scaffold SMILES (de-substituted) to define R-group "
                         "attachment points; defaults to FMCS auto-detection")
    ap.add_argument("--summarize", action="store_true",
                    help="after the RGA/3D-QSAR run, also write a narrative "
                         "sar_summary.md/json (sar.summarize). Uses the LLM "
                         "configured via SAR_LLM_BASE_URL/SAR_LLM_MODEL/"
                         "SAR_LLM_API_KEY, or a deterministic digest if no "
                         "API key is set")
    ap.add_argument("--summary-lang", choices=["zh", "en"], default="zh")
    a = ap.parse_args()
    outdir = Path(a.outdir)
    res = run(a.dataset, a.smiles_col, a.id_col, a.activity_col, a.target,
              outdir, a.spacing, a.n_comp, core=a.core)
    print(f"compounds={res['n']} | RGA positions={res['positions']} | "
          f"boxplots={res['boxplots']} | -> {res['outdir']}/sar_report.md")
    if a.summarize:
        from . import summarize as _sum
        sres = _sum.summarize(a.dataset, a.target, str(outdir),
                              rgroup_table=str(outdir / "rgroup_table.csv"),
                              activity_col=a.activity_col, id_col=a.id_col,
                              lang=a.summary_lang)
        print(f"summary: {len(sres['claims'])} claims verified, "
              f"{len(sres['quarantined'])} quarantined | -> "
              f"{outdir}/sar_summary.{{json,md}}")


if __name__ == "__main__":
    main()
