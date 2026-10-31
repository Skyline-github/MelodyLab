"""Field-based 3D-QSAR (open-source CoMFA-style analog of Flare's FIFST).

Pipeline:
  1. build a low-energy 3D conformer per molecule (ETKDG + MMFF94)
  2. align all molecules by their common scaffold (FMCS core atom map)
  3. compute steric (Lennard-Jones) + electrostatic (Coulomb/Gasteiger) fields
     on a 3D grid around the aligned set
  4. PLS regression vs pEC50, LOO cross-validation -> R2, Q2, RMSE
  5. report the strongest steric/electrostatic coefficient regions

This is an educational analog: it uses a single aligned conformer and a coarse
grid, so it will not match Flare's force-field/FEP-based models. It does
demonstrate the full field-based 3D-QSAR workflow end to end on CPU.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdFMCS, rdMolAlign, rdPartialCharges
from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler

RDLogger.DisableLog("rdApp.*")


@dataclass
class Model3D:
    cids: list[str]
    y: np.ndarray            # pEC50
    X: np.ndarray            # n x features (steric then electrostatic)
    grid: np.ndarray         # npts x 3
    n_steric: int
    stats: dict


def _conformer(smiles: str, seed: int = 0xF00D):
    mol0 = Chem.MolFromSmiles(smiles)
    if mol0 is None:
        return None
    mol = Chem.AddHs(mol0)
    if mol.GetNumAtoms() == 0:
        return None
    ps = AllChem.ETKDGv3(); ps.randomSeed = seed
    if AllChem.EmbedMolecule(mol, ps) != 0:
        if AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=seed) != 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        pass
    return mol


def _core(mols):
    return Chem.MolFromSmarts(rdFMCS.FindMCS(
        [m for m in mols], maximizeBonds=True, ringMatchesRingOnly=True,
        completeRingsOnly=True).smartsString)


def build_aligned(smiles_list):
    """Return (aligned mols, kept_indices) aligned to the first molecule's frame.

    kept_indices map back into smiles_list so callers keep cid/activity aligned.
    """
    mols, kept = [], []
    for i, s in enumerate(smiles_list):
        m = _conformer(s)
        if m is not None:
            mols.append(m); kept.append(i)
    if len(mols) < 2:
        return mols, kept
    heavy = [Chem.RemoveHs(Chem.Mol(m)) for m in mols]
    core = _core(heavy)
    ref = mols[0]
    ref_map = ref.GetSubstructMatch(core) if core else ()
    if not ref_map:
        return mols, kept
    aligned = [ref]
    for m, h in zip(mols[1:], heavy[1:]):
        hmatch = h.GetSubstructMatch(core)
        if not hmatch or len(hmatch) != len(ref_map):
            aligned.append(m)
            continue
        atom_map = [(hmatch[p], ref_map[p]) for p in range(len(ref_map))]
        try:
            rdMolAlign.AlignMol(m, ref, atomMap=atom_map)
        except Exception:
            pass
        aligned.append(m)
    return aligned, kept


def _gasteiger(mol):
    cp = Chem.Mol(mol)
    try:
        rdPartialCharges.ComputeGasteigerCharges(cp)
        return np.array([float(a.GetProp("_GasteigerCharge")) for a in cp.GetAtoms()])
    except Exception:
        return np.zeros(cp.GetNumAtoms())


def compute_fields(mols, spacing=2.0, margin=2.0, cutoff=10.0):
    coords = np.vstack([m.GetConformer().GetPositions() for m in mols])
    lo = coords.min(0) - margin
    hi = coords.max(0) + margin
    gx = np.arange(lo[0], hi[0] + spacing, spacing)
    gy = np.arange(lo[1], hi[1] + spacing, spacing)
    gz = np.arange(lo[2], hi[2] + spacing, spacing)
    grid = np.stack(np.meshgrid(gx, gy, gz, indexing="ij"), -1).reshape(-1, 3)
    npts = len(grid)
    steric = np.zeros((len(mols), npts))
    elec = np.zeros((len(mols), npts))
    for mi, m in enumerate(mols):
        pos = m.GetConformer().GetPositions()
        q = _gasteiger(m)
        d = np.linalg.norm(pos[:, None, :] - grid[None, :, :], axis=2)  # atoms x pts
        d = np.clip(d, 0.8, None)
        with np.errstate(divide="ignore", over="ignore"):
            s = 0.36 / d**12 - 0.15 / d**10
            e = 332.0 * q[:, None] / d
        s = np.clip(s, 0, cutoff)
        e = np.clip(e, -cutoff, cutoff)
        steric[mi] = s.sum(0)
        elec[mi] = e.sum(0)
    # drop all-zero columns
    keep = (np.abs(steric).sum(0) + np.abs(elec).sum(0)) > 0
    return grid[keep], steric[:, keep], elec[:, keep]


def fit_pls(X, y, n_comp=2):
    n = len(y)
    # keep grid columns that are (a) non-trivially variable and (b) populated in
    # most molecules; sparse/constant cells make LOO PLS extrapolate wildly.
    occ = (np.abs(X) > 1e-6).mean(axis=0) >= 0.5
    keep = occ & (X.std(axis=0) > 1e-2)
    Xf = X[:, keep]
    if Xf.shape[1] == 0:
        return {"n": n, "n_comp": 0, "R2": 0.0, "Q2_loo": 0.0, "RMSE_fit": 0.0,
                "RMSE_loo": 0.0, "tau_loo": 0.0}, np.zeros(X.shape[1]), np.zeros(n)
    n_comp = max(1, min(n_comp, n - 1, Xf.shape[1]))
    # LOO cross-validation for Q2 (PRESS)
    press = 0.0
    yhat_cv = np.zeros(n)
    for i in range(n):
        tr = [j for j in range(n) if j != i]
        sc = StandardScaler().fit(Xf[tr])
        p = PLSRegression(n_components=n_comp).fit(sc.transform(Xf[tr]), y[tr])
        yhat_cv[i] = float(np.asarray(p.predict(sc.transform(Xf[i:i+1]))).ravel()[0])
    press = float(((y - yhat_cv) ** 2).sum())
    ss = float(((y - y.mean()) ** 2).sum())
    q2 = 1 - press / ss if ss > 0 else 0.0
    rmse_loo = math.sqrt(press / n)
    # full-data model for R2 and coefficients
    sc = StandardScaler().fit(Xf)
    pls = PLSRegression(n_components=n_comp).fit(sc.transform(Xf), y)
    yfit = np.asarray(pls.predict(sc.transform(Xf))).ravel()
    r2 = 1 - float(((y - yfit) ** 2).sum()) / ss if ss > 0 else 0.0
    rmse = float(np.sqrt(((y - yfit) ** 2).mean()))
    # coefficients back on original (unscaled) features, expanded to full width
    coef_f = pls.coef_.ravel() / sc.scale_
    coef = np.zeros(X.shape[1])
    coef[keep] = coef_f
    tau = float(_kendall_tau(y, yhat_cv))
    return {"n": n, "n_comp": n_comp, "R2": round(r2, 3), "Q2_loo": round(q2, 3),
            "RMSE_fit": round(rmse, 3), "RMSE_loo": round(rmse_loo, 3),
            "tau_loo": round(tau, 3)}, coef, yhat_cv


def _kendall_tau(a, b):
    n = len(a)
    conc = disc = 0
    for i in range(n):
        for j in range(i + 1, n):
            s = (a[i] - a[j]) * (b[i] - b[j])
            if s > 0:
                conc += 1
            elif s < 0:
                disc += 1
    tot = n * (n - 1) / 2
    return (conc - disc) / tot if tot else 0.0


def top_coefficient_regions(model: Model3D, k=8):
    """Return the grid points with the largest |field coefficients|."""
    n = model.n_steric
    cs = model.stats.get("_coef_steric", np.zeros(n))
    ce = model.stats.get("_coef_elec", np.zeros(n))
    pts = []
    for idx in range(n):
        pts.append((abs(cs[idx]) + abs(ce[idx]), model.grid[idx], cs[idx], ce[idx]))
    pts.sort(key=lambda t: -t[0])
    return pts[:k]


def build_model(cids, smiles, pec50, spacing=2.0, n_comp=2) -> Model3D:
    aligned, kept = build_aligned(smiles)
    if len(aligned) < 3:
        raise ValueError("need >=3 aligned 3D structures")
    cids = [cids[i] for i in kept]
    pec50 = [pec50[i] for i in kept]
    grid, steric, elec = compute_fields(aligned, spacing=spacing)
    X = np.hstack([steric, elec])
    y = np.array(pec50, dtype=float)
    stats, coef, yhat_cv = fit_pls(X, y, n_comp=n_comp)
    stats["_yhat_loo"] = yhat_cv
    stats["_coef_steric"] = coef[:steric.shape[1]]
    stats["_coef_elec"] = coef[steric.shape[1]:]
    return Model3D(cids, y, X, grid, steric.shape[1], stats)


def free_wilson(rows: list[dict], positions: list[int], n_comp=None):
    """Additive Free-Wilson QSAR from the R-group table (dummy-coded).

    rows: from rgroup.rgroup_table (each has 'groups' {pos->Rsmiles} and
    'activity' in nM). Returns (stats dict, coefficient list of (pos, R, coef)).
    This is the Free-Wilson / AIR-style analysis the blog references; it is
    robust on congeneric series where 3D fields are noisy.
    """
    import math
    from sklearn.linear_model import LinearRegression
    from sklearn.model_selection import LeaveOneOut, cross_val_predict

    # collect distinct R-groups per position
    cats: dict[int, list[str]] = {}
    for p in positions:
        cats[p] = sorted({r["groups"].get(p) for r in rows if r["groups"].get(p)})
    # design matrix: one column per (position, R) except drop first per position (ref)
    cols = []
    for p in positions:
        for g in cats[p][1:]:
            cols.append((p, g))
    X, y = [], []
    for r in rows:
        pa = _pact(r["activity"])
        if pa is None:
            continue
        vec = [1.0 if r["groups"].get(p) == g else 0.0 for (p, g) in cols]
        X.append(vec)
        y.append(pa)
    X = np.array(X); y = np.array(y)
    n = len(y)
    if n < len(cols) + 2:
        return {"n": n, "note": "too few compounds for Free-Wilson"}, []
    lr = LinearRegression()
    yfit = lr.fit(X, y).predict(X)
    ss = float(((y - y.mean()) ** 2).sum())
    r2 = 1 - float(((y - yfit) ** 2).sum()) / ss if ss > 0 else 0.0
    ycv = cross_val_predict(LinearRegression(), X, y, cv=LeaveOneOut())
    press = float(((y - ycv) ** 2).sum())
    q2 = 1 - press / ss if ss > 0 else 0.0
    stats = {"n": n, "n_terms": len(cols), "R2": round(r2, 3), "Q2_loo": round(q2, 3),
             "RMSE_fit": round(math.sqrt(((y - yfit) ** 2).mean()), 3),
             "RMSE_loo": round(math.sqrt(press / n), 3)}
    coefs = [((p, g), round(float(lr.coef_[i]), 3)) for i, (p, g) in enumerate(cols)]
    return stats, coefs


def _pact(ec50_nm):
    try:
        v = float(ec50_nm)
    except (TypeError, ValueError):
        return None
    return -math.log10(v / 1e9) if v > 0 else None  # nM -> M -> pEC50
