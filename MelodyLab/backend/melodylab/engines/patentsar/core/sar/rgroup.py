"""R-group decomposition via RDKit FMCS (open-source analog of Flare RGA).

Given a congeneric series we find the maximum common substructure (FMCS) = the
core scaffold, then for every molecule split off the substituents attached to
each core position, expressed as R-group SMILES with a single '*' attachment.
The result is an R-group table (compound x position) plus a per-position
inventory, ready for activity correlation (see sar.plots).
"""
from __future__ import annotations

import csv
from dataclasses import dataclass

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFMCS, rdmolops

RDLogger.DisableLog("rdApp.*")


@dataclass
class Compound:
    cid: str
    smiles: str
    activity: float | None = None


def load_dataset(path: str, smiles_col: str, id_col: str | None = None,
                 activity_col: str | None = None) -> list[Compound]:
    rows = list(csv.DictReader(open(path)))
    out = []
    for i, r in enumerate(rows):
        smi = (r.get(smiles_col) or "").strip()
        if not smi:
            continue
        cid = (r.get(id_col) if id_col else None) or f"c{i+1}"
        act = None
        if activity_col and r.get(activity_col) not in (None, "", "ND"):
            try:
                act = float(r[activity_col])
            except ValueError:
                act = None
        out.append(Compound(str(cid), smi, act))
    return out


def find_core(comps: list[Compound], core_smiles: str | None = None) -> Chem.Mol | None:
    """Return the core query mol. If core_smiles (a scaffold with '*' attachment
    points) is given, use it directly (preferred, matches the blog's predefined
    Markush core); otherwise fall back to FMCS over the series."""
    if core_smiles:
        core = Chem.MolFromSmiles(core_smiles) or Chem.MolFromSmarts(core_smiles)
        return core
    mols = [Chem.MolFromSmiles(c.smiles) for c in comps]
    mols = [m for m in mols if m is not None]
    if len(mols) < 2:
        return None
    res = rdFMCS.FindMCS(mols, maximizeBonds=True, ringMatchesRingOnly=True,
                         completeRingsOnly=True)
    return Chem.MolFromSmarts(res.smartsString)


def _build_r_group(mol: Chem.Mol, frag: list[int], attach_atom: int) -> str | None:
    """Build an R-group SMILES: fragment atoms + one '*' at the attachment."""
    rw = Chem.RWMol()
    dummy = rw.AddAtom(Chem.Atom(0))
    mapping: dict[int, int] = {}
    for x in frag:
        a = mol.GetAtomWithIdx(x)
        na = Chem.Atom(a.GetAtomicNum())
        na.SetIsAromatic(a.GetIsAromatic())
        na.SetNoImplicit(True)
        na.SetNumExplicitHs(a.GetTotalNumHs())
        mapping[x] = rw.AddAtom(na)
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in mapping and j in mapping:
            rw.AddBond(mapping[i], mapping[j], b.GetBondType())
    rw.AddBond(dummy, mapping[attach_atom], Chem.BondType.SINGLE)
    try:
        m = rw.GetMol()
        Chem.SanitizeMol(m)
        return Chem.MolToSmiles(m)
    except Exception:
        return None


def decompose(mol: Chem.Mol, core: Chem.Mol) -> dict[int, str]:
    """Return {core_position_index -> R-group SMILES with '*' attachment}."""
    match = mol.GetSubstructMatch(core)
    if not match:
        return {}
    cset = set(match)
    inv = {mol_idx: pos for pos, mol_idx in enumerate(match)}  # mol atom -> core pos

    # connected substituent fragments over non-core atoms
    seen: set[int] = set()
    frags: list[list[int]] = []
    for a in mol.GetAtoms():
        i = a.GetIdx()
        if i in cset or i in seen:
            continue
        stack, comp = [i], []
        while stack:
            x = stack.pop()
            if x in seen or x in cset:
                continue
            seen.add(x)
            comp.append(x)
            for nb in mol.GetAtomWithIdx(x).GetNeighbors():
                if nb.GetIdx() not in cset and nb.GetIdx() not in seen:
                    stack.append(nb.GetIdx())
        frags.append(comp)

    out: dict[int, str] = {}
    for frag in frags:
        # attachment atom = the fragment atom bonded to a core atom; position = core idx
        attach = pos = None
        for x in frag:
            for nb in mol.GetAtomWithIdx(x).GetNeighbors():
                if nb.GetIdx() in inv:
                    attach, pos = x, inv[nb.GetIdx()]
                    break
            if attach is not None:
                break
        if attach is None:
            continue
        rsmi = _build_r_group(mol, frag, attach)
        if rsmi:
            out[pos] = rsmi
    return out


def rgroup_table(comps: list[Compound], core_smiles: str | None = None) -> tuple[list[dict], list[int]]:
    """Return (records with R-group dict per position, ordered positions)."""
    core = find_core(comps, core_smiles)
    if core is None:
        return [], []
    rows = []
    positions: set[int] = set()
    for c in comps:
        mol = Chem.MolFromSmiles(c.smiles)
        groups = decompose(mol, core) if mol else {}
        positions.update(groups.keys())
        rows.append({"cid": c.cid, "smiles": c.smiles, "activity": c.activity,
                     "groups": groups})
    return rows, sorted(positions)


def position_names(positions: list[int]) -> dict[int, str]:
    return {p: f"R{i+1}" for i, p in enumerate(positions)}


def inventory(rows: list[dict], positions: list[int]) -> dict[str, dict[str, list[str]]]:
    """{position_name -> {r_group_smiles: [cids...]}}."""
    names = position_names(positions)
    inv: dict[str, dict[str, list[str]]] = {}
    for p in positions:
        d: dict[str, list[str]] = {}
        for r in rows:
            g = r["groups"].get(p)
            if g:
                d.setdefault(g, []).append(r["cid"])
        inv[names[p]] = d
    return inv


if __name__ == "__main__":
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else "output/fixture_patent_dataset.csv"
    comps = load_dataset(path, smiles_col="canonical_smiles", id_col="example_no",
                         activity_col="ec50_amy3r")
    rows, positions = rgroup_table(comps)
    names = position_names(positions)
    print(f"{len(rows)} compounds | core positions: {[names[p] for p in positions]}")
    for r in rows:
        print(f"  {r['cid']:>4} act={r['activity']}  " +
              "  ".join(f"{names[p]}={r['groups'].get(p,'-')}" for p in positions))
