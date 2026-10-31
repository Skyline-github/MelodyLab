"""R-group activity correlation plots (open-source analog of Flare RGA charts).

boxplot_by_position  -> per-position R-group vs activity (like Flare RGA Boxplot)
combination_heatmap  -> two-position R-group pair activity matrix (like Heatmap)

Activity is EC50 in nM; we convert to pEC50 (-log10 M) so higher = more potent,
which makes the box/heat ordering intuitive.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def to_pactivity(ec50_nm: float | None) -> float | None:
    if ec50_nm is None:
        return None
    try:
        v = float(ec50_nm)
    except (TypeError, ValueError):
        return None
    if v <= 0:
        return None
    return -math.log10(v / 1e9)  # nM -> M -> pEC50


def _trunc(s: str, n: int = 12) -> str:
    return s if len(s) <= n else s[: n - 1] + "…"


def boxplot_by_position(rows: list[dict], positions: list[int], names: dict[int, str],
                        out_dir: str | Path) -> list[Path]:
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    made = []
    for p in positions:
        groups: dict[str, list[float]] = {}
        for r in rows:
            g = r["groups"].get(p)
            pa = to_pactivity(r["activity"])
            if g and pa is not None:
                groups.setdefault(g, []).append(pa)
        groups = {g: v for g, v in groups.items() if len(v) >= 1}
        if len(groups) < 2:
            continue
        # keep the most populated groups (then order by potency) so the plot stays
        # readable on positions with many singletons
        groups = dict(sorted(groups.items(), key=lambda kv: -len(kv[1]))[:12])
        labels = sorted(groups, key=lambda g: -np.mean(groups[g]))
        data = [groups[g] for g in labels]
        fig, ax = plt.subplots(figsize=(max(8, 0.9 * len(labels)), 5.2))
        bp = ax.boxplot(data, tick_labels=[f"{_trunc(g, 16)}\n(n={len(groups[g])})" for g in labels],
                        patch_artist=True, widths=0.6)
        for pc in bp["boxes"]:
            pc.set_facecolor("#cfe3ef")
        ax.set_ylabel("pEC50 (higher = more potent)", fontsize=11)
        ax.set_title(f"RGA boxplot - position {names[p]}", fontsize=13)
        ax.tick_params(axis="x", labelsize=9)
        ax.tick_params(axis="y", labelsize=10)
        ax.grid(axis="y", alpha=0.3)
        fig.autofmt_xdate(rotation=35)
        fig.tight_layout()
        out = out_dir / f"rga_box_{names[p]}.png"
        fig.savefig(out, dpi=120); plt.close(fig)
        made.append(out)
    return made


def combination_heatmap(rows: list[dict], pos_a: int, pos_b: int, names: dict[int, str],
                        out_dir: str | Path) -> Path | None:
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    cells: dict[tuple[str, str], list[float]] = {}
    for r in rows:
        a, b = r["groups"].get(pos_a), r["groups"].get(pos_b)
        pa = to_pactivity(r["activity"])
        if a and b and pa is not None:
            cells.setdefault((a, b), []).append(pa)
    if len(cells) < 2:
        return None
    xs = sorted({k[0] for k in cells})
    ys = sorted({k[1] for k in cells})
    grid = np.full((len(ys), len(xs)), np.nan)
    for i, y in enumerate(ys):
        for j, x in enumerate(xs):
            if (x, y) in cells:
                grid[i, j] = np.mean(cells[(x, y)])
    fig, ax = plt.subplots(figsize=(max(5, 0.7 * len(xs)), max(4, 0.6 * len(ys))))
    im = ax.imshow(grid, aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(xs)), [_trunc(x) for x in xs], rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(ys)), [_trunc(y) for y in ys], fontsize=7)
    ax.set_xlabel(names[pos_a]); ax.set_ylabel(names[pos_b])
    ax.set_title(f"RGA heatmap - {names[pos_a]} x {names[pos_b]} (mean pEC50)")
    fig.colorbar(im, ax=ax, label="pEC50")
    fig.tight_layout()
    out = out_dir / f"rga_heat_{names[pos_a]}x{names[pos_b]}.png"
    fig.savefig(out, dpi=120); plt.close(fig)
    return out


def _render_rgroup(rsmi: str, out_png: Path, size: int = 90) -> bool:
    """Draw an R-group (SMILES with a single '*' attachment) to a small PNG."""
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D
    mol = Chem.MolFromSmiles(rsmi)
    if mol is None:
        return False
    d = rdMolDraw2D.MolDraw2DCairo(size, size)
    d.drawOptions().addStereoAnnotation = False
    d.drawOptions().padding = 0.02
    try:
        d.DrawMolecule(mol); d.FinishDrawing()
        out_png.write_bytes(d.GetDrawingText())
        return True
    except Exception:
        return False


def combination_heatmap_structural(rows: list[dict], pos_a: int, pos_b: int,
                                   names: dict[int, str], out_dir: str | Path,
                                   max_axis: int = 8) -> Path | None:
    """Like combination_heatmap but labels the axes with drawn R-group structures
    instead of SMILES text (much more intuitive)."""
    from matplotlib.offsetbox import OffsetImage, AnnotationBbox
    from PIL import Image
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    tmp = out_dir / "_rimgs"; tmp.mkdir(exist_ok=True)

    cells: dict[tuple[str, str], list[float]] = {}
    for r in rows:
        a, b = r["groups"].get(pos_a), r["groups"].get(pos_b)
        pa = to_pactivity(r["activity"])
        if a and b and pa is not None:
            cells.setdefault((a, b), []).append(pa)
    if len(cells) < 2:
        return None
    # rank axis groups by mean potency so the informative ones are shown
    def order(idx):
        agg: dict[str, list[float]] = {}
        for (a, b), v in cells.items():
            agg.setdefault(a if idx == 0 else b, []).extend(v)
        keys = sorted(agg, key=lambda k: -np.mean(agg[k]))[:max_axis]
        return keys
    xs, ys = order(0), order(1)
    grid = np.full((len(ys), len(xs)), np.nan)
    for i, y in enumerate(ys):
        for j, x in enumerate(xs):
            if (x, y) in cells:
                grid[i, j] = np.mean(cells[(x, y)])

    fig, ax = plt.subplots(figsize=(1.05 * len(xs) + 2.2, 1.05 * len(ys) + 2.2))
    im = ax.imshow(grid, aspect="auto", cmap="viridis", vmin=np.nanmin(grid), vmax=np.nanmax(grid))
    ax.set_xticks(range(len(xs))); ax.set_yticks(range(len(ys)))
    ax.set_xticklabels([]); ax.set_yticklabels([])
    for i, y in enumerate(ys):
        for j, x in enumerate(xs):
            if not np.isnan(grid[i, j]):
                ax.text(j, i, f"{grid[i, j]:.1f}", ha="center", va="center",
                        color="white", fontsize=8)

    def place(rsmi, xy, xybox, ha):
        p = tmp / f"rg{abs(hash(rsmi)) % 10**9}.png"
        if not _render_rgroup(rsmi, p):
            return
        img = OffsetImage(Image.open(p), zoom=0.9)
        ab = AnnotationBbox(img, xy, xybox=xybox, frameon=False,
                           box_alignment=(0.5, 0.5) if ha == "x" else (1.0, 0.5),
                           pad=0.2)
        ax.add_artist(ab)

    for j, x in enumerate(xs):
        place(x, (j, -0.5), (j, -1.1), "x")
    for i, y in enumerate(ys):
        place(y, (-0.5, i), (-1.1, i), "y")

    ax.set_xlim(-0.5, len(xs) - 0.5); ax.set_ylim(len(ys) - 0.5, -0.5)
    ax.set_title(f"RGA heatmap — {names[pos_a]} (x) × {names[pos_b]} (y), mean pEC50")
    fig.colorbar(im, ax=ax, label="pEC50", fraction=0.046, pad=0.15)
    fig.subplots_adjust(left=0.28, bottom=0.22, top=0.9, right=0.86)
    out = out_dir / f"rga_heat_{names[pos_a]}x{names[pos_b]}.png"
    fig.savefig(out, dpi=120); plt.close(fig)
    return out

