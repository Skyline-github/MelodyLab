"""CB-Dock2 风格的几何口袋检测（纯 numpy 实现，无外部服务依赖）。

思路：在蛋白表面附近采样网格点，保留「凹腔」点——
距 ≥N 个原子的探针距离内、但不陷入任何原子内部的点；
对候选点做密度聚类，每个簇给出中心、估计半径与埋藏度评分。
"""
from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np

VDW = {  # 简化范德华半径 (Å)
    "C": 1.70, "A": 1.70, "N": 1.55, "O": 1.52, "S": 1.80,
    "P": 1.80, "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98,
    "ZN": 1.25, "FE": 1.25, "MG": 1.25, "CA": 1.25, "MN": 1.25,
    "H": 0.00, "D": 0.00,
}


@dataclass
class Pocket:
    id: str
    score: float
    center_x: float
    center_y: float
    center_z: float
    radius: float
    residues: str = ""
    n_points: int = 0

    def as_json(self) -> dict:
        return asdict(self)


def parse_pdb_atoms(path: Path, keep_hetatm: bool = False, remove_waters: bool = True):
    atoms = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            if line.startswith("HETATM") and not keep_hetatm:
                continue
            record = line[:6].strip()
            altloc = line[16]
            if altloc not in (" ", "A"):
                continue
            resname = line[17:20].strip()
            if remove_waters and resname in ("HOH", "WAT", "TIP3", "TIP3P", "SOL", "SPC"):
                continue
            element = line[76:78].strip().upper() or line[12:14].strip().strip("0123456789").upper()
            if element in ("H", "D"):
                continue  # 重原子模型
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            chain = line[21].strip()
            resid = line[22:26].strip()
            atoms.append({"element": element, "xyz": (x, y, z),
                          "residue": f"{chain}{resid}:{resname}" if chain else f"{resid}:{resname}"})
    return atoms


def detect_pockets(
    protein_path: Path,
    *,
    grid_spacing: float = 2.0,
    probe_radius: float = 1.8,
    burial_radius: float = 6.0,
    burial_min_neighbors: int = 8,
    min_radius: float = 2.2,
    max_radius: float = 8.0,
    max_pockets: int = 5,
    min_pocket_separation: float = 10.0,
) -> list[Pocket]:
    atoms = parse_pdb_atoms(protein_path)
    if len(atoms) < 20:
        raise ValueError(f"蛋白原子过少（{len(atoms)}），文件可能无法解析: {protein_path.name}")
    coords = np.array([a["xyz"] for a in atoms])
    radii = np.array([VDW.get(a["element"], 1.7) for a in atoms])

    # 1. 在包围盒内布网格（留 padding）
    pad = probe_radius + 2.0
    lo, hi = coords.min(0) - pad, coords.max(0) + pad
    axes = [np.arange(lo[i], hi[i] + grid_spacing, grid_spacing) for i in range(3)]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, 3)

    # 2. 每个格点到最近重原子的距离 & 邻域内原子数
    from scipy.spatial import cKDTree

    atom_tree = cKDTree(coords)
    d_near, i_near = atom_tree.query(grid, k=1)
    inside = d_near < radii[i_near]  # 原子内部，丢弃
    neighbors = np.array([len(ix) for ix in atom_tree.query_ball_point(grid[~inside], burial_radius)]) if (~inside).any() else np.array([])
    free = grid[~inside]
    concave = neighbors >= burial_min_neighbors  # 埋藏度：周围原子多 → 凹腔
    pts = free[concave]
    if len(pts) == 0:
        return []

    # 3. 密度聚类（贪心：按局部密度从高到低收簇）
    order = np.argsort(-neighbors[concave])
    pts_sorted = pts[order]
    assigned = np.zeros(len(pts_sorted), dtype=bool)
    clusters = []
    sep2 = min_pocket_separation ** 2
    for idx in order:
        if assigned[idx]:
            continue
        p = pts[idx]
        d2 = ((pts_sorted - p) ** 2).sum(1)
        members = np.where((d2 < sep2) & (~assigned))[0]
        assigned[members] = True
        if len(members) >= 4:
            center = pts[members].mean(0)
            clusters.append((center, len(members)))
    if not clusters:
        return []

    # 4. 评分：簇大小 × 中心埋藏深度；半径按点云范围估计
    pockets: list[Pocket] = []
    for rank, (center, size) in enumerate(sorted(clusters, key=lambda c: -c[1])[: max_pockets * 3]):
        near_mask = (((pts - center) ** 2).sum(1) < sep2)
        cloud = pts[near_mask]
        if len(cloud) == 0:
            continue
        extent = np.sqrt(((cloud - center) ** 2).sum(1).mean())
        radius = float(np.clip(max(extent * 1.2, min_radius), min_radius, max_radius))
        d_center, _ = atom_tree.query(center, k=1)
        score = float(len(cloud) * (probe_radius - min(d_center, probe_radius) + 0.5))
        pockets.append(Pocket(
            id=f"pocket_{rank+1}", score=round(score, 2),
            center_x=round(float(center[0]), 2), center_y=round(float(center[1]), 2),
            center_z=round(float(center[2]), 2), radius=radius, n_points=int(len(cloud)),
        ))
        if len(pockets) >= max_pockets:
            break

    # 5. 相近口袋去重（中心距离过近保留高分）
    deduped: list[Pocket] = []
    for p in sorted(pockets, key=lambda x: -x.score):
        if all(((p.center_x - q.center_x) ** 2 + (p.center_y - q.center_y) ** 2 + (p.center_z - q.center_z) ** 2) ** 0.5 >= min_pocket_separation * 0.6 for q in deduped):
            deduped.append(p)
    for i, p in enumerate(deduped[:max_pockets]):
        p.id = f"pocket_{i+1}"
    return deduped[:max_pockets]


def pocket_residues(protein_path: Path, pocket: Pocket, cutoff: float | None = None) -> list[str]:
    cutoff = cutoff or pocket.radius + 3.0
    atoms = parse_pdb_atoms(protein_path)
    c = np.array([pocket.center_x, pocket.center_y, pocket.center_z])
    residues, seen = [], set()
    for a in atoms:
        d = np.linalg.norm(np.array(a["xyz"]) - c)
        if d <= cutoff and a["residue"] not in seen:
            seen.add(a["residue"])
            residues.append(a["residue"])
    return residues


def write_pockets_json(out_json: Path, protein_path: Path, pockets: list[Pocket]) -> None:
    payload = {
        "protein": str(protein_path.name),
        "method": "melodylab_grid_cavity",
        "pockets": [
            {**p.as_json(), "residues": ";".join(pocket_residues(protein_path, p))[:2000]}
            for p in pockets
        ],
    }
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
