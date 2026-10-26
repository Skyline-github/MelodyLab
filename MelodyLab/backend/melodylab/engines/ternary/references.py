"""DeepTernary 输入制备：参照晶体加载 + PROTAC→三元复合物输入装配。

装配协议（算法源自公开研究方法）：
  1. PROTAC 全分子 ETKDG 构象（seed 42）；
  2. POI 侧：MCS 对齐到 TD2 晶体姿态（5E89，galectin-3C 链 A）→ 配体进入 gal3 坐标系；
  3. E3 侧：晶体 E3 配体 MCS 对齐到 PROTAC 的 E3 warhead（排除与 POI warhead 重叠的
     错误匹配），E3 蛋白随之变换，锚点外推到连接子伸展距离，随机扭转采样至无碰撞；
  4. 写出 lig.pdb / p1.pdb / p2.pdb / unbound_lig1.pdb / unbound_lig2.pdb /
     lig1_mask.pdb / lig2_mask.pdb（predict.py 的原子替换需要 <1 Å 对齐，
     故 warhead 原子硬锚定晶体坐标）。
"""
from __future__ import annotations

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdFMCS
from rdkit.Geometry import Point3D

from . import config as cfg

RDLogger.DisableLog("rdApp.*")


# ---------------------------------------------------------------- 基础工具（原样移植）
def kabsch(P, Q):
    """Return R, t so that (R @ P.T).T + t ≈ Q. P, Q: (n,3)."""
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    t = Q.mean(0) - R @ P.mean(0)
    return R, t


def read_ligand_xtal(pdb_path, resname, chain=None, template=None):
    """提取 HETATM 残基为带晶体坐标的 RDKit mol；键由邻近重建，键级从 ideal SDF 模板转移。"""
    lines = []
    for line in open(pdb_path):
        if line.startswith("HETATM") and line[17:20].strip() == resname:
            if chain is None or line[21] == chain:
                lines.append(line)
    block = "REMARK   crystal pose\n" + "".join(lines) + "END\n"
    mol = Chem.MolFromPDBBlock(block, sanitize=True, removeHs=True,
                               proximityBonding=True)
    if template is not None and mol is not None:
        try:
            mol = AllChem.AssignBondOrdersFromTemplate(
                Chem.RemoveHs(Chem.Mol(template)), mol)
            Chem.SanitizeMol(mol)
        except Exception as e:  # noqa: BLE001 —— 键级转移失败时保留邻近建键结果
            print(f"  [warn] bond-order assignment failed for {resname}: {e}")
    return mol, len(lines)


def protein_chains_pdb(pdb_path, chains):
    out = []
    for line in open(pdb_path):
        if line.startswith("ATOM") and line[21] in chains:
            out.append(line)
        elif line.startswith("TER") and line[21] in chains:
            out.append(line)
    return out


def ca_coords(pdb_lines):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in pdb_lines
                     if l.startswith("ATOM") and l[12:16].strip() == "CA"])


def constrained_embed(mol, mask1, mask2, coords1, coords2, seed=cfg.SEED):
    """DeepTernary 的双端约束嵌构象：把两个 warhead 的内距约束到晶体姿态后再 ETKDG。"""
    from rdkit.Chem import rdDistGeom
    from rdkit.DistanceGeometry.DistGeom import DoTriangleSmoothing
    work = Chem.Mol(mol)
    ch1 = work.GetSubstructMatch(mask1)
    ch2 = work.GetSubstructMatch(mask2)
    if not (len(ch1) and len(ch2)):
        return None
    bounds = rdDistGeom.GetMoleculeBoundsMatrix(work)
    for mapping, crystal_coords in ((ch1, coords1), (ch2, coords2)):
        for i in range(len(mapping)):
            ai = mapping[i]
            for j in range(i + 1, len(mapping)):
                aj = mapping[j]
                a, b = min(ai, aj), max(ai, aj)
                d = float(np.linalg.norm(crystal_coords[i] - crystal_coords[j]))
                bounds[a, b] = d + 0.1
                bounds[b, a] = max(d - 0.1, 0.1)
    ok = DoTriangleSmoothing(bounds)
    params = rdDistGeom.ETKDGv3()
    params.randomSeed = seed
    params.SetBoundsMat(bounds)
    try:
        cids = rdDistGeom.EmbedMultipleConfs(work, 1, params)
    except Exception:  # noqa: BLE001
        return None
    if not len(cids):
        return None
    work = Chem.RemoveHs(work)
    return work, cids[0]


def transform_pdb_lines(lines, R, t, chain_remap=None):
    out = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            xt, yt, zt = R @ np.array([x, y, z]) + t
            line = line[:30] + f"{xt:8.3f}{yt:8.3f}{zt:8.3f}" + line[54:]
            if chain_remap is not None:
                line = line[:21] + chain_remap + line[22:]
        out.append(line)
    return out


def write_pdb(mol, path, name="LIG", chain="C", relink=False,
              anchor_poi=None, anchor_e3=None):
    """RDKit mol 写标准 PDB（去芳香标志避免 kekulize）；relink=True 时把连接子
    原子沿锚点轴重排（几何退化时的兜底）。"""
    m = Chem.RWMol(mol)
    for atom in m.GetAtoms():
        atom.SetIsAromatic(False)
        atom.SetNoImplicit(True)
    for bond in m.GetBonds():
        if bond.GetIsAromatic():
            bond.SetIsAromatic(False)
            bond.SetBondType(Chem.BondType.SINGLE)
    if relink and anchor_poi is not None:
        conf = m.GetConformer()
        pos = conf.GetPositions()
        pf = pos[anchor_poi]
        pt = pos[anchor_e3]
        axis = (pt - pf) / max(np.linalg.norm(pt - pf), 1e-6)
        tmp = np.array([1.0, 0.0, 0.0])
        perp = np.cross(axis, tmp)
        if np.linalg.norm(perp) < 0.3:
            perp = np.cross(axis, np.array([0.0, 1.0, 0.0]))
        perp = perp / max(np.linalg.norm(perp), 1e-6)
        span = float(np.linalg.norm(pt - pf))
        linker = [i for i in range(m.GetNumAtoms())
                  if i != anchor_poi and i != anchor_e3]
        for idx, l in enumerate(sorted(linker)):
            frac = 0.15 + 0.7 * (idx / max(len(linker), 1))
            p_ = pf + axis * (span * frac) + perp * 0.7 * ((idx % 5) - 2)
            conf.SetAtomPosition(l, Point3D(*p_))
    block = Chem.MolToPDBBlock(m.GetMol())
    with open(path, "w") as f:
        f.write(block)


def submol(mol, atom_indices):
    """提取原子子集（保坐标）。"""
    em = Chem.RWMol(mol)
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(atom_indices),
                      reverse=True):
        em.RemoveAtom(idx)
    return em.GetMol()


def mcs_pairs(mol_a, mol_b, ref_xtal_conf, mol_a_conf, max_pairs=6):
    """MCS 对应关系：枚举匹配组合，取 Kabsch 拟合最紧的一组。"""
    res = rdFMCS.FindMCS([mol_a, mol_b], timeout=60,
                         atomCompare=rdFMCS.AtomCompare.CompareElements,
                         bondCompare=rdFMCS.BondCompare.CompareAny,
                         ringMatchesRingOnly=True, matchValences=False)
    if res.numAtoms < cfg.MIN_POI_MCS:
        return None
    patt = Chem.MolFromSmarts(res.smartsString)
    a_matches = mol_a.GetSubstructMatches(patt, uniquify=False, maxMatches=2)
    b_matches = mol_b.GetSubstructMatches(patt, uniquify=False, maxMatches=max_pairs)
    best = None
    for ma in a_matches:
        for mb in b_matches:
            P = np.array(mol_a_conf.GetPositions())[list(ma)]
            Q = np.array(ref_xtal_conf.GetPositions())[list(mb)]
            R, t = kabsch(P, Q)
            rmsd = float(np.sqrt((((R @ P.T).T + t - Q) ** 2).sum(1).mean()))
            if best is None or rmsd < best[2]:
                best = (list(ma), list(mb), rmsd, R, t)
    return best


# ---------------------------------------------------------------- 参照包
class ReferencePack:
    """参照晶体一次加载，整个任务复用。"""

    def __init__(self):
        spec = cfg.POI_REF
        self.td2_ideal = Chem.MolFromMolFile(str(cfg.REFERENCES_DIR / spec["ideal"]),
                                             removeHs=True)
        self.td2_xtal, n_td2 = read_ligand_xtal(
            cfg.REFERENCES_DIR / spec["pdb"], spec["resname"], None,
            template=self.td2_ideal)
        if n_td2 <= 0 or self.td2_xtal is None:
            raise RuntimeError("TD2 参照提取失败（5E89.pdb）")
        self.poi_prot_lines = protein_chains_pdb(
            cfg.REFERENCES_DIR / spec["pdb"], spec["protein_chains"])

        self.e3_cache = {}
        for e3_name, s in cfg.E3_REF.items():
            ideal = Chem.MolFromMolFile(str(cfg.REFERENCES_DIR / s["ideal"]),
                                        removeHs=True)
            lig, n = read_ligand_xtal(cfg.REFERENCES_DIR / s["pdb"], s["resname"],
                                      s["lig_chain"], template=ideal)
            if n <= 0 or lig is None:
                raise RuntimeError(f"{s['resname']} 参照提取失败（{s['pdb']}）")
            prot = protein_chains_pdb(cfg.REFERENCES_DIR / s["pdb"],
                                      s["protein_chains"])
            self.e3_cache[e3_name] = {"ideal": ideal, "lig": lig, "prot": prot}

        self.poi_ca = ca_coords(self.poi_prot_lines)
        self.e3_ca = {k: ca_coords(v["prot"]) for k, v in self.e3_cache.items()}
        self.poi_centroid = np.mean(
            [[float(l[30:38]), float(l[38:46]), float(l[46:54])]
             for l in self.poi_prot_lines if l.startswith("ATOM")], axis=0)


_PACK: ReferencePack | None = None


def get_pack() -> ReferencePack:
    global _PACK
    if _PACK is None:
        _PACK = ReferencePack()
    return _PACK


# ---------------------------------------------------------------- 单分子装配
def build_inputs(smiles: str, e3: str, outdir, pack: ReferencePack,
                 cid: str = "LIG") -> dict:
    """装配一个 PROTAC 的 DeepTernary 输入（研究项目 main() 单行逻辑的原样移植）。"""
    e3 = (e3 or "").strip().upper()
    if e3 not in cfg.E3_REF:
        return {"status": "failed", "reason": f"不支持的 E3：{e3 or '未指定'}（可选 {', '.join(cfg.E3_REF)}）"}
    outdir = str(outdir)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"status": "failed", "reason": "SMILES 无法解析"}
    mol = Chem.AddHs(mol)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = cfg.SEED
    ps.useRandomCoords = True
    if AllChem.EmbedMolecule(mol, ps) != 0:
        ps.randomSeed = cfg.SEED + 1
        if AllChem.EmbedMolecule(mol, ps) != 0:
            return {"status": "failed", "reason": "ETKDG 构象生成失败"}
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    except Exception:  # noqa: BLE001 —— 保留未优化构象，DeepTernary 会重嵌
        pass
    mol = Chem.RemoveHs(mol)

    # --- POI 对齐：PROTAC → TD2 晶体坐标系 --------------------------------
    fit = mcs_pairs(mol, pack.td2_ideal, pack.td2_xtal.GetConformer(),
                    mol.GetConformer())
    if fit is None:
        return {"status": "failed", "reason": f"POI warhead MCS 过小（<{cfg.MIN_POI_MCS} 原子）"}
    a_idx, b_idx, rmsd_poi, R, t = fit
    conf = mol.GetConformer()
    td2_conf = pack.td2_xtal.GetConformer()
    td2_pos = td2_conf.GetPositions()
    for i in range(mol.GetNumAtoms()):
        p = np.array(list(conf.GetAtomPosition(i)))
        xt, yt, zt = R @ p + t
        conf.SetAtomPosition(i, Point3D(float(xt), float(yt), float(zt)))
    # graft：POI warhead 原子直接取晶体坐标（predict.py 的原子替换需要 <1 Å）
    for local_i, j in zip(a_idx, b_idx):
        conf.SetAtomPosition(local_i, Point3D(*td2_pos[j]))

    write_pdb(submol(mol, a_idx), f"{outdir}/lig2_mask.pdb", name="M2")
    write_pdb(pack.td2_xtal, f"{outdir}/unbound_lig2.pdb", name="U2")
    with open(f"{outdir}/p2.pdb", "w") as f:
        f.writelines(pack.poi_prot_lines)

    # --- E3 对齐：晶体 E3 配体 → PROTAC 的 E3 warhead ---------------------
    e3_spec = cfg.E3_REF[e3]
    e3_xtal = pack.e3_cache[e3]["lig"]
    ideal = pack.e3_cache[e3]["ideal"]
    e3_prot_lines = pack.e3_cache[e3]["prot"]
    res2 = rdFMCS.FindMCS([mol, ideal], timeout=60,
                          atomCompare=rdFMCS.AtomCompare.CompareElements,
                          bondCompare=rdFMCS.BondCompare.CompareAny,
                          ringMatchesRingOnly=True, matchValences=False)
    if res2.numAtoms < cfg.MIN_POI_MCS:
        return {"status": "failed", "reason": f"E3 warhead MCS 过小（<{cfg.MIN_POI_MCS} 原子）"}
    patt = Chem.MolFromSmarts(res2.smartsString)
    # 关键：K/J 系列有重复 motif，SMARTS 会匹配到 POI warhead 附近的错误位点，
    # 必须枚举全部匹配并剔除与 POI warhead 重叠者（研究项目实测教训）
    ma_all = [m for m in mol.GetSubstructMatches(patt, uniquify=False, maxMatches=50)
              if not (set(m) & set(a_idx))]
    mb = e3_xtal.GetSubstructMatches(patt, uniquify=False, maxMatches=6)
    if not ma_all or not mb:
        return {"status": "failed", "reason": "E3 warhead 匹配失败"}
    prot_pos_all = np.array(mol.GetConformer().GetPositions())
    best = None
    for a2 in ma_all:
        for b2 in mb:
            Q = np.array(e3_xtal.GetConformer().GetPositions())[list(b2)]
            P = prot_pos_all[list(a2)]
            R2, t2 = kabsch(Q, P)          # 晶体配体 → PROTAC
            rmsd = float(np.sqrt((((R2 @ Q.T).T + t2 - P) ** 2).sum(1).mean()))
            if best is None or rmsd < best[2]:
                best = (list(a2), list(b2), rmsd, R2, t2)
    a2, b2, rmsd_e3, R2, t2 = best

    # --- 摆放可行性：E3 蛋白外推到连接子伸展距离，随机扭转采样至无碰撞 -------
    linker_atoms = sorted(set(range(mol.GetNumAtoms())) - set(a_idx) - set(a2))
    n_link = max(len(linker_atoms), 1)
    d_target = float(np.clip(0.9 * 1.45 * (n_link + 1) + 8.0,
                             e3_spec["min_anchor"], 55.0))
    prot_pos = np.array([list(conf.GetAtomPosition(k))
                         for k in range(mol.GetNumAtoms())])
    link_c = prot_pos[linker_atoms].mean(0)
    anchor_poi = a_idx[int(np.argmin(((prot_pos[a_idx] - link_c) ** 2).sum(1)))]
    anchor_e3 = a2[int(np.argmin(((prot_pos[a2] - link_c) ** 2).sum(1)))]
    from scipy.spatial.transform import Rotation
    rng = np.random.default_rng(cfg.SEED)
    anchor_poi_xyz = prot_pos[anchor_poi]
    # E3 锚点沿 gal3 表面外法向摆放（折叠构象的 E3 warhead 常位于 gal3 包络内，
    # 不能用它定义方向）
    poi_out = anchor_poi_xyz - pack.poi_centroid
    poi_out /= max(np.linalg.norm(poi_out), 1e-6)
    anchor_target = anchor_poi_xyz + poi_out * d_target

    e3_conf = e3_xtal.GetConformer()
    e3_pos = e3_conf.GetPositions()
    e3_ca_coords = pack.e3_ca[e3]
    poi_ca = pack.poi_ca
    k = list(a2).index(anchor_e3)
    anchor_xtal = e3_pos[b2[k]]
    u_out = anchor_xtal - e3_ca_coords.mean(0)
    u_out /= max(np.linalg.norm(u_out), 1e-6)
    # 口袋底部朝外 → 口袋开口朝向 POI 锚点
    target_dir = poi_out
    R0, _ = Rotation.align_vectors([target_dir], [u_out])
    e3_all = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                       for l in e3_prot_lines if l.startswith("ATOM")])
    poi_all = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                        for l in pack.poi_prot_lines if l.startswith("ATOM")])
    best_clear, best_R = -1.0, None
    trial = 0
    for trial in range(400):
        psi = rng.uniform(0, 2 * np.pi)
        Rt = Rotation.from_rotvec(np.asarray(target_dir) * psi) * R0
        if trial > 0:
            Rt = Rotation.from_rotvec(rng.normal(0, 0.45, 3)) * Rt
        Rr = Rt.as_matrix()
        ca_moved = (Rr @ (e3_ca_coords - anchor_xtal).T).T + anchor_target
        ca_min = np.sqrt(((ca_moved[:, None, :] - poi_ca[None, :, :]) ** 2)
                         .sum(-1)).min()
        if ca_min < 7.0 and trial % 50 != 0:
            continue
        p1_moved = (Rr @ (e3_all - anchor_xtal).T).T + anchor_target
        dmin = np.sqrt(((p1_moved[:, None, :] - poi_all[None, :, :]) ** 2)
                       .sum(-1)).min()
        if dmin > best_clear:
            best_clear = dmin
            best_R = Rr
        if dmin > 2.8:
            break
    clear_warn = None
    if best_clear <= 2.8:
        clear_warn = f"未找到完全无碰撞朝向（最好 {best_clear:.2f} Å），使用重叠最小者"
    Rr = best_R
    R2 = Rr
    t2 = anchor_target - Rr @ anchor_xtal
    e3_pos_moved = (Rr @ (e3_pos - anchor_xtal).T).T + anchor_target
    for local_i, j in zip(a2, b2):
        conf.SetAtomPosition(local_i, Point3D(*e3_pos_moved[j]))

    e3_lig_moved = Chem.Mol(e3_xtal)
    c2 = e3_lig_moved.GetConformer()
    for i in range(e3_lig_moved.GetNumAtoms()):
        p = np.array(list(c2.GetAtomPosition(i)))
        xt, yt, zt = R2 @ p + t2
        c2.SetAtomPosition(i, Point3D(float(xt), float(yt), float(zt)))
    write_pdb(e3_lig_moved, f"{outdir}/unbound_lig1.pdb", name="U1")

    # --- 最终配体构象：双端约束嵌入 → 直接嫁接 → 连接子重排，逐级回退 -------
    pose1 = e3_pos_moved[b2]
    pose2 = td2_pos[b_idx]
    m1 = submol(mol, a2)
    c1 = m1.GetConformer()
    for i, p in enumerate(pose1):
        c1.SetAtomPosition(i, Point3D(*p))
    m2 = submol(mol, a_idx)
    cc2 = m2.GetConformer()
    for i, p in enumerate(pose2):
        cc2.SetAtomPosition(i, Point3D(*p))

    lig_path = f"{outdir}/lig.pdb"

    def lig_readback_ok():
        # 回读校验对齐真实消费者：DeepTernary 用 biopandas 读 lig.pdb（不解析键级）。
        # 新版 RDKit 的 MolFromPDBFile 会感知出错误键级（研究项目产物同样读不过），
        # 故改为文本解析 + 原子数 + 最小间距检查（等价于 biopandas 校验，零额外依赖）。
        try:
            atoms = [l for l in open(lig_path).read().splitlines()
                     if l.startswith("HETATM")]
            if len(atoms) != mol.GetNumAtoms():
                return False
            pos = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                            for l in atoms])
            d2 = ((pos[:, None, :] - pos[None, :, :]) ** 2).sum(-1)
            np.fill_diagonal(d2, np.inf)
            return float(d2.min()) > 0.04  # 无塌缩（>0.2 Å）
        except Exception:  # noqa: BLE001
            return False

    chosen = None
    res = constrained_embed(mol, m1, m2, pose1, pose2, seed=cfg.SEED)
    if res is not None:
        work, conf_id = res
        conf_new = work.GetConformer(conf_id)
        pos_new = conf_new.GetPositions()
        fit_R, fit_t = kabsch(pos_new[list(work.GetSubstructMatch(m1))], pose1)
        for i in range(work.GetNumAtoms()):
            p = fit_R @ pos_new[i] + fit_t
            conf_new.SetAtomPosition(i, Point3D(*p))
        write_pdb(work, lig_path, name=(cid or "LIG")[:3])
        if lig_readback_ok():
            chosen = "constrained"
    if chosen is None:
        write_pdb(mol, lig_path, name=(cid or "LIG")[:3])
        if lig_readback_ok():
            chosen = "grafted"
    if chosen is None:
        write_pdb(mol, lig_path, name=(cid or "LIG")[:3], relink=True,
                  anchor_poi=anchor_poi, anchor_e3=anchor_e3)
        if lig_readback_ok():
            chosen = "relinked"
    if chosen is None:
        return {"status": "failed", "reason": "无可读配体构象（三级回退均失败）"}
    write_pdb(submol(mol, a2), f"{outdir}/lig1_mask.pdb", name="M1")
    with open(f"{outdir}/p1.pdb", "w") as f:
        f.writelines(transform_pdb_lines(e3_prot_lines, R2, t2))

    return {"status": "ok", "e3": e3, "mcs_poi": len(a_idx), "mcs_e3": len(a2),
            "poi_rmsd": round(float(rmsd_poi), 2), "e3_rmsd": round(float(rmsd_e3), 2),
            "lig_kind": chosen, "best_clear": round(float(best_clear), 2),
            "warning": clear_warn, "input_dir": outdir}
