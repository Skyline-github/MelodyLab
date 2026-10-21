"""逆合成路线规划引擎：RDKit BRICS + RECAP 混合断键 → 递归分层路线树。

纯本地计算。置信度 = 键型先验（酰胺/酯/芳基键等易合成键加分）× 前体复杂度惩罚。
"""
from __future__ import annotations

from rdkit import Chem

# 断键类型先验（BRICS 键型 -1..16；常见可逆合成键加分）
BOND_PRIOR = {
    0: 0.75,   # amide 类
    1: 0.70,   # ester 类
    3: 0.65,   # amine
    5: 0.60,   # olefin
    9: 0.55,   # aromatic C-N
    10: 0.55,  # aromatic C-C
    11: 0.50,  # sulfonamide?/urea 类
    -1: 0.35,  # 默认单键
}


def _mol_from_smiles(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"无法解析 SMILES: {smiles}")
    return mol


def _complexity(mol) -> float:
    """前体复杂度 0..1（重原子数 + 环数 + 杂原子比例启发式）。"""
    n_heavy = mol.GetNumHeavyAtoms()
    n_rings = mol.GetRingInfo().NumRings()
    n_het = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() not in (6, 1))
    score = 0.02 * n_heavy + 0.06 * n_rings + 0.03 * n_het
    return min(1.0, score)


def _fragment(smiles: str) -> str | None:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    return Chem.MolToSmiles(m)


def one_step(smiles: str, engine: str) -> list[dict]:
    """单步断键：返回 [{bond_type, precursors: [smiles...], score}]。"""
    import re

    m = _mol_from_smiles(smiles)
    cuts: list[tuple[int, list[str]]] = []

    def _strip_dummies(smi: str) -> list[str]:
        plain = re.sub(r"\[\d+\*\]", "[*]", smi)
        frags = []
        for part in plain.split("."):
            mm = Chem.MolFromSmiles(part)
            if mm is not None and mm.GetNumAtoms() >= 1:
                frags.append(Chem.MolToSmiles(mm))
        return frags

    if engine in ("brics", "hybrid"):
        from rdkit.Chem import BRICS

        try:
            frags_all: list[str] = []
            btypes: list[int] = []
            for frag in BRICS.BRICSDecompose(m, minFragmentSize=2):
                smi = frag if isinstance(frag, str) else Chem.MolToSmiles(frag)
                btypes += [int(x) for x in re.findall(r"\[(\d+)\*\]", smi)]
                frags_all += _strip_dummies(smi)
            frags_all = sorted(set(frags_all))
            if len(frags_all) >= 2:
                cuts.append((btypes[0] if btypes else -1, frags_all))
        except Exception:
            pass
    if engine in ("recap", "hybrid"):
        from rdkit.Chem import Recap

        try:
            recap = Recap.RecapDecompose(m)
            if recap is not None:
                pre = []
                for k in recap.GetAllChildren().values():
                    mm = Chem.MolFromSmiles(k.mol_smiles)
                    if mm is not None:
                        pre.append(Chem.MolToSmiles(mm))
                pre = sorted(set(pre))
                if len(pre) >= 2:
                    has_n = any(re.search(r"[nN]", p) for p in pre)
                    cuts.append((9 if has_n else -1, pre))
        except Exception:
            pass
    # 去重（同一组前体只留一个）
    seen: set[tuple] = set()
    steps = []
    for btype, pre in cuts:
        key = tuple(sorted(pre))
        if key in seen:
            continue
        seen.add(key)
        comp = sum(_complexity(Chem.MolFromSmiles(p)) for p in pre) / len(pre)
        score = round(max(0.05, min(0.95, BOND_PRIOR.get(btype, 0.4) * (1.0 - 0.5 * comp))), 3)
        steps.append({"bond_type": int(btype), "precursors": pre, "score": score})
    steps.sort(key=lambda s: -s["score"])
    return steps


def plan_routes(smiles: str, engine: str = "hybrid", max_routes: int = 5,
                max_depth: int = 2) -> dict:
    """递归展开路线树（深度 max_depth；叶子 = 可购买/简单前体，这里只按深度截断）。"""
    root = _fragment(smiles) or smiles

    def expand(smi: str, depth: int) -> dict:
        node: dict = {"smiles": smi, "children": [], "score": 1.0}
        if depth <= 0:
            return node
        try:
            steps = one_step(smi, engine)
        except Exception:
            return node
        node["children"] = [
            {"bond_type": s["bond_type"], "score": s["score"],
             "precursors": [expand(p, depth - 1) for p in s["precursors"]]}
            for s in steps[:max_routes]
        ]
        return node

    tree = expand(root, max_depth)

    def flatten(node: dict, path: list, acc: list, top_score: float = 1.0) -> None:
        for ch in node.get("children", []):
            route = {"steps": path + [{"smiles": node["smiles"], "bond_type": ch["bond_type"],
                                       "score": ch["score"],
                                       "precursors": [p["smiles"] for p in ch["precursors"]]}],
                     "confidence": round(top_score * ch["score"], 3)}
            acc.append(route)

    routes: list[dict] = []

    def walk(node: dict, acc_path: list) -> None:
        if not node.get("children"):
            return
        for ch in node["children"]:
            step = {"smiles": node["smiles"], "bond_type": ch["bond_type"],
                    "score": ch["score"],
                    "precursors": [p["smiles"] for p in ch["precursors"]]}
            if ch["score"] > 0.45:
                routes.append({"steps": acc_path + [step], "confidence": ch["score"]})
            for sub in ch["precursors"]:
                if sub.get("children"):
                    walk(sub, acc_path + [step])

    walk(tree, [])
    routes.sort(key=lambda r: -r["confidence"])
    for i, r in enumerate(routes[:max_routes], 1):
        r["route_id"] = i
    return {"target": root, "engine": engine,
            "routes": routes[:max_routes],
            "note": "置信度=键型先验×前体复杂度（启发式，非反应模板预测）"}
