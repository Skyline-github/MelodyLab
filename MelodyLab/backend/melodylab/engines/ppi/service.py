"""PPI 蛋白互作网络引擎：STRING-DB v12 + networkx 图分析 + 多格式导出。

STRING 不可用（无网络）时降级为"仅节点"模式。图分析全套：
Degree / Betweenness / Closeness / Eigenvector / MCC(极大团) / MNC / DMNC /
k-core / 标签传播社区发现。
"""
from __future__ import annotations

import json
from pathlib import Path

import httpx

from ...config import settings

STRING_BASE = "https://version-12-0.string-db.org/api"


def _string_get(path: str, params: dict) -> dict | list:
    params = {**params, "species": params.get("species") or 9606}
    with httpx.Client(timeout=httpx.Timeout(60.0, connect=15.0),
                      follow_redirects=True) as c:
        r = c.get(f"{STRING_BASE}/{path}", params=params)
        r.raise_for_status()
        if params.get("output") == "tsv" or path.endswith((".tsv", ".tsv-no-header")):
            return _parse_tsv(r.text)
        return r.json()


def _parse_tsv(text: str) -> list[dict]:
    lines = text.strip().splitlines()
    if not lines:
        return []
    header = lines[0].split("\t")
    return [dict(zip(header, l.split("\t"))) for l in lines[1:]]


def fetch_network(identifiers: list[str], species: int = 9606, required_score: int = 400,
                  max_extra: int = 0, network_type: str = "functional") -> tuple[list[dict], list[dict]]:
    """返回 (nodes, edges)。include_partners 时追加一阶伙伴（受 max_extra 限制）。"""
    body = {"identifiers": "%0d".join(identifiers), "species": species,
            "required_score": required_score, "network_type": network_type,
            "caller_identity": "MelodyLab"}
    with httpx.Client(timeout=httpx.Timeout(90.0, connect=15.0), follow_redirects=True) as c:
        r = c.post(f"{STRING_BASE}/json/network", data=body)
        r.raise_for_status()
        raw_edges = r.json()
    nodes_seen = {}
    edges = []
    for e in raw_edges:
        for a, b in ((e["preferredName_A"], e["preferredName_B"]),):
            nodes_seen.setdefault(a, {"id": a, "string_id": e["stringId_A"]})
            nodes_seen.setdefault(b, {"id": b, "string_id": e["stringId_B"]})
            edges.append({"source": a, "target": b,
                          "score": float(e.get("score", 0))})
    if max_extra > 0:
        r2 = c_post_interactors(list(nodes_seen), species, max_extra)
        for it in r2:
            for a, b in ((it.get("preferredName_A"), it.get("preferredName_B")),):
                if a and b and (a not in nodes_seen or b not in nodes_seen):
                    nodes_seen.setdefault(a, {"id": a, "string_id": it.get("stringId_A")})
                    nodes_seen.setdefault(b, {"id": b, "string_id": it.get("stringId_B")})
                    edges.append({"source": a, "target": b, "score": float(it.get("score", 0)),
                                  "extra": True})
    return list(nodes_seen.values()), edges


def c_post_interactors(ids: list[str], species: int, limit: int) -> list[dict]:
    with httpx.Client(timeout=httpx.Timeout(90.0, connect=15.0), follow_redirects=True) as c:
        r = c.post(f"{STRING_BASE}/json/interactors", data={
            "identifiers": "%0d".join(ids), "species": species, "limit": str(limit),
            "caller_identity": "MelodyLab"})
        r.raise_for_status()
        return r.json()


def fetch_enrichment(identifiers: list[str], species: int = 9606) -> list[dict]:
    with httpx.Client(timeout=httpx.Timeout(90.0, connect=15.0), follow_redirects=True) as c:
        r = c.post(f"{STRING_BASE}/json/enrichment", data={
            "identifiers": "%0d".join(identifiers), "species": species,
            "caller_identity": "MelodyLab"})
        r.raise_for_status()
        rows = r.json()
    return [{"category": e.get("category"), "term": e.get("term"),
             "description": e.get("description"),
             "fdr": e.get("fdr"), "p_value": e.get("p_value"),
             "genes": e.get("preferredNames", "")} for e in rows[:40]]


# ---------------- 图分析 ----------------
def analyze_graph(nodes: list[dict], edges: list[dict]) -> dict:
    import networkx as nx

    G = nx.Graph()
    for n in nodes:
        G.add_node(n["id"])
    for e in edges:
        if G.has_node(e["source"]) and G.has_node(e["target"]):
            G.add_edge(e["source"], e["target"], weight=e.get("score", 1.0))
    if G.number_of_nodes() == 0:
        return {"centrality": {}, "communities": {}, "stats": {}}
    centrality = {}
    deg = dict(G.degree())
    bet = nx.betweenness_centrality(G) if G.number_of_nodes() <= 400 else {}
    clo = nx.closeness_centrality(G) if G.number_of_nodes() <= 400 else {}
    try:
        eig = nx.eigenvector_centrality_numpy(G) if G.number_of_nodes() <= 400 else {}
    except Exception:
        eig = {}
    mcc = _mcc(G)
    mnc = _mnc(G)
    dmnc = _dmnc(G)
    core = nx.core_number(G)
    for n in G.nodes():
        centrality[n] = {
            "degree": deg.get(n, 0),
            "betweenness": round(bet.get(n, 0.0), 6) if bet else None,
            "closeness": round(clo.get(n, 0.0), 6) if clo else None,
            "eigenvector": round(eig.get(n, 0.0), 6) if eig else None,
            "mcc": mcc.get(n), "mnc": mnc.get(n), "dmnc": dmnc.get(n),
            "kcore": core.get(n),
        }
    communities = _label_propagation(G)
    return {
        "centrality": centrality,
        "communities": communities,
        "stats": {"n_nodes": G.number_of_nodes(), "n_edges": G.number_of_edges(),
                  "density": round(nx.density(G), 4),
                  "n_communities": len(set(communities.values()))},
    }


def _mcc(G) -> dict:
    """Maximal Clique Centrality（CytoHubba 定义）。"""
    import networkx as nx

    cliques = list(nx.find_cliques(G))
    out = {n: 0 for n in G.nodes()}
    for cl in cliques:
        s = len(cl)
        for n in cl:
            out[n] += s * (s - 1) // 2
    return out


def _mnc(G) -> dict:
    import networkx as nx

    out = {}
    for n in G.nodes():
        best = 0
        for cl in nx.find_cliques(G.subgraph(list(G.neighbors(n)) + [n]) if n in G else []):
            best = max(best, len(cl))
        out[n] = best
    return out


def _dmnc(G) -> dict:
    out = {}
    for n in G.nodes():
        nb = list(G.neighbors(n))
        sub = G.subgraph(nb)
        e = sub.number_of_edges()
        k = len(nb)
        out[n] = round(e / max(1, k - 1) + k, 3) if k > 1 else k
    return out


def _label_propagation(G) -> dict:
    import networkx as nx

    if G.number_of_edges() == 0:
        return {n: 0 for n in G.nodes()}
    comms = nx.algorithms.community.label_propagation_communities(G)
    out = {}
    for i, c in enumerate(sorted(comms, key=len, reverse=True)):
        for n in c:
            out[n] = i
    return out


# ---------------- 导出 ----------------
def export_graph(nodes: list[dict], edges: list[dict], out_base: Path,
                 analysis: dict) -> list[Path]:
    import networkx as nx

    G = nx.Graph()
    for n in nodes:
        G.add_node(n["id"], **{k: v for k, v in n.items() if k != "id"})
    for e in edges:
        G.add_edge(e["source"], e["target"], weight=e.get("score", 1.0))
    files = []
    try:
        gpml = out_base.with_suffix(".graphml")
        nx.write_graphml(G, gpml)
        files.append(gpml)
    except Exception:
        pass
    cyto = {
        "data": {"name": "MelodyLab PPI network"},
        "elements": {
            "nodes": [{"data": {"id": n["id"],
                                "centrality": (analysis.get("centrality") or {}).get(n["id"], {})}}
                      for n in nodes],
            "edges": [{"data": {"source": e["source"], "target": e["target"],
                                "score": e.get("score", 1.0)}} for e in edges],
        },
    }
    cj = out_base.with_suffix(".cytoscape.json")
    cj.write_text(json.dumps(cyto, ensure_ascii=False, indent=1), encoding="utf-8", newline="\n")
    files.append(cj)
    return files
