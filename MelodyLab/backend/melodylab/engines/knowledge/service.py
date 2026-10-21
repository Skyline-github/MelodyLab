"""知识库检索引擎：ChEMBL / PubChem / UniProt / RCSB / AlphaFold / ZINC-22 / GEO。

全部走官方公开 API（httpx，线程内同步调用——检索量小，无需全局事件循环）。
每个结果附来源 URL；结构文件缓存到任务 work/discovery_cache/。
"""
from __future__ import annotations

import csv
import json
import time
import zipfile
from pathlib import Path

import httpx

from ...config import settings

UA = "MelodyLab/0.2 (local CADD workbench; open-source)"
TIMEOUT = httpx.Timeout(30.0, connect=10.0)


def _get(url: str, **kw) -> httpx.Response:
    headers = {"User-Agent": UA, **kw.pop("headers", {})}
    last = None
    for attempt in range(3):
        try:
            with httpx.Client(timeout=TIMEOUT, follow_redirects=True, headers=headers) as c:
                return c.get(url, **kw)
        except (httpx.TransportError, ConnectionError) as exc:
            last = exc
            time.sleep(1.5 * (attempt + 1))
    raise last


def _post_json(url: str, payload: dict, headers: dict | None = None) -> dict:
    with httpx.Client(timeout=TIMEOUT, follow_redirects=True) as c:
        r = c.post(url, json=payload, headers={"User-Agent": UA, **(headers or {})})
        r.raise_for_status()
        return r.json()


# ---------------- ChEMBL：靶点 / 分子 / 活性 ----------------
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"


def chembl_activity(query: str, query_type: str, activity_types: list[str], max_records: int) -> list[dict]:
    """按靶点（基因名/UniProt/ChEMBL ID）检索活性数据。"""
    rows: list[dict] = []
    if query_type in ("gene", "protein", "uniprot"):
        # 解析成 ChEMBL target id：基因名走 UniProt(reviewed)→xref 关联，失败回退 pref_name
        if query.upper().startswith("CHEMBL"):
            target_id = query.upper()
        else:
            targets = []
            if query_type == "gene":
                try:
                    # 基因名 → UniProt(reviewed) → 蛋白名 → ChEMBL pref_name 关联
                    up = _get("https://rest.uniprot.org/uniprotkb/search", params={
                        "query": f"gene:{query} AND organism_id:9606 AND reviewed:true",
                        "format": "json", "size": 1})
                    up.raise_for_status()
                    res = up.json().get("results", [])
                    if res:
                        desc = res[0].get("proteinDescription", {})
                        names = [((desc.get("recommendedName") or {}).get("fullName", {}) or {}).get("value", "")]
                        names += [a.get("fullName", {}).get("value", "")
                                  for a in desc.get("alternativeNames", [])]
                        for nm in [n for n in names if n]:
                            r = _get(f"{CHEMBL}/target.json", params={
                                "pref_name__icontains": nm.split("(")[0].strip(), "limit": 3})
                            if r.status_code == 200:
                                targets = r.json().get("targets", [])
                                if targets:
                                    break
                except Exception:
                    targets = []
            if not targets:
                r = _get(f"{CHEMBL}/target.json", params={"pref_name__icontains": query, "limit": 1})
                r.raise_for_status()
                targets = r.json().get("targets", [])
            if not targets:
                return []
            target_id = targets[0]["target_chembl_id"]
        acts = _get(f"{CHEMBL}/activity.json", params={
            "target_chembl_id": target_id, "limit": min(max_records, 200)})
        acts.raise_for_status()
        for a in acts.json().get("activities", []):
            if activity_types and a.get("standard_type") not in activity_types:
                continue
            rows.append({
                "source": "chembl", "target": target_id,
                "compound": a.get("molecule_chembl_id"),
                "smiles": None, "name": a.get("molecule_pref_name"),
                "activity_type": a.get("standard_type"),
                "value": a.get("standard_value"), "units": a.get("standard_units"),
                "assay": a.get("assay_chembl_id"),
                "url": f"https://www.ebi.ac.uk/chembl/compound_report_card/{a.get('molecule_chembl_id')}/",
            })
            if len(rows) >= max_records:
                break
    elif query_type == "inhibitor":
        # 分子名/SMILES 相似检索走 PubChem；此处按名字查 ChEMBL 分子
        mols = _get(f"{CHEMBL}/molecule.json", params={"pref_name__iexact": query, "limit": 10})
        mols.raise_for_status()
        for m in mols.json().get("molecules", []):
            rows.append({
                "source": "chembl", "compound": m.get("molecule_chembl_id"),
                "name": m.get("pref_name"),
                "smiles": ((m.get("molecule_structures") or {}).get("canonical_smiles")),
                "url": f"https://www.ebi.ac.uk/chembl/compound_report_card/{m.get('molecule_chembl_id')}/",
            })
    return rows


# ---------------- PubChem：生物测定 / 相似分子 ----------------
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def pubchem_assays(query: str, max_records: int) -> list[dict]:
    """PubChem 生物测定检索（PUG-REST 明文检索端点；不可用时空结果降级）。"""
    r = _get(f"{PUBCHEM}/assay/text/{query}/json", params={"listkey_count": min(max_records, 50)})
    if r.status_code != 200:
        return []
    out = []
    for c in (r.json().get("PC_AssayContainer") or [])[:max_records]:
        sub = c.get("PC_AssaySubmit") or {}
        aid = (sub.get("assay") or {}).get("id")
        descr = (sub.get("assay") or {}).get("descr") or {}
        out.append({"source": "pubchem", "assay": descr.get("name"), "aid": aid,
                    "url": f"https://pubchem.ncbi.nlm.nih.gov/bioassay/{aid}" if aid else None})
    return [x for x in out if x["aid"]]


def pubchem_similarity(smiles: str, max_records: int) -> list[dict]:
    # POST 相似性检索 → ListKey → 轮询属性
    with httpx.Client(timeout=TIMEOUT, follow_redirects=True) as c:
        r = c.post(f"{PUBCHEM}/compound/similarity/smiles/json",
                   data={"smiles": smiles, "Threshold": "85", "MaxRecords": str(min(max_records, 100))})
        r.raise_for_status()
        info = r.json().get("Waiting", {}).get("ListKey", {}).get("id")
        if not info:
            return []
        for _ in range(20):
            time.sleep(1.0)
            rr = c.get(f"{PUBCHEM}/compound/listkey/{info}/property/CanonicalSMILES,IsomericSMILES,MolecularFormula/json")
            if rr.status_code == 200:
                data = rr.json()
                props = ((data.get("PropertyTable") or {}).get("Properties")) or []
                return [{
                    "source": "pubchem", "compound": f"CID{p.get('CID')}",
                    "smiles": p.get("IsomericSMILES") or p.get("CanonicalSMILES"),
                    "formula": p.get("MolecularFormula"),
                    "url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{p.get('CID')}",
                } for p in props[:max_records]]
            if rr.status_code != 202:
                rr.raise_for_status()
    return []


# ---------------- UniProt ----------------
def uniprot_search(query: str, max_records: int) -> list[dict]:
    r = _get("https://rest.uniprot.org/uniprotkb/search", params={
        "query": query, "format": "json", "size": min(max_records, 50)})
    r.raise_for_status()
    out = []
    for e in r.json().get("results", []):
        names = e.get("proteinDescription", {}).get("recommendedName", {})
        seq = e.get("sequence", {})
        out.append({
            "source": "uniprot", "accession": e.get("primaryAccession"),
            "name": names.get("fullName", {}).get("value"),
            "organism": ((e.get("organism") or {}).get("scientificName")),
            "length": seq.get("length"),
            "url": f"https://www.uniprot.org/uniprotkb/{e.get('primaryAccession')}",
        })
    return out


# ---------------- RCSB + AlphaFold：结构检索 ----------------
def rcsb_search(query: str, max_records: int, resolution: float | None = None,
                require_ligand: bool = False, method: str | None = None) -> list[dict]:
    """RCSB v2 全文/序列检索 + 元数据过滤。query 可为基因名或 PDB id。"""
    node = {"type": "terminal", "service": "full_text",
            "parameters": {"value": query}}
    body = {
        "query": node,
        "return_type": "polymer_entity",
        "request_options": {"paginate": {"start": 0, "rows": min(max_records, 50)},
                            "results_content_type": ["experimental"]},
    }
    ids = _post_json("https://search.rcsb.org/rcsbsearch/v2/query", body)
    ids_list = [r["identifier"].split("_")[0] for r in ids.get("result_set", [])]
    out: list[dict] = []
    for pid in list(dict.fromkeys(ids_list))[:max_records]:
        meta = _get(f"https://data.rcsb.org/rest/v1/core/entry/{pid}")
        if meta.status_code != 200:
            continue
        d = meta.json()
        res = d.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0]
        exp = d.get("rcsb_entry_info", {}).get("experimental_method")
        if resolution and res and res > resolution:
            continue
        if method and exp and method.lower() not in str(exp).lower():
            continue
        if require_ligand and not d.get("rcsb_entry_info", {}).get("nonpolymer_entity_count"):
            continue
        out.append({
            "source": "rcsb", "pdb_id": pid, "resolution": res, "method": exp,
            "title": d.get("struct", {}).get("title"),
            "has_ligand": bool(d.get("rcsb_entry_info", {}).get("nonpolymer_entity_count")),
            "url": f"https://www.rcsb.org/structure/{pid}",
            "pdb_url": f"https://files.rcsb.org/download/{pid}.pdb",
        })
    # AlphaFold 兜底：查不到实验结构时给预测模型
    if not out and len(query) >= 3:
        up = uniprot_search(f'gene:{query} AND organism_id:9606', 1)
        if up:
            acc = up[0]["accession"]
            out.append({
                "source": "alphafold", "uniprot": acc,
                "url": f"https://alphafold.ebi.ac.uk/entry/{acc}",
                "pdb_url": f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v4.pdb",
                "title": f"AlphaFold model for {query} ({acc})", "method": "predicted",
            })
    return out


# ---------------- ZINC-22（CartBlanche） ----------------
def zinc_search(smiles: str, max_records: int) -> list[dict]:
    r = _get("https://cartblanche22.docking.org/search/result.txt", params={
        "smiles": smiles, "sim": "88", "count": min(max_records, 100)})
    r.raise_for_status()
    lines = r.text.strip().splitlines()
    if not lines:
        return []
    header = lines[0].split("\t")
    out = []
    for l in lines[1:max_records + 1]:
        parts = l.split("\t")
        row = dict(zip(header, parts))
        out.append({"source": "zinc22", "zinc_id": row.get("zmol_id") or row.get("zinc_id"),
                    "smiles": row.get("smiles"), "url": f"https://cartblanche22.docking.org/compound/{row.get('zmol_id') or row.get('zinc_id')}"})
    return out


# ---------------- GEO（NCBI eutils） ----------------
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def geo_search(query: str, max_records: int) -> list[dict]:
    r = _get(f"{EUTILS}/esearch.fcgi", params={
        "db": "gds", "term": f"{query} AND gse[ETYP]", "retmax": min(max_records, 50),
        "retmode": "json", "tool": "MelodyLab", "email": "local@melodylab"})
    r.raise_for_status()
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    out = []
    if ids:
        s = _get(f"{EUTILS}/esummary.fcgi", params={
            "db": "gds", "id": ",".join(ids), "retmode": "json",
            "tool": "MelodyLab", "email": "local@melodylab"})
        s.raise_for_status()
        doc = s.json().get("result", {})
        for uid in doc.get("uids", []):
            d = doc.get(uid, {})
            out.append({
                "source": "geo", "uid": uid, "accession": d.get("accession"),
                "title": d.get("title"), "summary": (d.get("summary") or "")[:300],
                "organism": d.get("taxon"), "n_samples": d.get("n_samples"),
                "url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={d.get('accession')}",
            })
    return out


# ---------------- 任务 worker：并发检索 + 缓存 + ZIP ----------------
def run_knowledge_job(record: dict, control) -> dict:
    """JobManager worker（kind=knowledge）。params: task + 各检索参数。"""
    from concurrent.futures import ThreadPoolExecutor

    params = record["params"]
    task = params["task"]
    job_dir = settings.jobs_dir / record["id"]
    cache = job_dir / "work" / "discovery_cache"
    cache.mkdir(parents=True, exist_ok=True)
    results: dict[str, list] = {}
    errors: dict[str, str] = {}

    def call(name, fn, *a):
        try:
            res = fn(*a)
            if res:
                (cache / f"{name}.json").write_text(
                    json.dumps(res, ensure_ascii=False, indent=1), encoding="utf-8")
            results[name] = res
            control.checkpoint(message=f"{name} 完成（{len(res)} 条）")
        except Exception as exc:
            errors[name] = str(exc)[:300]
            results[name] = []
            control.checkpoint(message=f"{name} 失败: {str(exc)[:80]}")

    jobs = []
    if task == "activity":
        q, qt = params["query"], params.get("query_type", "gene")
        ats = params.get("activity_types") or []
        mx = int(params.get("max_records", 50))
        jobs.append(("chembl_activity", lambda: chembl_activity(q, qt, ats, mx)))
        if qt == "gene":
            jobs.append(("pubchem_assays", lambda: pubchem_assays(q, mx)))
    elif task == "structures":
        jobs.append(("rcsb", lambda: rcsb_search(
            params["query"], int(params.get("max_records", 20)),
            params.get("resolution"), params.get("require_ligand", False), params.get("method"))))
    elif task == "zinc":
        jobs.append(("zinc", lambda: zinc_search(params["smiles"], int(params.get("max_records", 50)))))
    elif task == "geo":
        jobs.append(("geo", lambda: geo_search(params["query"], int(params.get("max_records", 20)))))
    else:
        raise ValueError(f"未知知识库任务类型: {task}")

    with ThreadPoolExecutor(max_workers=4) as ex:
        futures = [ex.submit(call, j[0], j[1]) for j in jobs]
        for f in futures:
            f.result()

    # 结构任务：下载 PDB 打包
    artifacts = []
    if task == "structures" and results.get("rcsb"):
        zip_path = job_dir / "outputs" / "structures.zip"
        manifest = []
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
            for item in results["rcsb"]:
                pu = item.get("pdb_url")
                if not pu:
                    continue
                try:
                    pdb = _get(pu)
                    pdb.raise_for_status()
                    name = (f"AF-{item.get('uniprot')}.pdb" if item["source"] == "alphafold"
                            else f"{item['pdb_id']}.pdb")
                    zf.writestr(name, pdb.text)
                    manifest.append({"file": name, "origin": pu})
                except Exception:
                    continue
            if manifest:
                zf.writestr("MANIFEST.json", json.dumps(manifest, indent=1))
        if manifest:
            artifacts.append({"name": "structures.zip", "path": str(zip_path), "kind": "zip"})

    flat = [{"task": k, **r} for k, v in results.items() for r in v]
    # 活性/相似分子检索：结果明细导出 CSV，供页面任务卡片下载
    if flat:
        csv_path = job_dir / "outputs" / "knowledge_results.csv"
        keys: list[str] = []
        for row in flat:
            for k2 in row:
                if k2 not in keys:
                    keys.append(k2)
        with open(csv_path, "w", encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
            w.writeheader()
            w.writerows(flat)
        artifacts.append({"name": "knowledge_results.csv", "path": str(csv_path), "kind": "csv"})
    return {"stats": {"counts": {k: len(v) for k, v in results.items()},
                      "errors": errors},
            "message": f"检索完成：{sum(len(v) for v in results.values())} 条结果"
                       + (f"；{len(errors)} 个源失败" if errors else ""),
            "artifacts": artifacts}

