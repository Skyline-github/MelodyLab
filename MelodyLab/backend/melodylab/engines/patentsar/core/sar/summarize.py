"""LLM-powered SAR narrative summarizer with citation enforcement.

The missing "环节⑤": turn the mining + R-group analysis outputs (dataset CSV,
R-group table) into a narrative SAR summary. The LLM only ever sees a compact,
deterministic "evidence pack" built from the data — never the raw patent — and
every claim it makes is then machine-checked:

  1. every cited compound id must exist in the pack;
  2. every data_ref (a number the claim relies on) must match the pack within
     tolerance;
  3. every bare numeric token in the statement must be covered by (2), by a
     pack-level statistic (counts, fold-change, per-group means), or by a
     fold / pEC50 delta computed between cited compounds.

Claims that fail the check are quarantined into an appendix instead of the
main narrative — the anti-hallucination counterpart of pipeline.chem_qc's
m/z cross-validation. Outputs (in --outdir):

  sar_summary.json   machine-readable: pack, claims, verification status
  sar_summary.md     human-readable narrative (verified claims only + appendix)

LLM: any OpenAI-compatible chat endpoint (Zhipu bigmodel, OpenAI, DeepSeek,
Ollama, ...). Configure with --base-url / --model / --api-key or env vars
SAR_LLM_BASE_URL / SAR_LLM_MODEL / SAR_LLM_API_KEY. With --llm none (or when
no API key is available) a deterministic digest-only summary is produced
instead — offline-safe and CI-friendly. stdlib only; RDKit is imported lazily
and only when R-groups must be recomputed from the dataset.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import statistics
import sys
import urllib.request
from datetime import date
from pathlib import Path

ACTIVITY_COL_CANDIDATES = ("activity_ec50_nm",)  # as written by sar.run
ID_COL_CANDIDATES = ("id", "example_no")


# --------------------------------------------------------------------------
# Evidence pack (deterministic, built from data only)
# --------------------------------------------------------------------------

def load_activity_rows(paths: list[Path], activity_col: str | None,
                       id_col: str | None) -> list[dict]:
    """Merge rows from dataset.csv / rgroup_table.csv into one id-keyed view."""
    rows: dict[str, dict] = {}
    for path in paths:
        if not path or not Path(path).exists():
            continue
        with open(path, newline="") as fh:
            for r in csv.DictReader(fh):
                rid = next((r[c].strip() for c in (id_col,) if c and r.get(c)),
                           None) or next((r[c].strip() for c in ID_COL_CANDIDATES
                                          if r.get(c)), None)
                if not rid:
                    continue
                rec = rows.setdefault(rid, {"id": rid})
                if not rec.get("smiles") and r.get("canonical_smiles", "").strip():
                    rec["smiles"] = r["canonical_smiles"].strip()
                elif not rec.get("smiles") and r.get("smiles", "").strip():
                    rec["smiles"] = r["smiles"].strip()
                act_src = activity_col or next(
                    (c for c in r if c and c.lower().startswith("ec50")),
                    None) or next(
                    (c for c in r if c and "activity" in c.lower()), None)
                if act_src and rec.get("ec50_nm") is None:
                    v = (r.get(act_src) or "").strip()
                    if v and v != "ND":
                        try:
                            rec["ec50_nm"] = float(v)
                        except ValueError:
                            pass
                for c in r:
                    if c and re.fullmatch(r"R\d+", c) and r[c].strip():
                        rec.setdefault("rgroups", {})[c] = r[c].strip()
                if r.get("flags", "").strip():
                    rec.setdefault("flags", set()).update(
                        f.strip() for f in r["flags"].split(";") if f.strip())
    for rec in rows.values():
        if "flags" in rec:
            rec["flags"] = sorted(rec["flags"])
    return [rows[k] for k in sorted(rows, key=lambda x: (len(x), x))]


def to_pec50(ec50_nm: float | None) -> float | None:
    """Same convention as sar.plots.to_pactivity: EC50 nM -> pEC50 (M)."""
    if ec50_nm is None or ec50_nm <= 0:
        return None
    return -math.log10(ec50_nm / 1e9)


def build_pack(rows: list[dict], target: str) -> dict:
    act = [r for r in rows if r.get("ec50_nm") is not None]
    pacts = {r["id"]: to_pec50(r["ec50_nm"]) for r in act}
    positions: dict[str, set[str]] = {}
    for r in rows:
        for p, g in r.get("rgroups", {}).items():
            positions.setdefault(p, set()).add(g)

    def entry(r):
        return {"id": r["id"], "ec50_nm": r.get("ec50_nm"),
                "pec50": pacts.get(r["id"]),
                "R": r.get("rgroups", {}), "smiles": r.get("smiles", "")}

    best = min(act, key=lambda r: r["ec50_nm"]) if act else None
    worst = max(act, key=lambda r: r["ec50_nm"]) if act else None
    pack: dict = {
        "target": target,
        "activity_metric": "EC50 in nM (lower = more potent); "
                           "pEC50 = -log10(EC50[M]), higher = more potent",
        "n_total": len(rows),
        "n_with_activity": len(act),
        "n_no_activity": len(rows) - len(act),
        "no_activity_ids": [r["id"] for r in rows if r.get("ec50_nm") is None],
        "flagged_ids": [r["id"] for r in rows if r.get("flags")],
        "compounds": [entry(r) for r in rows],
    }
    if act:
        vals = [r["ec50_nm"] for r in act]
        pack["activity_span"] = {"min_nm": min(vals), "max_nm": max(vals),
                                 "fold": max(vals) / min(vals)}
        pec = [pacts[r["id"]] for r in act]
        pack["activity_stats"] = {"pec50_median": round(statistics.median(pec), 3),
                                  "pec50_mean": round(statistics.mean(pec), 3)}
        pack["best"] = entry(best)
        pack["worst"] = entry(worst)
    pack["positions"] = {}
    ec50_by_id = {r["id"]: r["ec50_nm"] for r in act}
    for p, groups in sorted(positions.items()):
        entries = []
        for g in sorted(groups):
            ids = [r["id"] for r in rows if r.get("rgroups", {}).get(p) == g]
            gvals = [pacts[i] for i in ids if i in pacts]
            geval = [ec50_by_id[i] for i in ids if i in ec50_by_id]
            entries.append({
                "r": g, "ids": ids, "n": len(ids),
                "n_with_activity": len(gvals),
                "mean_pec50": round(statistics.mean(gvals), 3) if gvals else None,
                "ec50_nm_range": ([round(min(geval), 3), round(max(geval), 3)]
                                  if geval else None),
            })
        entries.sort(key=lambda e: -(e["mean_pec50"] if e["mean_pec50"] is not None else -99))
        pack["positions"][p] = {"n_distinct": len(groups), "groups": entries}
    return pack


# --------------------------------------------------------------------------
# Citation verification (the anti-hallucination gate)
# --------------------------------------------------------------------------

def _close(a: float, b: float, rel: float, abs_: float) -> bool:
    return abs(a - b) <= max(rel * abs(b), abs_)


NUM_CTX_RE = re.compile(r"(?<![\w.])\d+(?:\.\d+)?\s*(?:nM|µM|uM|pEC50|倍数?|fold|-fold)?")
ACT_SUFFIX = {"nm", "µm", "um"}
FOLD_SUFFIX = {"倍", "fold", "-fold"}


def _classify_numbers(statement: str) -> list[tuple[float, str]]:
    """Extract numeric tokens with their semantic suffix: (value, kind) where
    kind is 'act' (nM/µM), 'pec' (pEC50), 'fold' (倍/fold) or 'bare'."""
    out = []
    for m in NUM_CTX_RE.finditer(statement):
        tok = m.group(0)
        num_part = re.match(r"\d+(?:\.\d+)?", tok).group(0)
        suffix = tok[len(num_part):].strip().lower()
        if suffix in ACT_SUFFIX:
            kind = "act"
        elif suffix == "pec50":
            kind = "pec"
        elif suffix in FOLD_SUFFIX:
            kind = "fold"
        else:
            kind = "bare"
        out.append((float(num_part), kind))
    return out


def verify_claim(claim: dict, pack: dict) -> tuple[bool, list[str]]:
    """Return (verified, problems). A claim must cite existing ids and every
    numeric token in its statement must be traceable to the pack. Tokens are
    routed by their suffix (nM -> activity values, 倍/fold -> fold changes,
    pEC50 -> p-values, bare -> counts/ids/values) so e.g. a fabricated
    '12.3-fold' cannot slip through by matching a 12.0 nM EC50."""
    problems: list[str] = []
    ids = claim.get("evidence_ids") or []
    if not ids:
        return False, ["no evidence_ids"]
    known = {c["id"]: c for c in pack["compounds"]}
    for i in ids:
        if i not in known:
            problems.append(f"unknown compound id cited: {i!r}")
    if problems:
        return False, problems

    act = {c["id"]: c["ec50_nm"] for c in pack["compounds"] if c["ec50_nm"] is not None}
    pec = {c["id"]: c["pec50"] for c in pack["compounds"] if c["pec50"] is not None}

    # 1) data_refs must match the pack
    ref_vals: list[float] = []
    for ref in claim.get("data_refs") or []:
        rid, kind, val = ref.get("id"), ref.get("kind"), ref.get("value")
        if rid not in known:
            problems.append(f"data_ref id {rid!r} not in pack")
        elif kind == "ec50_nm":
            if rid not in act or not _close(float(val), act[rid], 0.02, 0.05):
                problems.append(f"data_ref ec50_nm mismatch for {rid}: "
                                f"claimed {val}, pack has {act.get(rid)}")
            else:
                ref_vals.append(float(val))
        elif kind == "pec50":
            if rid not in pec or not _close(float(val), pec[rid], 0.0, 0.055):
                problems.append(f"data_ref pec50 mismatch for {rid}: "
                                f"claimed {val}, pack has {pec.get(rid)}")
            else:
                ref_vals.append(float(val))
        else:
            problems.append(f"data_ref kind {kind!r} not recognised")
    if problems:
        return False, problems

    # 2) numeric tokens in the statement must be covered, routed by suffix
    idset = set(ids)
    cited_act = [act[i] for i in ids if i in act]
    cited_pec = [pec[i] for i in ids if i in pec]
    fold_reg = ([pack["activity_span"]["fold"]] if "activity_span" in pack else []) \
        + [a / b for a in cited_act for b in cited_act if b > 0] + list(ref_vals)
    act_reg = cited_act + list(ref_vals)
    if "activity_span" in pack:
        s = pack["activity_span"]
        act_reg += [s["min_nm"], s["max_nm"]]
    if "best" in pack:
        act_reg += [pack["best"]["ec50_nm"], pack["worst"]["ec50_nm"]]
    pec_reg = cited_pec + list(ref_vals)
    if "activity_stats" in pack:
        pec_reg += [pack["activity_stats"]["pec50_median"],
                    pack["activity_stats"]["pec50_mean"]]
    for pd in pack.get("positions", {}).values():
        pec_reg += [g["mean_pec50"] for g in pd["groups"] if g["mean_pec50"] is not None]
    counts = {float(pack["n_total"]), float(pack["n_with_activity"]),
              float(pack["n_no_activity"]), float(len(pack.get("positions", {})))}
    for pd in pack.get("positions", {}).values():
        counts.add(float(pd["n_distinct"]))
        counts.update(float(g["n"]) for g in pd["groups"])

    def hit(reg, num, rel, abs_):
        return any(_close(num, v, rel, abs_) for v in reg)

    for num, kind in _classify_numbers(claim.get("statement", "")):
        if num == int(num) and str(int(num)) in idset:
            continue  # a cited compound id
        if kind == "act":
            ok = hit(act_reg, num, 0.02, 0.05)
        elif kind == "pec":
            ok = hit(pec_reg, num, 0.0, 0.055)
        elif kind == "fold":
            ok = hit(fold_reg, num, 0.05, 0.3)
        else:  # bare: count, or a value written without unit (tight tolerance)
            ok = (num in counts
                  or hit(act_reg, num, 0.01, 0.01)
                  or hit(pec_reg, num, 0.0, 0.011)
                  or hit(fold_reg, num, 0.01, 0.05))
        if not ok:
            problems.append(f"unverifiable number in statement: "
                            f"{num:g}{'' if kind == 'bare' else ' (' + kind + ')'}")
    return (not problems), problems


# --------------------------------------------------------------------------
# LLM access (OpenAI-compatible chat endpoint, stdlib only)
# --------------------------------------------------------------------------

def _chat(messages: list[dict], base_url: str, api_key: str, model: str,
          temperature: float = 0.2, timeout: int = 240) -> str:
    url = base_url.rstrip("/") + "/chat/completions"
    body = json.dumps({"model": model, "messages": messages,
                       "temperature": temperature}).encode("utf-8")
    req = urllib.request.Request(url, data=body, headers={
        "Content-Type": "application/json",
        "Authorization": f"Bearer {api_key}",
    })
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    return data["choices"][0]["message"]["content"]


def _extract_json(text: str) -> dict | None:
    text = re.sub(r"```(?:json)?", "", text).strip()
    start = text.find("{")
    if start < 0:
        return None
    depth = 0
    for i in range(start, len(text)):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                try:
                    return json.loads(text[start:i + 1])
                except json.JSONDecodeError:
                    return None
    return None


SYSTEM_PROMPT = """你是资深药物化学家，基于从专利中提取的 SAR 数据集撰写构效关系总结。
你只会看到一个"证据包" JSON（化合物 id、R 基团、EC50/pEC50、统计量）。严格遵守：
1. 只能引用证据包中存在的化合物 id，绝不发明 id。
2. statement 中出现的每一个数字（EC50、pEC50、n、倍数、差值）都必须给出来源：要么列入
   data_refs（某化合物的 ec50_nm/pec50），要么是证据包顶层统计（n_total、fold 等），
   要么是已引用化合物之间的 EC50 倍数 / pEC50 差值。写不出依据的数字一律不要写。
3. kind=trend 的结论至少需要 2 个 evidence_ids；单化合物事实用 kind=fact；
   无活性/数据缺失用 kind=gap；质控相关用 kind=qc；超出数据的推测用 kind=design
   且 confidence=low。
4. 不要外推证据包之外的组合；不要输出 SMILES 之外的化学式缩写猜测；不要用 markdown。
5. 输出且仅输出一个 JSON 对象，schema：
{"headline": "...", "claims": [{"kind": "fact|trend|gap|qc|design",
  "statement": "...", "evidence_ids": ["..."], "confidence": "high|medium|low",
  "data_refs": [{"kind": "ec50_nm|pec50", "id": "...", "value": 0}]}],
 "suggested_next": ["..."]}
data_refs 必须逐个列出 statement 中出现的所有非统计量数字。"""


def call_llm(pack: dict, cfg: dict, lang: str) -> dict | None:
    """Ask the LLM for a claim list; returns parsed JSON or None."""
    user = (f"证据包（全部数据，真实可信，直接以此为唯一事实来源）\n"
            f"```json\n{json.dumps(pack, ensure_ascii=False)}\n```\n"
            f"请输出 SAR 总结 JSON（语言：{'中文' if lang == 'zh' else 'English'}）。"
            f"覆盖：总体活性概况、R 基团趋势（按位点）、活性悬崖或异常、数据缺口、"
            f"下一步设计建议（suggested_next 不含数字）。")
    messages = [{"role": "system", "content": SYSTEM_PROMPT},
                {"role": "user", "content": user}]
    content = _chat(messages, cfg["base_url"], cfg["api_key"], cfg["model"])
    parsed = _extract_json(content)
    if parsed is None:
        content = _chat(messages + [{"role": "user", "content":
                          "上一次输出无法解析为 JSON。请重新输出，仅含一个 JSON 对象。"}],
                        cfg["base_url"], cfg["api_key"], cfg["model"])
        parsed = _extract_json(content)
    return parsed


# --------------------------------------------------------------------------
# Deterministic digest (used as fallback and always included in the report)
# --------------------------------------------------------------------------

def digest_claims(pack: dict) -> list[dict]:
    def ref(kind, cid, val):
        return {"kind": kind, "id": cid, "value": val}

    out: list[dict] = []
    if "best" in pack:
        b, w = pack["best"], pack["worst"]
        ev = [b["id"]] if b["id"] == w["id"] else [b["id"], w["id"]]
        out.append({"kind": "fact", "label": "potency_range",
                    "statement": (f"活性最优为 Ex {b['id']}（EC50 {b['ec50_nm']:g} nM，"
                                  f"pEC50 {b['pec50']:.1f}），最弱为 Ex {w['id']}"
                                  f"（EC50 {w['ec50_nm']:g} nM），整体跨度约 "
                                  f"{pack['activity_span']['fold']:.1f} 倍。"),
                    "evidence_ids": ev,
                    "data_refs": [ref("ec50_nm", b["id"], b["ec50_nm"]),
                                  ref("pec50", b["id"], b["pec50"]),
                                  ref("ec50_nm", w["id"], w["ec50_nm"])],
                    "confidence": "high"})
    for p, pd in pack.get("positions", {}).items():
        top = next((g for g in pd["groups"] if g["mean_pec50"] is not None
                    and g["n_with_activity"] >= 1), None)
        if top and pd["n_distinct"] >= 2:
            out.append({"kind": "trend", "label": "r_group_best",
                        "statement": (f"{p} 位点中 `{top['r']}` 组平均 pEC50 最高"
                                      f"（{top['mean_pec50']:.2f}，n={top['n']}），"
                                      f"涉及 Ex {', '.join(top['ids'][:6])}。"),
                        "evidence_ids": top["ids"][:6],
                        "data_refs": [],
                        "confidence": "high"})
    if pack.get("n_no_activity"):
        out.append({"kind": "gap", "label": "no_activity",
                    "statement": (f"共 {pack['n_no_activity']} 个化合物无活性数据"
                                  f"（ND）：Ex {', '.join(pack['no_activity_ids'])}，"
                                  f"总结未包含它们。"),
                    "evidence_ids": pack["no_activity_ids"],
                    "data_refs": [], "confidence": "high"})
    if pack.get("flagged_ids"):
        out.append({"kind": "qc", "label": "flagged",
                    "statement": (f"质控异常标记化合物：Ex "
                                  f"{', '.join(pack['flagged_ids'])}，引用其数据前需人工复核。"),
                    "evidence_ids": pack["flagged_ids"],
                    "data_refs": [], "confidence": "high"})
    return out


# --------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------

LABELS_ZH = {"fact": "事实", "trend": "趋势", "gap": "数据缺口", "qc": "质控",
             "design": "设计假设", "potency_range": "活性范围",
             "r_group_best": "位点最优取代基", "no_activity": "无活性数据",
             "flagged": "质控异常"}


def render_md(pack: dict, headline: str | None, claims: list[dict],
              rejected: list[dict], suggested: list[str], meta: dict,
              lang: str) -> str:
    zh = lang == "zh"
    L = (lambda k: LABELS_ZH.get(k, k)) if zh else (lambda k: k)
    sec = (lambda a, b: a) if zh else (lambda a, b: b)
    out = [f"# {sec('SAR 总结', 'SAR summary')} — {pack['target']}", "",
           f"> {sec('生成方式', 'How this was produced')}: "
           f"{meta['llm_desc']}; {sec('所有结论经引用校验', 'all claims citation-checked')}"
           f" ({meta['n_verified']} {sec('条通过', 'passed')} / "
           f"{len(rejected)} {sec('条隔离', 'quarantined')}) · {meta['date']}", "",
           f"- dataset: `{meta['dataset']}`"
           + (f" | rgroup_table: `{meta['rgroup_table']}`" if meta["rgroup_table"] else ""),
           f"- {sec('化合物', 'compounds')}: {pack['n_total']} "
           f"({pack['n_with_activity']} {sec('有活性数据', 'with activity')})", ""]
    out += [f"## {sec('总览', 'Headline')}", "", headline or "-", ""]
    if claims:
        out += [f"## {sec('主要结论', 'Key claims')} "
                f"({sec('已通过引用校验', 'citation-verified')})", ""]
        for c in claims:
            ev = ", ".join(f"Ex {i}" for i in c["evidence_ids"])
            lab = L(c.get("label") or c.get("kind", ""))
            out.append(f"- **[{lab} / {c.get('confidence', '?')}]** "
                       f"{c['statement']} *({sec('证据', 'evidence')}: {ev})*")
        out.append("")
    if suggested:
        out += [f"## {sec('下一步建议', 'Suggested next steps')}", ""]
        out += [f"- {s}" for s in suggested] + [""]
    if rejected:
        out += [f"## {sec('未通过引用校验（已隔离，需人工确认）', 'Quarantined claims (failed verification)')}", ""]
        for c in rejected:
            out.append(f"- ~~{c['statement']}~~")
            out.append(f"  - {sec('原因', 'reason')}: {', '.join(c['problems'])}")
        out.append("")
    out += [f"## {sec('附录：R 基团位点统计', 'Appendix: R-group position stats')}", ""]
    for p, pd in pack.get("positions", {}).items():
        out.append(f"### {p} ({pd['n_distinct']} {sec('个取代基', 'distinct R-groups')})")
        out += ["| R | n | EC50 range (nM) | mean pEC50 | examples |", "|---|---|---|---|---|"]
        for g in pd["groups"][:12]:
            rng = (f"{g['ec50_nm_range'][0]:g}–{g['ec50_nm_range'][1]:g}"
                   if g.get("ec50_nm_range") else "-")
            out.append(f"| `{g['r']}` | {g['n']} | {rng} | "
                       f"{g['mean_pec50'] if g['mean_pec50'] is not None else '-'} "
                       f"| {', '.join(g['ids'][:5])} |")
        out.append("")
    out += [f"## {sec('附录：化合物活性表', 'Appendix: compound activity table')}", "",
            "| Ex | EC50 (nM) | pEC50 | R groups | flags |", "|---|---|---|---|---|"]
    for c in pack["compounds"]:
        rg = "; ".join(f"{k}=`{v}`" for k, v in c["R"].items()) or "-"
        ec = f"{c['ec50_nm']:g}" if c["ec50_nm"] is not None else "ND"
        pc = f"{c['pec50']:.2f}" if c["pec50"] is not None else "-"
        out.append(f"| {c['id']} | {ec} | {pc} "
                   f"| {rg} | {', '.join(c.get('flags', [])) or '-'} |")
    out += ["", "---", sec("*本报告由 patent-sar-mining sar.summarize 自动生成；"
                           "结论仅覆盖实施例层面（不含权利要求 Markush 空间），"
                           "引用前请对照 output/*_check.html 复核结构归属。*",
                           "*Auto-generated by patent-sar-mining sar.summarize; "
                           "claims cover patent examples only (no Markush claim space); "
                           "verify structure assignment against output/*_check.html.*")]
    return "\n".join(out)


# --------------------------------------------------------------------------
# Orchestration + CLI
# --------------------------------------------------------------------------

def summarize(dataset: str | None, target: str, outdir: str,
              rgroup_table: str | None = None, activity_col: str | None = None,
              id_col: str | None = None, core: str | None = None,
              llm: str | None = None, base_url: str | None = None,
              model: str | None = None, api_key: str | None = None,
              lang: str = "zh") -> dict:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    paths = [Path(p) for p in (dataset, rgroup_table) if p]
    # recompute R-groups from the dataset when no table was supplied (needs RDKit)
    computed = None
    if not rgroup_table and dataset:
        from . import rgroup as _rg  # lazy: needs RDKit
        comps = _rg.load_dataset(dataset,
                                 smiles_col="canonical_smiles",
                                 id_col=id_col or "example_no",
                                 activity_col=activity_col)
        rows0, positions = _rg.rgroup_table(comps, core_smiles=core)
        names = _rg.position_names(positions)
        computed = outdir / "rgroup_table.csv"
        _rg_rows_to_csv(rows0, names, computed)
        paths.append(computed)

    rows = load_activity_rows(paths, activity_col, id_col)
    pack = build_pack(rows, target)

    cfg = {"base_url": base_url or os.environ.get("SAR_LLM_BASE_URL",
           "https://open.bigmodel.cn/api/paas/v4"),
           "model": model or os.environ.get("SAR_LLM_MODEL", "glm-4.6"),
           "api_key": api_key or os.environ.get("SAR_LLM_API_KEY", "")}
    use_llm = llm != "none" and bool(cfg["api_key"])
    llm_desc, llm_json = None, None
    if use_llm:
        try:
            llm_json = call_llm(pack, cfg, lang)
            llm_desc = f"LLM {cfg['model']} via {cfg['base_url']}"
        except Exception as e:  # noqa: BLE001 - degrade to digest, never crash batch runs
            print(f"[summarize] LLM call failed ({type(e).__name__}: {e}); "
                  f"falling back to digest-only", file=sys.stderr)
            llm_json = None
    else:
        llm_desc = "deterministic digest (no LLM)"

    digest = digest_claims(pack)
    raw_claims = list(llm_json.get("claims", [])) if llm_json else []
    for c in raw_claims:
        c.setdefault("label", c.get("kind", ""))
    verified, rejected = [], []
    for c in raw_claims:
        ok, problems = verify_claim(c, pack)
        (verified if ok else rejected).append({**c, "problems": problems})
    # digest claims are pack-derived; verify anyway and keep only passing ones
    for c in digest:
        ok, problems = verify_claim(c, pack)
        if ok:
            verified.append({**c, "problems": [], "source": "digest"})
    # digest headline
    headline = (llm_json or {}).get("headline")
    if headline is not None:
        ok, problems = verify_claim({"statement": headline,
                                     "evidence_ids": list(
                                         {i for c in verified for i in c["evidence_ids"]} or
                                         [pack["compounds"][0]["id"]]),
                                     "data_refs": []}, pack)
        if not ok:
            headline = None
    if headline is None and "best" in pack:
        headline = digest_claims(pack)[0]["statement"]

    result = {"meta": {"date": date.today().isoformat(), "target": target,
                       "dataset": dataset, "rgroup_table": rgroup_table,
                       "llm_desc": llm_desc,
                       "n_verified": len(verified)},
              "evidence_pack": pack,
              "headline": headline,
              "claims": verified,
              "quarantined": rejected,
              "suggested_next": (llm_json or {}).get("suggested_next", [])}
    (outdir / "sar_summary.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=1))
    md = render_md(pack, headline, verified, rejected, result["suggested_next"],
                   {"date": result["meta"]["date"], "dataset": dataset,
                    "rgroup_table": rgroup_table, "llm_desc": llm_desc,
                    "n_verified": len(verified)}, lang)
    (outdir / "sar_summary.md").write_text(md)
    result["outdir"] = str(outdir)
    return result


def _rg_rows_to_csv(rows0, names, path: Path):
    """Write an R-group table in the same format as sar.run.write_rgroup_csv."""
    positions = sorted(names)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "activity_ec50_nm"] + [names[p] for p in positions] + ["smiles"])
        for r in rows0:
            w.writerow([r["cid"],
                        "" if r["activity"] is None else r["activity"]] +
                       [r["groups"].get(p, "") for p in positions] + [r["smiles"]])


def _selftest(pack: dict) -> int:
    """Offline check of the citation gate: good claims pass, bad ones don't."""
    act = {c["id"]: c["ec50_nm"] for c in pack["compounds"] if c["ec50_nm"] is not None}
    best_id = min(act, key=act.get)
    ids = sorted(act, key=lambda i: (len(i), i))
    pair = (ids[0], ids[-1])
    real_fold = act[pair[1]] / act[pair[0]]
    cases = [
        ({"kind": "fact", "statement":
            f"最优化合物为 Ex {best_id}，EC50 {act[best_id]:g} nM。",
          "evidence_ids": [best_id],
          "data_refs": [{"kind": "ec50_nm", "id": best_id, "value": act[best_id]}]},
         True, "good claim (id + EC50)"),
        ({"kind": "trend", "statement":
            f"Ex {pair[1]} 比 Ex {pair[0]} 弱约 {real_fold:.1f} 倍。",
          "evidence_ids": [pair[0], pair[1]], "data_refs": []},
         True, "good fold claim between cited pair"),
        ({"kind": "trend", "statement": "完全不存在的化合物说明活性更好。",
          "evidence_ids": ["NO_SUCH_ID"], "data_refs": []},
         False, "bad id"),
        ({"kind": "trend", "statement": "该取代基使活性提升 12.3 倍。",
          "evidence_ids": ids[:2], "data_refs": []},
         False, "fabricated fold number"),
        ({"kind": "fact", "statement":
            f"Ex {ids[0]} 的 EC50 是 {act[ids[0]] * 1.5:g} nM。",
          "evidence_ids": [ids[0]], "data_refs": []},
         False, "wrong EC50 value"),
    ]
    fails = 0
    for claim, expect, name in cases:
        ok, problems = verify_claim(claim, pack)
        status = "PASS" if ok == expect else "FAIL"
        if ok != expect:
            fails += 1
        print(f"  [{status}] selftest {name}: verified={ok} {problems}")
    # LLM response parsing (offline): fenced / prose-wrapped JSON
    json_cases = [
        ('```json\n{"headline": "ok", "claims": []}\n```', True, "fenced json"),
        ('好的，以下是总结：\n{"headline": "ok", "claims": []}\n如上。', True,
         "prose-wrapped json"),
        ('没有 JSON 在这里', False, "no json"),
    ]
    for text, expect, name in json_cases:
        parsed = _extract_json(text) is not None
        status = "PASS" if parsed == expect else "FAIL"
        if parsed != expect:
            fails += 1
        print(f"  [{status}] selftest {name}")
    return fails


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", default=None, help="mining dataset CSV")
    ap.add_argument("--rgroup-table", default=None,
                    help="R-group table CSV (sar.run output); if omitted and "
                         "--dataset is given, recomputed via RDKit")
    ap.add_argument("--activity-col", default=None,
                    help="activity column (auto-detected: first ec50*/activity* col)")
    ap.add_argument("--id-col", default=None)
    ap.add_argument("--core", default=None, help="scaffold SMILES for R-group "
                    "recomputation (same as sar.run --core)")
    ap.add_argument("--target", default="TARGET")
    ap.add_argument("--outdir", default="output/sar")
    ap.add_argument("--llm", choices=["api", "none"], default="api",
                    help="'api' (default, needs SAR_LLM_API_KEY) or 'none' for "
                         "deterministic digest only")
    ap.add_argument("--base-url", default=None,
                    help="OpenAI-compatible endpoint (env SAR_LLM_BASE_URL; "
                         "default Zhipu bigmodel)")
    ap.add_argument("--model", default=None, help="env SAR_LLM_MODEL")
    ap.add_argument("--api-key", default=None, help="env SAR_LLM_API_KEY")
    ap.add_argument("--lang", choices=["zh", "en"], default="zh")
    ap.add_argument("--selftest", action="store_true",
                    help="run the citation-gate unit test on the loaded pack "
                         "and exit (no LLM call)")
    a = ap.parse_args()

    if not a.dataset and not a.rgroup_table:
        ap.error("--dataset and/or --rgroup-table required")
    res = summarize(a.dataset, a.target, a.outdir, rgroup_table=a.rgroup_table,
                    activity_col=a.activity_col, id_col=a.id_col, core=a.core,
                    llm=a.llm, base_url=a.base_url, model=a.model,
                    api_key=a.api_key, lang=a.lang)
    if a.selftest:
        fails = _selftest(res["evidence_pack"])
        sys.exit(1 if fails else 0)
    print(f"claims={len(res['claims'])} verified | "
          f"{len(res['quarantined'])} quarantined | -> "
          f"{res['outdir']}/sar_summary.{{json,md}}")


if __name__ == "__main__":
    main()
