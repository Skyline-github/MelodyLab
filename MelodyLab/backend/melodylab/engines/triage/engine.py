"""立项五问引擎：并发调取各数据源，按 rules.yaml 阈值生成评估报告。

移植自 target-triage（app/engine/five_questions.py），仅调整包内导入。
"""
from __future__ import annotations

import asyncio
import uuid
from functools import lru_cache
from pathlib import Path

import yaml

from . import charts
from . import sources as services
from .config import FIVE_QUESTIONS
from .models import (AssessmentReport, Evidence, GREEN, QuestionReport, RED,
                     SourceLink, UNKNOWN, YELLOW)
from .sources import oncokb

RULES_PATH = Path(__file__).resolve().parent / "rules.yaml"


@lru_cache(maxsize=1)
def rules() -> dict:
    return yaml.safe_load(RULES_PATH.read_text())


# ---------------------------------------------------------------- Q1 正确的靶点
async def q_target(symbol: str, ensembl_id: str, indication: str) -> QuestionReport:
    r = rules()["target"]
    ev: list[Evidence] = []
    manual: list[Evidence] = [
        Evidence(source="DisGeNET", check="基因-疾病关联综合数据",
                 summary="需注册核查：学术免费申请，商业付费",
                 links=[SourceLink(label="disgenet.org", url="https://www.disgenet.org")],
                 manual=True),
        Evidence(source="TTD", check="已知和在研治疗靶点",
                 summary="免费开放数据库，交叉验证靶点-疾病关系",
                 links=[SourceLink(label="db.idrblab.net/ttd", url="https://db.idrblab.net/ttd/")],
                 manual=True),
    ]

    try:
        assoc = await services.opentargets.get_associations(ensembl_id)
    except Exception as exc:  # noqa: BLE001
        assoc = []
        ev.append(Evidence(source="Open Targets Platform", check="靶点-疾病关联",
                           status=UNKNOWN, summary=f"查询失败：{exc}"))

    if assoc:
        focus = None
        if indication:
            focus = next((a for a in assoc
                          if indication.lower() in a["disease"].lower()), None)
        top = focus or assoc[0]
        label = f"{top['disease']}（整体 {top['overall']}，遗传学 {top['genetic']}）"
        if top["overall"] >= r["association_strong"]:
            status, summary = GREEN, f"关联证据充分：{label}"
        elif top["overall"] >= r["association_weak"]:
            status, summary = YELLOW, f"关联证据中等：{label}"
        else:
            status, summary = RED, f"关联证据薄弱：{label}"
        if focus is None and indication:
            summary += f"；注意：Top 关联疾病中未直接匹配“{indication}”"
        ev.append(Evidence(
            source="Open Targets Platform", check="靶点-疾病关联强度",
            status=status, summary=summary,
            detail={"top_diseases": assoc[:10]},
            chart=charts.hbar(
                [a["disease"] for a in assoc[:10]],
                "整体关联分数", [a["overall"] for a in assoc[:10]],
                second_label="其中遗传学", second_values=[a["genetic"] for a in assoc[:10]],
                max_value=1.0),
            links=[SourceLink(
                label=f"OT 靶点页 ({ensembl_id})",
                url=f"https://platform.opentargets.org/target/{ensembl_id}")],
        ))

        # 人类遗传学信号
        max_genetic = max((a["genetic"] for a in assoc), default=0)
        gwas = await services.gwas.get_associations_for_gene(symbol)
        gwas_rows = gwas.get("associations", []) if gwas.get("available") else []
        if max_genetic >= r["genetic_strong"] or gwas_rows:
            ev.append(Evidence(
                source="GWAS Catalog / OT 遗传学", check="人类遗传学信号",
                status=GREEN,
                summary=f"存在人类遗传学证据（OT 遗传学分型最高 {max_genetic}；"
                        f"GWAS Catalog 命中 {len(gwas_rows)} 条）",
                detail={"gwas_catalog": gwas_rows[:8],
                        "gwas_note": gwas.get("note", "")},
                links=[SourceLink(label="GWAS Catalog",
                                  url=f"https://www.ebi.ac.uk/gwas/genes/{symbol}")],
            ))
        else:
            ev.append(Evidence(
                source="GWAS Catalog / OT 遗传学", check="人类遗传学信号",
                status=UNKNOWN,
                summary="未发现人类遗传学信号——创新靶点缺人体数据属常态，"
                        "不能直接判死刑；应转向动物模型/机制研究等其他证据类别",
                detail={"gwas_note": gwas.get("note", "")},
            ))
    return QuestionReport(key="target", title=FIVE_QUESTIONS[0]["title"],
                          question=FIVE_QUESTIONS[0]["question"],
                          evidences=ev, manual_checklist=manual)


# ---------------------------------------------------------------- Q2 正确的组织
async def q_tissue(symbol: str, ensembl_id: str = "") -> QuestionReport:
    r = rules()["tissue"]
    ev: list[Evidence] = []
    manual: list[Evidence] = [
        Evidence(source="KEGG / Reactome", check="通路信息",
                 summary="Reactome 开源可用；KEGG 学术免费、商业需许可",
                 links=[SourceLink(label="Reactome", url="https://reactome.org"),
                        SourceLink(label="KEGG", url="https://www.genome.jp/kegg/")],
                 manual=True),
        Evidence(source="STRING", check="蛋白互作网络全景",
                 summary="查互作伙伴以判断通路级联影响",
                 links=[SourceLink(label=f"STRING ({symbol})",
                                   url=f"https://string-db.org/network/9606.{symbol}")],
                 manual=True),
    ]

    try:
        hpa = await services.protein_atlas.get_tissue_expression(symbol, ensembl_id=ensembl_id)
    except Exception as exc:  # noqa: BLE001
        hpa = {"found": False, "note": f"HPA 数据加载失败：{exc}"}

    if hpa.get("found"):
        table = hpa.get("tissue_table") or []
        detected = [t["tissue"] for t in table if t.get("protein_level", "Not detected")
                    not in ("", "Not detected")]
        high = [t["tissue"] for t in table
                if t.get("protein_level") in ("High", "Medium")]
        kws = r["cns_tissue_keywords"]
        cns_only = bool(detected) and all(
            any(k in t.lower() for k in kws) for t in detected)
        status = YELLOW if cns_only else GREEN
        parts = [f"RNA 表达：{hpa['rna_tissue_specificity'] or '—'} / "
                 f"{hpa['rna_tissue_distribution'] or '—'}",
                 f"蛋白表达：{hpa['protein_tissue_specificity'] or '—'} / "
                 f"{hpa['protein_tissue_distribution'] or '—'}"]
        summary = "；".join(parts)
        if table:
            summary += f"；{len(detected)} 种组织检出蛋白（中高表达 {len(high)} 种）"
        if hpa.get("top_rna_tissues"):
            tops = ", ".join(t["tissue"] for t in hpa["top_rna_tissues"][:4])
            summary += f"；RNA 最富集组织：{tops}"
        if cns_only:
            summary += "；表达局限于中枢神经系统——分子过血脑屏障是递送前提，基本锁定小分子路线"
        tissue_chart = None
        if table:
            shown = table[:12]
            tissue_chart = charts.hbar(
                [t["tissue"] for t in shown],
                "蛋白表达强度（IHC 等级）",
                [{"High": 3, "Medium": 2, "Low": 1}.get(t["protein_level"], 0)
                 for t in shown],
                colors=[charts.LEVEL_COLORS.get(t["protein_level"], "#d8dedb")
                        for t in shown],
                max_value=3)
        ev.append(Evidence(
            source="Human Protein Atlas", check="组织分布与表达特异性",
            status=status, summary=summary,
            detail={"tissue_table": table[:40] if table else [],
                    "table_note": ("组织级蛋白 IHC Level（HPA v22，按等级排序）"
                                   if table else
                                   "完整组织表未就绪：v22 归档文件正在后台下载"
                                   "（约 6MB，慢网络下可能需要多跑几次评估才完成），"
                                   "本次以 HPA 单基因端点数据为准，下载完成后自动带全表"),
                    "top_rna_tissues": hpa.get("top_rna_tissues", []),
                    "single_cell_specificity": hpa.get("single_cell_specificity", ""),
                    "subcellular_location": hpa.get("subcellular_location", []),
                    "secretome": hpa.get("secretome_location", "")},
            chart=tissue_chart,
            links=[SourceLink(label=f"HPA ({symbol})",
                              url=f"https://www.proteinatlas.org/{ensembl_id}/tissue")],
        ))
    else:
        ev.append(Evidence(source="Human Protein Atlas", check="组织分布与表达特异性",
                           status=UNKNOWN, summary=hpa.get("note", "无数据")))

    try:
        uni = await services.uniprot.get_protein(symbol)
    except Exception as exc:  # noqa: BLE001
        uni = None
        ev.append(Evidence(source="UniProt", check="蛋白功能与亚细胞定位",
                           status=UNKNOWN, summary=f"查询失败：{exc}"))
    if uni:
        loc = "；".join(uni["subcellular_location"])[:200] or "未注释"
        ev.append(Evidence(
            source="UniProt", check="蛋白功能与亚细胞定位",
            status=GREEN,
            summary=f"{uni['protein_name']}：{(uni['function'] or [''])[0][:150]}",
            detail={"subcellular_location": loc, "uniprot_id": uni["uniprot_id"]},
            links=[SourceLink(label=f"UniProt {uni['uniprot_id']}",
                              url=f"https://www.uniprot.org/uniprotkb/{uni['uniprot_id']}")],
        ))
    return QuestionReport(key="tissue", title=FIVE_QUESTIONS[1]["title"],
                          question=FIVE_QUESTIONS[1]["question"],
                          evidences=ev, manual_checklist=manual)


# ---------------------------------------------------------------- Q3 正确的安全性
async def q_safety(symbol: str, ensembl_id: str) -> QuestionReport:
    r = rules()["safety"]
    ev: list[Evidence] = []

    try:
        ko = await services.opentargets.get_mouse_knockout(ensembl_id)
    except Exception as exc:  # noqa: BLE001
        ko = {"has_knockout_data": False, "lethal_phenotype": False,
              "models": [], "error": str(exc)}
    if ko["has_knockout_data"]:
        if ko["lethal_phenotype"]:
            status, summary = RED, ("小鼠敲除记录到致死类表型"
                                    f"（如：{', '.join(ko['lethal_examples'][:2])}）——"
                                    "提示该基因在正常生理中可能有不可替代的功能，"
                                    "抑制它的安全窗需要重点评估")
        else:
            status, summary = GREEN, (f"小鼠敲除共 {ko['phenotype_count']} 条表型记录，"
                                      "无致死类表型")
        ev.append(Evidence(
            source="Open Targets（小鼠表型）", check="必需基因判定（敲除致死性）",
            status=status, summary=summary, detail={"phenotypes": ko["phenotypes"][:30]},
        ))
    else:
        ev.append(Evidence(source="Open Targets（小鼠表型）", check="必需基因判定（敲除致死性）",
                           status=UNKNOWN, summary="无小鼠敲除数据，请用 MGI 手动核查"))

    try:
        overview = await services.opentargets.get_target_overview(ensembl_id)
        lof_z = overview["lof_z"]
        is_essential = overview["is_essential"]
    except Exception:  # noqa: BLE001
        lof_z, is_essential = None, None
    if is_essential:
        ev.append(Evidence(source="Open Targets", check="必需基因判定（isEssential）",
                           status=YELLOW,
                           summary="Open Targets 标记为必需基因——结合敲除数据与 DepMap 证据，安全窗可能很窄"))
    if lof_z is not None:
        if lof_z >= r["lof_z_strong"]:
            status = RED
            summary = (f"LOF 约束 z = {lof_z:.2f}（≥{r['lof_z_strong']}，强约束）——"
                       "人群中天然 LOF 携带者可能有严重表型，抑制该靶点的治疗窗可能非常窄")
        elif lof_z >= r["lof_z_moderate"]:
            status = YELLOW
            summary = f"LOF 约束 z = {lof_z:.2f}，中等约束，建议到 gnomAD 核查 LOF 携带者表型"
        else:
            status = GREEN
            summary = (f"LOF 约束 z = {lof_z:.2f}，弱约束——天然功能缺失携带者表型大概率良性，"
                       "是最强的安全性背书")
        ev.append(Evidence(
            source="gnomAD（经 Open Targets）", check="人群 LOF 约束（z 分数）",
            status=status, summary=summary,
            links=[SourceLink(label=f"gnomAD gene ({symbol})",
                              url=f"https://gnomad.broadinstitute.org/gene/{ensembl_id}")],
        ))
    else:
        ev.append(Evidence(source="gnomAD（经 Open Targets）", check="人群 LOF 约束（z 分数）",
                           status=UNKNOWN, summary="无遗传约束数据"))

    try:
        liabilities = await services.opentargets.get_safety(ensembl_id)
    except Exception as exc:  # noqa: BLE001
        liabilities = []
        ev.append(Evidence(source="Open Targets Safety", check="已知安全事件",
                           status=UNKNOWN, summary=f"查询失败：{exc}"))
    n = len(liabilities)
    if n >= r["liabilities_red"]:
        status = RED
    elif n >= r["liabilities_yellow"]:
        status = YELLOW
    else:
        status = GREEN
    events = sorted({l["event"] for l in liabilities})
    by_ds: dict[str, int] = {}
    for l in liabilities:
        ds = l.get("datasource") or "其他"
        by_ds[ds] = by_ds.get(ds, 0) + 1
    ev.append(Evidence(
        source="Open Targets Safety", check="已知安全事件清单",
        status=status,
        summary=(f"{n} 项已报告安全事件（来自 PharmGKB、ToxCast、文献 curated）——"
                 "多为成药性相关信号，需结合给药方式与适应症判断"
                 if n else "无已报告安全事件——但创新靶点也可能只是还没人查过"),
        detail={"events": events},
        chart=charts.doughnut(list(by_ds), list(by_ds.values())) if by_ds else None,
        links=[SourceLink(label=f"OT Safety ({symbol})",
                          url=f"https://platform.opentargets.org/target/{ensembl_id}/safety")],
    ))

    # 同靶点已上市药物的安全记录线索
    try:
        drugs = await services.chembl.get_drugs(symbol)
    except Exception as exc:  # noqa: BLE001
        drugs = {"resolved": False, "drugs": [], "note": str(exc)}
    approved = drugs.get("approved_count", 0) if drugs.get("resolved") else 0
    if drugs.get("resolved") and approved:
        ev.append(Evidence(
            source="ChEMBL", check="同靶点已上市药物安全记录",
            status=GREEN,
            summary=f"该靶点已有 {approved} 个上市药物、{drugs.get('total_count', 0)} 个药物记录，"
                    "安全记录可从上市药物追溯（严格说这是同靶点，家族级需手动扩展）",
            detail={"drugs": [
                {"name": d["name"], "phase": d["max_phase"], "moa": d["moa"]}
                for d in drugs["drugs"][:15]]},
            links=[SourceLink(label="ChEMBL target",
                              url=f"https://www.ebi.ac.uk/chembl/target_report_card/{drugs.get('chembl_id', '')}/")],
        ))
    elif drugs.get("resolved"):
        ev.append(Evidence(source="ChEMBL", check="同靶点已上市药物安全记录",
                           status=UNKNOWN,
                           summary=f"该靶点无上市药物（共 {drugs.get('total_count', 0)} 条药物记录），"
                                   "同家族药物安全记录请按蛋白家族手动核查"))
    else:
        ev.append(Evidence(source="ChEMBL", check="同靶点已上市药物安全记录",
                           status=UNKNOWN, summary=drugs.get("note", "未解析到 ChEMBL 靶点")))
    return QuestionReport(key="safety", title=FIVE_QUESTIONS[2]["title"],
                          question=FIVE_QUESTIONS[2]["question"],
                          evidences=ev)


# ---------------------------------------------------------------- Q4 正确的患者
async def q_patient(symbol: str, indication: str) -> QuestionReport:
    ev: list[Evidence] = []
    manual: list[Evidence] = [
        Evidence(source="FDA 伴随诊断列表", check="有无配套伴随诊断",
                 summary="FDA 已批准/cleared 的伴随诊断器械清单（含对应药物与适应症），"
                         "标志物已有获批检测说明可检测性经过监管验证",
                 links=[SourceLink(label="FDA List of Cleared or Approved CDx",
                                   url="https://www.fda.gov/medical-devices/in-vitro-diagnostics/"
                                       "list-cleared-or-approved-companion-diagnostic-devices-"
                                       "in-vitro-and-imaging-tools")],
                 manual=True),
        Evidence(source="PubMed", check="人群规模（标志物亚组患病率）",
                 summary="检索式：disease + prevalence + biomarker + subgroup；"
                         "标志物把获益人群从 100% 缩到 5% 时，罕见病可能撑不起一个项目",
                 links=[SourceLink(
                     label="PubMed 检索（自动拼接）",
                     url=("https://pubmed.ncbi.nlm.nih.gov/?term="
                          + (f"{indication}+prevalence+biomarker" if indication
                             else "prevalence+biomarker+subgroup")))],
                 manual=True),
        Evidence(source="Clarivate Epidemiology Intelligence", check="人群规模（商业工具）",
                 summary="订阅制；预算有限时用上面的 PubMed 检索替代",
                 links=[SourceLink(label="clarivate.com", url="https://clarivate.com")],
                 manual=True),
    ]

    oncokb_data = await oncokb.get_gene_actionability(symbol)
    if oncokb_data.get("available"):
        n_var = len(oncokb_data.get("variants", []))
        ev.append(Evidence(
            source="OncoKB", check="伴随诊断/精准肿瘤学证据",
            status=GREEN if n_var else UNKNOWN,
            summary=f"OncoKB 收录该基因 {n_var} 个注释变异（含治疗级别信息）",
            detail={"variants": oncokb_data.get("variants", [])[:10]},
            links=[SourceLink(label="OncoKB", url="https://www.oncokb.org")],
        ))
    else:
        ev.append(Evidence(
            source="OncoKB", check="伴随诊断/精准肿瘤学证据",
            status=UNKNOWN, summary=oncokb_data.get("note", "未启用"),
            manual=False,
        ))

    try:
        fda = await services.openfda.get_label_mentions(symbol)
    except Exception as exc:  # noqa: BLE001
        fda = {"total": 0, "labels": [], "note": str(exc)}
    if fda.get("total"):
        ev.append(Evidence(
            source="openFDA", check="药物标签中的靶点/突变提及",
            status=GREEN,
            summary=f"FDA 药物标签中有 {fda['total']} 条提及 {symbol}（适应症/机制层面存在靶向生态）",
            detail={"labels": fda.get("labels", [])[:8]},
        ))
    else:
        ev.append(Evidence(
            source="openFDA", check="药物标签中的靶点/突变提及",
            status=UNKNOWN,
            summary="药物标签中未检索到提及——非肿瘤靶点常见，肿瘤靶点建议重点走 OncoKB/FDA CDx 清单",
            detail={"note": fda.get("note", "")},
        ))
    return QuestionReport(key="patient", title=FIVE_QUESTIONS[3]["title"],
                          question=FIVE_QUESTIONS[3]["question"],
                          evidences=ev, manual_checklist=manual)


# ---------------------------------------------------------------- Q5 正确的商业
async def q_commercial(symbol: str, indication: str) -> QuestionReport:
    r = rules()["commercial"]
    ev: list[Evidence] = []
    manual: list[Evidence] = [
        Evidence(source="FTO 专利分析", check="专利自由空间（Freedom to Operate）",
                 summary="临床前阶段发现的专利可以用有限资源挑战或绕开；"
                         "三期才发现就只能高价谈许可或重复临床——立项前必查",
                 links=[SourceLink(label="Google Patents",
                                   url=f"https://patents.google.com/?q={symbol}")],
                 manual=True),
        Evidence(source="Cortellis / Evaluate Pharma", check="市场预测与交易动态",
                 summary="订阅制商业情报：同赛道药物的市场预测、交易格局",
                 links=[SourceLink(label="cortellis.com", url="https://www.cortellis.com"),
                        SourceLink(label="evaluate.com", url="https://www.evaluate.com")],
                 manual=True),
    ]

    try:
        trials = await services.clinicaltrials.get_trials(symbol, indication)
    except Exception as exc:  # noqa: BLE001
        trials = {"total": 0, "status_dist": {}, "terminated": [],
                  "top_conditions": [], "note": str(exc)}
    total = trials.get("total", 0)
    scope = f"（限定“{indication}”）" if indication else "（不限适应症）"
    if total > r["trials_crowded"]:
        status = YELLOW
        summary = f"ClinicalTrials.gov 已有 {total} 项相关试验{scope}——赛道拥挤，立项需回答差异化问题"
    elif total <= r["trials_empty"]:
        status = YELLOW
        summary = (f"ClinicalTrials.gov 无相关试验{scope}——别人推进了可以借力；"
                   "没人推进可能是空白，也可能是坑：把不确定性拆成可逐步验证的问题")
    else:
        status = GREEN
        summary = f"ClinicalTrials.gov 有 {total} 项相关试验{scope}——赛道有验证但不拥挤"
    trial_charts = []
    sd = trials.get("status_dist") or {}
    if sd:
        trial_charts.append(charts.doughnut(list(sd), list(sd.values())))
    conds = trials.get("top_conditions") or []
    if conds:
        trial_charts.append(charts.hbar(
            [c["condition"] for c in conds[:8]], "试验数",
            [c["count"] for c in conds[:8]]))
    ev.append(Evidence(
        source="ClinicalTrials.gov", check="竞争格局（临床试验数量与终止原因）",
        status=status, summary=summary,
        detail={"status_dist": sd,
                "terminated_samples": trials.get("terminated", [])[:5],
                "top_conditions": conds,
                "note": trials.get("note", "")},
        chart={"panels": trial_charts} if trial_charts else None,
        links=[SourceLink(
            label="ClinicalTrials.gov 检索",
            url=(f"https://clinicaltrials.gov/search?intr={symbol}"
                 + (f"&cond={indication.replace(' ', '+')}" if indication else "")))],
    ))

    try:
        fda = await services.openfda.get_label_mentions(symbol)
    except Exception:  # noqa: BLE001
        fda = {"total": 0, "labels": []}
    if fda.get("total"):
        ev.append(Evidence(
            source="openFDA", check="已上市产品线索（未满足需求参照）",
            status=GREEN,
            summary=f"已有 {fda['total']} 条 FDA 药物标签提及该靶点——评估未满足需求时以此为竞品参照",
            detail={"labels": fda.get("labels", [])[:8]},
        ))
    else:
        ev.append(Evidence(source="openFDA", check="已上市产品线索（未满足需求参照）",
                           status=UNKNOWN,
                           summary="无已上市产品直接提及——若因果证据扎实，FIC 潜力本身就是机遇"))
    return QuestionReport(key="commercial", title=FIVE_QUESTIONS[4]["title"],
                          question=FIVE_QUESTIONS[4]["question"],
                          evidences=ev, manual_checklist=manual)


# ---------------------------------------------------------------- 编排入口
# DATA_SOURCES 定义在 .config（供 availability 与前端展示共用）


async def assess(gene_symbol: str, indication: str = "") -> AssessmentReport:
    """执行五问评估：并发取数 → 组装报告 → 聚合状态与一票否决项。"""
    report = AssessmentReport(
        id=uuid.uuid4().hex[:12],
        gene_symbol=gene_symbol.strip().upper(),
        indication=indication.strip(),
    )

    resolved = await services.opentargets.search_target(report.gene_symbol)
    if not resolved:
        raise ValueError(f"Open Targets 未收录基因符号 “{report.gene_symbol}”，请检查拼写")
    report.ensembl_id = resolved["ensembl_id"]
    report.gene_name = resolved["approved_name"]
    try:
        uni = await services.uniprot.get_protein(report.gene_symbol)
        report.uniprot_id = uni["uniprot_id"] if uni else ""
    except Exception:  # noqa: BLE001
        report.uniprot_id = ""

    results = await asyncio.gather(
        q_target(report.gene_symbol, report.ensembl_id, report.indication),
        q_tissue(report.gene_symbol, report.ensembl_id),
        q_safety(report.gene_symbol, report.ensembl_id),
        q_patient(report.gene_symbol, report.indication),
        q_commercial(report.gene_symbol, report.indication),
        return_exceptions=True,
    )
    for i, res in enumerate(results):
        if isinstance(res, BaseException):
            results[i] = QuestionReport(
                key=FIVE_QUESTIONS[i]["key"], title=FIVE_QUESTIONS[i]["title"],
                question=FIVE_QUESTIONS[i]["question"],
                evidences=[Evidence(source="engine", check="模块执行",
                                    status=UNKNOWN, summary=f"执行异常：{res}")],
            )
    report.questions = list(results)
    report.finalize()
    return report
