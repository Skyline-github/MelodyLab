"""靶点立项评估（立项五问）配置：外部 API 端点、超时、缓存与报告目录。

从 target-triage 项目（https://github.com/Skyline-github/target-triage）移植。
目录挂到 MelodyLab 数据根（settings.data_root），换数据根时行为一致。
"""
from __future__ import annotations

import os

from ...config import settings as _app_settings

# ---------------------------------------------------------------- 目录
CACHE_DIR = _app_settings.data_root / "cache" / "triage"
REPORTS_DIR = _app_settings.data_root / "triage_reports"
for _d in (CACHE_DIR, REPORTS_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- 外部 API
OPEN_TARGETS_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
GWAS_CATALOG_V2 = "https://api.gwas.scientific-data.io/api/v2"
HUMAN_PROTEIN_ATLAS_TSV = (
    "https://www.proteinatlas.org/download/normal_tissue.tsv.zip"
)
STRING_API = "https://string-db.org/api"
CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"
CLINICALTRIALS_API = "https://clinicaltrials.gov/api/v2/studies"
OPENFDA_DRUG_LABEL = "https://api.fda.gov/drug/label.json"
ONCOKB_API = "https://www.oncokb.org/api/v1"

# OncoKB token：环境变量 > 数据根 oncokb_token.txt（设置页/手动放置均可）
ONCOKB_TOKEN = os.getenv("ONCOKB_TOKEN", "").strip()
if not ONCOKB_TOKEN:
    _tok_file = _app_settings.data_root / "oncokb_token.txt"
    if _tok_file.is_file():
        ONCOKB_TOKEN = _tok_file.read_text(encoding="utf-8").strip()

# ---------------------------------------------------------------- 请求行为
REQUEST_TIMEOUT = 30.0          # 单个外部请求超时（秒）
MAX_RETRIES = 2                 # 网络错误重试次数
CACHE_TTL = 7 * 24 * 3600       # 成功结果缓存：7 天
CACHE_TTL_FAILURE = 3600        # 失败结果缓存：1 小时（避免一次网络抖动挂一周）
USER_AGENT = "MelodyLab-triage/0.1 (research use; https://github.com/Skyline-github/target-triage)"

# 五问定义（标题 + 问题原文，与 rules.yaml 的 key 对应）
FIVE_QUESTIONS = [
    {"key": "target", "title": "正确的靶点", "question": "它和疾病之间有没有扎实的因果证据？"},
    {"key": "tissue", "title": "正确的组织", "question": "药物能不能到达、并在正确的组织里起作用？"},
    {"key": "safety", "title": "正确的安全性", "question": "打这个靶点会不会出事？"},
    {"key": "patient", "title": "正确的患者", "question": "你能不能识别出“对的那群病人”？"},
    {"key": "commercial", "title": "正确的商业", "question": "值不值得做？"},
]

DATA_SOURCES = [
    {"name": "Open Targets Platform", "url": "https://platform.opentargets.org", "license": "开源可用"},
    {"name": "GWAS Catalog", "url": "https://www.ebi.ac.uk/gwas/", "license": "开源可用"},
    {"name": "Human Protein Atlas", "url": "https://www.proteinatlas.org", "license": "免费开放"},
    {"name": "UniProt", "url": "https://www.uniprot.org", "license": "开源可用"},
    {"name": "STRING", "url": "https://string-db.org", "license": "学术免费"},
    {"name": "ChEMBL", "url": "https://www.ebi.ac.uk/chembl/", "license": "开源免费"},
    {"name": "ClinicalTrials.gov", "url": "https://clinicaltrials.gov", "license": "免费"},
    {"name": "openFDA", "url": "https://open.fda.gov", "license": "免费（限流）"},
    {"name": "OncoKB（可选）", "url": "https://www.oncokb.org", "license": "需免费学术 token"},
]
