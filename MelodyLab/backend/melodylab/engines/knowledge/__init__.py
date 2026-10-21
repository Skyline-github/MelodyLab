"""二期 · 数据库知识检索引擎骨架。

设计参考 DockScope knowledge_service（已逆向），数据源与端点：
  ChEMBL   https://www.ebi.ac.uk/chembl/api/data/{target,molecule,activity}
  PubChem  https://pubchem.ncbi.nlm.nih.gov/rest/pug/{assay,compound}
  UniProt  https://rest.uniprot.org/uniprotkb/search
  RCSB     https://search.rcsb.org/rcsbsearch/v2/query + data.rcsb.org 元数据
  AlphaFold https://alphafold.ebi.ac.uk/files/AF-<uniprot>-F1-model_v4.pdb
  ZINC-22  https://cartblanche22.docking.org/{smiles.txt,search/result}
  GEO      https://eutils.ncbi.nlm.nih.gov/entrez/eutils (esearch → 下载 series matrix)
全部官方公开 API；需要网络。结果缓存 + 结构 ZIP 打包。
"""
from ...config import settings


def availability() -> dict:
    import importlib.util

    has_http = importlib.util.find_spec("httpx") is not None or importlib.util.find_spec("aiohttp") is not None
    return {
        "engine": "knowledge",
        "available": has_http,  # 纯 API 检索，无本地资产
        "sources": ["chembl", "pubchem", "uniprot", "rcsb", "alphafold", "zinc22", "geo"],
        "hint": None if has_http else "二期实现时补充 httpx 依赖",
    }
