"""服务包：显式导入各数据源子模块，便于 `services.xxx` 直接访问。"""
from . import (chembl, clinicaltrials, gwas, oncokb, openfda, opentargets,
               protein_atlas, stringdb, uniprot)

__all__ = ["chembl", "clinicaltrials", "gwas", "oncokb", "openfda",
           "opentargets", "protein_atlas", "stringdb", "uniprot"]
