"""二期 · 蛋白互作网络引擎骨架（STRING-DB）。

设计参考 DockScope ppi_service（已逆向）：
  https://version-12-0.string-db.org/api  (network / enrichment / interactors)
  离线兜底；图分析：Degree/Betweenness/Closeness/Eigenvector/MCC/MNC/DMNC/k-core/社区
  导出：GraphML / Cytoscape JSON / SVG。
"""
def availability() -> dict:
    return {"engine": "string_ppi", "available": True, "source": "string-db v12",
            "hint": "需网络；离线兜底模式实现时再评估"}
