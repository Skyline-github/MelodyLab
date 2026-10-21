"""二期 · GEO 表达差异分析引擎骨架。

设计参考 DockScope geo_array_pipeline v6.2（已逆向，功能完整）：
  NCBI FTP/HTTPS 下载 series matrix → GPL 探针注释映射 → 分组设计（自动/手动）
  → limma / t-test → QC（PCA/相关性）→ 火山图/MA/热图 → HTML 报告 + Excel。
复刻时优先移植 limma 统计（Python 侧无 limma，用 statsmodels/scipy 等效实现或内嵌 R）。
"""
def availability() -> dict:
    import importlib.util

    ok = all(importlib.util.find_spec(m) for m in ("pandas", "numpy", "scipy"))
    return {"engine": "geo_diff", "available": ok,
            "hint": None if ok else "缺少科学计算依赖", "needs": ["statsmodels(评估)", "GEOparse(评估)"]}
