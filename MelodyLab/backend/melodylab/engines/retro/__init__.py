"""二期 · 逆合成路线规划引擎骨架。

设计参考 DockScope workflow_service.synthesis_plan（已逆向）：
  RDKit BRICS + RECAP 混合断键 → 逐层前体构建 → 置信度启发式评分
  → 交互式路线图（前端分层布局）→ 导出 SVG/PNG/PDF。
纯本地计算，无外部模型依赖。
"""
def availability() -> dict:
    try:
        from rdkit import Chem  # noqa: F401
        from rdkit.Chem import BRICS  # noqa: F401

        return {"engine": "retrosynthesis", "available": True,
                "engines": ["brics", "recap", "hybrid"]}
    except Exception:
        return {"engine": "retrosynthesis", "available": False, "hint": "需要 rdkit"}
