"""三元复合物预测报告数据模型（pydantic）。"""
from __future__ import annotations

from pydantic import BaseModel, Field


class PrepInfo(BaseModel):
    """输入制备指标（装配质量，决定下游可信度）。"""

    mcs_poi: int | None = None          # POI warhead MCS 原子数
    mcs_e3: int | None = None
    poi_rmsd: float | None = None       # 对齐到晶体姿态的 RMSD（Å）
    e3_rmsd: float | None = None
    lig_kind: str = ""                  # constrained / grafted / relinked
    best_clear: float | None = None     # 装配时蛋白间最小原子距（Å）
    prep_warning: str = ""


class TernaryRow(BaseModel):
    compound_id: str
    smiles: str
    e3: str = ""
    status: str = "ok"                  # ok | prep_failed | infer_failed
    reason: str = ""                    # 失败原因（status != ok 时）
    prep: PrepInfo = Field(default_factory=PrepInfo)
    pred_p2_rmsd: float | None = None   # DeepTernary 姿态质量代理指标（越小越好）
    clash: float | None = None          # 重叠原子占比
    bsa: float | None = None            # 蛋白-蛋白埋藏表面积（Å²，只有相对意义）
    complex_file: str = ""              # 最优 seed 复合物 PDB（相对报告目录）


class TernaryReport(BaseModel):
    id: str
    created_at: str
    name: str = ""
    seeds: int = 2
    n_query: int = 0
    n_ok: int = 0
    rows: list[TernaryRow] = Field(default_factory=list)
    notes: list[str] = Field(default_factory=list)
