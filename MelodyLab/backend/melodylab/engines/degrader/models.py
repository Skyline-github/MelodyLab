"""降解剂评估报告数据模型（pydantic）。"""
from __future__ import annotations

from pydantic import BaseModel, Field


class CvMetric(BaseModel):
    r2: float | None = None
    spearman: float | None = None


class Calibration(BaseModel):
    """骨架 5 折交叉验证的诚实技能指标 + conformal 区间参数。"""

    n_train: int
    n_folds: int
    q90_residual_log: float
    cv: dict[str, CvMetric] = Field(default_factory=dict)
    trained_at: str = ""
    db_snapshot: str = ""


class RankedRow(BaseModel):
    compound_id: str
    smiles: str
    e3: str = ""
    status: str = "ok"                       # ok | invalid_smiles
    source_line: int = 0
    pred_log_dc50: float | None = None       # log nM
    pred_dc50_nm: float | None = None
    lo90_nm: float | None = None             # 90% 区间（nM）
    hi90_nm: float | None = None
    rank_consensus: float | None = None      # 三族排名均值（越小越优）
    rank_mlp: float | None = None
    rank_et: float | None = None
    rank_hgb: float | None = None
    max_tanimoto: float | None = None        # 与训练集最近邻相似度（适用域）
    in_db: bool = False                      # PROTAC-DB 已收录该 SMILES
    db_dc50: str = ""                        # 库中实测 DC50（原始字符串）
    warnings: list[str] = Field(default_factory=list)


class DegraderReport(BaseModel):
    id: str
    created_at: str
    name: str = ""
    model_status: str = "cached"             # cached | trained
    n_query: int = 0
    n_invalid: int = 0
    n_duplicates: int = 0
    calibration: Calibration | None = None
    rows: list[RankedRow] = Field(default_factory=list)
    notes: list[str] = Field(default_factory=list)
