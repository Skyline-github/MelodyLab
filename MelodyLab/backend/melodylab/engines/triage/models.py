"""数据模型：评估报告的证据与五问结构。

红/黄/绿只是数据驱动的筛查信号，不是结论——红色项列为"重点关注"，
最终判断权在评估者。
"""
from __future__ import annotations

from datetime import datetime
from typing import Any, Optional

from pydantic import BaseModel, Field, model_validator

# 检查点状态
GREEN, YELLOW, RED, UNKNOWN = "green", "yellow", "red", "unknown"

STATUS_LABEL = {GREEN: "证据充分", YELLOW: "需注意", RED: "风险信号", UNKNOWN: "信息不足"}

# 问题状态 → 数值（雷达图用）
STATUS_SCORE = {GREEN: 3, YELLOW: 2, RED: 1, UNKNOWN: 1.5}


class SourceLink(BaseModel):
    label: str
    url: str


class Evidence(BaseModel):
    """单个检查点的证据条目。chart 为可选的 Chart.js 配置。"""

    source: str                       # 数据库名，如 "Open Targets Platform"
    check: str                        # 检查点，如 "靶点-疾病关联强度"
    status: str = UNKNOWN             # green / yellow / red / unknown
    summary: str = ""                 # 一句话结论
    detail: Any = None                # 结构化证据（表格数据等）
    chart: Optional[dict] = None      # Chart.js 配置（None 表示无图）
    links: list[SourceLink] = Field(default_factory=list)
    manual: bool = False              # True 表示手动核查项（无自动化数据）


class QuestionReport(BaseModel):
    """一问的报告：多个证据条目 + 汇总。"""

    key: str
    title: str
    question: str
    evidences: list[Evidence] = Field(default_factory=list)
    manual_checklist: list[Evidence] = Field(default_factory=list)
    status: str = UNKNOWN

    def aggregate(self) -> None:
        """聚合提示级别：任一 red → red；否则任一 yellow → yellow；全 unknown → unknown。"""
        statuses = [e.status for e in self.evidences]
        if RED in statuses:
            self.status = RED
        elif YELLOW in statuses:
            self.status = YELLOW
        elif statuses and all(s == GREEN for s in statuses):
            self.status = GREEN
        else:
            self.status = UNKNOWN


class AssessmentReport(BaseModel):
    """一次完整评估。"""

    id: str
    gene_symbol: str
    gene_name: str = ""
    ensembl_id: str = ""
    uniprot_id: str = ""
    indication: str = ""              # 可选适应症
    created_at: str = Field(default_factory=lambda: datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    questions: list[QuestionReport] = Field(default_factory=list)
    overall_status: str = UNKNOWN
    risk_flags: list[str] = Field(default_factory=list)   # 红色项：需人工判断的风险信号
    watch_flags: list[str] = Field(default_factory=list)  # 黄色项：建议进一步核查

    @model_validator(mode="before")
    @classmethod
    def _accept_legacy_veto(cls, data: Any) -> Any:
        """兼容 v0.1 归档报告里的 veto_flags 字段。"""
        if isinstance(data, dict) and "veto_flags" in data and "risk_flags" not in data:
            data = dict(data)
            data["risk_flags"] = data.pop("veto_flags")
        return data

    def finalize(self) -> None:
        for q in self.questions:
            q.aggregate()
        self.risk_flags, self.watch_flags = [], []
        for q in self.questions:
            for e in q.evidences:
                if e.status == RED:
                    self.risk_flags.append(f"{q.title} · {e.check}：{e.summary}")
                elif e.status == YELLOW:
                    self.watch_flags.append(f"{q.title} · {e.check}：{e.summary}")
        statuses = [q.status for q in self.questions]
        if RED in statuses:
            self.overall_status = RED
        elif UNKNOWN in statuses:
            self.overall_status = UNKNOWN
        elif YELLOW in statuses:
            self.overall_status = YELLOW
        elif statuses:
            self.overall_status = GREEN

    @property
    def radar_data(self) -> dict:
        """五问雷达图数据：维度 = 五问，值 = 状态分（3 绿 / 2 黄 / 1.5 未知 / 1 红）。"""
        return {
            "labels": [q.title for q in self.questions],
            "values": [STATUS_SCORE.get(q.status, 1.5) for q in self.questions],
            "statuses": [q.status for q in self.questions],
        }
