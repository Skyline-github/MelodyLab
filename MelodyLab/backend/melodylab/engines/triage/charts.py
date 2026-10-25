"""Chart.js 配置生成器：引擎侧把证据转成图表配置，模板只负责渲染。"""
from __future__ import annotations

PALETTE = ["#2e6b5e", "#6aa38f", "#a3c9bd", "#d9b25f", "#c97b5e",
           "#7e8fa6", "#b5a264", "#8f6f8e", "#5f8fa8", "#9aa37a"]
LEVEL_COLORS = {"High": "#c0504d", "Medium": "#d9b25f", "Low": "#a3c9bd",
                "Not detected": "#d8dedb"}


def hbar(labels: list[str], dataset_label: str, values: list[float],
         second_label: str | None = None, second_values: list[float] | None = None,
         colors: list[str] | None = None, max_value: float | None = None) -> dict:
    """水平条形图（可叠第二组数据，如 遗传学分数）。"""
    datasets: list[dict] = [{
        "label": dataset_label, "data": values,
        "backgroundColor": colors or PALETTE[0], "borderRadius": 4,
    }]
    if second_label and second_values is not None:
        datasets.append({
            "label": second_label, "data": second_values,
            "backgroundColor": PALETTE[3], "borderRadius": 4,
        })
    return {
        "type": "bar",
        "data": {"labels": labels, "datasets": datasets},
        "options": {
            "indexAxis": "y",
            "responsive": True, "maintainAspectRatio": False,
            "plugins": {"legend": {"display": second_label is not None}},
            "scales": {
                "x": {"beginAtZero": True,
                      "max": max_value,
                      "grid": {"color": "#eef2f0"}},
                "y": {"grid": {"display": False}},
            },
        },
    }


def doughnut(labels: list[str], values: list[float],
             colors: list[str] | None = None) -> dict:
    return {
        "type": "doughnut",
        "data": {"labels": labels,
                 "datasets": [{"data": values,
                               "backgroundColor": colors or PALETTE}]},
        "options": {
            "responsive": True, "maintainAspectRatio": False,
            "plugins": {"legend": {"position": "right"}},
            "cutout": "55%",
        },
    }
