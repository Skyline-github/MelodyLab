"""三元复合物预测（DeepTernary）配置：vendor 路径、解释器解析、参照与报告目录。

协议源自 Test_009_protac 研究项目（galectin-3 POI + CRBN/VHL 参照晶体装配 →
DeepTernary 推理 → BSA/clash 结构健康评估）。DeepTernary 依赖 dgl，与主环境
torch 版本不兼容，一律以**子进程**调用独立 conda 环境执行推理。
"""
from __future__ import annotations

import os
from pathlib import Path

from ...config import settings as _app_settings

# ---------------------------------------------------------------- 目录
_DATA_DIR = Path(__file__).resolve().parent / "data"
REFERENCES_DIR = _DATA_DIR / "references"

CACHE_DIR = _app_settings.data_root / "cache" / "ternary"
REPORTS_DIR = _app_settings.data_root / "ternary_reports"
for _d in (CACHE_DIR, REPORTS_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- vendor
VENDOR_DIR = _app_settings.assets / "vendor" / "DeepTernary"
CHECKPOINT_REL = "output/checkpoints/PROTAC/checkpoint.pth"
CONFIG_REL = "deepternary/configs/protac.py"


def checkpoint_path() -> Path | None:
    p = VENDOR_DIR / CHECKPOINT_REL
    return p if p.is_file() else None


def interpreter_path() -> Path | None:
    """DeepTernary 子进程 python 解析：env 覆盖 > 数据根标记文件 > 常见 conda env 目录。"""
    env = os.environ.get("MELODYLAB_DEEPTERNARY_PYTHON", "").strip()
    if env and Path(env).is_file():
        return Path(env)
    marker = _app_settings.data_root / "deepternary_python.txt"
    if marker.is_file():
        p = Path(marker.read_text(encoding="utf-8").strip())
        if str(p) and p.is_file():
            return p
    candidates = []
    for root_env in (os.environ.get("CONDA_PREFIX"), os.environ.get("CONDA_EXE")):
        if root_env:
            candidates.append(Path(root_env))
    for base in (Path("D:/Anaconda"), Path("C:/Anaconda"),
                 Path("C:/ProgramData/Anaconda3"), Path.home() / "Anaconda3",
                 Path.home() / "miniconda3", Path.home() / "anaconda3"):
        candidates.append(base / "envs" / "deepternary" / "python.exe")
        candidates.append(base / "envs" / "deepternary" / "bin" / "python")
    for c in candidates:
        if c.is_file():
            return c
    return None


# ---------------------------------------------------------------- 参照晶格
# 与研究项目 prepare_inputs.py 完全一致：
# POI 参照：galectin-3C + TD2（5E89 链 A）；E3 参照：CRBN/LVY（4CI2 链 B）、
# VHL/3JF（4W9H 链 C）。min_anchor：CRBN 口袋内凹，锚点距离需外推到 >=34 Å。
POI_REF = {"resname": "TD2", "pdb": "5E89.pdb", "lig_chain": None,
           "protein_chains": ["A"], "ideal": "TD2_ideal.sdf"}
E3_REF = {
    "CRBN": {"resname": "LVY", "pdb": "4CI2.pdb", "lig_chain": "B",
             "protein_chains": ["B"], "ideal": "LVY_ideal.sdf", "min_anchor": 34.0},
    "VHL": {"resname": "3JF", "pdb": "4W9H.pdb", "lig_chain": "C",
            "protein_chains": ["C"], "ideal": "3JF_ideal.sdf", "min_anchor": 22.0},
}

SEED = 42                # 构象/采样随机种子（与研究项目一致）
DEFAULT_SEEDS = 2        # DeepTernary 推理 seed 数（研究项目用 4；2 兼顾速度）
MAX_QUERY = 100          # 单次任务分子数上限（制备+推理较重）
MIN_POI_MCS = 8          # MCS 原子数下限（低于则该分子跳过）
COMPLEX_CHAINS = ("A", "B")   # 预测复合物中 POI/E3 的链名
DB_NOTE = ("参照晶体：galectin-3C+TD2 (5E89)、CRBN+lenalidomide (4CI2)、"
           "VHL+VH032 (4W9H)；协议移植自 Test_009_protac 研究项目")
