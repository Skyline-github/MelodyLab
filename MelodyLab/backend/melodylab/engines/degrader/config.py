"""降解剂评估（DC50 排序）配置：数据文件、缓存与报告目录、模型配方。

方法论源自 Test_009_protac 研究项目（PROTAC-DB + 树模型/深度集成的多族共识排序，
conformal 不确定性区间 + 适用域检查）。定位是**排序参考工具**，不是绝对值预测器。
"""
from __future__ import annotations

import hashlib
from pathlib import Path

from ...config import settings as _app_settings

# ---------------------------------------------------------------- 目录
_DATA_DIR = Path(__file__).resolve().parent / "data"
DB_FILE = _DATA_DIR / "protacdb_snapshot_2024-10-06.csv"
EXAMPLE_FILE = _DATA_DIR / "example_galectin3_protacs.csv"

CACHE_DIR = _app_settings.data_root / "cache" / "degrader"
REPORTS_DIR = _app_settings.data_root / "degrader_reports"
for _d in (CACHE_DIR, REPORTS_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- 数据来源
DB_SNAPSHOT = "PROTAC-DB 2024-10-06 快照（protacdb.weizmann.ac.il，开源学术资源）"
DB_SOURCE_URL = "https://protacdb.weizmann.ac.il"

# ---------------------------------------------------------------- 特征
ECFP_BITS = 2048
ECFP_RADIUS = 4          # 与研究项目一致（deepchem CircularFingerprint radius=4, chiral）
TANIMOTO_RADIUS = 2      # 适用域检查用标准 ECFP4
DESC_NAMES = ["MW", "CrippenLogP", "TPSA", "HBA", "HBD",
              "RotB", "RingCount", "HeavyAtoms", "Fsp3", "AmideBonds"]

# ---------------------------------------------------------------- 训练配方
N_FOLDS = 5              # 骨架（Murcko）分组交叉验证
SEED = 42
MLP_HIDDEN = (256, 32)   # 与教程 MultitaskRegressor 相同的三层结构
MLP_DROPOUT = 0.2
MLP_BATCH = 10
MLP_FINAL_SEEDS = 8      # 最终集成成员数
MLP_FINAL_EPOCHS = 60    # 研究项目实测 20 epochs 后收益趋平（60ep 复核过），60 兼顾质量与时长
MLP_CV_SEEDS = 3         # CV 折内轻量配置（控制首次训练时长）
MLP_CV_EPOCHS = 40
TORCH_THREADS = 4        # 小批量矩阵乘，多线程仅小幅提速（默认全核反而有争用）
CONFORMAL_Q = 0.90       # 折外残差分位数 → 90% 区间的下限
Z90 = 1.645              # 高斯近似 90% 半宽
MAX_QUERY = 500          # 单次任务分子数上限
MIN_TRAIN_ROWS = 300     # 低于此值视为数据异常，拒绝训练

# ---------------------------------------------------------------- 缓存
MODEL_VERSION = "v2"     # 配方变更时递增，旧缓存自动失效


def db_signature() -> str:
    """训练数据文件指纹（内容 hash），用于缓存失效判断。"""
    h = hashlib.sha256()
    h.update(DB_FILE.read_bytes())
    return h.hexdigest()[:12]


def cache_path() -> Path:
    return CACHE_DIR / f"models_{MODEL_VERSION}_{db_signature()}.joblib"
