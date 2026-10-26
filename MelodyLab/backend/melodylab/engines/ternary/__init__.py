"""三元复合物预测引擎（DeepTernary 子进程封装）。

参照晶体装配（主环境 rdkit）→ DeepTernary 推理（deepternary 独立环境子进程，
dgl 与主环境 torch 不兼容）→ BSA/clash 结构健康评估（主环境 biopython）。
POI 参照内置 galectin-3C；E3 支持 CRBN/VHL。
"""
from __future__ import annotations

from ...config import settings
from . import config as cfg


def availability() -> dict:
    vendor_ok = cfg.VENDOR_DIR.is_dir()
    ckpt = cfg.checkpoint_path()
    refs_ok = all((cfg.REFERENCES_DIR / f).is_file() for f in
                  ("5E89.pdb", "4CI2.pdb", "4W9H.pdb",
                   "TD2_ideal.sdf", "LVY_ideal.sdf", "3JF_ideal.sdf"))
    import importlib.util
    libs = {name: importlib.util.find_spec(name) is not None
            for name in ("rdkit", "scipy", "Bio", "pandas")}
    interp = cfg.interpreter_path()
    missing = [k for k, ok in libs.items() if not ok]
    available = vendor_ok and bool(ckpt) and refs_ok and not missing and bool(interp)
    hint = None
    if missing:
        hint = f"主环境缺少依赖：{', '.join(missing)}"
    elif not vendor_ok or not ckpt:
        hint = ("DeepTernary 资产未就绪：请运行 scripts/setup_assets.py 放置 "
                "assets/vendor/DeepTernary（代码 + output/checkpoints/PROTAC 权重，仅限研究用途）")
    elif not refs_ok:
        hint = f"参照晶体缺失：{cfg.REFERENCES_DIR}"
    elif not interp:
        hint = ("未找到 deepternary 独立环境 python：创建 conda env（python=3.10, "
                "torch==2.2.1, dgl==2.2.1, torchdata==0.7.1），并把其 python 路径写入 "
                f"{settings.data_root / 'deepternary_python.txt'} 或设置环境变量 "
                "MELODYLAB_DEEPTERNARY_PYTHON")
    return {
        "engine": "ternary",
        "available": available,
        "hint": hint,
        "vendor_dir": str(cfg.VENDOR_DIR) if vendor_ok else "",
        "checkpoint": bool(ckpt),
        "references_ready": refs_ok,
        "interpreter": str(interp) if interp else "",
        "poi_reference": "galectin-3C（5E89 + TD2）",
        "e3_supported": list(cfg.E3_REF),
        "db_note": cfg.DB_NOTE,
    }
