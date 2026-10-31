"""开发启动脚本：python run.py --port 8300"""
import os
import sys
from pathlib import Path

# 模型权重缓存重定向到工作区（DECIMER/easyocr 等 pystow 系默认落 C 盘 ~/.data）
_MODEL_CACHE = Path(__file__).resolve().parent.parent / ".model_cache"
_MODEL_CACHE.mkdir(exist_ok=True)
os.environ.setdefault("PYSTOW_HOME", str(_MODEL_CACHE))

# 供 `python -m melodylab` 与直接执行两种方式共用
sys.path.insert(0, str(Path(__file__).resolve().parent))

from melodylab.main import main  # noqa: E402

if __name__ == "__main__":
    main()
