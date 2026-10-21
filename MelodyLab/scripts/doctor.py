#!/usr/bin/env python3
"""MelodyLab 环境检索器（命令行入口）。

只依赖 Python 标准库 —— 在还没安装任何第三方包的电脑上也能直接运行。

用法（在 MelodyLab 目录）：
    python scripts/doctor.py                 # 扫描全部 conda/系统环境，给出推荐
    python scripts/doctor.py --install       # 向推荐环境自动补装缺失依赖
    python scripts/doctor.py --env NAME      # 指定环境名（如 base / melodylab）
    python scripts/doctor.py --json          # 机器可读输出
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from melodylab.core import envscan  # noqa: E402  (纯标准库模块，安全)


def main() -> int:
    ap = argparse.ArgumentParser(description="MelodyLab local environment doctor")
    ap.add_argument("--env", help="只看指定环境（名称或 python 路径）")
    ap.add_argument("--install", action="store_true", help="向选定/推荐环境补装缺失依赖")
    ap.add_argument("--json", action="store_true", help="JSON 输出")
    args = ap.parse_args()

    result = envscan.scan_all()
    envs = result["envs"]
    if args.env:
        envs = [e for e in envs if e["name"] == args.env or e["python"] == args.env]
        if not envs:
            print(f"[!] 未找到环境: {args.env}")
            return 1

    if args.json:
        import json

        print(json.dumps({"envs": envs, "recommended": result["recommended"]},
                         ensure_ascii=False, indent=2))
        return 0

    print("=" * 74)
    print(" MelodyLab 本地环境检索")
    print(f" 需要: Python >= {envscan.MIN_PY[0]}.{envscan.MIN_PY[1]} + {len(envscan.REQUIRED)} 个依赖包")
    print("=" * 74)
    for e in envs:
        mark = "✓" if e["ready"] else ("~" if e["python_ok"] else "✗")
        print(f"\n[{mark}] {e['name']}  (Python {e['python_version']})")
        print(f"     {e['python']}")
        if not e["python_ok"]:
            print("     [!] Python 版本过低，跳过")
            continue
        have = [p for p, info in e["packages"].items() if info["ok"]]
        miss = e["missing"]
        print(f"     已装 {len(have)}/{len(envscan.REQUIRED)}: {', '.join(have) if have else '无'}")
        if miss:
            print(f"     缺失: {', '.join(miss)}")
        if e["ready"]:
            print("     [OK] 可直接运行 MelodyLab")

    ready_names = [e["name"] for e in envs if e["ready"]]
    if ready_names:
        print("\n>>> 推荐直接使用环境: " + ", ".join(ready_names))
        print(">>> 启动方式（Anaconda Prompt）:")
        for n in ready_names:
            print(f"    conda activate {n}")
        print("    cd /d <MelodyLab>\\backend  &&  python run.py --port 8300")
        return 0

    target = envs[0] if envs and envs[0]["python_ok"] else None
    if not target:
        print("\n[!] 没有 Python >= 3.10 的环境。请先安装: https://www.python.org/downloads/")
        return 1

    print(f"\n>>> 推荐: 向环境 [{target['name']}] 补装缺失包:")
    print("    " + target["install_cmd"])
    if args.install:
        print("\n正在安装 ...")
        r = envscan.install_missing(target["python"])
        print(("\n[OK] 完成，已安装: " + ", ".join(r.get("installed", []))) if r.get("ok")
              else ("\n[FAIL] " + r.get("error", "")))
        return 0 if r.get("ok") else 1
    print("\n(自动补装: python scripts/doctor.py --install)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
