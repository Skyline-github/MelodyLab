#!/usr/bin/env python3
"""MelodyLab 资产准备脚本（在本机开发机上运行一次）。

从本地已有的 FlashLab 发行目录与 DockScope 解包目录复制：
  1. 对接引擎二进制  → assets/tools/          （Windows 优先；Linux/macOS 需自备同名文件）
  2. 模型 vendor     → assets/vendor/         （FlashDiff / DiffPROTACs / PLANET / Uni-Mol / p2rank）

用法：
  python scripts/setup_assets.py [--flashlab PATH] [--dockscope PATH] [--skip-vendor]

默认探测路径：
  ../FlashLab-Win            （相对本仓库）
  ../_analysis/dockscope_extract/resources
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ASSETS = REPO / "assets"


def copy_file(src: Path, dst: Path) -> bool:
    if not src.is_file():
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    print(f"  ✓ {dst.relative_to(ASSETS)}  ({src.stat().st_size/1e6:.1f} MB)")
    return True


def copy_tree(src: Path, dst: Path, exclude=("__pycache__", ".DS_Store", "*.pyc")) -> bool:
    if not src.is_dir():
        return False

    def ignore(dir_, names):
        bad = set()
        for n in names:
            if n in exclude or n.endswith(".pyc"):
                bad.add(n)
        return bad

    dst.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, dst, dirs_exist_ok=True, ignore=ignore)
    print(f"  ✓ {dst.relative_to(ASSETS)}/  ({sum(1 for _ in dst.rglob('*'))} 项)")
    return True


def setup_tools(flashlab: Path | None, dockscope: Path | None, with_gromacs: bool = False) -> None:
    print("\n[1/2] 对接引擎二进制 → assets/tools/")
    tools = ASSETS / "tools"
    tools.mkdir(parents=True, exist_ok=True)
    ds_tools = (dockscope / "tools") if dockscope else None

    got = False
    for src in [ds_tools / "vina.exe" if ds_tools else None,
                flashlab / "core" / "vina_1.2.7_win.exe" if flashlab else None]:
        if src and copy_file(src, tools / "vina_1.2.7_win.exe"):
            got = True
            break
    if not got:
        print("  ! 未找到 vina.exe（可从 https://github.com/ccsb-scripps/AutoDock-Vina/releases 下载 v1.2.7）")

    if ds_tools:
        if copy_file(ds_tools / "Vina-GPU-CUDA.exe", tools / "Vina-GPU-CUDA.exe"):
            copy_file(ds_tools / "Kernel1_Opt.bin", tools / "Kernel1_Opt.bin")
            copy_file(ds_tools / "Kernel2_Opt.bin", tools / "Kernel2_Opt.bin")
        else:
            print("  ! Vina-GPU-CUDA.exe 未找到（无 GPU 或不需要 GPU 对接可忽略）")
        if (ds_tools / "smina-runtime").is_dir():
            copy_tree(ds_tools / "smina-runtime", tools / "smina-runtime")
        else:
            print("  ! smina-runtime 未找到（可选）")
        if (ds_tools / "gromacs").is_dir():
            if with_gromacs:
                copy_tree(ds_tools / "gromacs", tools / "gromacs")
            else:
                print("  i GROMACS 属二期 MD 功能，未复制（运行 --with-gromacs 添加，约 500MB）")


def setup_vendor(flashlab: Path | None) -> None:
    print("\n[2/2] 模型 vendor → assets/vendor/")
    if not flashlab:
        print("  ! 未提供 --flashlab，跳过模型资产（生成/Linker/亲和力功能将不可用）")
        return
    vendor = ASSETS / "vendor"
    vendor.mkdir(parents=True, exist_ok=True)
    others = flashlab / "others"

    if not copy_tree(others / "FlashDiff_app", vendor / "FlashDiff_app") if others.is_dir() else False:
        print("  ! FlashDiff_app 未找到")
    if not copy_tree(others / "linker_design" / "DiffPROTACs", vendor / "DiffPROTACs"):
        print("  ! DiffPROTACs 未找到")
    if not copy_tree(others / "PLANET", vendor / "PLANET"):
        print("  ! PLANET 未找到")
    if not copy_tree(others / "Uni-Mol" / "unimol_docking_v2", vendor / "unimol_docking_v2"):
        print("  ! Uni-Mol docking 未找到")
    if not copy_tree(others / "p2rank_2.5", vendor / "p2rank_2.5"):
        print("  ! p2rank_2.5 未找到（可选，需 Java 17）")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--flashlab", type=Path, default=REPO.parent / "FlashLab-Win")
    ap.add_argument("--dockscope", type=Path,
                    default=REPO.parent / "_analysis" / "dockscope_extract" / "resources")
    ap.add_argument("--skip-vendor", action="store_true", help="只复制引擎二进制，跳过大体积模型")
    ap.add_argument("--with-gromacs", action="store_true", help="复制内置 GROMACS（二期 MD 功能，约 500MB）")
    args = ap.parse_args()

    flashlab = args.flashlab if args.flashlab and args.flashlab.is_dir() else None
    dockscope = args.dockscope if args.dockscope and args.dockscope.is_dir() else None
    print(f"MelodyLab 资产准备 → {ASSETS}")
    print(f"  FlashLab 来源: {flashlab or '未找到'}")
    print(f"  DockScope 来源: {dockscope or '未找到'}")

    setup_tools(flashlab, dockscope, with_gromacs=args.with_gromacs)
    if not args.skip_vendor:
        setup_vendor(flashlab)

    print("\n完成。重启后端后在首页查看环境状态。")
    print("提示：DiffPROTACs / FlashDiff / PLANET / Uni-Mol 均为开源研究许可，商用前请核对各项目 License。")


if __name__ == "__main__":
    sys.exit(main())
