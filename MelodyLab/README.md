# MelodyLab

**开源 AI 药物发现工作台** —— 整合 FlashLab（AI 分子生成）与 DockScope（批量 GPU 对接与筛选流水线）功能的本地一体化平台。

> 形态：本地服务（FastAPI + Vue3），Windows 打开即用（打包方案见下）；无用户体系，面向研究与开源社区。

## 功能总览

| 模块 | 说明 | 状态 |
|---|---|---|
| 分子准入筛选 | RDKit 标准化、Ro5/QED/SA、PAINS/BRENK 警报、多样性选择 | ✅ 一期 |
| 批量分子对接 | 多蛋白×多配体、CB-Dock2 风格口袋检测、Vina-GPU/Vina/Smina 自动回退、断点续跑 | ✅ 一期 |
| 3D 分子生成 | FlashDiff 扩散模型口袋内从头生成（复用权重） | ✅ 一期（需资产） |
| Linker 生成 | DiffPROTACs PROTAC 连接子设计（复用权重） | ✅ 一期（需资产） |
| 亲和力预测 | PLANET pKd（挂接对接流水线） | 🟡 适配器就绪 |
| MD 模拟 | 内置 GROMACS（GPU）、OpenFF 柔性配体、NVT 预检门控、RMSD/RMSF/Rg 分析、3Dmol 轨迹动画 | ✅ 0.2.0 |
| 知识库检索 | ChEMBL/PubChem/RCSB/AlphaFold/ZINC/GEO | ✅ 0.2.0 |
| PPI 网络 | STRING v12 + 图中心性/社区发现 + Cytoscape 导出 | ✅ 0.2.0 |
| GEO 差异分析 | series matrix 解析 + t-test/BH-FDR + QC/火山图/热图 + HTML 报告 | ✅ 0.2.0 |
| 逆合成规划 | BRICS+RECAP 混合路线 + 递归分层 | ✅ 0.2.0 |

## 目录结构

```
MelodyLab/
├── backend/melodylab/        # FastAPI 后端
│   ├── core/jobs.py          #   任务中心（job.json 断点恢复）
│   ├── core/sysinfo.py       #   环境自检
│   ├── engines/              #   引擎插件层（screening / docking / generation / linker / affinity）
│   └── routers/              #   API（system / jobs / workloads）
├── frontend/                 # Vue3 + Vite 前端
├── assets/                   # 引擎二进制与模型权重（setup_assets.py 生成，git 忽略）
└── scripts/                  # setup_assets.py 资产收编 / doctor.py 环境自检
```

## 快速开始（开发模式）

```bash
# 1. 后端（Python 3.11+）
cd backend
pip install -r requirements.txt
python run.py --port 8300

# 2. 前端（开发热更）
cd frontend
npm install && npm run dev     # http://localhost:5173

# 3. 生产构建（后端直接托管前端）
cd frontend && npm run build   # 产物 dist/，重启后端后 http://127.0.0.1:8300 即是完整应用
```

## 资产准备（对接引擎与模型权重）

```bash
python scripts/setup_assets.py --flashlab ../FlashLab-Win --dockscope ../_analysis/dockscope_extract/resources
```

- `assets/tools/`：vina（必需）、smina / Vina-GPU-CUDA（可选）
  - 若无本地来源：从 [AutoDock Vina releases](https://github.com/ccsb-scripps/AutoDock-Vina/releases) 下载 `vina_1.2.7_win.exe` 放入 `assets/tools/`
- `assets/vendor/`：FlashDiff / DiffPROTACs / PLANET / Uni-Mol / p2rank（生成与亲和力功能需要）
- MD 模拟：GROMACS 2020.6 Windows 版已收编于 `assets/tools/gromacs/`；柔性配体参数化需要独立 conda 环境
  `mdff`（openff-toolkit 0.19 + torch≥2.6，与生产环境隔离）：`conda create -n mdff python=3.11 -y && conda run -n mdff pip install openff-toolkit torch rdkit`
- 未就绪的功能会在首页与对应页面给出明确提示，不影响其他模块使用。

## Windows 打开即用（打包方案）

1. 后端 PyInstaller onedir：`melodylab-backend.exe`（`build/` 提供打包脚本，产出单目录含 `_internal`）
2. 前端 `npm run build` → `frontend/dist` 由后端托管
3. 可选 Electron 薄壳：拉起后端（随机端口）+ 加载工作台（对齐 DockScope 体验，`build/electron/`）
4. NSIS 打包安装器（内嵌 `app-64.7z`）

详细步骤见 `build/README.md`（一期收尾时补全自动化脚本）。

## 设计要点（继承自逆向调研）

- **job.json 断点恢复**：任务目录内 `job.json` 是唯一恢复源，进程崩溃/重启后自动识别；
- **GPU 优先自动回退**：vina_gpu → vina_cpu → smina；
- **插件式引擎**：每个能力一个引擎模块（可用性探测 + 统一 worker 接口），新功能零侵入接入；
- **协作式暂停/停止**：工作单元之间检查点，长任务可暂停续跑；
- **资产与代码分离**：开源代码入库，大体积权重/二进制由脚本收编。

## 许可注意

复用的模型均为开源研究项目：FlashDiff、DiffPROTACs（**仅限研究用途**）、PLANET、Uni-Mol、p2rank、AutoDock Vina（Apache-2.0）、Smina（GPL-2.0）、Vina-GPU（GPL-3.0）。开源发布时请在发行说明中完整保留各自 LICENSE 与引用信息。
