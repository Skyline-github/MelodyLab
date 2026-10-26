# MelodyLab

**开源 AI 药物发现工作台** —— 一站式覆盖小分子药物研发全流程：**靶点立项评估（选靶）→ 分子筛选 → 批量对接 → 3D 分子生成 → Linker 设计 → 降解剂评估 → 三元复合物预测 → 分子动力学验证**，并内置知识库检索、PPI 网络分析、GEO 差异分析与逆合成规划。

- **本地运行，数据不出机器**：全部计算在工作站本地完成，适合保密性课题；
- **GPU 加速**：Vina-GPU 批量对接 + GROMACS CUDA 动力学模拟；
- **断点续跑**：每个任务落盘 `job.json`，进程崩溃或重启后从断点继续，长任务可暂停/恢复；
- **插件式引擎**：各能力以「可用性探测 + 统一 worker 接口」接入，互不干扰、按需启用。

> 形态：本地服务（FastAPI + Vue3），Windows 打开即用；无用户体系，面向研究与开源社区。当前版本 **v0.2.0**。

## 功能总览

| 模块 | 说明 | 状态 |
|---|---|---|
| 靶点立项评估 | 立项五问（AstraZeneca 5R）：Open Targets / GWAS Catalog / Human Protein Atlas / UniProt / ChEMBL / ClinicalTrials.gov / openFDA 等 9 大公开数据库自动取证据，红/黄/绿风险信号 + 五问雷达图，报告导出 JSON/Markdown/Excel，阈值可在 `rules.yaml` 自调 | ✅ 新 |
| 分子准入筛选 | RDKit 标准化、Ro5/QED/SA、PAINS/BRENK 警报、多样性选择 | ✅ |
| 批量分子对接 | 多蛋白×多配体、口袋检测、Vina-GPU/Vina/Smina 自动回退、断点续跑 | ✅ |
| 3D 分子生成 | FlashDiff 扩散模型口袋内从头生成（复用权重） | ✅ 需资产 |
| Linker 生成 | DiffPROTACs PROTAC 连接子设计（复用权重） | ✅ 需资产 |
| 降解剂评估 | PROTAC-DB 训练的 MLP 集成 / ExtraTrees / HistGB 三族共识 DC50 排序 + conformal 90% 区间 + 适用域检查（排序参考工具，非绝对值预测） | ✅ 新 |
| 三元复合物预测 | DeepTernary：参照晶体装配（Galectin-3C × CRBN/VHL）→ 三元复合物预测 → BSA/clash 结构健康，3Dmol 交互查看 | ✅ 需资产 |
| 亲和力预测 | PLANET pKd（挂接对接流水线） | 🟡 适配器就绪 |
| MD 模拟 | 内置 GROMACS（GPU）、OpenFF 柔性配体、NVT 预检门控、RMSD/RMSF/Rg 分析、3Dmol 轨迹动画 | ✅ 0.2.0 |
| 三元复合物 MD | 三元复合物预测结果一键动力学验证：界面姿态配体、几何解穿插、三体接触曲线（LIG–POI / LIG–E3 / POI–E3） | ✅ 新 |
| 知识库检索 | ChEMBL/PubChem/RCSB/AlphaFold/ZINC/GEO | ✅ 0.2.0 |
| PPI 网络 | STRING v12 + 图中心性/社区发现 + Cytoscape 导出 | ✅ 0.2.0 |
| GEO 差异分析 | series matrix 解析 + t-test/BH-FDR + QC/火山图/热图 + HTML 报告 | ✅ 0.2.0 |
| 逆合成规划 | BRICS+RECAP 混合路线 + 递归分层 | ✅ 0.2.0 |

## 目录结构

```
.
├── MelodyLab/                # 项目主体
│   ├── backend/melodylab/    #   FastAPI 后端
│   │   ├── core/jobs.py      #     任务中心（job.json 断点恢复）
│   │   ├── engines/          #     引擎插件层（triage/screening/docking/generation/linker/affinity/md/knowledge/ppi/geo/retro）
│   │   └── routers/          #     API 路由
│   ├── frontend/             #   Vue3 + Vite 前端（npm run build 后 dist 由后端托管）
│   ├── assets/               #   引擎二进制与模型权重（setup_assets.py 本地收编，不入库）
│   ├── scripts/              #   setup_assets.py 资产收编 / doctor.py 环境自检
│   └── build/                #   PyInstaller/Electron 打包方案
├── example_data/             # 示例数据：EGFR/Mpro 蛋白 + 激酶抑制剂分子，开箱试玩
└── Start MelodyLab.bat       # Windows 双击启动脚本（需按本机 Python 路径微调）
```

## 快速开始

### 方式一：Windows 双击启动

双击 `Start MelodyLab.bat`（首次使用请用编辑器打开，把其中的 Python 路径改成本机 Anaconda 环境路径）。

### 方式二：命令行（开发模式）

```bash
# 1. 后端（Python 3.11+）
cd MelodyLab/backend
pip install -r requirements.txt
python run.py --port 8300

# 2. 前端（开发热更，可选）
cd MelodyLab/frontend
npm install && npm run dev     # http://localhost:5173

# 3. 生产构建（后端直接托管前端 dist）
cd MelodyLab/frontend && npm run build
```

启动后访问 **http://127.0.0.1:8300** 即是完整应用。在「靶点立项评估」输入基因符号（如 EGFR）即可做立项五问筛查；上传 `example_data/` 里的蛋白与分子即可体验 筛选 → 对接 → MD 全流程。

## 资产准备（对接引擎与模型权重）

大体积第三方资产不入库，从本地发行目录收编：

```bash
cd MelodyLab
python scripts/setup_assets.py --flashlab ../FlashLab-Win --dockscope ../_analysis/dockscope_extract/resources
```

- `assets/tools/`：vina（必需）、smina / Vina-GPU-CUDA（可选）
  - 若无本地来源：从 [AutoDock Vina releases](https://github.com/ccsb-scripps/AutoDock-Vina/releases) 下载 `vina_1.2.7_win.exe` 放入 `assets/tools/`
- `assets/vendor/`：FlashDiff / DiffPROTACs / PLANET / Uni-Mol / p2rank（生成与亲和力功能需要）
  - 三元复合物预测：`assets/vendor/DeepTernary/`（代码 + `output/checkpoints/PROTAC/` 权重，[DeepTernary v1.0.1](https://github.com/youqingxiaozhua/DeepTernary/releases)，仅限研究用途），并需独立 conda 环境
    `deepternary`（dgl 与主环境 torch 不兼容，必须隔离）：`conda create -n deepternary python=3.10 -y`，
    `pip install torch==2.2.1 "dgl==2.2.1" "numpy==1.26.4" "torchdata==0.7.1" "mmengine==0.10.3" "mmcv-lite==2.2.0" "rdkit==2023.9.3" biopandas "biopython==1.84" POT importlib_metadata prettytable`
    （注意 torchdata 必须 0.7.1、numpy 必须 1.x；环境路径写入 `数据根/deepternary_python.txt` 或设 `MELODYLAB_DEEPTERNARY_PYTHON`）
- MD 模拟：GROMACS 2020.6 Windows 版收编于 `assets/tools/gromacs/`；柔性配体参数化需要独立 conda 环境
  `mdff`（openff-toolkit 0.19 + torch≥2.6，与生产环境隔离）：`conda create -n mdff python=3.11 -y && conda run -n mdff pip install openff-toolkit torch rdkit`
- 未就绪的功能会在首页与对应页面给出明确提示，不影响其他模块使用。

## Windows 打包分发（可选）

1. 后端 PyInstaller onedir：`melodylab-backend.exe`（`build/` 提供打包脚本）
2. 前端 `npm run build` → `frontend/dist` 由后端托管
3. 可选 Electron 薄壳 + NSIS 安装器，详见 `MelodyLab/build/README.md`

## 设计要点

- **job.json 断点恢复**：任务目录内 `job.json` 是唯一恢复源，进程崩溃/重启后自动识别；
- **GPU 优先自动回退**：vina_gpu → vina_cpu → smina；
- **插件式引擎**：每个能力一个引擎模块（可用性探测 + 统一 worker 接口），新功能零侵入接入；
- **协作式暂停/停止**：工作单元之间检查点，长任务可暂停续跑；
- **资产与代码分离**：代码入库，大体积权重/二进制由脚本本地收编。

## 许可证

本项目代码以 [MIT](LICENSE) 许可证发布。

## 第三方许可注意

`assets/` 中收编的模型与工具（不入库）均为开源研究项目：FlashDiff、DiffPROTACs（**仅限研究用途**）、PLANET、Uni-Mol、p2rank、AutoDock Vina（Apache-2.0）、Smina（GPL-2.0）、Vina-GPU（GPL-3.0）。分发时请完整保留各自 LICENSE 与引用信息。
