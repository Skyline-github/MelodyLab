# Windows「打开即用」打包方案（一期收尾）

目标产物：一个 Windows 安装包（NSIS），安装后双击 MelodyLab 即进入工作台，
对齐 DockScope 的桌面体验：Electron 薄壳 + PyInstaller 后端 + 静态前端。

## 分层产物

```
build/
├── backend.spec          # PyInstaller：melodylab-backend.exe（onedir）
├── electron/             # 薄壳（复刻 DockScope main.cjs 结构，自研实现）
│   ├── main.cjs          #   随机端口拉起后端 → /api/health 轮询 → 加载 #/home
│   ├── package.json      #   electron + electron-builder 配置
│   └── (无授权门 —— 开源版不需要)
└── installer.nsi         # NSIS 脚本：内嵌 app-64.7z
```

## 步骤

1. **前端构建**：`cd frontend && npm run build`（产物 frontend/dist）
2. **后端打包**（Windows 上执行）：
   ```powershell
   pip install pyinstaller
   pyinstaller build/backend.spec
   # 产出 build/dist/melodylab-backend/{melodylab-backend.exe, _internal/}
   ```
   spec 关键点：
   - `--collect-all rdkit`（数据文件：BaseFeatures/FragmentDescriptors 等）
   - 隐式导入：melodylab.engines.* 适配器（延迟导入不被静态发现时需显式声明）
   - 不打包 assets/（体积原因，安装器单独携带）
3. **目录组装**（.resources 布局，与 DockScope 一致）：
   ```
   resources/
   ├── backend/{melodylab-backend.exe, _internal/}
   ├── frontend/dist/
   ├── assets/tools/          # vina / smina / Vina-GPU（安装时可选组件）
   └── assets/vendor/         # 模型权重（可选下载器或单独包）
   ```
4. **Electron 壳**：`cd build/electron && npx electron-builder --win nsis`
   - `extraResources` 映射上述 resources/
   - main.cjs 要点：`findPort()` → `spawn(backend, ['--port', port], {env: 注入
     MELODYLAB_DATA_ROOT / MELODYLAB_ASSETS / MELODYLAB_FRONTEND_DIST})` →
     轮询 `/api/health` → `loadURL('http://127.0.0.1:port/#/home')`
   - 数据根：`C:\MelodyLabData`（短路径规避 MAX_PATH）
5. **安装包**：electron-builder 直接产出 NSIS 安装器；或手写 installer.nsi 嵌 `app-64.7z`。

## 验收清单（Win10 全新机）

- [ ] 双击安装 → 桌面图标 → 启动 ≤ 30s 出现工作台
- [ ] 首页环境自检显示 vina 引擎可用
- [ ] 筛选：粘贴 SMILES → 出结果 CSV 下载
- [ ] 对接：1 蛋白 × 5 SMILES → 出结果与 TOP 包
- [ ] 生成/Linker（带资产包时）：任务完成产出 SDF
