"""MelodyLab 后端入口：FastAPI 应用工厂。

本地服务形态：uvicorn 启动后由前端（或 Electron 壳）访问 http://127.0.0.1:{port}。
若存在 frontend/dist 则直接托管前端静态资源。
"""
from __future__ import annotations

import argparse
from pathlib import Path

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from .config import frontend_dist, settings
from .routers import jobs, system, workloads
from .routers import md as md_router
from .routers import phase2
from .routers import triage as triage_router


def create_app() -> FastAPI:
    app = FastAPI(title="MelodyLab", version="0.2.0",
                  description="开源 AI 药物发现工作台（FlashLab + DockScope 功能整合）")
    app.add_middleware(
        CORSMiddleware,
        allow_origins=settings.cors_origins or ["*"],
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )
    app.include_router(system.router)
    app.include_router(jobs.router)
    app.include_router(workloads.router)
    app.include_router(md_router.router)
    app.include_router(phase2.router)
    app.include_router(triage_router.router)

    dist = frontend_dist()
    if dist:
        app.mount("/assets", StaticFiles(directory=dist / "assets"), name="assets")

        @app.get("/{full_path:path}", include_in_schema=False)
        def spa(full_path: str):
            candidate = dist / full_path
            if full_path and candidate.is_file():
                return FileResponse(candidate)
            # index.html 禁止缓存：前端重新构建（资产哈希变化）后浏览器必须立即拿到新入口，
            # 否则用户会一直用旧 bundle，"改了没生效"
            return FileResponse(dist / "index.html", headers={"Cache-Control": "no-cache"})
    return app


app = create_app()


def main() -> None:
    parser = argparse.ArgumentParser(description="MelodyLab backend")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, required=True)
    parser.add_argument("--log-level", default="info")
    args = parser.parse_args()
    import uvicorn

    print("=" * 60)
    print(f"MelodyLab Server Starting...  data root: {settings.data_root}")
    print(f"Serving on http://{args.host}:{args.port}")
    print("Press Ctrl+C to stop")
    print("=" * 60)
    uvicorn.run(app, host=args.host, port=args.port, log_level=args.log_level, access_log=False)


if __name__ == "__main__":
    main()
