"""统一任务中心：job.json 文件即恢复源（借鉴 DockScope 的可靠性设计）。

每个任务一个目录：
    <jobs_dir>/<job_id>/
        job.json     # 状态记录（原子写入）
        logs/run.log
        inputs/      # 用户上传
        work/        # 引擎工作区（含断点文件）
        outputs/     # 结果与打包产物

状态机：queued → running ⇄ paused → stopping → stopped
                       ↘ completed / failed / cancelled
暂停/停止都是协作式：worker 在每个工作单元之间调用 control.checkpoint()。
"""
from __future__ import annotations

import json
import os
import shutil
import threading
import time
import traceback
import uuid
from pathlib import Path
from typing import Any, Callable

from ..config import settings

WORKER_FN = Callable[[dict, "JobControl"], dict]


class JobPaused(Exception):
    pass


class JobInterrupted(Exception):
    pass


def _now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


class JobControl:
    """交给 worker 的控制柄：进度上报 + 暂停/停止检查点。"""

    def __init__(self, manager: "JobManager", job_id: str) -> None:
        self.manager = manager
        self.job_id = job_id
        self._stop = threading.Event()
        self._pause = threading.Event()
        self.cancel_event = self._stop  # 引擎子进程取消回调别名
        self.heartbeat = time.time()

    # ---- 引擎轮询的标志 ----
    @property
    def stop_requested(self) -> bool:
        return self._stop.is_set()

    @property
    def pause_requested(self) -> bool:
        return self._pause.is_set()

    def request_stop(self) -> None:
        self._stop.set()

    def request_pause(self) -> None:
        self._pause.set()

    # ---- worker 在工作单元之间调用 ----
    def checkpoint(
        self,
        progress: int | None = None,
        phase: str | None = None,
        message: str | None = None,
        **stats: Any,
    ) -> None:
        self.heartbeat = time.time()
        self.manager.update(self.job_id, progress=progress, phase=phase, message=message, stats=stats)
        if self._stop.is_set():
            raise JobInterrupted()
        if self._pause.is_set():
            raise JobPaused()

    def subprocess_callback(self, proc):  # 与引擎的 process_callback 约定对齐
        if self._stop.is_set() or self._pause.is_set():
            try:
                proc.terminate()
            except Exception:
                pass


class JobManager:
    def __init__(self, jobs_dir: Path | None = None) -> None:
        self.jobs_dir = Path(jobs_dir) if jobs_dir else settings.jobs_dir
        self.jobs_dir.mkdir(parents=True, exist_ok=True)
        self._threads: dict[str, threading.Thread] = {}
        self._controls: dict[str, JobControl] = {}
        self._wake = threading.Event()  # 暂停挂起时等待恢复信号

    # ---------- 路径 ----------
    def job_dir(self, job_id: str) -> Path:
        return self.jobs_dir / job_id

    def _job_file(self, job_id: str) -> Path:
        return self.job_dir(job_id) / "job.json"

    # ---------- 基础读写 ----------
    def _read(self, job_id: str) -> dict | None:
        f = self._job_file(job_id)
        if not f.is_file():
            return None
        try:
            return json.loads(f.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            return None

    def _write(self, record: dict) -> None:
        record["updated_at"] = _now()
        f = self._job_file(record["id"])
        tmp = f.with_name(f.name + f".{os.getpid()}.tmp")
        tmp.write_text(json.dumps(record, ensure_ascii=False, indent=2), encoding="utf-8")
        tmp.replace(f)

    # ---------- 生命周期 ----------
    def create(self, kind: str, name: str, params: dict) -> dict:
        job_id = f"{kind}_{time.strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:6]}"
        d = self.job_dir(job_id)
        for sub in ("logs", "inputs", "work", "outputs"):
            (d / sub).mkdir(parents=True, exist_ok=True)
        record = {
            "id": job_id,
            "kind": kind,
            "name": name or job_id,
            "status": "queued",
            "progress": 0,
            "phase": "",
            "message": "",
            "params": params,
            "stats": {},
            "artifacts": [],
            "error": None,
            "created_at": _now(),
            "started_at": None,
            "finished_at": None,
        }
        self._write(record)
        return record

    def update(
        self,
        job_id: str,
        *,
        status: str | None = None,
        progress: int | None = None,
        phase: str | None = None,
        message: str | None = None,
        error: str | None = None,
        artifacts: list | None = None,
        stats: dict | None = None,
    ) -> dict | None:
        record = self._read(job_id)
        if record is None:
            return None
        if status:
            record["status"] = status
            if status == "running" and not record["started_at"]:
                record["started_at"] = _now()
            if status in ("completed", "failed", "stopped", "cancelled", "paused"):
                record["finished_at"] = _now() if status != "paused" else None
        if progress is not None:
            record["progress"] = max(0, min(100, int(progress)))
        if phase is not None:
            record["phase"] = phase
        if message is not None:
            record["message"] = str(message)[:500]
        if error is not None:
            record["error"] = str(error)[:2000]
        if artifacts is not None:
            record["artifacts"] = artifacts
        if stats:
            record["stats"].update({k: v for k, v in stats.items() if v is not None})
        self._write(record)
        return record

    def add_artifact(self, job_id: str, name: str, path: Path, kind: str = "file") -> None:
        record = self._read(job_id)
        if record is None:
            return
        artifacts = [a for a in record["artifacts"] if a["path"] != str(path)]
        artifacts.append({"name": name, "path": str(path), "kind": kind})
        record["artifacts"] = artifacts
        self._write(record)

    def get(self, job_id: str) -> dict | None:
        record = self._read(job_id)
        if record is None:
            return None
        thread = self._threads.get(job_id)
        # 线程已死但状态未落盘的兜底（崩溃恢复）
        if record["status"] == "running" and (thread is None or not thread.is_alive()):
            record["status"] = "failed"
            record["error"] = record.get("error") or "工作线程意外退出（进程可能重启），可尝试断点恢复"
            self._write(record)
        return record

    def list(self, kind: str | None = None, limit: int = 200) -> list[dict]:
        out = []
        for d in sorted(self.jobs_dir.glob("*/job.json"), key=lambda p: p.stat().st_mtime, reverse=True):
            try:
                record = json.loads(d.read_text(encoding="utf-8"))
            except Exception:
                continue
            if kind and record.get("kind") != kind:
                continue
            if record.get("status") == "running":
                record = self.get(record["id"]) or record
            out.append(record)
            if len(out) >= limit:
                break
        return out

    def delete(self, job_id: str) -> bool:
        control = self._controls.get(job_id)
        if control:
            control.request_stop()
        d = self.job_dir(job_id)
        if not d.is_dir():
            return False
        shutil.rmtree(d, ignore_errors=True)
        self._threads.pop(job_id, None)
        self._controls.pop(job_id, None)
        return True

    # ---------- 执行 ----------
    def submit(self, job_id: str, worker: WORKER_FN, resume: bool = False) -> dict:
        record = self.get(job_id)
        if record is None:
            raise KeyError(job_id)
        old = self._threads.get(job_id)
        if old and old.is_alive():
            record["message"] = "任务已在运行中"
            self._write(record)
            return record
        control = JobControl(self, job_id)
        self._controls[job_id] = control
        if not resume:
            self.update(job_id, status="running", progress=0, error=None)
        else:
            self.update(job_id, status="running", error=None)

        def runner() -> None:
            log_path = self.job_dir(job_id) / "logs" / "run.log"
            try:
                result = worker(record, control) or {}
                artifacts = result.pop("artifacts", [])
                self.update(
                    job_id,
                    status="completed",
                    progress=100,
                    phase="done",
                    message=result.pop("message", "完成"),
                    stats=result.get("stats"),
                )
                for a in artifacts:
                    self.add_artifact(job_id, a.get("name", Path(a["path"]).name), Path(a["path"]), a.get("kind", "file"))
                with open(log_path, "a", encoding="utf-8") as f:
                    f.write(f"[{_now()}] completed: {json.dumps(result.get('stats', {}), ensure_ascii=False)}\n")
            except JobPaused:
                self.update(job_id, status="paused", phase="paused", message="已暂停，可从断点恢复")
            except JobInterrupted:
                self.update(job_id, status="stopped", phase="stopped", message="已停止，现有结果已保留")
            except Exception as exc:  # noqa: BLE001
                with open(log_path, "a", encoding="utf-8") as f:
                    f.write(f"[{_now()}] FAILED:\n{traceback.format_exc()}\n")
                self.update(job_id, status="failed", phase="failed", message=str(exc), error=traceback.format_exc()[-2000:])
            finally:
                self._threads.pop(job_id, None)

        t = threading.Thread(target=runner, name=f"job-{job_id}", daemon=True)
        self._threads[job_id] = t
        t.start()
        return self.get(job_id) or record

    def pause(self, job_id: str) -> dict | None:
        control = self._controls.get(job_id)
        if control:
            control.request_pause()
        return self.update(job_id, message="已请求暂停，当前工作单元完成后写入断点")

    def resume(self, job_id: str, worker: WORKER_FN) -> dict | None:
        record = self.get(job_id)
        if record is None:
            return None
        if record["status"] != "paused":
            record["message"] = "只有已暂停的任务可以恢复"
            self._write(record)
            return record
        return self.submit(job_id, worker, resume=True)

    def stop(self, job_id: str) -> dict | None:
        control = self._controls.get(job_id)
        if control:
            control.request_stop()
        record = self.get(job_id)
        if record and record["status"] in ("queued", "running", "paused", "stopping"):
            self.update(job_id, status="stopping", message="已请求停止，现有结果将保留")
        return record


manager = JobManager()
