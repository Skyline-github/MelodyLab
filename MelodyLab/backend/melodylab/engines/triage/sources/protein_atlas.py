"""Human Protein Atlas —— 组织表达谱。

HPA 无正式 REST API，策略：
1. 主数据：单基因 JSON 端点 https://www.proteinatlas.org/{ensembl}.json
   （当前版本：组织特异/分布分类、top 特异组织 nTPM、亚细胞定位等）
2. 补充：v22 归档批量文件 normal_tissue.tsv.zip（44 种组织的 IHC 蛋白
   表达 Level 完整表格，压缩包约 6MB，缓存于 <data_root>/cache/triage/，
   由后台线程下载、绝不阻塞评估，可离线复用）
"""
from __future__ import annotations

import io
import json
import threading
import zipfile
from pathlib import Path

import httpx

from ..cache import FAILURE_TTL, cache_get, cache_set
from ..config import CACHE_DIR, REQUEST_TIMEOUT, USER_AGENT
from .http_client import request_json

_SERVICE = "hpa"
HPA_GENE_JSON = "https://www.proteinatlas.org/{id}.json"
HPA_TISSUE_TSV_V22 = "https://v22.proteinatlas.org/download/normal_tissue.tsv.zip"

_TSV_PATH = CACHE_DIR / "normal_tissue_v22.tsv"
_TSV_DONE = CACHE_DIR / "normal_tissue_v22.done"
_tsv_index: dict[str, list[dict]] | None = None


# ---------------------------------------------------------------- 主数据
async def get_tissue_expression(symbol: str, ensembl_id: str = "") -> dict:
    """组织表达概要：特异性分类 + top 特异组织 +（若可用）完整组织表。"""
    key = f"{symbol}|{ensembl_id}"
    cached = cache_get(_SERVICE, key)
    if cached is not None:
        return cached

    ident = ensembl_id or symbol
    try:
        data = await request_json("GET", HPA_GENE_JSON.format(id=ident))
    except Exception as exc:  # noqa: BLE001
        out = {"found": False, "note": f"HPA 查询失败（可到 proteinatlas.org 手动确认）：{exc}"}
        cache_set(_SERVICE, key, out, ttl=FAILURE_TTL)
        return out

    if data.get("Gene", "").upper() != symbol.upper():
        # 用符号再试一次
        data = await request_json("GET", HPA_GENE_JSON.format(id=symbol))
        if data.get("Gene", "").upper() != symbol.upper():
            out = {"found": False,
                   "note": "HPA 未收录该基因，请到 proteinatlas.org 手动确认"}
            cache_set(_SERVICE, key, out)
            return out

    tissue_table = await _bulk_tissue_table(symbol)

    out = {
        "found": True,
        "gene": data.get("Gene"),
        "description": data.get("Gene description"),
        "rna_tissue_specificity": data.get("RNA tissue specificity", ""),
        "rna_tissue_distribution": data.get("RNA tissue distribution", ""),
        "protein_tissue_specificity": data.get("Protein tissue specificity", ""),
        "protein_tissue_distribution": data.get("Protein tissue distribution", ""),
        "top_rna_tissues": _top_items(data.get("RNA tissue specific nTPM"), "nTPM"),
        "single_cell_specificity": data.get("RNA single cell type specificity", ""),
        "subcellular_location": (data.get("Subcellular location") or [])[:8],
        "secretome_location": data.get("Secretome location", ""),
        "blood_concentration_ms": data.get("Blood concentration - Conc. blood MS [pg/L]"),
        "tissue_table": tissue_table,
        "tissue_table_version": "v22" if tissue_table else None,
    }
    cache_set(_SERVICE, key, out)
    return out


def _top_items(mapping: dict | None, unit: str, n: int = 8) -> list[dict]:
    if not mapping:
        return []
    items = [{"tissue": k, unit: float(v)} for k, v in mapping.items()]
    return sorted(items, key=lambda x: -x[unit])[:n]


# ---------------------------------------------------------------- 批量表格
# v22 归档下载在慢网络下可能极慢（实测有被限速到 ~500 B/s 的环境），
# 因此绝不阻塞评估：未就绪时立即返回 []，由后台守护线程下载，
# 完成后写 .done 标记，之后的评估离线复用。
_TSV_LOCK = threading.Lock()
_BULK_TIMEOUT = 1800.0  # 后台下载超时（秒）


async def _bulk_tissue_table(symbol: str) -> list[dict]:
    """从 v22 批量文件取完整组织表；未就绪返回 []（不阻塞报告）。"""
    idx = _load_if_ready()
    if idx is None:
        _start_background_download()
        return []
    return idx.get(symbol.upper(), [])[:60]


def _load_if_ready() -> dict[str, list[dict]] | None:
    global _tsv_index
    if _tsv_index is not None:
        return _tsv_index
    if _TSV_DONE.exists() and _TSV_PATH.exists():
        import pandas as pd
        _tsv_index = _load_bulk(_TSV_PATH, pd)
        return _tsv_index
    return None


def _start_background_download() -> None:
    if not _TSV_LOCK.acquire(blocking=False):
        return  # 已有下载在跑

    def _worker():
        global _tsv_index
        try:
            import asyncio as _asyncio

            async def _fetch() -> bytes:
                async with httpx.AsyncClient(timeout=_BULK_TIMEOUT,
                                             follow_redirects=True) as client:
                    resp = await client.get(HPA_TISSUE_TSV_V22,
                                            headers={"User-Agent": USER_AGENT})
                    resp.raise_for_status()
                    return resp.content

            content = _asyncio.run(_fetch())
            with zipfile.ZipFile(io.BytesIO(content)) as zf:
                name = zf.namelist()[0]
                with zf.open(name) as fh:
                    _TSV_PATH.write_bytes(fh.read())
            _TSV_DONE.write_text("ok")
            import pandas as pd
            _tsv_index = _load_bulk(_TSV_PATH, pd)
        except Exception:  # noqa: BLE001 —— 保持未就绪，下次评估再试
            pass
        finally:
            _TSV_LOCK.release()

    threading.Thread(target=_worker, name="hpa-bulk-download", daemon=True).start()


def _load_bulk(path: Path, pd) -> dict[str, list[dict]]:
    """把 v22 normal_tissue.tsv 建成 基因→组织行 的索引。

    v22 文件为 组织×细胞类型 粒度（无 NX 列），按组织聚合取最高 IHC Level。
    """
    df = pd.read_csv(path, sep="\t")
    rank = {"Not detected": 0, "Low": 1, "Medium": 2, "High": 3}
    df["_rank"] = df["Level"].map(rank).fillna(0)
    idx: dict[str, list[dict]] = {}
    for gene, grp in df.groupby("Gene name"):
        agg = grp.groupby("Tissue")["_rank"].max()
        inv = {v: k for k, v in rank.items()}
        rows = [{"tissue": t, "protein_level": inv[int(v)], "rna_nx": None}
                for t, v in agg.items()]
        rows.sort(key=lambda x: -rank.get(x["protein_level"], 0))
        idx[str(gene).upper()] = rows
    return idx
