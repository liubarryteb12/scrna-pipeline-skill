#!/usr/bin/env python3
"""
00_fetch.py — 取数、读入、硬校验

这一步是整条流水线的**硬门禁**：数据必须真的是计数矩阵（非负整数），
细胞数和基因数必须落在可分析范围内，否则直接抛异常，不进入任何分析。

**为什么要校验"是不是整数计数"。** 单细胞流程里最常见的静默错误是把
已经 normalize 过（或 log 过）的矩阵再当原始计数喂一遍：QC 指标（total_counts）
失去意义、HVG 的 seurat_v3 口味会报错或给出错的结果、pseudobulk 的
DESeq2 负二项模型前提被破坏 —— 而**图看起来完全正常**。

输出：
  data/<id>/raw.h5ad          原始计数（未过滤）
  data/<id>/dataset_info.json 数据来源、形状、计数性质、校验结论
"""

from __future__ import annotations

import sys
import tarfile
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))

from common import (ensure_dirs, load_config, log_info, log_warn,  # noqa: E402
                    parse_args, record_step, set_seed, write_json)

# 硬门禁：低于/高于这些值的输入不做分析
MIN_CELLS = 50
MIN_GENES = 50
MAX_CELLS_WARN = 200_000   # 超过只警告（内存与时间），不拒绝


def load_registry(repo_root: Path) -> dict:
    import yaml
    p = repo_root / "assets" / "datasets.yml"
    if not p.exists():
        return {}
    with open(p, encoding="utf-8") as fh:
        return (yaml.safe_load(fh) or {}).get("datasets", {}) or {}


def download(url: str, dest: Path, retries: int = 3) -> Path:
    """
    带重试的下载。已存在且非空则跳过（CI 里 data/ 不持久化，但本地会命中）。

    **必须带 User-Agent。** 10x 的 CDN（cf.10xgenomics.com）对
    `Python-urllib/3.x` 直接返回 403，而同样的 URL 用浏览器或 PowerShell
    请求是 200 —— 实测踩过，报错信息只有 "HTTP Error 403: Forbidden"，
    看不出是 UA 的问题，很容易误判成"链接失效了"。
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        log_info(f"已缓存，跳过下载: {dest.name} ({dest.stat().st_size/1e6:.1f} MB)")
        return dest
    last = None
    for attempt in range(1, retries + 1):
        try:
            log_info(f"下载 ({attempt}/{retries}): {url}")
            tmp = dest.with_suffix(dest.suffix + ".part")
            req = urllib.request.Request(url, headers={
                "User-Agent": ("Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
                               "(KHTML, like Gecko) Chrome/120.0 Safari/537.36"),
                "Accept": "*/*",
            })
            with urllib.request.urlopen(req, timeout=180) as r, open(tmp, "wb") as fh:
                total = int(r.headers.get("Content-Length") or 0)
                got = 0
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
                    got += len(chunk)
            if total and got != total:
                raise IOError(f"下载不完整: {got}/{total} 字节")
            if got == 0:
                raise IOError("下载到 0 字节")
            tmp.replace(dest)
            log_info(f"下载完成: {dest.name} ({got/1e6:.1f} MB)")
            return dest
        except Exception as e:  # noqa: BLE001
            last = e
            log_warn(f"下载失败 ({attempt}/{retries}): {e}")
    raise RuntimeError(f"下载 {url} 失败（重试 {retries} 次）: {last}")


def read_10x_tar(tar_path: Path, cache_dir: Path):
    """
    解压 10x tar.gz 并读入。

    **10x 的 tar 里目录名带版本后缀**（filtered_gene_bc_matrices/hg19/），
    不能写死路径 —— 不同版本的目录结构不同。这里找到含 matrix.mtx 的目录。
    """
    import scanpy as sc

    extract_to = cache_dir / "10x_extracted"
    marker = extract_to / ".done"
    if not marker.exists():
        extract_to.mkdir(parents=True, exist_ok=True)
        log_info(f"解压 {tar_path.name} -> {extract_to}")
        with tarfile.open(tar_path, "r:gz") as tf:
            # Python 3.12+ 的 tarfile 有 data 过滤；老版本没有这个参数
            try:
                tf.extractall(extract_to, filter="data")
            except TypeError:
                tf.extractall(extract_to)
        marker.write_text("ok", encoding="utf-8")
    else:
        log_info(f"已解压，跳过: {extract_to}")

    # 找含 matrix.mtx(.gz) 的目录
    hits = [p.parent for p in extract_to.rglob("matrix.mtx*")]
    if not hits:
        raise RuntimeError(f"{tar_path.name} 里找不到 matrix.mtx —— 不是 10x 计数矩阵格式")
    # 有多个时取路径最浅的（避免命中子目录里的备份）
    mtx_dir = sorted(hits, key=lambda p: (len(p.parts), str(p)))[0]
    log_info(f"10x 矩阵目录: {mtx_dir.relative_to(extract_to)}")
    adata = sc.read_10x_mtx(mtx_dir, var_names="gene_symbols", cache=False)
    return adata


def read_h5ad(path: Path):
    import anndata as ad
    return ad.read_h5ad(path)


def read_10x_h5(path: Path):
    import scanpy as sc
    return sc.read_10x_h5(path)


def validate_counts(adata, cfg: dict) -> dict:
    """
    判断 X 是不是原始整数计数。

    **这是本步骤存在的理由。** 把已 normalize 的矩阵当计数喂进去，
    QC 指标失真、seurat_v3 HVG 前提被破坏、pseudobulk 的负二项模型不成立,
    而所有图和数字看起来都正常。
    """
    import numpy as np
    from scipy import sparse

    X = adata.X
    sample = X[:min(200, X.shape[0]), :]
    vals = sample.data if sparse.issparse(sample) else np.asarray(sample).ravel()
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return {"is_counts": False, "reason": "矩阵里没有有限值"}

    is_int = bool(np.allclose(vals, np.round(vals), atol=1e-8))
    min_v = float(vals.min())
    max_v = float(vals.max())
    frac_int = float(np.mean(np.isclose(vals, np.round(vals), atol=1e-8)))

    ok = is_int and min_v >= 0
    reason = ("非负整数 -> 判定为原始计数"
              if ok else
              f"**不是整数计数**（最小 {min_v:.4g}，最大 {max_v:.4g}，"
              f"整数值占比 {frac_int:.1%}）—— 可能已经 normalize/log 过")
    return {"is_counts": ok, "reason": reason, "min": min_v, "max": max_v,
            "frac_integer": round(frac_int, 6),
            "sampled_values": int(vals.size)}


def run_00_fetch(cfg: dict) -> dict:
    import scanpy as sc  # noqa: F401  (确认依赖存在)
    import anndata as ad

    ensure_dirs(cfg)
    set_seed(cfg)
    repo_root = Path(__file__).resolve().parent.parent
    data_dir = Path(cfg["output"]["data_dir"])
    cache_dir = data_dir / "cache"

    src = cfg.get("source") or {}
    if src.get("dataset"):
        reg = load_registry(repo_root)
        key = src["dataset"]
        if key not in reg:
            raise KeyError(f"source.dataset='{key}' 不在 assets/datasets.yml 里；"
                           f"可选: {', '.join(sorted(reg))}")
        entry = dict(reg[key])
        entry.update({k: v for k, v in src.items() if k != "dataset"})
        log_info(f"数据集 '{key}': {entry.get('note', '')}")
    elif src.get("url"):
        entry = dict(src)
    else:
        raise ValueError("配置缺少 source.dataset 或 source.url")

    kind = entry.get("kind")
    url = entry.get("url")
    if not kind or not url:
        raise ValueError(f"数据源条目缺少 kind/url: {entry}")

    # ---- 取数 --------------------------------------------------------------
    if kind == "10x_tar":
        tar_path = download(url, cache_dir / Path(url).name)
        adata = read_10x_tar(tar_path, cache_dir)
    elif kind == "10x_h5":
        h5_path = download(url, cache_dir / Path(url).name)
        adata = read_10x_h5(h5_path)
    elif kind == "h5ad":
        h5_path = download(url, cache_dir / Path(url).name)
        adata = read_h5ad(h5_path)
    else:
        raise ValueError(f"不支持的 kind: {kind}（可选 10x_tar / 10x_h5 / h5ad）")

    adata.var_names_make_unique()
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise RuntimeError(f"读入后为空: {adata.n_obs} 细胞 x {adata.n_vars} 基因")

    # ---- 硬校验 ------------------------------------------------------------
    counts = validate_counts(adata, cfg)
    log_info(f"计数性质: {counts['reason']}")
    if not counts["is_counts"]:
        raise RuntimeError(
            "输入矩阵不是原始整数计数，拒绝继续。\n"
            f"  依据: {counts['reason']}\n"
            "  为什么必须停: QC 指标会失真、seurat_v3 HVG 前提被破坏、\n"
            "  pseudobulk 的负二项模型不成立 —— 而所有图和数字看起来都正常。\n"
            "  修法: 用原始计数矩阵（10x 的 matrix.mtx / filtered_feature_bc_matrix）。")

    if adata.n_obs < MIN_CELLS:
        raise RuntimeError(f"细胞数 {adata.n_obs} < {MIN_CELLS}，单细胞分析无从谈起")
    if adata.n_vars < MIN_GENES:
        raise RuntimeError(f"基因数 {adata.n_vars} < {MIN_GENES}")
    if adata.n_obs > MAX_CELLS_WARN:
        log_warn(f"细胞数 {adata.n_obs} 较大，CI 上可能超出内存/时限")

    # ---- 落盘 --------------------------------------------------------------
    raw_path = data_dir / "raw.h5ad"
    adata.write_h5ad(raw_path)
    log_info(f"已写出 {raw_path}（{adata.n_obs} 细胞 x {adata.n_vars} 基因）")

    info = {
        "dataset_id": cfg["dataset_id"],
        "source_kind": kind,
        "source_url": url,
        "source_note": entry.get("note", ""),
        "organism": entry.get("organism"),
        "tissue": entry.get("tissue"),
        "condition": entry.get("condition"),
        "n_cells_raw": int(adata.n_obs),
        "n_genes_raw": int(adata.n_vars),
        "counts_check": counts,
        "min_cells_gate": MIN_CELLS,
        "min_genes_gate": MIN_GENES,
    }
    write_json(data_dir / "dataset_info.json", info)
    return info


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    t0 = __import__("time").time()
    try:
        run_00_fetch(cfg)
        record_step(cfg, "fetch", "ok", __import__("time").time() - t0)
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "fetch", "failed", __import__("time").time() - t0,
                    message=str(e))
        raise
