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
                    parse_args, record_step, result_status_of, set_seed, write_json)

# 硬门禁：低于/高于这些值的输入不做分析
MIN_CELLS = 50
MIN_GENES = 50
MAX_CELLS_WARN = 200_000   # 超过只警告（内存与时间），不拒绝
# 为什么是 20 万（L18）：本地实测 ~2 万细胞用 ~1.5 GB、CI runner 7 GB 内存，
# 20 万约 15 GB —— 已经**超过 runner 内存**，只告警不拒绝是因为有些分析
# （只做 QC + 聚类）确实跑得完。这个依据原先只存在于作者脑子里，
# 现在写下来，否则没人知道这个数字能不能改。


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

    返回 `(adata, extract_filter)`。第二个值是 `"data"` 或
    `"none_fallback"`（M14：老 Python 上 `filter=` 不可用时**必须留下痕迹**，
    不能"用了"和"没用"在产物里长得一样）。
    """
    import scanpy as sc

    extract_to = cache_dir / "10x_extracted"
    marker = extract_to / ".done"
    extract_filter = "cached"
    if not marker.exists():
        extract_to.mkdir(parents=True, exist_ok=True)
        log_info(f"解压 {tar_path.name} -> {extract_to}")
        with tarfile.open(tar_path, "r:gz") as tf:
            # Python 3.12+ 的 tarfile 有 `filter=` 参数（防路径穿越）；老版本
            # 没有，传了会 `TypeError`。
            #
            # **M14（R-03 裁决）**：原来的回退分支 `except TypeError:
            # tf.extractall(extract_to)` **静默放弃安全过滤且不打任何日志** ——
            # 于是"这份 tar 被检查过"和"这份 tar 被无条件解压"在日志里
            # 长得一样。修法：回退时显式 WARN + 记进 status，把降级说出来。
            # （不改成 raise：老 Python 上会直接不可用，而这是数据获取步骤。）
            try:
                tf.extractall(extract_to, filter="data")
                extract_filter = "data"
            except TypeError:
                extract_filter = "none_fallback"
                log_warn(
                    f"本机 Python {sys.version_info.major}.{sys.version_info.minor} "
                    f"的 tarfile 不支持 filter= 参数（需 3.12+）—— "
                    f"**本次解压未做路径穿越过滤**（10x 官方 tar 包，风险低，"
                    f"但这是一次降级，不是正常路径）")
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
    return adata, extract_filter


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

    **抽查的是"跨全表的非零元素"，不是"前 200 个细胞"**（审计 S10）。
    旧实现取 `X[:200, :]` —— 若数据按样本排序、前 200 个恰好来自某个
    未归一化的样本，校验会给出 `is_counts: True` 而整体不是计数。
    稀疏矩阵的零天然是整数，所以用非零值判整数性是**合理的**；
    要修的是**抽样方式**：按行等距抽两批，两批结论一致才算数。
    """
    import numpy as np
    from scipy import sparse

    X = adata.X
    n_obs = int(X.shape[0])

    def _vals_of(rows):
        if rows.size == 0:
            return np.array([])
        sub = X[rows, :]
        v = sub.data if sparse.issparse(sub) else np.asarray(sub).ravel()
        return v[np.isfinite(v)]

    def _strided(k: int = 200):
        """跨全表等距抽 k 行 —— **主判据**。"""
        if n_obs == 0:
            return np.array([]), 0
        step = max(1, n_obs // k)
        rows = np.arange(0, n_obs, step)[:k]
        return _vals_of(rows), int(rows.size)

    def _contiguous(start: int, k: int = 200):
        """从 `start` 起抽**连续** k 行 —— **交叉核对**。

        两批必须抽法不同才有意义：若两批都是等距抽样，它们覆盖的是
        同一批细胞，占比天然接近，`sampling_consistent` 就成了恒真的
        摆设。连续块才能暴露"数据按样本拼接、各样本处理方式不同"——
        那正是单看前 200 行会误判的场景。
        """
        if n_obs == 0:
            return np.array([]), 0
        rows = np.arange(start, min(start + k, n_obs))
        return _vals_of(rows), int(rows.size)

    vals, n_rows1 = _strided()
    vals2, n_rows2 = _contiguous(n_obs // 2)
    if vals.size == 0:
        return {"is_counts": False, "reason": "矩阵里没有有限值",
                "n_cells_checked": 0}

    is_int = bool(np.allclose(vals, np.round(vals), atol=1e-8))
    min_v = float(vals.min())
    max_v = float(vals.max())
    frac_int = float(np.mean(np.isclose(vals, np.round(vals), atol=1e-8)))

    # 两批结论必须一致 —— 不一致说明数据内部不均匀（例如按样本拼接），
    # 此时"是/不是计数"这个二值结论本身不可靠，要如实说出来。
    frac_int2 = (float(np.mean(np.isclose(vals2, np.round(vals2), atol=1e-8)))
                 if vals2.size else None)
    consistent = (frac_int2 is None) or (abs(frac_int - frac_int2) < 0.01)

    ok = is_int and min_v >= 0 and consistent
    if not consistent:
        reason = (f"**两批抽查结论不一致**（等距抽样整数值占比 {frac_int:.1%} vs "
                  f"中段连续抽样 {frac_int2:.1%}）—— 数据可能按样本拼接且各样本"
                  "处理方式不同，此时『是/不是计数』的二值判定不可靠")
    elif ok:
        reason = (f"非负整数 -> 判定为原始计数（跨 {n_obs} 个细胞等距抽 "
                  f"{n_rows1} 行 + 中段连续抽 {n_rows2} 行，两批一致）")
    else:
        reason = (f"**不是整数计数**（最小 {min_v:.4g}，最大 {max_v:.4g}，"
                  f"整数值占比 {frac_int:.1%}）—— 可能已经 normalize/log 过")
    return {"is_counts": ok, "reason": reason, "min": min_v, "max": max_v,
            "frac_integer": round(frac_int, 6),
            "frac_integer_second_batch": (round(frac_int2, 6)
                                          if frac_int2 is not None else None),
            "sampling_consistent": consistent,
            "sampling_note": ("第一批 = 跨全表等距抽样；第二批 = 从中段起连续抽样。"
                              "两批抽法不同，占比若差 >1 个百分点即判为不一致"),
            "n_cells_total": n_obs,
            "n_cells_checked": int(n_rows1 + n_rows2),
            "sampled_values": int(vals.size + vals2.size)}


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
    extract_filter = None
    if kind == "10x_tar":
        tar_path = download(url, cache_dir / Path(url).name)
        adata, extract_filter = read_10x_tar(tar_path, cache_dir)
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
        # **M14**：解压时有没有做路径穿越过滤。`none_fallback` 是降级路径。
        "extract_filter": extract_filter,
        # **L18（R-03 裁决）**：这个阈值原来只 WARN、且不落盘 ——
        # 于是"细胞数 21 万，已经越过提示线"这件事在产物里查不到，
        # 事后想解释 CI 为什么慢/为什么内存吃紧时没有任何依据。
        "max_cells_warn_gate": MAX_CELLS_WARN,
        "n_cells_over_warn_gate": bool(adata.n_obs > MAX_CELLS_WARN),
    }
    write_json(data_dir / "dataset_info.json", info)
    return info


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    t0 = __import__("time").time()
    try:
        res = run_00_fetch(cfg)
        # E-56：接住返回值。本步骤的返回 dict 没有 `status` 键（它写的是
        # `dataset_info.json`），`result_status_of` 会给出 "ok" —— 与其他步骤
        # 保持同一种记法，将来加了 status 键也自动生效。
        record_step(cfg, "fetch", "ok", __import__("time").time() - t0,
                    result_status=result_status_of(res))
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "fetch", "failed", __import__("time").time() - t0,
                    message=str(e), result_status="failed")
        raise
