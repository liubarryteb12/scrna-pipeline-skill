#!/usr/bin/env python3
"""
01_qc.py — 质量控制与过滤

做四件事，每一件都留记录：
  1. QC 指标（n_genes / total_counts / pct_mt / pct_ribo / pct_hb）
  2. 硬阈值过滤（可配）
  3. 双细胞检测（scrublet）—— **失败不等于跳过，要写明为什么**
  4. ambient RNA（胞外游离 RNA）—— 默认不做，**并写明不做**

**为什么每一条都要落盘。** 用户文档把双细胞与 ambient RNA 列为单细胞 QC
必做项。不做是可以的（有些数据确实做不了），但"没做"必须和"做了没问题"
在产物里长得不一样 —— 否则读报告的人以为该查的都查了。
Part 1 的教训（规则 24）：可选步骤失败时 job 仍然是绿的。
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))

import numpy as np  # noqa: E402
import scanpy as sc  # noqa: E402

from common import (ensure_dirs, load_config, log_info, log_warn,  # noqa: E402
                    parse_args, record_step, save_fig, set_seed, write_json,
                    W_DOUBLE, W_SINGLE, mm, PAL,)

# 血红蛋白基因（红细胞污染）与核糖体基因的前缀
HB_PREFIXES = ("HBA", "HBB", "HBD", "HBE", "HBG", "HBM", "HBQ", "HBZ")
RIBO_PREFIXES = ("RPS", "RPL")


def add_qc_metrics(adata, organism: str = "Homo sapiens") -> list:
    """算 QC 指标。线粒体前缀按物种选 —— 小鼠是 mt-，人是 MT-。"""
    if organism.lower().startswith("mus") or organism.lower().startswith("mouse"):
        mt_pref, ribo_pref = ("mt-",), ("Rps", "Rpl")
    else:
        mt_pref, ribo_pref = ("MT-",), RIBO_PREFIXES

    adata.var["mt"] = adata.var_names.str.startswith(mt_pref)
    adata.var["ribo"] = adata.var_names.str.startswith(ribo_pref)
    adata.var["hb"] = adata.var_names.str.startswith(HB_PREFIXES)

    qc_vars = [v for v in ("mt", "ribo", "hb") if bool(adata.var[v].any())]
    sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, percent_top=None,
                               log1p=False, inplace=True)
    return qc_vars


def run_scrublet(adata, cfg: dict) -> dict:
    """
    双细胞检测。

    **scrublet 依赖 scikit-image**（做模拟双细胞的邻居图）。缺它时
    scanpy 抛 ImportError —— 这里捕获并记录 `package_missing`，
    而不是让整个步骤失败，也不是静默跳过。

    返回的记录会进 qc_status.json：读者能看出"跑了、没跑、还是跑失败了"。
    """
    if not cfg["qc"].get("scrublet", True):
        return {"status": "disabled", "reason": "配置 qc.scrublet=false"}

    try:
        import skimage  # noqa: F401
    except ImportError:
        return {"status": "package_missing",
                "reason": "scikit-image 未安装（scrublet 依赖它做模拟双细胞的邻居图）",
                "fix": "pip install scikit-image"}

    try:
        # scrublet 的期望双细胞率：10x 数据常见 5%-10%。用细胞数反推，
        # 但不低于 0.05 —— 小数据上用默认 0.05 会让阈值不稳。
        n = adata.n_obs
        rate = float(np.clip(5000 / max(n, 1) * 0.01, 0.05, 0.10))
        sc.pp.scrublet(adata, expected_doublet_rate=rate, random_state=cfg["analysis"]["seed"])
        n_db = int(adata.obs["predicted_doublet"].sum())
        return {"status": "ok",
                "expected_doublet_rate": round(rate, 4),
                "n_predicted_doublets": n_db,
                "frac_predicted_doublets": round(n_db / max(n, 1), 5)}
    except Exception as e:  # noqa: BLE001
        return {"status": "failed", "reason": f"{type(e).__name__}: {e}"}


def assess_ambient_rna(adata, cfg: dict) -> dict:
    """
    ambient RNA（胞外游离 RNA）评估。

    **默认不做，而且如实写明。** 正规做法是 SoupX（R，需要空液滴）或
    CellBender（Python，需要 GPU 或极长的 CPU 训练）。本流水线在
    GitHub 托管 runner 上跑，没有 GPU，CellBender 不现实；SoupX 是 R 包，
    混进 Python 流水线会让依赖变得很脆。

    **不编一个代理指标冒充。** 能做的诚实事是：报一个下界信号 ——
    高丰度基因在"不该表达它的细胞"里的背景水平（例如血红蛋白基因在
    非红细胞里的表达）。这**不是** ambient RNA 校正，只是一个提示。
    """
    mode = cfg["qc"].get("ambient_rna", "none")
    if mode == "none":
        return {"status": "not_done",
                "reason": ("配置 qc.ambient_rna=none。正规做法需要 SoupX（R，需空液滴）"
                           "或 CellBender（需 GPU）；GitHub 托管 runner 无 GPU，"
                           "故未执行，也未用代理指标冒充校正结果")}

    if mode == "simple":
        # 粗略背景提示：血红蛋白基因在所有细胞里的中位占比。
        # 红细胞污染高时这个值会明显抬高，但它分不清"污染"与"真的有红细胞"。
        out = {"status": "heuristic_only",
               "note": ("**这不是 ambient RNA 校正**，只是背景水平提示；"
                        "要真正校正需 SoupX/CellBender")}
        if "pct_counts_hb" in adata.obs.columns:
            v = adata.obs["pct_counts_hb"].astype(float)
            out["pct_counts_hb_median"] = round(float(v.median()), 4)
            out["pct_counts_hb_p95"] = round(float(v.quantile(0.95)), 4)
        return out

    return {"status": "unknown_mode", "reason": f"qc.ambient_rna='{mode}' 不是 none/simple"}


def run_01_qc(cfg: dict) -> dict:
    ensure_dirs(cfg)
    set_seed(cfg)
    data_dir = Path(cfg["output"]["data_dir"])
    res_dir = Path(cfg["output"]["results_dir"])

    info = None
    try:
        import json
        with open(data_dir / "dataset_info.json", encoding="utf-8") as fh:
            info = json.load(fh)
    except Exception:  # noqa: BLE001
        pass
    organism = (info or {}).get("organism") or "Homo sapiens"

    adata = sc.read_h5ad(data_dir / "raw.h5ad")
    n0, g0 = adata.n_obs, adata.n_vars
    log_info(f"读入 {n0} 细胞 x {g0} 基因")

    # ---- 1. 指标 ------------------------------------------------------------
    qc_vars = add_qc_metrics(adata, organism)
    log_info(f"QC 指标已算（{', '.join(qc_vars) if qc_vars else '无线粒体/核糖体基因命中'}）")

    # ---- 2. 过滤前的图 ------------------------------------------------------
    # 先画再滤 —— 滤完再画就看不到"滤掉了什么"，而那正是要判断的东西。
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    keys = [k for k in ("n_genes_by_counts", "total_counts", "pct_counts_mt",
                        "pct_counts_ribo", "pct_counts_hb") if k in adata.obs.columns]
    fig, axes = plt.subplots(1, len(keys), figsize=(W_DOUBLE, mm(58)))
    axes = np.atleast_1d(axes)
    for ax, k in zip(axes, keys):
        ax.violinplot(adata.obs[k].astype(float).values, showmedians=True)
        ax.set_title(k)
        ax.set_xticks([])
    fig.suptitle(f"QC metrics before filtering (n={n0})")
    save_fig(cfg, "qc_violin_before", fig)

    # 阈值线：基因数 vs 线粒体比例 —— 双细胞和死细胞在这张图上是两个角
    fig, ax = plt.subplots(figsize=(W_SINGLE, mm(64)))
    sc_ = ax.scatter(adata.obs["total_counts"], adata.obs["n_genes_by_counts"],
                     c=adata.obs["pct_counts_mt"] if "pct_counts_mt" in adata.obs else None,
                     s=3, cmap="viridis", alpha=0.6)
    if sc_ is not None:
        fig.colorbar(sc_, ax=ax, label="pct_counts_mt")
    q = cfg["qc"]
    ax.axhline(q["min_genes"], color=PAL["highlight"], lw=1, ls="--")
    ax.axhline(q["max_genes"], color=PAL["highlight"], lw=1, ls="--")
    ax.set_xlabel("total_counts"); ax.set_ylabel("n_genes_by_counts")
    ax.set_title("QC thresholds (red = cut-offs)")
    save_fig(cfg, "qc_scatter_thresholds", fig)

    # ---- 3. 过滤 ------------------------------------------------------------
    q = cfg["qc"]
    n_before = adata.n_obs
    adata.var["mt"] = adata.var_names.str.startswith(
        ("mt-",) if organism.lower().startswith("mus") else ("MT-",))
    sc.pp.filter_cells(adata, min_genes=int(q["min_genes"]))
    sc.pp.filter_genes(adata, min_cells=int(q["min_cells"]))
    n_after_gene = adata.n_obs
    if q.get("max_genes"):
        adata = adata[adata.obs["n_genes_by_counts"] < int(q["max_genes"])].copy()
    n_after_maxg = adata.n_obs
    if "pct_counts_mt" in adata.obs.columns and q.get("max_pct_mt") is not None:
        adata = adata[adata.obs["pct_counts_mt"] < float(q["max_pct_mt"])].copy()
    n_after_mt = adata.n_obs

    log_info(f"过滤: {n_before} -> {n_after_gene} (min_genes/min_cells) "
             f"-> {n_after_maxg} (max_genes) -> {n_after_mt} (max_pct_mt)")

    # ---- 4. 双细胞 ----------------------------------------------------------
    db = run_scrublet(adata, cfg)
    n_after_db = adata.n_obs
    if db["status"] == "ok" and db["n_predicted_doublets"] > 0:
        adata = adata[~adata.obs["predicted_doublet"]].copy()
        n_after_db = adata.n_obs
        log_info(f"去除双细胞: {n_after_db + db['n_predicted_doublets']} -> {n_after_db}")
    elif db["status"] == "ok":
        log_info("scrublet 未标出双细胞")
    else:
        log_warn(f"双细胞检测未执行: {db['status']} —— {db.get('reason', '')}")

    # ---- 5. ambient RNA -----------------------------------------------------
    amb = assess_ambient_rna(adata, cfg)
    if amb["status"] != "ok":
        log_warn(f"ambient RNA: {amb['status']} —— {amb.get('reason', amb.get('note', ''))}")

    if adata.n_obs < 50:
        raise RuntimeError(f"过滤后只剩 {adata.n_obs} 个细胞，无法继续分析")

    # ---- 6. 落盘 ------------------------------------------------------------
    out = data_dir / "qc_filtered.h5ad"
    adata.write_h5ad(out)
    log_info(f"已写出 {out}（{adata.n_obs} 细胞 x {adata.n_vars} 基因）")

    # QC 汇总表：每个细胞一行，便于事后复查"到底滤掉了谁"
    cols = [c for c in ("n_genes_by_counts", "total_counts", "pct_counts_mt",
                        "pct_counts_ribo", "pct_counts_hb", "predicted_doublet",
                        "doublet_score") if c in adata.obs.columns]
    adata.obs[cols].to_csv(res_dir / "qc_cells.csv")

    status = {
        "dataset_id": cfg["dataset_id"],
        "organism": organism,
        "n_cells_raw": int(n0),
        "n_genes_raw": int(g0),
        "filtering": {
            "min_genes": int(q["min_genes"]),
            "max_genes": q.get("max_genes"),
            "max_pct_mt": q.get("max_pct_mt"),
            "min_cells": int(q["min_cells"]),
            "n_after_min_genes": int(n_after_gene),
            "n_after_max_genes": int(n_after_maxg),
            "n_after_mt": int(n_after_mt),
            "n_after_doublets": int(n_after_db),
            "n_removed_total": int(n0 - adata.n_obs),
            "frac_removed": round((n0 - adata.n_obs) / max(n0, 1), 5),
        },
        "doublet_detection": db,
        "ambient_rna": amb,
        "n_cells_final": int(adata.n_obs),
        "n_genes_final": int(adata.n_vars),
        "status": "ok",
    }
    write_json(res_dir / "qc_status.json", status)
    return status


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    t0 = time.time()
    try:
        run_01_qc(cfg)
        record_step(cfg, "qc", "ok", time.time() - t0)
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "qc", "failed", time.time() - t0, message=str(e))
        raise
