#!/usr/bin/env python3
"""
05_trajectory.py — 轨迹推断（PAGA + 扩散拟时序）

**轨迹分析的结论比聚类弱得多，这一点必须写进产物。**
  - PAGA 给出的是**簇之间的连通图**（谁和谁相连、连接强度），
    不是"分化方向"。箭头方向来自 DPT，而 DPT 需要**根**，
    根选错则整条轨迹反向。
  - 拟时序是**一维坐标**，把高维的分化过程压成一条线。分支过程
    （一个祖先进两种细胞）在一维坐标里会被压成"先后"，而实际是"并列"。
  - 没有 RNA 速率（spliced/unspliced）时，拟时序**不能**说明方向性，
    只能说明"相似度排序"。本流水线只有计数矩阵，**没有速率**。

所以这里报的是：PAGA 连通性 + DPT 排序 + **根是怎么选的**。
不报"分化轨迹"这种说法。
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402

from common import (df_to_records, ensure_dirs, load_config, log_info,  # noqa: E402
                    log_warn, parse_args, record_step, save_fig, set_seed, write_json)


def run_05_trajectory(cfg: dict) -> dict:
    ensure_dirs(cfg)
    set_seed(cfg)
    data_dir = Path(cfg["output"]["data_dir"])
    res_dir = Path(cfg["output"]["results_dir"])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    traj = cfg.get("trajectory") or {}
    if not traj.get("enabled", True):
        status = {"dataset_id": cfg["dataset_id"], "status": "disabled",
                  "reason": "配置 trajectory.enabled=false"}
        write_json(res_dir / "trajectory_status.json", status)
        log_info("轨迹分析已按配置关闭")
        return status

    adata = sc.read_h5ad(data_dir / "clustered.h5ad")
    if "leiden" not in adata.obs.columns:
        status = {"dataset_id": cfg["dataset_id"], "status": "missing_clusters",
                  "reason": "obs 里没有 leiden（step 03 未产出）"}
        write_json(res_dir / "trajectory_status.json", status)
        log_warn(status["reason"])
        return status

    n_clusters = int(adata.obs["leiden"].nunique())
    if n_clusters < 2:
        status = {"dataset_id": cfg["dataset_id"], "status": "not_applicable",
                  "reason": f"只有 {n_clusters} 个簇，轨迹分析无从谈起"}
        write_json(res_dir / "trajectory_status.json", status)
        log_warn(status["reason"])
        return status

    # ---- 1. PAGA ------------------------------------------------------------
    try:
        sc.tl.paga(adata, groups="leiden")
    except Exception as e:  # noqa: BLE001
        status = {"dataset_id": cfg["dataset_id"], "status": "failed",
                  "reason": f"PAGA 失败: {type(e).__name__}: {e}"}
        write_json(res_dir / "trajectory_status.json", status)
        log_warn(status["reason"])
        return status

    conn = adata.uns["paga"]["connectivities"]
    conn = np.asarray(conn.todense()) if hasattr(conn, "todense") else np.asarray(conn)
    cats = [str(c) for c in adata.obs["leiden"].cat.categories]
    conn_df = pd.DataFrame(conn, index=cats, columns=cats)
    conn_df.to_csv(res_dir / "paga_connectivities.csv")

    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    sc.pl.paga(adata, show=False, ax=ax)
    ax.set_title("PAGA graph (edge width = connectivity)", fontsize=10)
    save_fig(cfg, "paga_graph", fig)

    # ---- 2. 选根 ------------------------------------------------------------
    root = traj.get("root_cluster")
    root_record = {"configured": root}
    if root is not None:
        root = str(root)
        if root not in cats:
            status = {"dataset_id": cfg["dataset_id"], "status": "bad_root",
                      "reason": f"trajectory.root_cluster='{root}' 不在簇列表里: {cats}"}
            write_json(res_dir / "trajectory_status.json", status)
            log_warn(status["reason"])
            return status
        root_record["method"] = "configured"
        root_record["reason"] = "来自配置 trajectory.root_cluster"
    else:
        # 自动选根：用"连通度最高的簇"作起点。
        # **这是一个启发式，不是生物学判断。** 连通度高只说明它在图上居中，
        # 而居中既可能是祖细胞，也可能是被各种中间态包围的终末态。
        deg = conn_df.sum(axis=1).sort_values(ascending=False)
        root = str(deg.index[0])
        root_record.update({
            "method": "auto_max_connectivity",
            "reason": ("未配置 root_cluster，取 PAGA 连通度最高的簇。"
                       "**这是启发式**：连通度高只说明它在图上居中，"
                       "居中既可能是祖细胞也可能是终末态。要下方向性结论必须人工定根"),
            "connectivity_ranking": {str(k): round(float(v), 4) for k, v in deg.items()},
        })
        log_warn(f"根为自动选取（簇 {root}）—— 拟时序方向依赖这个选择，见 trajectory_status.json")

    adata.uns["iroot"] = int(np.flatnonzero(adata.obs["leiden"].astype(str).values == root)[0])

    # ---- 3. DPT -------------------------------------------------------------
    sc.tl.diffmap(adata, random_state=cfg["analysis"]["seed"])
    sc.tl.dpt(adata)
    dpt = adata.obs["dpt_pseudotime"].astype(float)
    log_info(f"DPT 完成（根=簇 {root}），拟时序范围 {dpt.min():.3f}-{dpt.max():.3f}")

    # UMAP 上色拟时序
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2))
    xy = adata.obsm["X_umap"]
    s0 = axes[0].scatter(xy[:, 0], xy[:, 1], c=dpt.values, s=4, cmap="viridis")
    axes[0].set_title(f"DPT pseudotime (root = cluster {root})", fontsize=9)
    fig.colorbar(s0, ax=axes[0], label="pseudotime")
    s1 = axes[1].scatter(xy[:, 0], xy[:, 1],
                         c=adata.obs["leiden"].astype(str).astype("category").cat.codes,
                         s=4, cmap="tab20")
    axes[1].set_title("Leiden clusters", fontsize=9)
    for ax in axes:
        ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    save_fig(cfg, "pseudotime_umap", fig)

    # 每簇的拟时序分布
    per_cluster = (adata.obs.groupby("leiden", observed=True)["dpt_pseudotime"]
                   .agg(["mean", "median", "std", "min", "max", "count"])
                   .reset_index().rename(columns={"leiden": "cluster"}))
    per_cluster.to_csv(res_dir / "pseudotime_by_cluster.csv", index=False)

    fig, ax = plt.subplots(figsize=(max(5.5, 0.5 * n_clusters + 2), 3.8))
    order = per_cluster.sort_values("median")["cluster"].astype(str).tolist()
    data = [dpt[adata.obs["leiden"].astype(str) == c].values for c in order]
    # matplotlib 3.9 起 boxplot 的 `labels` 改名为 `tick_labels`；老名字在
    # 3.11 直接 TypeError。用 tick_labels 并给旧版本留回退。
    try:
        ax.boxplot(data, tick_labels=order, showfliers=False)
    except TypeError:
        ax.boxplot(data, labels=order, showfliers=False)
    ax.set_xlabel("cluster (ordered by median pseudotime)")
    ax.set_ylabel("DPT pseudotime")
    ax.set_title("Pseudotime distribution per cluster", fontsize=10)
    save_fig(cfg, "pseudotime_by_cluster", fig)

    status = {
        "dataset_id": cfg["dataset_id"],
        "status": "ok",
        "n_clusters": n_clusters,
        "root_cluster": root,
        "root_selection": root_record,
        "dpt_range": [round(float(dpt.min()), 4), round(float(dpt.max()), 4)],
        "pseudotime_by_cluster": df_to_records(per_cluster),
        # **方法学限定必须落盘。** 这些不是免责声明，是结论的适用范围。
        "limitations": [
            "PAGA 给出的是簇间连通性，不是分化方向",
            "DPT 方向完全依赖根的选择；根选错则轨迹反向",
            "拟时序是一维坐标，分支过程会被压成先后关系",
            "**没有 RNA 速率（spliced/unspliced）时不能下方向性结论** —— "
            "本流水线只有计数矩阵，没有速率，所以这里报的是相似度排序而非分化方向",
        ],
    }
    write_json(res_dir / "trajectory_status.json", status)
    return status


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    t0 = time.time()
    try:
        run_05_trajectory(cfg)
        record_step(cfg, "trajectory", "ok", time.time() - t0)
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "trajectory", "failed", time.time() - t0, message=str(e))
        raise
