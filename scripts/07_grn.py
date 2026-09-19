#!/usr/bin/env python3
"""
07_grn.py — 转录因子调控网络（共表达推断 + 调控子活性）

**这里做的是共表达推断，不是 SCENIC。** 差别必须说清楚：

  SCENIC 的流程是 GRNBoost2/GENIE3（梯度提升/随机森林）推断 TF-靶基因
  关系，**然后用 cisTarget 做 motif 富集剪枝** —— 只有靶基因启动子区真的
  有该 TF 的结合 motif 才保留。剪枝是 SCENIC 的关键一步，它把
  "共表达" 收紧成 "可能有直接调控"。

  本流水线**没有 motif 剪枝**（需要 cisTarget 的排名数据库，几百 MB）。
  所以这里报的是**共表达模块**，用 TF 命名而已。

后果：一个 TF 和一组基因共表达，可能是它调控它们，也可能是：
  - 它们被同一个上游因子调控
  - 它们只是同一种细胞类型的标志物（**最常见**）
  - 拷贝数变异导致的共表达

所以每个调控子都报 `n_cells_expressing_tf` 与 `cluster_specificity` ——
让"这个调控子是不是只是细胞类型的代理"这件事可以被看出来。

调控子活性用 AUCell 式打分（scanpy 的 score_genes），比直接取靶基因
平均表达稳 —— 它校正了基因的表达水平与测序深度。
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import yaml  # noqa: E402
from scipy import sparse  # noqa: E402

from common import (df_to_records, ensure_dirs, load_config, log_info,  # noqa: E402
                    log_warn, parse_args, record_step, save_fig, set_seed, write_json)

# 每个调控子保留多少个共表达靶基因
N_TARGETS = 30
# 推断用的高变基因数上限（控制耗时）
MAX_GENES_FOR_INFERENCE = 3000


def load_tfs() -> list:
    p = Path(__file__).resolve().parent.parent / "assets" / "tf_list.yml"
    with open(p, encoding="utf-8") as fh:
        doc = yaml.safe_load(fh)
    return doc.get("tfs", []) or []


def get_full_expression(adata):
    """取全基因集表达（同 06 —— 转录因子大多是低表达的，不在 HVG 里）。"""
    if adata.raw is None:
        raise RuntimeError("adata.raw 为空 —— 无法取全基因集表达做 GRN 推断")
    src = adata.raw
    X = src.X
    X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    return X, list(src.var_names)


def run_07_grn(cfg: dict) -> dict:
    ensure_dirs(cfg)
    set_seed(cfg)
    data_dir = Path(cfg["output"]["data_dir"])
    res_dir = Path(cfg["output"]["results_dir"])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    grn = cfg.get("grn") or {}
    if not grn.get("enabled", True):
        status = {"dataset_id": cfg["dataset_id"], "status": "disabled",
                  "reason": "配置 grn.enabled=false"}
        write_json(res_dir / "grn_status.json", status)
        log_info("GRN 分析已按配置关闭")
        return status

    adata = sc.read_h5ad(data_dir / "clustered.h5ad")
    X, var_names = get_full_expression(adata)
    lookup = {g: i for i, g in enumerate(var_names)}
    log_info(f"表达矩阵（全基因集）: {X.shape[0]} 细胞 x {X.shape[1]} 基因")

    tfs = [t for t in load_tfs() if t in lookup]
    if not tfs:
        status = {"dataset_id": cfg["dataset_id"], "status": "no_tfs_in_data",
                  "reason": "TF 列表里的基因一个都不在数据里"}
        write_json(res_dir / "grn_status.json", status)
        log_warn(status["reason"])
        return status
    log_info(f"转录因子: {len(tfs)}/{len(load_tfs())} 在数据里存在")

    # ---- 1. 选推断用的基因集 ------------------------------------------------
    # 用表达方差最高的基因（而不是 HVG 标记），因为需要全基因集上的方差
    var = X.var(axis=0)
    order = np.argsort(var)[::-1]
    tf_idx = {t: lookup[t] for t in tfs}
    # 保证所有 TF 都在候选里
    cand = list(dict.fromkeys(list(tf_idx.values()) +
                              [int(i) for i in order[:MAX_GENES_FOR_INFERENCE]]))
    cand = sorted(set(cand))
    log_info(f"推断用基因: {len(cand)} 个（含全部 {len(tfs)} 个 TF）")

    Xc = X[:, cand]
    # 中心化 + 标准化，相关系数就等于内积
    Xc = Xc - Xc.mean(axis=0, keepdims=True)
    sd = Xc.std(axis=0, keepdims=True)
    sd[sd == 0] = 1.0
    Xc /= sd
    cand_names = [var_names[i] for i in cand]
    pos = {g: i for i, g in enumerate(cand_names)}

    # ---- 2. 推断调控子 ------------------------------------------------------
    rows, activities = [], {}
    cluster = adata.obs["leiden"].astype(str).values
    clusters = sorted(set(cluster), key=lambda x: int(x) if x.isdigit() else x)

    for tf in tfs:
        j = pos[tf]
        # 与所有候选基因的相关（已标准化，内积 / n）
        corr = (Xc.T @ Xc[:, j]) / Xc.shape[0]
        corr[j] = -np.inf          # 排除自己
        top = np.argsort(corr)[::-1][:N_TARGETS]
        targets = [cand_names[k] for k in top if corr[k] > 0]
        if len(targets) < 5:
            continue

        # 调控子活性：AUCell 式打分（scanpy score_genes，校正表达水平与深度）
        #
        # **必须立刻 pop 走。** `score_genes` 是往 adata.obs 里插列；循环里
        # 插 200 多列会让 DataFrame 严重碎片化 —— pandas 刷一屏
        # PerformanceWarning，而且插入越来越慢。pop 掉后 obs 宽度不增长。
        name = f"_act_{tf}"
        sc.tl.score_genes(adata, targets, score_name=name,
                          random_state=42, use_raw=True)
        act = adata.obs.pop(name).astype(float).values
        activities[tf] = act

        per_cluster = {c: float(act[cluster == c].mean()) for c in clusters}
        best = max(per_cluster, key=per_cluster.get)
        vals = np.array(list(per_cluster.values()))
        # 簇特异性：最高簇与其余簇均值的差，除以整体标准差。
        # **这是用来识别"这个调控子是不是只是细胞类型的代理"的。**
        spec = float((vals.max() - np.median(vals)) / (vals.std() + 1e-9))
        rows.append({
            "tf": tf,
            "n_targets": len(targets),
            "top_targets": ",".join(targets[:12]),
            "mean_corr": round(float(corr[top[:len(targets)]].mean()), 4),
            "best_cluster": best,
            "best_cluster_activity": round(per_cluster[best], 4),
            "cluster_specificity": round(spec, 3),
            "n_cells_expressing_tf": int((X[:, lookup[tf]] > 0).sum()),
            "frac_cells_expressing_tf": round(float((X[:, lookup[tf]] > 0).mean()), 4),
        })

    if not rows:
        status = {"dataset_id": cfg["dataset_id"], "status": "no_regulons",
                  "reason": f"{len(tfs)} 个 TF 里没有一个找到 >=5 个正相关靶基因"}
        write_json(res_dir / "grn_status.json", status)
        log_warn(status["reason"])
        return status

    reg = pd.DataFrame(rows).sort_values("cluster_specificity", ascending=False)
    reg.to_csv(res_dir / "tf_regulons.csv", index=False)
    log_info(f"推断出 {len(reg)} 个调控子；"
             f"最高簇特异性 {reg['cluster_specificity'].iloc[0]:.2f}（{reg['tf'].iloc[0]}）")

    # 活性矩阵（簇 x TF）
    act_mat = pd.DataFrame(
        {tf: [float(activities[tf][cluster == c].mean()) for c in clusters]
         for tf in activities},
        index=[f"cluster_{c}" for c in clusters])
    act_mat.to_csv(res_dir / "tf_activity_by_cluster.csv")

    # ---- 3. 出图 ------------------------------------------------------------
    top_tfs = reg.head(int(grn.get("top_tfs", 10)) * 2)["tf"].tolist()
    if top_tfs:
        sub = act_mat[top_tfs]
        fig, ax = plt.subplots(figsize=(max(6.0, 0.42 * len(top_tfs) + 2.4),
                                        max(3.4, 0.34 * len(sub) + 1.8)))
        im = ax.imshow(sub.values, aspect="auto", cmap="RdBu_r",
                       vmin=-np.abs(sub.values).max(), vmax=np.abs(sub.values).max())
        ax.set_xticks(range(len(top_tfs)))
        ax.set_xticklabels(top_tfs, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels(sub.index, fontsize=7)
        ax.set_title("TF regulon activity by cluster (top by specificity)", fontsize=10)
        fig.colorbar(im, ax=ax, label="mean AUCell-style score")
        save_fig(cfg, "tf_activity_heatmap", fig)

    # 特异性 vs 表达细胞比例：识别"只是细胞类型代理"的调控子
    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    ax.scatter(reg["frac_cells_expressing_tf"], reg["cluster_specificity"],
               s=22, color="#2C7FB8", alpha=0.75)
    for r in reg.head(8).itertuples():
        ax.annotate(r.tf, (r.frac_cells_expressing_tf, r.cluster_specificity),
                    fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("fraction of cells expressing TF")
    ax.set_ylabel("cluster specificity")
    ax.set_title("Regulon specificity vs TF detection", fontsize=10)
    save_fig(cfg, "tf_specificity_scatter", fig)

    status = {
        "dataset_id": cfg["dataset_id"],
        "status": "ok",
        "n_tfs_in_list": len(load_tfs()),
        "n_tfs_present": len(tfs),
        "n_regulons": int(len(reg)),
        "n_targets_per_regulon": N_TARGETS,
        "n_genes_used_for_inference": len(cand),
        "top_regulons": df_to_records(reg.head(15)),
        "method": ("共表达推断（Pearson 相关取 top 靶基因）+ AUCell 式调控子活性打分。"
                   "**不是 SCENIC**"),
        "limitations": [
            "**没有 motif 剪枝**：SCENIC 用 cisTarget 做 motif 富集把共表达收紧成"
            "可能有直接调控；本流水线没有这一步，所以报的是共表达模块",
            "共表达不等于调控：可能是共同上游、或只是同一细胞类型的标志物",
            "cluster_specificity 高且 frac_cells_expressing_tf 高的调控子，"
            "很可能是细胞类型的代理而非真实调控程序",
            "靶基因取自表达方差最高的基因，低表达靶基因会被漏掉",
        ],
    }
    write_json(res_dir / "grn_status.json", status)
    return status


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    t0 = time.time()
    try:
        run_07_grn(cfg)
        record_step(cfg, "grn", "ok", time.time() - t0)
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "grn", "failed", time.time() - t0, message=str(e))
        raise
