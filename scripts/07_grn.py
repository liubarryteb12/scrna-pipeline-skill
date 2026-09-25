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
                    log_warn, parse_args, place_labels, record_step, save_fig,
                    set_seed, write_json, W_DOUBLE, W_ONE_HALF, W_SINGLE, mm, PAL,)

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
    rows, activities, edges = [], {}, []
    cluster = adata.obs["leiden"].astype(str).values
    clusters = sorted(set(cluster), key=lambda x: int(x) if x.isdigit() else x)

    for tf in tfs:
        j = pos[tf]
        # 与所有候选基因的相关（已标准化，内积 / n）
        corr = (Xc.T @ Xc[:, j]) / Xc.shape[0]
        corr[j] = -np.inf          # 排除自己
        top = np.argsort(corr)[::-1][:N_TARGETS]
        # **靶基因和权重必须成对取。** 分开写（一个列表取名字、另一个按
        # 位置取 corr）在过滤 `corr > 0` 之后会错位 —— 而错位不报错，
        # 只会给每个靶基因配上一个别人的相关系数。
        pairs = [(cand_names[k], float(corr[k])) for k in top if corr[k] > 0]
        targets = [p[0] for p in pairs]
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
        # 完整边表（带权重）。tf_regulons.csv 只存了 top 12 的字符串，
        # 下游的虚拟扰动需要**全部**靶基因及其相关系数。
        for tgt, w in pairs:
            edges.append({"tf": tf, "target": tgt, "corr": round(w, 4)})
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

    # 完整边表：下游虚拟扰动（§1.7/§1.8）要用全部靶基因与权重，
    # 而 tf_regulons.csv 里的 top_targets 只有前 12 个、且没有权重。
    edge_df = pd.DataFrame(edges)
    edge_df.to_csv(res_dir / "tf_regulon_edges.csv", index=False)
    log_info(f"调控子边表: {len(edge_df)} 条边，"
             f"{edge_df['tf'].nunique()} 个 TF，{edge_df['target'].nunique()} 个靶基因")

    # 活性矩阵（簇 x TF）
    act_mat = pd.DataFrame(
        {tf: [float(activities[tf][cluster == c].mean()) for c in clusters]
         for tf in activities},
        index=[f"cluster_{c}" for c in clusters])
    act_mat.to_csv(res_dir / "tf_activity_by_cluster.csv")

    # ---- 2b. regulon 活性 × 拟时序 ------------------------------------------
    #
    # 把调控子活性投到 step 05 的共识拟时序上，找"沿轨迹动态变化"的调控子。
    # 这一步依赖 05 的产物，所以放在 07（顺序：05 在 07 之前）。
    #
    # **p 值必须校正。** 200 多个调控子同时检验，不做 BH 的话按 alpha=0.05
    # 会有 10 个左右纯靠运气"显著"。
    pt_df = None
    pt_path = res_dir / "pseudotime_per_cell.csv"
    if pt_path.exists():
        try:
            pt_df = pd.read_csv(pt_path)
        except Exception as e:  # noqa: BLE001
            log_warn(f"读取 {pt_path.name} 失败: {type(e).__name__}: {e}")

    traj_rows, traj_status = [], {"status": "not_available",
                                  "reason": f"{pt_path.name} 不存在（step 05 未产出）"}
    if pt_df is not None and "consensus_pseudotime" in pt_df.columns:
        try:
            from scipy.stats import spearmanr

            pt = pt_df["consensus_pseudotime"].astype(float).values
            # 细胞顺序必须与 activities 对齐：两者都来自同一个 clustered.h5ad，
            # 但保险起见按 cell 名重排，而不是假设顺序一致。
            if "cell" in pt_df.columns and len(pt_df) == adata.n_obs:
                idx = pd.Index(pt_df["cell"].astype(str)).get_indexer(
                    adata.obs_names.astype(str))
                if (idx >= 0).all():
                    pt = pt[idx]
                else:
                    log_warn(f"pseudotime_per_cell.csv 有 {(idx < 0).sum()} 个细胞"
                             "在 adata 里找不到，按原顺序使用（**可能错位**）")
            if len(pt) != adata.n_obs:
                raise ValueError(f"拟时序长度 {len(pt)} != 细胞数 {adata.n_obs}")

            for tf in activities:
                r, p = spearmanr(activities[tf], pt)
                if not np.isfinite(r):
                    continue
                traj_rows.append({"tf": tf,
                                  "rho_with_pseudotime": round(float(r), 4),
                                  "p_value": float(p),
                                  "direction": "increases" if r > 0 else "decreases"})
            if traj_rows:
                tdf = pd.DataFrame(traj_rows)
                # Benjamini-Hochberg
                p = tdf["p_value"].values
                order = np.argsort(p)
                n = len(p)
                q = np.empty(n)
                prev = 1.0
                for rank, i in enumerate(order[::-1]):
                    k = n - rank
                    val = min(prev, p[i] * n / k)
                    q[i] = val
                    prev = val
                tdf["q_value_BH"] = q
                tdf = tdf.sort_values("q_value_BH")
                tdf.to_csv(res_dir / "tf_activity_vs_pseudotime.csv", index=False)

                n_sig = int((tdf["q_value_BH"] < 0.05).sum())
                traj_status = {
                    "status": "ok",
                    "n_regulons_tested": int(len(tdf)),
                    "n_significant_BH05": n_sig,
                    "top": df_to_records(tdf.head(15)),
                    "method": "调控子活性（AUCell 式）vs 共识拟时序的 Spearman 相关，BH 校正",
                    "limitations": [
                        "相关不等于沿轨迹的因果驱动；调控子活性本身是共表达推断的产物",
                        "拟时序方向若整体反掉，这里的 rho 符号会全部反过来"
                        "（见 trajectory_status.json 的 direction_source）",
                        "拟时序是一维坐标，分支上的反向变化会被压掉",
                        # **这条必须写：n 大时显著性很廉价。**
                        f"**{n_sig}/{len(tdf)} 个调控子 BH<0.05，但这个数字本身信息量很低** ——"
                        f"细胞数 n={adata.n_obs}，|rho| 只要约 0.06 就能过 BH<0.05。"
                        "要看的是效应量：本数据 |rho| 中位数 "
                        f"{np.median(np.abs(tdf['rho_with_pseudotime'].values)):.3f}，"
                        f"|rho|>0.3 的 {int((tdf['rho_with_pseudotime'].abs() > 0.3).sum())} 个，"
                        f">0.5 的 {int((tdf['rho_with_pseudotime'].abs() > 0.5).sum())} 个。"
                        "**报显著个数而不报效应量分布，等于把弱关联说成发现**",
                    ],
                }
                log_info(f"regulon×拟时序：{len(tdf)} 个调控子，"
                         f"BH<0.05 的 {n_sig} 个")

                # 图：前 20 个显著调控子的活性沿拟时序分箱
                show = tdf.head(20)["tf"].tolist()
                if show:
                    nb = 20
                    bins = pd.qcut(pt, nb, labels=False, duplicates="drop")
                    nb = int(np.nanmax(bins)) + 1
                    mat = np.zeros((len(show), nb), dtype=float)
                    for i, tf in enumerate(show):
                        v = activities[tf]
                        for b in range(nb):
                            m = bins == b
                            mat[i, b] = float(np.mean(v[m])) if m.sum() else np.nan
                    mu = np.nanmean(mat, axis=1, keepdims=True)
                    sd = np.nanstd(mat, axis=1, keepdims=True)
                    sd[sd == 0] = 1.0
                    matz = np.nan_to_num((mat - mu) / sd)

                    # **`figsize` 是两个元素的元组，第三个元素是历史事故。**
                    # 原写作 `figsize=(W_ONE_HALF, W_SINGLE, max(mm(56), 0.26*len(show)+1.6))`
                    # —— 三元素元组会被 matplotlib 当成 `(nrows, ncols, figsize)` 的
                    # 参数串位，抛 `ValueError: Invalid unit 6.8 in 'figsize'`。
                    # 该行在 `try:` 内、被下面的 `except` 吞成
                    # `grn_status.regulon_vs_pseudotime.status = failed`，
                    # 于是热图与紧随其后的 4 张单 TF 面板**从未产出过**
                    # （commit a3a86fb 引入，早于 K-01b），而顶层 `status` 仍是
                    # `ok`、`acceptance.json` 全绿 —— 图没了却没人发现。
                    # 第三个值随 `len(show)`（行数）变化，是**高度**，所以
                    # 宽度取 W_ONE_HALF、高度取它。
                    fig, ax = plt.subplots(
                        figsize=(W_ONE_HALF, max(mm(56), 0.26 * len(show) + 1.6)))
                    im = ax.imshow(matz, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
                    ax.set_yticks(range(len(show)))
                    ax.set_yticklabels(show, fontsize=7)
                    ax.set_xticks(range(nb))
                    ax.set_xticklabels([str(b) for b in range(nb)], fontsize=7)
                    ax.set_xlabel("consensus pseudotime bin (higher = later)")
                    ax.set_title("Regulon activity along pseudotime (z-scored)")
                    fig.colorbar(im, ax=ax, label="z-scored activity")
                    save_fig(cfg, "02-07-01-unit1-tf-activity-vs-pseudotime", fig)
                    DYNAMIC_FIG_BASES = {"01": 4}
                    TF_FIG_BASE = "-".join(["02", "07", "01", "unit"])
                    # **单 TF 置信带面板**（差距清单 #18，文献范式）：
                    # 热图看得到"哪些 TF 相关"，看不到单个 TF 的趋势形状
                    # 与不确定性。每个 top TF 一张小面板：均值线 + ±1SD ribbon。
                    xs = np.arange(nb)
                    for ui, tf in enumerate(show[:4], start=1):
                        v = activities[tf]
                        mu_b = np.array([float(np.nanmean(v[bins == b])) if np.sum(bins == b) else np.nan
                                         for b in range(nb)])
                        sd_b = np.array([float(np.nanstd(v[bins == b])) if np.sum(bins == b) else np.nan
                                         for b in range(nb)])
                        fig_t, ax_t = plt.subplots(figsize=(W_SINGLE, mm(46)))
                        ax_t.plot(xs, mu_b, "-o", ms=2.5, lw=1.0, color=PAL["primary"])
                        ax_t.fill_between(xs, mu_b - sd_b, mu_b + sd_b,
                                          color=PAL["primary"], alpha=0.18,
                                          label="±1 SD")
                        ax_t.set_xlabel("pseudotime bin (higher = later)")
                        ax_t.set_ylabel("regulon activity")
                        ax_t.set_title(f"{tf} activity along pseudotime (mean ± 1 SD per bin)")
                        ax_t.legend(fontsize=6)
                        save_fig(cfg, TF_FIG_BASE + str(ui + 1) + "-tf-" + tf.lower() + "-trend", fig_t)
        except Exception as e:  # noqa: BLE001
            traj_status = {"status": "failed",
                           "reason": f"{type(e).__name__}: {e}"}
            log_warn(f"regulon×拟时序失败: {traj_status['reason']}")
    else:
        log_info(f"regulon×拟时序：{traj_status['reason']}")

    # ---- 3. 出图 ------------------------------------------------------------
    # **列必须按 TF 去重**（评审 3.8）：`reg.head(n)` 是行级 top，
    # 同一 TF 的多行（不同簇/来源）会让 x 轴出现 TBX21/TBX21 这样的重复列，
    # 评审实测最后 7 列里 6 列无法区分。先按 tf 去重再取前 N。
    top_tfs = (reg.drop_duplicates(subset=["tf"])
                  .head(int(grn.get("top_tfs", 10)) * 2)["tf"].tolist())
    if top_tfs:
        sub = act_mat[top_tfs]
        # 宽度夹在 [单栏半, 双栏]：类别少时不至于太空，类别多时也不会
        # 画出装不进一页的图
        fig, ax = plt.subplots(figsize=(min(W_DOUBLE, max(W_ONE_HALF, W_SINGLE, 0.42 * len(top_tfs) + 2.4)),
                                        max(3.4, 0.34 * len(sub) + 1.8)))
        im = ax.imshow(sub.values, aspect="auto", cmap="RdBu_r",
                       vmin=-np.abs(sub.values).max(), vmax=np.abs(sub.values).max())
        ax.set_xticks(range(len(top_tfs)))
        ax.set_xticklabels(top_tfs, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels(sub.index, fontsize=7)
        ax.set_title("TF regulon activity by cluster (top by specificity)")
        fig.colorbar(im, ax=ax, label="mean AUCell-style score")
        save_fig(cfg, "02-07-02-unit1-tf-activity-heatmap", fig)

    # 特异性 vs 表达细胞比例：识别"只是细胞类型代理"的调控子
    #
    # **标签必须真正防撞。** 原先用"上下交替 + 偶尔左右"的纯参数法 ——
    # 实测失效：top 8 的 TF 的 x 全挤在 0–0.05 的竖列里，交替只把标签分成
    # 两层，**同层内仍然互压**（用户 2026-09-24 反馈"基因标签有重叠"）。
    # 改用 `place_labels()`：贪心选第一个与已放标签不重叠的方位，
    # 重叠判据用**渲染后的真实包围盒**（不是估算字宽），几何上保证不重叠。
    #
    # 同时修颜色：原写死 `#2C7FB8`（规则 13 禁止就地写色值字面量，
    # 且它不在 PAL 里 —— 与 PAL["primary"] #0072B2 色相只差约 3°，
    # 属于"同一图上两个蓝"）。改走 PAL。
    fig, ax = plt.subplots(figsize=(W_ONE_HALF, mm(70)))
    ax.scatter(reg["frac_cells_expressing_tf"], reg["cluster_specificity"],
               s=22, color=PAL["primary"], alpha=0.75)
    _lab = reg.head(8).sort_values("cluster_specificity", ascending=False)
    place_labels(ax,
                 _lab["frac_cells_expressing_tf"].to_numpy(),
                 _lab["cluster_specificity"].to_numpy(),
                 _lab["tf"].tolist(), fontsize=7)
    ax.set_xlabel("fraction of cells expressing TF")
    ax.set_ylabel("cluster specificity")
    ax.set_title("Regulon specificity vs TF detection\n"
                 "labels = top 8 by cluster specificity (auto-placed, non-overlapping)")
    save_fig(cfg, "02-07-03-unit1-tf-specificity-scatter", fig)

    status = {
        "dataset_id": cfg["dataset_id"],
        "status": "ok",
        "n_tfs_in_list": len(load_tfs()),
        "n_tfs_present": len(tfs),
        "n_regulons": int(len(reg)),
        "n_targets_per_regulon": N_TARGETS,
        "n_genes_used_for_inference": len(cand),
        "n_regulon_edges": int(len(edge_df)),
        "n_regulon_target_genes": int(edge_df["target"].nunique()),
        "regulon_edges_file": "tf_regulon_edges.csv",
        "top_regulons": df_to_records(reg.head(15)),
        "regulon_vs_pseudotime": traj_status,
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
