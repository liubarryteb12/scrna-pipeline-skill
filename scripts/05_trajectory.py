#!/usr/bin/env python3
"""
05_trajectory.py — 轨迹推断（PAGA + 四种拟时序方法交叉验证）

**轨迹分析的结论比聚类弱得多，这一点必须写进产物。**

为什么至少两种方法：单一方法的拟时序是一个**一维坐标**，把高维状态
压成一条线；不同算法压出来的线可以完全不同。实测（PBMC3k）三种方法
的原始拟时序两两 Spearman 相关从 **-0.47 到 +0.25** —— 符号都不一样。
不做交叉验证就报"轨迹"，报的是某个算法的一次输出，不是数据里的结构。

为什么方向必须显式校正：拟时序的**符号是任意的**。DPT 从根出发，
根选在分化早期还是晚期，整条轴就反过来。所以本脚本：
  1. 每个方法各自算出拟时序（原始符号，不假设谁对）
  2. 用一个**方向参考**判断每个方法是否需要翻转
  3. 统一约定成「**值越大越晚**」后再互相比较

方向参考的来源（优先级）：
  a. 配置 `trajectory.early_markers` / `late_markers` —— 生物学判据
  b. 都没有时退回 CytoTRACE 分化潜能分（**记进 direction_source**）

**没有 RNA 速率（spliced/unspliced）时拟时序不能下方向性结论。**
本流水线只有计数矩阵，scVelo 因此记 not_done。见 status 的 limitations。
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
from scipy.stats import ks_2samp, spearmanr  # noqa: E402

from common import (df_to_records, ensure_dirs, load_config, log_info,  # noqa: E402
                    log_warn, parse_args, record_step, save_fig, set_seed, write_json,
                    W_DOUBLE, W_ONE_HALF, W_SINGLE, mm, PAL, PAL_CYCLE,)

# 校正后的统一方向：**值越大越晚**
N_MODULES = 6
MIN_GENES_FOR_MODULES = 200


# ============================================================================
# 方法 1：CytoTRACE 风格的 GCS（Gene Counts Signature）
# ============================================================================
def compute_cytotrace_gcs(adata, log=log_info) -> np.ndarray:
    """CytoTRACE 的核心统计量：转录组复杂度 → 分化潜能。

    **这是自行实现，不是 CytoTRACE/CytoTRACE2 包。** 本仓库的惯例是
    第三方包不可靠或不可得时自己实现统计量并写明（Part 3 的 Moran's I、
    Part 1 的 Harrell C 同理）。CytoTRACE2 不在 PyPI 上。

    算法（原论文的核心）：
      1. 每个细胞表达基因数 = 转录组复杂度
      2. 每个基因与"表达基因数"的 Spearman 相关
      3. 取相关性最高的一批基因（默认 1%，至少 50 个）
      4. 细胞得分 = 这批基因表达秩的均值，归一化到 [0,1]

    **约定：GCS 高 = 分化潜能高 = 早。** 与 CytoTRACE 原论文一致。
    """
    src = adata.raw.to_adata() if adata.raw is not None else adata
    X = src.X
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    X = X.astype(np.float32)

    n_genes_per_cell = (X > 0).sum(axis=1).astype(float)
    if n_genes_per_cell.std() == 0:
        raise ValueError("所有细胞的表达基因数相同，GCS 无从计算")

    # 逐基因与复杂度的相关。numba 不可用时这个循环是主要耗时（约 10 秒/万基因）。
    corr = np.zeros(X.shape[1], dtype=float)
    for j in range(X.shape[1]):
        col = X[:, j]
        if col.std() == 0:
            continue
        r = spearmanr(col, n_genes_per_cell).correlation
        corr[j] = 0.0 if not np.isfinite(r) else r

    top_n = max(50, int(0.01 * X.shape[1]))
    top_n = min(top_n, X.shape[1])
    top_idx = np.argsort(-corr)[:top_n]

    sub = X[:, top_idx]
    ranks = np.apply_along_axis(lambda c: pd.Series(c).rank().values, 0, sub)
    gcs = ranks.mean(axis=1)
    rng = gcs.max() - gcs.min()
    gcs = (gcs - gcs.min()) / (rng if rng > 0 else 1.0)
    log(f"CytoTRACE GCS（自实现）: top_genes={top_n}，范围 {gcs.min():.3f}-{gcs.max():.3f}")
    return gcs


# ============================================================================
# 方法 2：DPT（scanpy）
# ============================================================================
def compute_dpt(adata, root_idx: int, seed: int):
    """扩散拟时序。约定：DPT 0 = 根 = 早，所以**越小越早**。"""
    sc.tl.diffmap(adata, random_state=seed)
    adata.uns["iroot"] = int(root_idx)
    sc.tl.dpt(adata)
    return adata.obs["dpt_pseudotime"].astype(float).values


# ============================================================================
# 方法 3：Palantir（概率性命运建模）
# ============================================================================
def compute_palantir(adata, root_idx: int, seed: int, log=log_info):
    """Palantir：扩散图 + Markov 链，给出拟时序与**命运概率**。

    返回 (pseudotime, fate_probabilities_or_None)。约定：0 = 起始细胞 = 早。
    """
    import palantir

    p = sc.AnnData(X=adata.X.copy(), obs=adata.obs.copy(), var=adata.var.copy())
    for k in ("X_pca", "X_umap"):
        if k in adata.obsm:
            p.obsm[k] = adata.obsm[k].copy()

    palantir.utils.run_diffusion_maps(p, n_components=10, knn=30, seed=seed)
    palantir.utils.determine_multiscale_space(p)
    early = str(p.obs_names[int(root_idx)])
    palantir.core.run_palantir(p, early, num_waypoints=500, seed=seed)
    pt = p.obs["palantir_pseudotime"].astype(float).values

    fate = None
    if "palantir_fate_probabilities" in p.obsm:
        f = p.obsm["palantir_fate_probabilities"]
        fate = np.asarray(f.todense()) if hasattr(f, "todense") else np.asarray(f)
    log(f"Palantir: 起始细胞 {early}，范围 {pt.min():.3f}-{pt.max():.3f}")
    return pt, fate


# ============================================================================
# 方法 4：scFates（Slingshot 的 Python 移植，主曲线树）
# ============================================================================
def compute_scfates(adata, gcs: np.ndarray, seed: int, log=log_info):
    """scFates：主曲线树 + 沿树距离作为拟时序。

    **官方流程缺一不可：diffusion → tree → cleanup → root → pseudotime。**
    实测踩过：漏掉 `cleanup()` 时 `pseudotime()` 必定崩在
    `pd.Series(uns['graph']['milestones']) == t][0]`，IndexError ——
    因为 map_cells 读 milestones 时它还没被写入。

    `root()` 支持传 obs 里的列名做**自动选根**（文档原话：
    "a key (str) from obs/X (such as CytoTRACE) for automatic selection"），
    所以这里把 GCS 写进 obs 再传列名。

    返回 (pseudotime, segments, milestones, tips) 或抛异常由调用方处理。
    """
    import scFates as scf
    from scFates import tl as sct

    f = sc.AnnData(X=adata.X.copy(), obs=adata.obs.copy(), var=adata.var.copy())
    for k in ("X_pca", "X_umap"):
        if k in adata.obsm:
            f.obsm[k] = adata.obsm[k].copy()
    f.obs["CytoTRACE"] = gcs

    scf.pp.diffusion(f, n_components=10, knn=30, alpha=0, multiscale=True)
    sct.tree(f, use_rep="X_diffusion_multiscale", ndims_rep=5,
             method="ppt", Nodes=50, seed=seed)
    sct.cleanup(f)                      # **必需**，见 docstring
    sct.root(f, "CytoTRACE")
    sct.pseudotime(f, n_jobs=1, seed=seed)

    pt = f.obs["t"].astype(float).values
    seg = f.obs["seg"].astype(str).values
    mil = f.obs["milestones"].astype(str).values
    tips = np.asarray(f.uns["graph"]["tips"]).ravel().astype(int).tolist()
    log(f"scFates: 范围 {pt.min():.2f}-{pt.max():.2f}，"
        f"{len(set(seg))} 个片段，{len(set(mil))} 个 milestone")
    return pt, seg, mil, tips


# ============================================================================
# 方向校正
# ============================================================================
def marker_direction_score(adata, early, late, log=log_info):
    """用 marker 基因算一个"越晚值越大"的参考轴：mean(late) - mean(early)。

    返回 (score, 实际用到的基因) 或 (None, {})。
    """
    src = adata.raw.to_adata() if adata.raw is not None else adata
    present = set(src.var_names)
    e = [g for g in early if g in present]
    l = [g for g in late if g in present]
    if not e or not l:
        log(f"marker 方向参考不可用：early 命中 {e}，late 命中 {l}")
        return None, {"early": e, "late": l}

    X = src[:, e + l].X
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    X = X.astype(np.float32)
    score = X[:, len(e):].mean(axis=1) - X[:, :len(e)].mean(axis=1)
    log(f"marker 方向参考：early={e} late={l}")
    return score, {"early": e, "late": l}


def orient(raw: dict, reference: np.ndarray, convention: dict):
    """把每个方法的拟时序统一成「值越大越晚」。

    `convention[name]` 说明该方法**原始**输出里"大"代表什么：
      "later"  —— 大就是晚，直接用
      "earlier"— 大是早，取负
    然后与 reference 求相关；仍为负就翻转，并记录翻转与否。
    """
    out, rows = {}, []
    for name, v in raw.items():
        vv = np.asarray(v, dtype=float).copy()
        if convention.get(name) == "earlier":
            vv = -vv
        r_before = spearmanr(vv, reference).correlation
        flipped = False
        if r_before < 0:
            vv = -vv
            flipped = True
        r_after = spearmanr(vv, reference).correlation
        out[name] = vv
        rows.append({
            "method": name,
            "raw_convention": convention.get(name, "?"),
            "rho_vs_reference_before_flip": round(float(r_before), 4),
            "flipped": bool(flipped),
            "rho_vs_reference_after_flip": round(float(r_after), 4),
        })
    return out, rows


# ============================================================================
# 主流程
# ============================================================================
def run_05_trajectory(cfg: dict) -> dict:
    ensure_dirs(cfg)
    set_seed(cfg)
    data_dir = Path(cfg["output"]["data_dir"])
    res_dir = Path(cfg["output"]["results_dir"])
    seed = cfg["analysis"]["seed"]

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

    # ---- 1. PAGA 连通图 -----------------------------------------------------
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

    fig, ax = plt.subplots(figsize=(W_ONE_HALF, mm(84)))
    # **节点必须画在边之上**（评审 3.3）：默认把粗黑边画在节点上层，
    # 小节点（簇 9）被 5 条边直接切断 "9" 字形。分两层：先无标注骨架画边，
    # 再同一坐标轴上高 zorder 重画节点与标签。
    sc.pl.paga(adata, show=False, ax=ax, labels=False, node_size_scale=0.6)
    sc.pl.paga(adata, show=False, ax=ax, edges=False,
               labels=[str(c) for c in adata.obs["leiden"].cat.categories])
    ax.set_title("PAGA graph (edge width = connectivity)")
    save_fig(cfg, "02-05-01-unit1-paga-graph", fig)

    # ---- 2. CytoTRACE GCS：既是方法也是选根依据 ------------------------------
    gcs = compute_cytotrace_gcs(adata)
    root_idx = int(np.argmax(gcs))          # GCS 最高 = 分化潜能最高 = 早
    root_cluster = str(adata.obs["leiden"].astype(str).values[root_idx])
    root_record = {
        "method": "auto_max_cytotrace_gcs",
        "root_cluster": root_cluster,
        "root_cell": str(adata.obs_names[root_idx]),
        "reason": ("取 CytoTRACE GCS 最高的细胞作根（GCS 高 = 分化潜能高 = 早）。"
                   "**这是统计判据不是生物学判据** —— 要下方向性结论仍需人工"
                   "用已知早期/晚期 marker 复核。"),
    }
    configured_root = traj.get("root_cluster")
    if configured_root is not None:
        cr = str(configured_root)
        if cr not in cats:
            status = {"dataset_id": cfg["dataset_id"], "status": "bad_root",
                      "reason": f"trajectory.root_cluster='{cr}' 不在簇列表里: {cats}"}
            write_json(res_dir / "trajectory_status.json", status)
            log_warn(status["reason"])
            return status
        mask = adata.obs["leiden"].astype(str).values == cr
        # 该簇里挑 GCS 最高的细胞，避免任意取第一个
        idx_in = np.flatnonzero(mask)
        root_idx = int(idx_in[np.argmax(gcs[idx_in])])
        root_cluster = cr
        root_record = {
            "method": "configured",
            "root_cluster": cr,
            "root_cell": str(adata.obs_names[root_idx]),
            "reason": "来自配置 trajectory.root_cluster（簇内取 GCS 最高的细胞）",
        }
    log_info(f"根 = 簇 {root_cluster}，细胞 {adata.obs_names[root_idx]}")

    # ---- 3. 四种方法 --------------------------------------------------------
    raw_pt, methods_ok, methods_failed = {}, {}, {}

    try:
        raw_pt["dpt"] = compute_dpt(adata, root_idx, seed)
        methods_ok["dpt"] = "scanpy DPT（扩散拟时序）"
    except Exception as e:  # noqa: BLE001
        methods_failed["dpt"] = f"{type(e).__name__}: {e}"
        log_warn(f"DPT 失败: {methods_failed['dpt']}")

    fate_probs = None
    try:
        raw_pt["palantir"], fate_probs = compute_palantir(adata, root_idx, seed)
        methods_ok["palantir"] = "Palantir（扩散图 + Markov 链命运概率）"
    except Exception as e:  # noqa: BLE001
        methods_failed["palantir"] = f"{type(e).__name__}: {e}"
        log_warn(f"Palantir 失败: {methods_failed['palantir']}")

    seg = mil = None
    sf_tips = []
    try:
        raw_pt["scfates"], seg, mil, sf_tips = compute_scfates(adata, gcs, seed)
        methods_ok["scfates"] = "scFates（主曲线树，Slingshot 的 Python 移植）"
    except Exception as e:  # noqa: BLE001
        methods_failed["scfates"] = f"{type(e).__name__}: {e}"
        log_warn(f"scFates 失败: {methods_failed['scfates']}")

    # CytoTRACE 自己也是一种方法
    raw_pt["cytotrace"] = gcs
    methods_ok["cytotrace"] = "CytoTRACE 风格 GCS（本仓库自行实现，非 CytoTRACE 包）"

    if len(methods_ok) < 2:
        status = {
            "dataset_id": cfg["dataset_id"], "status": "insufficient_methods",
            "reason": (f"只有 {len(methods_ok)} 种方法成功，少于交叉验证所需的 2 种"),
            "methods_ok": methods_ok, "methods_failed": methods_failed,
        }
        write_json(res_dir / "trajectory_status.json", status)
        log_warn(status["reason"])
        return status

    # 原始符号约定：dpt/palantir/scfates 的 0 = 根 = 早（小=早）；
    # cytotrace 的 1 = 分化潜能高 = 早（大=早）。
    convention = {"dpt": "earlier", "palantir": "earlier",
                  "scfates": "earlier", "cytotrace": "later"}

    # ---- 4. 方向参考 --------------------------------------------------------
    early_m = list(traj.get("early_markers") or [])
    late_m = list(traj.get("late_markers") or [])
    marker_ref, marker_used = (None, {})
    if early_m and late_m:
        marker_ref, marker_used = marker_direction_score(adata, early_m, late_m)

    if marker_ref is not None:
        reference = marker_ref
        direction_source = "markers"
        direction_note = ("方向参考 = 配置的 early/late marker 基因"
                          "（mean(late) - mean(early)）")
    else:
        # **退回 CytoTRACE 并写明。** GCS 高 = 早，取负使其"大 = 晚"。
        reference = -gcs
        direction_source = "cytotrace_fallback"
        direction_note = (
            "**未配置 trajectory.early_markers/late_markers**，退回用 CytoTRACE "
            "分化潜能分作方向参考（取负使「大 = 晚」）。这是统计判据不是生物学"
            "判据：若该数据集里分化潜能与成熟度不同向，方向会整体反掉。"
            "要下方向性结论请在配置里给出已知的早期/晚期 marker。")
        log_warn(direction_note)

    corrected, direction_rows = orient(raw_pt, reference, convention)
    pd.DataFrame(direction_rows).to_csv(res_dir / "trajectory_direction.csv", index=False)
    for r in direction_rows:
        log_info(f"  方向 {r['method']:10} rho={r['rho_vs_reference_after_flip']:+.4f}"
                 f"{'（已翻转）' if r['flipped'] else ''}")

    # ---- 5. 交叉验证矩阵 ----------------------------------------------------
    names = [n for n in ("dpt", "palantir", "scfates", "cytotrace") if n in corrected]
    cmat = pd.DataFrame(index=names, columns=names, dtype=float)
    for a_ in names:
        for b_ in names:
            cmat.loc[a_, b_] = spearmanr(corrected[a_], corrected[b_]).correlation
    cmat.to_csv(res_dir / "trajectory_method_correlation.csv")

    # **方向参考不能算进"交叉验证一致性"。**
    # 退回模式下参考就是 CytoTRACE 本身，它与自己的相关恒为 ±1 ——
    # 那是定义，不是证据。把它算进去会把一致性整体抬高。
    reference_method = "cytotrace" if direction_source == "cytotrace_fallback" else None
    cv_names = [n for n in names if n != reference_method]
    off = [float(cmat.loc[a_, b_]) for i, a_ in enumerate(cv_names)
           for b_ in cv_names[i + 1:]]
    mean_rho = float(np.mean(off)) if off else float("nan")
    min_rho = float(np.min(off)) if off else float("nan")
    if reference_method:
        log_info(f"一致性统计已排除方向参考 {reference_method}"
                 f"（它与参考的相关是定义上的，不是证据）")

    fig, ax = plt.subplots(figsize=(W_SINGLE, mm(66)))
    im = ax.imshow(cmat.values.astype(float), cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(names))); ax.set_xticklabels(names, rotation=40, ha="right")
    ax.set_yticks(range(len(names))); ax.set_yticklabels(names)
    for i in range(len(names)):
        for j in range(len(names)):
            ax.text(j, i, f"{cmat.values[i, j]:.2f}", ha="center", va="center",
                    fontsize=8, color="black")
    ax.set_title("Pseudotime agreement (Spearman, direction-corrected)")
    fig.colorbar(im, ax=ax, shrink=0.8)
    save_fig(cfg, "02-05-02-unit1-trajectory-method-correlation", fig)

    # 共识拟时序：各方法 z-score 后取均值（方向已统一）
    Z = np.column_stack([
        (corrected[n] - np.nanmean(corrected[n])) / (np.nanstd(corrected[n]) + 1e-12)
        for n in names
    ])
    consensus = np.nanmean(Z, axis=1)
    log_info(f"共识拟时序完成；方法间平均 rho={mean_rho:+.4f}，最低 {min_rho:+.4f}")

    # ---- 6. 沿轨迹变化的基因 ------------------------------------------------
    genes_rows, gene_mat, gene_names, gene_rho = [], None, [], None
    try:
        src = adata.raw.to_adata() if adata.raw is not None else adata
        X = src.X
        X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
        X = X.astype(np.float32)
        # 只留有表达的基因，否则一堆 0 会拿到无意义的 rho
        expressed = (X > 0).sum(axis=0) >= max(10, int(0.01 * X.shape[0]))
        idx = np.flatnonzero(expressed)
        Xs = X[:, idx]
        # 与共识拟时序的 Spearman
        rho = np.array([spearmanr(Xs[:, k], consensus).correlation
                        for k in range(Xs.shape[1])])
        rho = np.nan_to_num(rho)
        names_g = [src.var_names[i] for i in idx]
        order = np.argsort(-np.abs(rho))
        top_n = min(400, len(order))
        sel = order[:top_n]
        gene_names = [names_g[i] for i in sel]
        gene_mat = Xs[:, sel]
        gene_rho = rho[sel]          # 与 gene_names 一一对应，供模块用
        genes_rows = [{"gene": names_g[i], "rho_with_pseudotime": round(float(rho[i]), 4),
                       "direction": "increases" if rho[i] > 0 else "decreases"}
                      for i in order[:200]]
        pd.DataFrame(genes_rows).to_csv(res_dir / "trajectory_genes.csv", index=False)
        log_info(f"沿轨迹变化基因：表达基因 {len(idx)} 个，报前 {len(genes_rows)} 个")
    except Exception as e:  # noqa: BLE001
        log_warn(f"沿轨迹基因分析失败: {type(e).__name__}: {e}")

    modules_rows = []
    if gene_mat is not None and gene_mat.shape[1] >= MIN_GENES_FOR_MODULES:
        try:
            # 按拟时序分箱 → 每个基因的平滑表达曲线 → k-means 聚成模块
            n_bins = 20
            bins = pd.qcut(consensus, n_bins, labels=False, duplicates="drop")
            nb = int(np.nanmax(bins)) + 1
            prof = np.zeros((gene_mat.shape[1], nb), dtype=np.float32)
            for b in range(nb):
                m = bins == b
                if m.sum() == 0:
                    prof[:, b] = np.nan
                    continue
                prof[:, b] = gene_mat[m].mean(axis=0)
            # 行 z-score，让模块反映"形状"而不是"表达量"
            mu = np.nanmean(prof, axis=1, keepdims=True)
            sd = np.nanstd(prof, axis=1, keepdims=True)
            sd[sd == 0] = 1.0
            profz = (prof - mu) / sd
            profz = np.nan_to_num(profz)

            from sklearn.cluster import KMeans

            km = KMeans(n_clusters=N_MODULES, n_init=10, random_state=seed)
            lab = km.fit_predict(profz)
            for m in range(N_MODULES):
                sel_m = np.flatnonzero(lab == m)
                if len(sel_m) == 0:
                    continue
                members = [gene_names[i] for i in sel_m]
                # 模块的平均 |rho|：用与 gene_names 对齐的 gene_rho，
                # 不要拿只存了前 200 个的 genes_rows 去比 —— 那样模块里
                # 大部分基因查不到，均值会变成 nan（踩过 RuntimeWarning）。
                modules_rows.append({
                    "module": m, "n_genes": len(members),
                    "mean_abs_rho": round(float(np.mean(np.abs(gene_rho[sel_m]))), 4),
                    "mean_rho": round(float(np.mean(gene_rho[sel_m])), 4),
                    "top_genes": ",".join(members[:20]),
                })
            pd.DataFrame(modules_rows).to_csv(res_dir / "trajectory_modules.csv",
                                              index=False)
            log_info(f"基因模块：{len(modules_rows)} 个（k-means on 分箱平滑曲线）")

            fig, axes = plt.subplots(1, 2, figsize=(W_DOUBLE, mm(64)))
            im = axes[0].imshow(profz[np.argsort(lab)], aspect="auto",
                                cmap="RdBu_r", vmin=-2, vmax=2)
            axes[0].set_xlabel("pseudotime bin"); axes[0].set_ylabel("gene (grouped by module)")
            axes[0].set_title(f"Genes along pseudotime, {N_MODULES} modules")
            fig.colorbar(im, ax=axes[0], label="z-scored mean expression")
            for m in range(N_MODULES):
                sel_m = np.flatnonzero(lab == m)
                if len(sel_m) == 0:
                    continue
                axes[1].plot(np.nanmean(profz[sel_m], axis=0), label=f"M{m} (n={len(sel_m)})")
            axes[1].set_xlabel("pseudotime bin"); axes[1].set_ylabel("z-scored mean")
            axes[1].set_title("Module profiles")
            axes[1].legend(fontsize=7)
            save_fig(cfg, "02-05-03-unit1-trajectory-modules", fig)
        except Exception as e:  # noqa: BLE001
            log_warn(f"基因模块分析失败: {type(e).__name__}: {e}")

    # ---- 7. 多条件拟时序分布比较（Kolmogorov-Smirnov）-----------------------
    ks_rows = []
    design = cfg.get("design") or {}
    group_key = design.get("group_key")
    if group_key and group_key in adata.obs.columns:
        groups = adata.obs[group_key].astype(str).values
        uniq = sorted(set(groups))
        if len(uniq) >= 2:
            for i in range(len(uniq)):
                for j in range(i + 1, len(uniq)):
                    a_ = consensus[groups == uniq[i]]
                    b_ = consensus[groups == uniq[j]]
                    if len(a_) < 5 or len(b_) < 5:
                        continue
                    st, p = ks_2samp(a_, b_)
                    ks_rows.append({
                        "group_a": uniq[i], "group_b": uniq[j],
                        "n_a": int(len(a_)), "n_b": int(len(b_)),
                        "ks_statistic": round(float(st), 4), "p_value": float(p),
                        "median_a": round(float(np.median(a_)), 4),
                        "median_b": round(float(np.median(b_)), 4),
                    })
            if ks_rows:
                pd.DataFrame(ks_rows).to_csv(res_dir / "trajectory_ks_by_group.csv",
                                             index=False)
                log_info(f"多条件 KS 检验：{len(ks_rows)} 对比较")
    else:
        log_info("未配置 design.group_key，跳过多条件拟时序分布比较")

    # ---- 8. 分支点 ----------------------------------------------------------
    branch_rows = []
    if seg is not None:
        per_seg = pd.DataFrame({"seg": seg, "consensus": consensus}) \
            .groupby("seg", observed=True)["consensus"] \
            .agg(["mean", "median", "count"]).reset_index()
        for _, r in per_seg.iterrows():
            branch_rows.append({"segment": r["seg"], "n_cells": int(r["count"]),
                                "median_pseudotime": round(float(r["median"]), 4)})
        pd.DataFrame(branch_rows).to_csv(res_dir / "trajectory_segments.csv", index=False)

    # ---- 9. 图：拟时序 UMAP + 每簇分布 --------------------------------------
    # **第 3 面板的簇色必须与 umap_clusters 同源**（评审 3.8：原先 tab20，
    # 同一 cluster 在两张图里颜色不同，跨图无法对照）；且 10 个簇仅颜色
    # 编码而无图例（评审 3.6）。逐簇按 PAL_CYCLE 上色 + 显式图例。
    fig, axes = plt.subplots(1, 3, figsize=(W_DOUBLE, mm(58)))
    xy = adata.obsm["X_umap"]
    s0 = axes[0].scatter(xy[:, 0], xy[:, 1], c=consensus, s=4, cmap="viridis")
    axes[0].set_title(f"Consensus pseudotime (root = cluster {root_cluster})")
    fig.colorbar(s0, ax=axes[0], label="pseudotime (higher = later)")
    if "dpt" in corrected:
        s1 = axes[1].scatter(xy[:, 0], xy[:, 1], c=corrected["dpt"], s=4, cmap="viridis")
        axes[1].set_title("DPT pseudotime (direction-corrected)")
        fig.colorbar(s1, ax=axes[1], label="pseudotime")
    leiden_str = adata.obs["leiden"].astype(str).values
    # **与 umap_clusters 完全同一排序与配色**：那边按数值序逐簇 scatter、
    # 颜色吃 axes.prop_cycle（=PAL_CYCLE）。这里用别的排序或 tab20 都会让
    # 同一簇跨图变色（评审 3.8 的原始问题），所以逐字对齐。
    leiden_cats = sorted(set(leiden_str), key=lambda x: int(x) if x.isdigit() else x)
    for ax_i in (axes[2],):
        from matplotlib.lines import Line2D
        for ci, cat in enumerate(leiden_cats):
            m = leiden_str == str(cat)
            ax_i.scatter(xy[m, 0], xy[m, 1], s=4,
                         color=PAL_CYCLE[ci % len(PAL_CYCLE)])
        ax_i.set_title("Leiden clusters")
        handles = [Line2D([0], [0], marker="o", ls="", markersize=4,
                          color=PAL_CYCLE[ci % len(PAL_CYCLE)], label=str(cat))
                   for ci, cat in enumerate(leiden_cats)]
        ax_i.legend(handles=handles, fontsize=5, ncol=2, loc="best",
                    framealpha=0.7)
    for ax in axes:
        ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    save_fig(cfg, "02-05-04-unit1-pseudotime-umap", fig)

    per_cluster = (pd.DataFrame({"cluster": adata.obs["leiden"].astype(str).values,
                                 "consensus": consensus,
                                 "dpt": corrected.get("dpt", np.full(len(consensus), np.nan))})
                   .groupby("cluster", observed=True)
                   .agg(n_cells=("consensus", "size"),
                        consensus_mean=("consensus", "mean"),
                        consensus_median=("consensus", "median"),
                        consensus_std=("consensus", "std"),
                        dpt_median=("dpt", "median"))
                   .reset_index())
    per_cluster.to_csv(res_dir / "pseudotime_by_cluster.csv", index=False)

    fig, ax = plt.subplots(figsize=(max(W_SINGLE, 0.5 * n_clusters + 2), mm(58)))
    order = per_cluster.sort_values("consensus_median")["cluster"].tolist()
    data = [consensus[adata.obs["leiden"].astype(str).values == c] for c in order]
    # matplotlib 3.9 起 boxplot 的 `labels` 改名 `tick_labels`；老名字在 3.11 直接 TypeError
    try:
        ax.boxplot(data, tick_labels=order, showfliers=False)
    except TypeError:
        ax.boxplot(data, labels=order, showfliers=False)
    ax.set_xlabel("cluster (ordered by median consensus pseudotime)")
    ax.set_ylabel("consensus pseudotime (higher = later)")
    ax.set_title("Pseudotime distribution per cluster")
    save_fig(cfg, "02-05-05-unit1-pseudotime-by-cluster", fig)

    # ---- 10. 每细胞拟时序落盘（供 07_grn 做 regulon×拟时序）-----------------
    cell_df = pd.DataFrame({"cell": adata.obs_names.astype(str),
                            "leiden": adata.obs["leiden"].astype(str).values,
                            "consensus_pseudotime": consensus})
    for n in names:
        cell_df[f"pseudotime_{n}"] = corrected[n]
    cell_df.to_csv(res_dir / "pseudotime_per_cell.csv", index=False)

    # ---- 11. 状态 -----------------------------------------------------------
    limitations = [
        "PAGA 给出的是簇间连通性，不是分化方向",
        "拟时序是一维坐标，分支过程会被压成先后关系",
        f"**拟时序的符号是任意的**，本脚本按「值越大越晚」统一了方向；"
        f"方向参考 = {direction_note}",
        "**没有 RNA 速率（spliced/unspliced）时不能下方向性结论** —— "
        "本流水线只有计数矩阵，scVelo 记 not_done，所以这里报的是相似度排序",
        f"方法间一致性只是**内部一致性**：{len(cv_names)} 种被交叉验证的方法"
        "都错向同一个伪轨迹时，它们依然彼此高度相关。一致不等于正确",
    ]
    if direction_source == "cytotrace_fallback":
        limitations.append(
            "**方向参考是 CytoTRACE 而非 marker 基因**：若该数据集中分化潜能"
            "与成熟度不同向，整条轴会反掉。要下方向性结论请在配置里给出"
            "trajectory.early_markers / late_markers")
    if mean_rho < 0.3:
        limitations.append(
            f"**方法间一致性偏低（平均 rho={mean_rho:+.3f}）**："
            "各算法对同一数据给出了差异较大的排序，此时不该报单一「轨迹」")

    status = {
        "dataset_id": cfg["dataset_id"],
        "status": "ok",
        "n_cells": int(adata.shape[0]),
        "n_clusters": n_clusters,
        "methods_ok": methods_ok,
        "methods_failed": methods_failed,
        "n_methods": len(methods_ok),
        "root_selection": root_record,
        "direction_source": direction_source,
        "direction_note": direction_note,
        "marker_genes_used": marker_used,
        "direction_table": direction_rows,
        "method_correlation": df_to_records(cmat.reset_index().rename(
            columns={"index": "method"})),
        "direction_reference_method": reference_method,
        "cross_validated_methods": cv_names,
        "method_correlation_mean_offdiag": round(mean_rho, 4),
        "method_correlation_min_offdiag": round(min_rho, 4),
        "method_correlation_note": (
            "一致性统计**不含方向参考方法**"
            + (f"（{reference_method}）：退回模式下参考就是它本身，"
               "它与参考的相关恒为 ±1，那是定义不是证据"
               if reference_method else "")),
        "pseudotime_range": [round(float(np.nanmin(consensus)), 4),
                             round(float(np.nanmax(consensus)), 4)],
        "pseudotime_by_cluster": df_to_records(per_cluster),
        "n_genes_along_trajectory": len(genes_rows),
        "n_modules": len(modules_rows),
        "scfates_segments": int(len(set(seg))) if seg is not None else 0,
        "scfates_milestones": int(len(set(mil))) if mil is not None else 0,
        "palantir_fate_probabilities": (
            {"n_terminal_states": int(fate_probs.shape[1]),
             "shape": list(fate_probs.shape)} if fate_probs is not None else None),
        "ks_by_group": ks_rows,
        "scvelo": {
            "status": "not_done",
            "reason": ("RNA 速率需要 spliced/unspliced 两套计数矩阵。本流水线的输入是"
                       "10x **filtered 表达矩阵**，只有一套计数，没有内含子/外显子"
                       "的区分，无法计算速率。要跑 scVelo 必须从 Cell Ranger 的 "
                       "`velocyto` 或 `--include-introns` 输出重新开始。"),
        },
        "limitations": limitations,
        # ---- 哪些数可以当确定值报，哪些不能（硬性规则 20）------------------
        #
        # **"跑通了"不等于"这个数可复现"。** 实测三轮 CI（同一 Python 3.12、
        # 同一批包版本）里，`dpt` 与 `palantir` 逐位相同，而 `scfates`
        # 从 +0.5644 变到 +0.5296 再到 +0.5328，**拓扑本身也变过**
        # （4 片段/6 milestone ↔ 6 片段/8 milestone）。
        #
        # 根因不在种子，在 scFates 的扩散图那一步：
        # `scf.pp.diffusion` 签名里**没有 seed**，内部走
        # `palantir.run_diffusion_maps` → `compute_kernel(backend="scanpy")`
        # → scanpy 默认 `method="umap"` → **pynndescent 近似 kNN**（Numba 并行，
        # 给了 random_state 也不保证逐位可复现）；再被 simpleppt 的 PPT
        # 主曲线树放大成不同的拓扑。
        #
        # 所以这一段跟着产物走，而不是只写在 AGENTS.md 里 ——
        # 拿到这份 JSON 的人必须能直接看到哪个数不能当定值用。
        "reproducibility": {
            "stable_methods": ["dpt", "palantir", "cytotrace"],
            "unstable_methods": ["scfates"],
            "evidence": ("六轮 CI（8f3f56c0 / 9af13b42 / 5b241a3 / db778bb / "
                         "f99eab1 / e549477，同一 Python 3.12、同一批包版本）："
                         "dpt +0.5856 六轮相同，palantir +0.4674 六轮相同，"
                         "scfates +0.5644 -> +0.5296 -> +0.5328 -> +0.5328 -> "
                         "+0.5328 -> +0.5646"),
            "scfates_rho_observed_range": [0.5296, 0.5646],
            "mean_rho_observed_range": [0.6254, 0.6363],
            "range_is_from": ("**历史观测值，不是本轮的** —— 本轮的值见 "
                              "`method_agreement`。这里记的是"
                              "『同一个 commit 重跑能漂多少』，用来判断"
                              "本轮的值该按定值报还是按范围报"),
            "cause": ("scFates 的扩散图没有 seed 可传（`scf.pp.diffusion` "
                      "签名里没有），内部经 palantir → scanpy 默认 "
                      "`method=\"umap\"` 的近似 kNN；simpleppt 的 PPT "
                      "主曲线树把上游的末位差异放大成可见的 ρ 变化。"
                      "**残留随机源尚未定位** —— 两次归因都被实测否证"),
            "falsified_hypotheses": [
                ("多线程 BLAS 归约顺序 —— 否证：钉住 OMP/OPENBLAS/MKL + "
                 "OPENBLAS_CORETYPE 之后 scfates 仍从 +0.5296 变成 +0.5328"),
                ("pynndescent 的 Numba 并行 —— 否证：补上 "
                 "NUMBA_NUM_THREADS=1 之后三轮仍给 +0.5328 / +0.5328 / +0.5646"),
            ],
            "what_pinning_did_fix": ("**离散的拓扑稳住了**：A 轮是 4 片段/6 "
                                     "milestone，之后六轮全部是 6 片段/8 "
                                     "milestone。所以钉并行度不是白做，"
                                     "**但它不足以让 ρ 稳定**"),
            "mitigation": ("workflow 在 job 级钉了 OMP/OPENBLAS/MKL/"
                           "NUMBA_NUM_THREADS=1 与 OPENBLAS_CORETYPE=Haswell。"
                           "**实测这五个变量不足以让 scFates 的 ρ 可复现**，"
                           "所以不要以为设了就可复现 —— 按范围报"),
            "how_to_report": ("`dpt`/`palantir`/`cytotrace` 可按确定值报；"
                              "**`scfates` 与「四条轨迹平均 rho」必须带范围报** —— "
                              "只报一个数会把方法间的不一致藏起来"),
        },
    }
    write_json(res_dir / "trajectory_status.json", status)
    log_info(f"轨迹分析完成：{len(methods_ok)} 种方法，"
             f"平均一致性 rho={mean_rho:+.4f}")
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
