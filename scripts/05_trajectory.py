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
                    log_warn, parse_args, record_step, result_status_of,
                    save_fig, set_seed, write_json,
                    W_DOUBLE, W_ONE_HALF, W_SINGLE, mm, PAL, PAL_CYCLE,)

# 校正后的统一方向：**值越大越晚**
N_MODULES = 6
MIN_GENES_FOR_MODULES = 200
# **M7（R-03 裁决）**：分支段的最小细胞数。低于这个数的段，其中位拟时序
# 由极少数细胞决定（实测有 11 个细胞的段），不该和 700+ 细胞的段并列展示。
# 取 30 是常见的小样本下限：再少的话中位数对单个细胞的进出就敏感了。
MIN_CELLS_PER_SEGMENT = 30
# **M5（R-03 裁决）**：山脊图至少要有这么多个簇真的画出来才落盘。
# 空坐标系**有墨**（有轴线、有刻度），所以 `check_figures.mjs` 的
# 空白/糊死/贴边三条判据全绿 —— 空图能一路混过所有门禁。
MIN_RIDGE_GROUPS = 2


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
        # **L6（R-03 裁决）**：原来只判 `r_before < 0` 才翻转，于是
        # `r_before == 0`（**与参考完全无关**）时 `flipped=False` —— 等价于
        # "随机选一个方向"然后当成有方向的结果用。`spearmanr` 在输入是
        # 常数列（某个方法退化）时正好返回 `nan`/`0`，而 `nan < 0` 为假，
        # 也会走到"不翻转"这条路上。
        #
        # 修法：把"方向无法判定"显式标出来（`flipped=None`），而不是让它
        # 伪装成"判断过、结论是不用翻"。**下游必须真的按
        # `direction_decided` 把这类方法剔除** —— `select_cv_methods()`
        # 是唯一的消费者，不要在这里只标不用（L6 第二半）。
        _decidable = bool(np.isfinite(r_before) and abs(r_before) > 0)
        flipped = None
        if _decidable and r_before < 0:
            vv = -vv
            flipped = True
        elif _decidable:
            flipped = False
        r_after = spearmanr(vv, reference).correlation
        out[name] = vv
        rows.append({
            "method": name,
            "raw_convention": convention.get(name, "?"),
            "rho_vs_reference_before_flip": round(float(r_before), 4),
            "flipped": flipped,
            "direction_decided": _decidable,
            "direction_note": (
                None if _decidable else
                f"与参考的相关 rho={r_before!r} —— **方向无法判定**，"
                f"该方法的拟时序方向未做校正，不应被当成有方向的结果"),
            "rho_vs_reference_after_flip": round(float(r_after), 4),
        })
    return out, rows


def select_cv_methods(names, direction_rows, reference_method):
    """从候选方法里挑出**能进交叉验证与共识**的那些（L6 第二半）。

    两道剔除，**理由不同所以分开记**：

    · `direction_reference` —— 退回模式下方向参考就是 CytoTRACE 本身，
      它与自己的相关恒为 ±1，那是**定义不是证据**；
    · `direction_undecided` —— 该方法的 `rho_vs_reference` 是 nan 或 0
      （退化成常数列），**方向未知**，与别人相关多少是随机的。

    旧实现只有第一道，第二类照样进 `cv_names`，于是**一个随机方向被
    算进一致性、又被算进共识**，而下游所有"沿轨迹变化的基因""模块"
    都建在那个随机方向上 —— 产物里完全看不出来。

    抽成纯函数是为了能被标定脚本直接调（AGENTS 规则 29.1：自检必须调
    真代码，不能自己重写一遍逻辑）。

    返回 `(cv_names, cv_excluded)`。
    """
    undecided = [r["method"] for r in direction_rows
                 if not r.get("direction_decided", True)]
    cv = [n for n in names if n != reference_method and n not in undecided]
    excluded = {
        "direction_reference": ([reference_method] if reference_method else []),
        # 显式去掉一次：两组本应互斥（参考与自己的 rho 恒为 1），
        # 但哪天参考换了之后两组都数它会把剔除数报大。
        "direction_undecided": [n for n in undecided if n != reference_method],
    }
    return cv, excluded


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

    fig, ax = plt.subplots(figsize=(W_ONE_HALF, mm(66)))   # 84->66：原空白图幅过大（评审 v2）
    # **边不能盖住节点标签**（评审 3.3：簇 9 的 "9" 被 5 条粗边切断）。
    # 不做双层重画的 hack —— 只调 scanpy 自己的参数：边宽减半（0.5）、
    # 节点加大（scale 0.8 / power 0.6）、字号提到 9，标签因此在最上层可见。
    sc.pl.paga(adata, show=False, ax=ax,
               edge_width_scale=0.5, node_size_scale=0.8, node_size_power=0.6,
               fontsize=9, fontoutline=2)
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
        if r["flipped"] is None:
            log_warn(f"  方向 {r['method']:10} "
                     f"rho={r['rho_vs_reference_before_flip']:+.4f} —— "
                     f"**无法判定方向**（L6），该方法的拟时序不做方向校正")
        else:
            log_info(f"  方向 {r['method']:10} rho={r['rho_vs_reference_after_flip']:+.4f}"
                     f"{'（已翻转）' if r['flipped'] else ''}")

    # **L6 的第二半（自查 2026-09-26）。** 上面只把"方向无法判定"标出来、
    # 打了 WARN —— 但 `direction_decided` **一个消费者都没有**，而
    # `orient()` 的注释写的是"下游按 `direction_decided` 决定是否把这个
    # 方法算进共识 —— 见 `consensus` 段"。**注释陈述的行为与代码实际行为
    # 对不上**（E-58 防复发② 的又一实例）。
    #
    # 后果不是"少一条告警"，而是**一个随机方向被当成结论用**：
    # `r_before` 是 nan（该方法退化成常数列）或恰为 0 时，`flipped=None`
    # 表示方向**没被校正**，可它照样进 `cv_names` → 进一致性统计 →
    # 进共识拟时序。而共识一旦带上一个随机方向的方法，下游**所有**
    # "沿轨迹变化的基因""模块"都建在一个随机方向上，且产物里看不出来
    # —— 这正是 E-58② 要防的"看起来在算、其实没在算"。
    #
    # 剔除动作在下面的 `select_cv_methods()`（纯函数，可被标定脚本直接调）。

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
    #
    # **方向无法判定的方法也不能算进去**（L6 第二半）—— 两道剔除都
    # 落在 `select_cv_methods()` 里，理由分开记。`direction_decided`
    # 在此之前**一个消费者都没有**，而 `orient()` 的注释却写着"下游按它
    # 决定是否算进共识"：**注释陈述的行为与代码实际行为对不上**
    # （E-58 防复发② 的又一实例），后果是一个随机方向被当成结论用。
    reference_method = "cytotrace" if direction_source == "cytotrace_fallback" else None
    cv_names, cv_excluded = select_cv_methods(names, direction_rows,
                                              reference_method)
    if cv_excluded["direction_undecided"]:
        log_warn(f"**{len(cv_excluded['direction_undecided'])} 个方法的方向无法判定"
                 f"（{', '.join(cv_excluded['direction_undecided'])}），"
                 f"已从交叉验证与共识中剔除** —— 不剔除等于给共识掺进一个随机方向")
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
    #
    # **共识不能包含方向参考本身**（审计 S7）。退回模式下方向参考就是
    # CytoTRACE，它既定了所有方法的符号、又与自己的相关恒为 ±1 ——
    # 把它算进共识，等于让"定义了方向的那个方法"再投一次票。
    # 上面算交叉验证一致性时已经排除了它（`cv_names`），共识必须同口径，
    # 否则两个数讲的是两件不同的事。
    #
    # **`cv_names if cv_names else names` 这个回退已删除（L6 第二半）。**
    # 它原本保证"至少有一个方法" —— 但在**所有**方法方向都未判定时，
    # `cv_names` 为空，回退恰好把刚剔除的方法**全部放回共识**，
    # 而且日志上什么都不会显示（`mean_rho` 是 nan，`nan < 0.3` 为假，
    # 连"一致性偏低"那条限制都不会加）。**一个静默放回排除项的兜底，
    # 与不排除等价。** 现在这种情况直接判 `no_consensus`。
    if not cv_names:
        status = {
            "dataset_id": cfg["dataset_id"],
            "status": "no_consensus",
            "reason": ("没有任何方法的方向可判定"
                       f"（方向未判定: {cv_excluded['direction_undecided']}；"
                       f"方向参考: {reference_method}）"
                       "—— 方向未知的拟时序不能取共识，否则共识方向是随机的"),
            "methods_ok": methods_ok, "methods_failed": methods_failed,
            "direction_table": direction_rows,
            "cv_excluded": cv_excluded,
        }
        write_json(res_dir / "trajectory_status.json", status)
        log_warn(f"拟时序共识不可用: {status['reason']}")
        return status
    consensus_names = cv_names
    Z = np.column_stack([
        (corrected[n] - np.nanmean(corrected[n])) / (np.nanstd(corrected[n]) + 1e-12)
        for n in consensus_names
    ])
    consensus = np.nanmean(Z, axis=1)
    # 同时算一个**含参考方法**的共识，供对比：两者差多少本身就是信息
    if len(consensus_names) < len(names):
        Z_all = np.column_stack([
            (corrected[n] - np.nanmean(corrected[n])) / (np.nanstd(corrected[n]) + 1e-12)
            for n in names
        ])
        consensus_with_ref = np.nanmean(Z_all, axis=1)
        rho_ref = float(np.nan_to_num(spearmanr(consensus, consensus_with_ref).correlation))
    else:
        consensus_with_ref, rho_ref = None, None
    log_info(f"共识拟时序完成（用 {len(consensus_names)} 个方法: "
             f"{', '.join(consensus_names)}"
             + (f"；已排除方向参考 {reference_method}" if reference_method else "")
             + "）；方法间平均 rho="
             f"{mean_rho:+.4f}，最低 {min_rho:+.4f}"
             + (f"；含参考方法的共识与它 rho={rho_ref:+.4f}" if rho_ref is not None else ""))

    # ---- 6. 沿轨迹变化的基因 ------------------------------------------------
    #
    # **M4（R-03 裁决）**：整段包在一个 `except Exception` 里只 `log_warn`。
    # 后果不是"少一个文件"，而是**一个看起来正常的 0**：
    # `n_genes_along_trajectory: len(genes_rows)` 静默变成 0、顶层 `status`
    # 仍是 `ok`、`trajectory_genes.csv` 不存在也没有任何检查 ——
    # 读者无法区分"这批数据里没有沿轨迹变化的基因"和"这一段崩了"。
    # 这正是 E-48 的形态：`except` 把异常降级成一个**没人读的字段**。
    # 修法：写**嵌套** `status["along_trajectory"]`（顶层保持 `ok`，
    # 因为后面还有模块/图等步骤能出结果），让崩了这件事在产物里可见。
    genes_rows, gene_mat, gene_names, gene_rho = [], None, [], None
    along_status = {"status": "not_run", "reason": None}
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
        along_status = {"status": "ok", "reason": None,
                        "n_expressed": int(len(idx)),
                        "n_reported": int(len(genes_rows))}
        log_info(f"沿轨迹变化基因：表达基因 {len(idx)} 个，报前 {len(genes_rows)} 个")
    except Exception as e:  # noqa: BLE001
        along_status = {"status": "failed",
                        "reason": f"{type(e).__name__}: {e}",
                        "n_expressed": None, "n_reported": 0}
        log_warn(f"沿轨迹基因分析失败: {type(e).__name__}: {e}"
                 f"（已记进 trajectory_status.json 的 along_trajectory.status）")

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

            # **单图原则拆分（D-006）**：原来是 `subplots(1, 2)` 一张图两个面板
            # （用户 2026-09-24 反馈"是双图，不符合单图准则"）。拆成两张单图，
            # 各自独立达标图幅与图例；共享色标 0->p99 的口径写进各图标题。
            # 组级叙事靠"同一图号 02-05-03 的 unit1/unit2"保留关联。
            #
            # unit1：基因 x 拟时序分箱的热图（看哪些基因在哪个阶段高）
            fig1, ax1 = plt.subplots(figsize=(W_ONE_HALF, mm(70)))
            im = ax1.imshow(profz[np.argsort(lab)], aspect="auto",
                            cmap="RdBu_r", vmin=-2, vmax=2)
            ax1.set_xlabel("pseudotime bin")
            ax1.set_ylabel("gene (grouped by module)")
            ax1.set_title(f"Genes along pseudotime, {N_MODULES} modules\n"
                          "colour = z-scored mean expression (shared scale -2..2)")
            fig1.colorbar(im, ax=ax1, label="z-scored mean expression")
            save_fig(cfg, "02-05-03-unit1-trajectory-modules-heatmap", fig1)

            # unit2：各模块的平均表达曲线（看每个模块随拟时序的走向）
            # **图例必须 fig.legend + ncol=1**（约定 v2）：原 `axes.legend()`
            # 是框内图例，会压住曲线（用户反馈"图例空间布局有问题"）。
            fig2, ax2 = plt.subplots(figsize=(W_ONE_HALF, mm(70)))
            for m in range(N_MODULES):
                sel_m = np.flatnonzero(lab == m)
                if len(sel_m) == 0:
                    continue
                ax2.plot(np.nanmean(profz[sel_m], axis=0), label=f"M{m} (n={len(sel_m)})")
            ax2.set_xlabel("pseudotime bin")
            ax2.set_ylabel("z-scored mean")
            ax2.set_title("Module profiles along pseudotime")
            fig2.legend(fontsize=7, ncol=1, loc="outside right center")
            save_fig(cfg, "02-05-03-unit2-trajectory-module-profiles", fig2)
        except Exception as e:  # noqa: BLE001
            log_warn(f"基因模块分析失败: {type(e).__name__}: {e}")

    # ---- 7. 多条件拟时序分布比较（Kolmogorov-Smirnov）-----------------------
    #
    # **M2（R-03 裁决）**：原来是裸 `ks_2samp` 落盘，**没有任何多重检验校正**
    # —— 而同一个仓库的 `06_communication.py:350` 与 `07_grn.py:243-254`
    # 都做了 BH。三处口径不一致，读者无法判断哪些 `p_value` 可比。
    # 两两比较的检验家庭是 `C(k,2)`，k 个条件时增长很快（k=6 → 15 对），
    # 不校正就是把"做了 15 次检验"当成"做了 1 次"。
    ks_rows = []
    design = cfg.get("design") or {}
    group_key = design.get("group_key")
    # **M2 修正（自查，2026-09-26）**：第一版只用 `p_adj_available` 一个布尔量，
    # 于是三种完全不同的处境（没配分组 / 配了但凑不出比较对 / 校正真失败）
    # 落盘成同一句话「statsmodels 不可用」。实测 CI artifact
    # `scrna-results-68` 的 `trajectory_status.json` 就是这形态：
    # `n_pairs: 0` 而 note 说「statsmodels 不可用」—— 而
    # `run_manifest.json` 的 `versions.statsmodels = 0.15.0`，**它装着呢**。
    # 这正是 E-58 防复发②「注释陈述的行为必须与代码实际行为对得上」的
    # 又一实例：**一个真实但错误的原因，比"没有原因"更糟** —— 下一个人会去
    # 查依赖装没装，而真正要改的是配置。所以改成显式的原因码 + 按码给话。
    ks_reason = "not_configured"   # not_configured / single_group / no_pairs / import_failed / ok
    ks_correction = None
    if group_key and group_key in adata.obs.columns:
        groups = adata.obs[group_key].astype(str).values
        uniq = sorted(set(groups))
        if len(uniq) < 2:
            # **与"没配分组"分开记**：配了但只有一个取值，说明这列在**本数据集上**
            # 没有可比的分组（例如单样本数据集的 sample 列）—— 与"根本没配"是
            # 两件事，混成一句话会让读者以为配置漏了。
            ks_reason = "single_group"
            log_info(f"分组列 {group_key} 只有 1 个取值（{uniq[0] if uniq else '空'}），"
                     f"没有可比较的两组，跳过 KS 比较")
        else:
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
            # **BH 校正**（与 06/07 同口径）。`multipletests` 不可用时退化成
            # 不校正但**显式记原因码** —— 静默不校正会让读者以为 `p_value`
            # 就是可用的判据。
            if not ks_rows:
                ks_reason = "no_pairs"
                log_info(f"分组列 {group_key} 有 {len(uniq)} 个取值，"
                         f"但没有任何一对满足最小细胞数（各 >=5），跳过 KS 比较")
            else:
                try:
                    from statsmodels.stats.multitest import multipletests
                    _pv = [r["p_value"] for r in ks_rows]
                    _rej, _padj, _, _ = multipletests(_pv, method="fdr_bh")
                    for r, padj, rej in zip(ks_rows, _padj, _rej):
                        r["p_adj_bh"] = float(padj)
                        r["significant_bh"] = bool(rej)
                    ks_reason = "ok"
                    ks_correction = "fdr_bh"
                except Exception as e:  # noqa: BLE001
                    ks_reason = "import_failed"
                    log_warn(f"KS 的 BH 校正失败（{type(e).__name__}: {e}）—— "
                             f"只报原始 p_value，状态里会记 ks_by_group_note.reason="
                             f"import_failed")
                pd.DataFrame(ks_rows).to_csv(res_dir / "trajectory_ks_by_group.csv",
                                             index=False)
                n_sig = sum(1 for r in ks_rows if r.get("significant_bh"))
                log_info(f"多条件 KS 检验：{len(ks_rows)} 对比较，"
                         + (f"BH 校正后显著 {n_sig} 对" if ks_reason == "ok" else
                            "**未做 BH 校正**（statsmodels.stats.multitest 导入失败）"))
    else:
        log_info("未配置 design.group_key，跳过多条件拟时序分布比较")
    p_adj_available = ks_reason == "ok"

    # ---- 8. 分支点 ----------------------------------------------------------
    #
    # **M7（R-03 裁决）**：原来每段都无条件写进 `trajectory_segments.csv`。
    # 产物证据：某轮的段 2 只有 **11 个细胞**，`median_pseudotime = -0.0918`
    # 与其余 5 段（|median| >= 0.196）明显不同 —— 11 个细胞的"中位拟时序"
    # 没有意义，但它和 748 个细胞的那段在 CSV 里长得**一模一样**，读者会
    # 把它当同等可信的结论。这和 spatial 规则 8 的"NNLS 给一个看起来像
    # 答案、但不含信息的解"是同一形态。
    # 修法：加最小细胞数守卫，并把 `is_reliable` 显式写进表 —— **保留
    # 小段而不是丢掉**（它可能是真的小分支，删掉等于隐藏信息），但让
    # 不可靠这件事在数据里可见，而不是只活在读者脑补里。
    branch_rows = []
    n_unreliable_segments = 0
    if seg is not None:
        per_seg = pd.DataFrame({"seg": seg, "consensus": consensus}) \
            .groupby("seg", observed=True)["consensus"] \
            .agg(["mean", "median", "count"]).reset_index()
        for _, r in per_seg.iterrows():
            n_c = int(r["count"])
            reliable = n_c >= MIN_CELLS_PER_SEGMENT
            if not reliable:
                n_unreliable_segments += 1
            branch_rows.append({"segment": r["seg"], "n_cells": n_c,
                                "median_pseudotime": round(float(r["median"]), 4),
                                "is_reliable": bool(reliable),
                                "min_cells_note": ("" if reliable else
                                                   f"n<{MIN_CELLS_PER_SEGMENT}，"
                                                   "中位数不可靠")})
        pd.DataFrame(branch_rows).to_csv(res_dir / "trajectory_segments.csv", index=False)
        if n_unreliable_segments:
            log_warn(f"分支段里有 {n_unreliable_segments} 段细胞数 < "
                     f"{MIN_CELLS_PER_SEGMENT}，已在 trajectory_segments.csv 里"
                     f"标 is_reliable=false")

    # ---- 9. 图：拟时序三联 -> **单图原则拆分（D-006）** ----------------------
    # 拆成三张独立单图（S2 拟时序方法对照链）：
    #   unit1 = 共识拟时序（方法主视图）
    #   unit2 = DPT 交叉验证（两种独立方法方向一致是可信度证据）
    #   unit3 = 细胞类型着色（**裁决 2**：原第 3 面板是"再换一套簇色"，
    #            与 umap_clusters 重复、无新语义；换 celltype 后组叙事变为
    #            "拟时序 -> 交叉验证 -> 拟时序与细胞类型的关系"）
    xy = adata.obsm["X_umap"]
    fig, ax = plt.subplots(figsize=(W_ONE_HALF, mm(58)))
    s0 = ax.scatter(xy[:, 0], xy[:, 1], c=consensus, s=4, cmap="viridis")
    ax.set_title(f"Consensus pseudotime (root = cluster {root_cluster})")
    fig.colorbar(s0, ax=ax, label="pseudotime (higher = later)")
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    save_fig(cfg, "02-05-04-unit1-pseudotime-consensus", fig)
    # **PAGA 连通骨架 + root 标记**（差距清单 #17，文献范式）：
    # 只有颜色渐变时读者看不出"轨迹从哪来到哪去"。叠加**有统计依据的
    # PAGA 连通边**（宽度正比于 connectivity），边按拟时序方向加箭头，
    # root 簇用白圈标出。
    # **不用"簇质心按拟时序排序直连"** —— 实测在分散嵌入（PBMC）上
    # 会画出横跨全图的假折线（三角形 + 长横线），看起来像推断出的轨迹，
    # 而实际只是排序连线的假象：**误导性的图不如不出**（坑 0.1）。
    fig_curve, ax_c = plt.subplots(figsize=(W_ONE_HALF, mm(58)))
    ax_c.scatter(xy[:, 0], xy[:, 1], c=consensus, s=4, cmap="viridis")
    try:
        cl = adata.obs["leiden"].astype(str).values
        cats_sorted = sorted(set(cl), key=lambda x: int(x) if x.isdigit() else x)
        cents = np.array([[xy[cl == c, 0].mean(), xy[cl == c, 1].mean()]
                          for c in cats_sorted])
        med = np.array([np.median(consensus[cl == c]) for c in cats_sorted])
        # PAGA connectivity：只画非平凡边（阈值过滤统计噪声）
        conn_m = np.asarray(conn)
        n_edges = 0
        for a in range(len(cats_sorted)):
            for b in range(a + 1, len(cats_sorted)):
                w = float(conn_m[a, b])
                if w < 0.15:  # 只画强连通边（PBMC 的 PAGA 普遍偏高，低阈值会画出跨簇长边）
                    continue
                n_edges += 1
                src, dst = (a, b) if med[a] <= med[b] else (b, a)
                ax_c.annotate("",
                              xy=tuple(cents[dst]), xytext=tuple(cents[src]),
                              arrowprops=dict(arrowstyle="->", color=PAL["black"],
                                              lw=0.5 + 2.0 * w, alpha=0.65,
                                              shrinkA=6, shrinkB=6),
                              zorder=4)
        ri = cats_sorted.index(root_cluster) if root_cluster in cats_sorted else -1
        if ri >= 0:
            ax_c.scatter(cents[ri, 0], cents[ri, 1], s=110,
                         facecolors="none", edgecolors=PAL["black"],
                         linewidths=1.4, zorder=6)
            ax_c.text(cents[ri, 0], cents[ri, 1] + 0.8, " root " + str(root_cluster),
                      fontsize=7, color=PAL["black"], va="bottom", ha="center", zorder=7)
        subtitle_curve = ("PAGA connectivity skeleton (" + str(n_edges)
                          + " strong edges only, width proportional to connectivity > 0.15); "
                          + "arrows point along increasing pseudotime; open circle = root cluster "
                          + str(root_cluster))
    except Exception as e:  # noqa: BLE001
        log_warn(f"PAGA 骨架叠加失败: {type(e).__name__}: {e}")
        subtitle_curve = "PAGA skeleton unavailable"
    ax_c.set_title("Consensus pseudotime with PAGA skeleton")
    ax_c.set_xlabel("UMAP1"); ax_c.set_ylabel("UMAP2")
    fig_curve.colorbar(s0, ax=ax_c, label="pseudotime (higher = later)")
    save_fig(cfg, "02-05-04-unit4-pseudotime-principal-path", fig_curve)
    if "dpt" in corrected:
        fig, ax = plt.subplots(figsize=(W_ONE_HALF, mm(58)))
        s1 = ax.scatter(xy[:, 0], xy[:, 1], c=corrected["dpt"], s=4, cmap="viridis")
        ax.set_title("DPT pseudotime (direction-corrected)")
        fig.colorbar(s1, ax=ax, label="pseudotime")
        ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
        save_fig(cfg, "02-05-04-unit2-pseudotime-dpt", fig)
    # unit3：按 celltype（若存在）或 leiden 着色。与 umap_clusters 的
    # 配色纪律一致（PAL_CYCLE、数值序）；图例显式。
    celltype_key = "celltype" if "celltype" in adata.obs.columns else "leiden"
    cat_vals = adata.obs[celltype_key].astype(str).values
    cats = sorted(set(cat_vals), key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x))
    fig, ax = plt.subplots(figsize=(W_ONE_HALF, mm(58)))
    from matplotlib.lines import Line2D
    for ci, cat in enumerate(cats):
        m = cat_vals == str(cat)
        ax.scatter(xy[m, 0], xy[m, 1], s=4,
                   color=PAL_CYCLE[ci % len(PAL_CYCLE)])
    handles = [Line2D([0], [0], marker="o", ls="", markersize=4,
                      color=PAL_CYCLE[ci % len(PAL_CYCLE)], label=str(cat))
               for ci, cat in enumerate(cats)]
    fig.legend(handles=handles, fontsize=4.5, ncol=1, loc="outside right center",
              framealpha=0.7)
    ax.set_title(f"{celltype_key} on the same UMAP (pseudotime context)")
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    if celltype_key != "celltype":
        log_warn("celltype 列缺失 —— 02-05-04-unit3 跳过（图名账目不含 leiden 回退）")
    else:
        save_fig(cfg, "02-05-04-unit3-celltype-on-umap", fig)
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

    # 宽度直接取标准档位：夹到区间只会得到区间内的任意值（n=5 时 114.3mm，非标），
    # 档位判据（89/136/183±1.5mm）照样判红。簇多信息多 → 一栏半；簇少 → 单栏。
    fig_w = W_ONE_HALF if n_clusters > 4 else W_SINGLE
    fig, ax = plt.subplots(figsize=(fig_w, mm(58)))
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
    # **ns 必须写出来**（差距清单 #22）：只标显著的星号会把"大多数簇
    # 彼此无差异"这件事藏起来。全局 Kruskal-Wallis 检验，不显著就写 ns。
    try:
        from scipy.stats import kruskal
        groups_ok = [d[np.isfinite(d)] for d in data if np.sum(np.isfinite(d)) >= 5]
        if len(groups_ok) >= 2:
            H, p_kr = kruskal(*groups_ok)
            mark = "p < 0.001" if p_kr < 1e-3 else (
                f"p = {p_kr:.3g}" if p_kr < 0.05 else f"ns (p = {p_kr:.3g})")
            ax.set_title(f"Pseudotime distribution per cluster (Kruskal-Wallis: {mark})")
    except Exception as e:  # noqa: BLE001
        log_warn(f"Kruskal 检验失败: {type(e).__name__}: {e}")
    save_fig(cfg, "02-05-05-unit1-pseudotime-by-cluster", fig)
    # **山脊图（ridgeline）**（差距清单 #20，文献范式）：箱线只给分位数，
    # 山脊图给每个簇的**分布形状**（双峰=该簇跨两个状态）。
    # 手写 KDE + 垂直错开（不引 ggridges/seaborn 新依赖）。
    #
    # **M5（R-03 裁决）**：`save_fig` 原来在循环**之后**、`try` 之内 ——
    # 循环体里的 `continue`（某簇细胞数 < 10）不会阻止落盘，所以**所有簇
    # 都被跳过时照样写出一张空图**。空坐标轴**有墨**（轴线 + 刻度），
    # `check_figures.mjs` 的空白/糊死/贴边三条判据全绿 ⇒ 空图混过所有门禁。
    # 修法：统计真的画了几条曲线，不够就**不落盘**并把原因写进状态 ——
    # 让"这张图为什么没有"在 `trajectory_status.json` 里看得见。
    n_ridge_drawn = 0
    ridge_skip_reason = None
    try:
        from scipy.stats import gaussian_kde
        fig_r, ax_r = plt.subplots(figsize=(W_ONE_HALF, mm(72)))
        grid = np.linspace(float(np.nanmin(consensus)),
                           float(np.nanmax(consensus)), 200)
        for ri, c in enumerate(order):
            v = consensus[adata.obs["leiden"].astype(str).values == c]
            v = v[np.isfinite(v)]
            if len(v) < 10:
                continue
            kde = gaussian_kde(v)
            dens = kde(grid)
            # 每条曲线缩放到固定高度后按簇错开
            dens = dens / float(np.nanmax(dens)) * 0.8
            ax_r.fill_between(grid, ri + dens, ri, color=PAL["primary"],
                             alpha=0.45)
            ax_r.plot(grid, ri + dens, "-", lw=0.7, color=PAL["primary"])
            n_ridge_drawn += 1
        if n_ridge_drawn < MIN_RIDGE_GROUPS:
            ridge_skip_reason = (
                f"只有 {n_ridge_drawn} 个簇的细胞数 >= 10（需要 >= "
                f"{MIN_RIDGE_GROUPS} 个）—— 落盘会是空图，故不落盘")
            log_warn(f"山脊图跳过：{ridge_skip_reason}")
            plt.close(fig_r)
        else:
            ax_r.set_yticks(range(len(order)))
            ax_r.set_yticklabels(order, fontsize=6)
            ax_r.set_xlabel("consensus pseudotime (higher = later)")
            ax_r.set_ylabel("cluster")
            ax_r.set_title("Pseudotime density per cluster (ridgeline) - KDE per cluster, offset vertically")
            save_fig(cfg, "02-05-05-unit2-pseudotime-ridgeline", fig_r)
    except Exception as e:  # noqa: BLE001
        ridge_skip_reason = f"{type(e).__name__}: {e}"
        log_warn(f"山脊图失败: {ridge_skip_reason}")

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
    limitations.append(
        f"**共识拟时序用的是 {len(consensus_names)} 个方法（{', '.join(consensus_names)}）"
        + (f"，已排除方向参考 {reference_method}" if reference_method else "")
        + "。** 共识是各方法 z-score 后的等权均值 —— z-score 只把尺度归一化到 1，"
          "**不改变分布形状**，所以分辨率高（分布更极端）的方法在共识里权重更大，"
          "这不是严格的等权平均"
        + (f"；含参考方法的共识与它 rho={rho_ref:+.4f}（两者差多少本身是信息）"
           if rho_ref is not None else ""))

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
        # 共识**用了哪些方法**（审计 S7）：不含方向参考，与 cv_names 同口径。
        "consensus_methods": consensus_names,
        "consensus_excludes_direction_reference": bool(reference_method),
        "consensus_vs_with_reference_rho": (round(rho_ref, 4)
                                            if rho_ref is not None else None),
        "method_correlation_mean_offdiag": round(mean_rho, 4),
        "method_correlation_min_offdiag": round(min_rho, 4),
        # **两组剔除分开记**（L6 第二半）：方向参考是"定义上相关"、
        # 方向未判定是"方向未知"，理由不同、处理也不同（前者无法修，
        # 后者应去看那个方法为什么退化）。合并成一个数会让读者以为
        # 剔除的都是同一类东西。
        "cv_excluded": cv_excluded,
        "method_correlation_note": (
            "一致性统计**不含方向参考方法**"
            + (f"（{reference_method}）：退回模式下参考就是它本身，"
               "它与参考的相关恒为 ±1，那是定义不是证据"
               if reference_method else "")),
        "pseudotime_range": [round(float(np.nanmin(consensus)), 4),
                             round(float(np.nanmax(consensus)), 4)],
        "pseudotime_by_cluster": df_to_records(per_cluster),
        "n_genes_along_trajectory": len(genes_rows),
        # **M4（R-03 裁决）**：`n_genes_along_trajectory` 为 0 时，读者
        # 分不清"数据里没有"和"这段崩了"。下面这个嵌套 status 把两者分开。
        "along_trajectory": along_status,
        # **M5**：山脊图没落盘时，这里必须说明为什么。
        "ridgeline": ({"status": "ok", "n_groups_drawn": n_ridge_drawn}
                      if ridge_skip_reason is None else
                      {"status": "skipped", "reason": ridge_skip_reason,
                       "n_groups_drawn": n_ridge_drawn}),
        # **M7**：分支段里有多少段细胞数不足。
        "n_unreliable_segments": n_unreliable_segments,
        "min_cells_per_segment": MIN_CELLS_PER_SEGMENT,
        "n_modules": len(modules_rows),
        "scfates_segments": int(len(set(seg))) if seg is not None else 0,
        "scfates_milestones": int(len(set(mil))) if mil is not None else 0,
        "palantir_fate_probabilities": (
            {"n_terminal_states": int(fate_probs.shape[1]),
             "shape": list(fate_probs.shape)} if fate_probs is not None else None),
        "ks_by_group": ks_rows,
        # **M2**：KS 的 p 值分辨率与校正状态必须显式落盘 —— 只给一串
        # `p_value` 而不说"校没校正、下限是多少"，读者会把它当成可用的
        # 判据（spatial 侧同问题记在 `p_value_resolution`，口径对齐）。
        "ks_by_group_note": {
            "p_adj_available": p_adj_available,
            "correction": ks_correction,
            "reason": ks_reason,
            "n_pairs": len(ks_rows),
            "note": {
                "ok": ("两两 KS 检验的家庭是全部条件对；"
                       "`p_adj_bh` 是 BH 校正后的值，判显著性看它"),
                "not_configured": ("未配置 `design.group_key`（或该列不在 obs 里），"
                                   "本轮**没有做**多条件拟时序分布比较 —— "
                                   "这是「没做」，不是「做了没问题」"),
                "single_group": ("`design.group_key` 配了，但该列在本数据集里"
                                 "只有 1 个取值，没有可比较的两组 —— "
                                 "这是「数据里没有分组」，不是「配置漏了」"),
                "no_pairs": ("分组列存在，但没有任何一对条件同时满足最小细胞数"
                             "（各 >=5），本轮**没有做** KS 比较 —— "
                             "这是「样本不够」，不是「分布没有差异」"),
                "import_failed": ("`statsmodels.stats.multitest` 导入失败，"
                                  "`p_value` 是**未校正**值，不能直接当判据"),
            }[ks_reason],
        },
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
        res = run_05_trajectory(cfg)
        # E-56：本步骤有 6 条早退路径，其中 `bad_root` / `missing_clusters` /
        # `insufficient_methods` / `failed` 都是**语义上的失败**却不抛异常 ——
        # 不接住返回值就会记成 ok。
        record_step(cfg, "trajectory", "ok", time.time() - t0,
                    result_status=result_status_of(res))
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "trajectory", "failed", time.time() - t0,
                    message=str(e), result_status="failed")
        raise
