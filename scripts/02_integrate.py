#!/usr/bin/env python3
"""
02_integrate.py — 标准化、高变基因、PCA、批次整合

**顺序很重要，而且 seurat_v3 与另外两种口味的顺序不同：**
  - seurat_v3 吃**原始计数**，必须在 normalize/log1p **之前**跑
  - seurat / cell_ranger 吃 log 后的数据，必须在之后
实测踩过的坑：顺序反了 seurat_v3 不报错，只是给出错的 HVG ——
而 HVG 错了下游全部跟着错，图上完全看不出来。

批次整合默认关闭。**单样本数据没有批次，关掉是正确的**，不是偷懒；
多供体数据必须开，否则聚类被供体差异主导，而那是技术差异不是细胞类型差异。
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))

import numpy as np  # noqa: E402
import scanpy as sc  # noqa: E402

from common import (ensure_dirs, load_config, log_info, log_warn,  # noqa: E402
                    parse_args, record_step, result_status_of, save_fig,
                    set_seed, write_json,
                    PAL, W_DOUBLE, W_ONE_HALF, mm,)


def run_02_integrate(cfg: dict) -> dict:
    ensure_dirs(cfg)
    set_seed(cfg)
    data_dir = Path(cfg["output"]["data_dir"])
    res_dir = Path(cfg["output"]["results_dir"])

    import matplotlib
    matplotlib.use("Agg")

    adata = sc.read_h5ad(data_dir / "qc_filtered.h5ad")
    log_info(f"读入 {adata.n_obs} 细胞 x {adata.n_vars} 基因")

    norm = cfg["norm"]
    flavor = norm.get("hvg_flavor", "seurat_v3")
    n_top = int(norm["n_top_genes"])
    target_sum = float(norm["target_sum"])
    batch_key = (cfg.get("design") or {}).get("batch_key")

    # 原始计数存进 layer —— 下游 pseudobulk 与 GRN 都要用计数，不能用 log 后的
    adata.layers["counts"] = adata.X.copy()

    # ---- 1. HVG（顺序按口味分）---------------------------------------------
    #
    # **seurat_v3 需要 scikit-misc（skmisc.loess）。** 缺它时 scanpy 直接
    # ModuleNotFoundError。这里降级到 `seurat` 口味，但**必须记录口味变了** ——
    # 两种口味选出的 HVG 集合不同，下游 PCA/聚类/所有 marker 都跟着变。
    # 静默降级会让"结果和上次不一样"变成无法解释的事。
    hvg_flavor_used = flavor
    hvg_fallback = None
    if flavor == "seurat_v3":
        try:
            import skmisc  # noqa: F401
        except ImportError:
            hvg_flavor_used = "seurat"
            hvg_fallback = {
                "requested": "seurat_v3",
                "used": "seurat",
                "reason": "scikit-misc (skmisc) 未安装，seurat_v3 的 loess 拟合不可用",
                "impact": "两种口味选出的 HVG 集合不同，下游 PCA/聚类/marker 全部跟着变",
                "fix": "pip install scikit-misc",
            }
            log_warn(f"HVG 口味降级: seurat_v3 -> seurat —— {hvg_fallback['reason']}")

    if hvg_flavor_used == "seurat_v3":
        # 必须在 normalize 之前：seurat_v3 内部对计数做负二项拟合
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor="seurat_v3",
                                    batch_key=batch_key, layer="counts")
        log_info("HVG (seurat_v3, 用原始计数) 先算，再 normalize")
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
    else:
        sc.pp.normalize_total(adata, target_sum=target_sum)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top,
                                    flavor=hvg_flavor_used, batch_key=batch_key)
        log_info(f"HVG ({hvg_flavor_used}, 用 log 后数据) 后算")

    n_hvg = int(adata.var["highly_variable"].sum())
    log_info(f"高变基因: {n_hvg}")
    adata.raw = adata

    # HVG 图：均值-离散度，标出被选中的。
    #
    # **不用 `sc.pl.highly_variable_genes`**：它是高层 API，自己建 figure、
    # figsize 不受控（实测 177.8mm，不在三档标准栏宽里），且样式门禁
    # （fig_dpi / legend_convention）对它也无从约束。`sc.pp.highly_variable_genes`
    # 已经把 `means` / 离散度指标 / `highly_variable` 写进 `adata.var`，
    # 直接用那些列自己画（不重算，口径与 scanpy 一致）。
    # 离散度列随 flavor 而异：seurat_v3 是 variances_norm，其余是 dispersions_norm。
    import matplotlib.pyplot as plt
    import matplotlib.patches

    disp_col = "variances_norm" if hvg_flavor_used == "seurat_v3" else "dispersions_norm"
    hv = adata.var["highly_variable"].to_numpy()
    means = adata.var["means"].to_numpy()
    disp = adata.var[disp_col].to_numpy()
    colors = np.where(hv, PAL["highlight"], PAL["muted"])
    fig, ax = plt.subplots(figsize=(W_ONE_HALF, mm(70)))
    ax.scatter(means, disp, c=colors, s=3, linewidths=0, alpha=0.7)
    ax.set_xlabel("mean expression of genes")
    ax.set_ylabel(f"{disp_col} of genes")
    ax.set_title(f"HVG selection ({hvg_flavor_used}, n={n_hvg})\n"
                 "orange = highly variable; grey = other genes")
    fig.legend(handles=[
        matplotlib.patches.Patch(color=PAL["highlight"], label="highly variable genes"),
        matplotlib.patches.Patch(color=PAL["muted"], label="other genes"),
    ], fontsize=7, ncol=1, loc="outside right center")
    save_fig(cfg, "02-02-01-unit1-hvg-selection", fig)

    # ---- 2. 子集到 HVG ------------------------------------------------------
    work = adata[:, adata.var["highly_variable"]].copy()

    # ---- 3. PCA -------------------------------------------------------------
    # **M13（R-03 裁决）**：原来没传 `n_comps`，于是算的是 scanpy 的默认 50 个，
    # 而状态文件里报的是 `cfg["reduce"]["n_pcs"]`（配置 40）—— **报的数和
    # 算的数不是一个数**。第 41–50 个 PC 算了却从不被报告或使用，白算；
    # 更糟的是读者拿 `n_pcs=40` 去核对 `pca_variance_ratio_top10` 时会以为
    # 两者同源。修法：把配置值真正传进去，两者从此必然一致。
    # 配置值可能超过数据维度（小数据集 / 少 HVG）—— scanpy 会直接抛错。
    # 夹到 `min(n_obs, n_vars)` 并把**实际算的个数**记进状态（不是配置值）。
    n_comps_max = int(min(work.n_obs, work.n_vars))
    n_comps = min(int(cfg["reduce"]["n_pcs"]), n_comps_max)
    n_comps_clamped = n_comps != int(cfg["reduce"]["n_pcs"])
    sc.tl.pca(work, n_comps=n_comps, svd_solver="arpack",
              random_state=cfg["analysis"]["seed"])
    import matplotlib.pyplot as plt
    # **不用 `sc.pl.pca_variance_ratio`。** 实测它出的图有三个问题：
    #
    #   1. **它给每一个 PC 都打一个标注**，而 `sc.tl.pca` 默认算 50 个 PC ——
    #      在 89 mm 宽的画布上 PC13 之后完全叠成一团，PC10/11/12 也挤在一起。
    #      标注沿着曲线排，看起来像刻度但不是刻度，改不了。
    #   2. **PC1 的标注顶到标题上**，把 "PCA variance ratio (elbow)" 横穿。
    #   3. **不给 y 轴标题**，读者只看到 -2.5 ~ -6.0 的数字，
    #      不知道那是 log10(方差解释比)。
    #
    # 自己画还顺带**去掉一个版本依赖** —— 仓库规则 6：绘图 API 会随版本变
    # （`sc.pl.highly_variable_genes(return_fig=True)` 就在 scanpy 1.12 被移除）。
    var_ratio = work.uns["pca"]["variance_ratio"]
    n_pc = len(var_ratio)
    fig, ax = plt.subplots(figsize=(W_ONE_HALF, mm(66)))
    xs = np.arange(1, n_pc + 1)
    ax.plot(xs, np.log10(var_ratio), marker="o", ms=2.5, lw=0.9,
            color=PAL["primary"], markeredgewidth=0)
    ax.set_xlabel("Principal component (rank)")
    ax.set_ylabel("log10(variance ratio)")
    ax.set_title("PCA variance ratio (elbow)")
    # **刻度每 5 个 PC 一个。** 逐个标就是上面那个叠字问题 ——
    # 而"多少个刻度放得下"取决于画布宽度，所以按宽度算而不是写死。
    step = 5 if n_pc > 12 else 1
    ticks = list(range(1, n_pc + 1, step))
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"PC{i}" for i in ticks])
    ax.grid(axis="y", alpha=0.25, linewidth=0.5)
    save_fig(cfg, "02-02-02-unit1-pca-variance-ratio", fig)
    log_info(f"PCA 方差比图: 标出 {len(ticks)} 个刻度（共 {n_pc} 个 PC）")

    # **L1（R-03 裁决）**：这里原来又写了一遍
    # `var_ratio = work.uns["pca"]["variance_ratio"]` —— 与上面那次赋值完全
    # 相同，是死代码。删掉（不删的代价不是性能，是读者会以为"这里重新读了
    # 一次，说明中间可能变过"）。
    log_info(f"PC1-{min(10, len(var_ratio))} 方差解释: "
             f"{', '.join(f'{v:.3f}' for v in var_ratio[:10])}")

    # ---- 4. 批次整合 --------------------------------------------------------
    integ = cfg.get("integration") or {}
    method = integ.get("method", "none")
    use_rep = "X_pca"
    integ_record = {"method": method, "batch_key": batch_key}

    if method == "none":
        integ_record["status"] = "not_applied"
        if batch_key:
            integ_record["reason"] = ("配置 method=none 但设了 batch_key —— "
                                      "多批次数据不做整合会让聚类被批次差异主导")
            log_warn(integ_record["reason"])
        else:
            integ_record["reason"] = "单样本数据没有批次，不做整合是正确的选择"
            log_info(integ_record["reason"])
    elif not batch_key:
        integ_record["status"] = "not_configured"
        integ_record["reason"] = "integration.method 不是 none 但 design.batch_key 为空"
        log_warn(integ_record["reason"])
    elif method == "harmony":
        try:
            sc.external.pp.harmony_integrate(work, batch_key,
                                             random_state=cfg["analysis"]["seed"])
            use_rep = "X_pca_harmony"
            integ_record["status"] = "ok"
            integ_record["use_rep"] = use_rep
            log_info(f"harmony 整合完成（batch_key={batch_key}）")
        except ImportError as e:
            integ_record["status"] = "package_missing"
            integ_record["reason"] = f"harmonypy 未安装: {e}"
            log_warn(integ_record["reason"])
    elif method == "combat":
        sc.pp.combat(work, key=batch_key)
        # M13：重跑 PCA 时**必须带同样的 `n_comps`** —— 否则 ComBat 分支算 50 个、
        # 非 ComBat 分支算 40 个，两条路径的 `n_pcs` 语义不同而状态文件
        # 报的是同一个配置值。
        sc.tl.pca(work, n_comps=n_comps, svd_solver="arpack",
                  random_state=cfg["analysis"]["seed"])
        var_ratio = work.uns["pca"]["variance_ratio"]
        integ_record["status"] = "ok"
        integ_record["note"] = "ComBat 后重跑了 PCA（ComBat 直接改 X，X_pca 已失效）"
        log_info(integ_record["note"])
    else:
        raise ValueError(f"不支持的 integration.method: {method}（none/harmony/combat）")

    # 批次效应可视化（有 batch_key 才有意义）
    if batch_key and batch_key in work.obs.columns:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(W_DOUBLE, mm(64)))
        for ax, rep in zip(axes, ["X_pca", use_rep]):
            if rep not in work.obsm:
                ax.axis("off"); continue
            xy = work.obsm[rep][:, :2]
            cats = work.obs[batch_key].astype(str).values
            for i, c in enumerate(sorted(set(cats))):
                m = cats == c
                ax.scatter(xy[m, 0], xy[m, 1], s=2, alpha=0.5, label=str(c))
            ax.set_title(f"{rep} by {batch_key}")
            ax.set_xlabel("dim 1"); ax.set_ylabel("dim 2")
            if len(set(cats)) <= 12:
                fig.legend(fontsize=6, markerscale=3, ncol=1, loc="outside right center")
        fig.suptitle("Batch mixing before/after integration")
        save_fig(cfg, "02-02-03-unit1-batch-mixing", fig)

    # ---- 5. 落盘 ------------------------------------------------------------
    out = data_dir / "integrated.h5ad"
    work.write_h5ad(out)
    log_info(f"已写出 {out}（{work.n_obs} 细胞 x {work.n_vars} HVG）")

    status = {
        "dataset_id": cfg["dataset_id"],
        "n_cells": int(work.n_obs),
        "n_hvg": n_hvg,
        "hvg_flavor_requested": flavor,
        "hvg_flavor_used": hvg_flavor_used,
        "hvg_fallback": hvg_fallback,
        "n_top_genes_requested": n_top,
        "target_sum": target_sum,
        # **M13（R-03 裁决）**：原来这里写的是配置值，而实际算的是 scanpy
        # 默认 50 —— **报的数和算的数不是一个数**。现在报实际值，并在被
        # 数据维度夹住时显式说明（"我按你说的做了"和"我改小了"必须能区分）。
        "n_pcs": int(n_comps),
        "n_pcs_requested": int(cfg["reduce"]["n_pcs"]),
        "n_pcs_clamped": bool(n_comps_clamped),
        "pca_variance_ratio_top10": [round(float(v), 5) for v in var_ratio[:10]],
        "pca_cumvar_top10": [round(float(v), 5) for v in np.cumsum(var_ratio[:10])],
        "integration": integ_record,
        "use_rep": use_rep,
        # **M11 配套**：批次混合图的产出与否取决于有没有 `batch_key`，
        # 而 `CONDITIONAL_FIGURES` 的判据需要一个能读的字段。以前这个
        # 条件只写在代码的 `if batch_key and ...` 里，验收层看不见 ——
        # 于是没有批次的数据集上这张图"声明了但没产出"会被判红。
        # 落盘成显式布尔值，判据就自愈了：哪天配了 batch_key，这张图
        # 立刻自动变成必需。
        "has_batch_key": bool(batch_key and batch_key in work.obs.columns),
        "status": "ok",
    }
    write_json(res_dir / "integration_status.json", status)
    return status


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    t0 = time.time()
    try:
        res = run_02_integrate(cfg)
        # E-56：接住返回值，别把"跑完了但没做成"记成 ok。
        record_step(cfg, "integrate", "ok", time.time() - t0,
                    result_status=result_status_of(res))
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "integrate", "failed", time.time() - t0,
                    message=str(e), result_status="failed")
        raise
