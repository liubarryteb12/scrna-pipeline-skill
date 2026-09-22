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
                    parse_args, record_step, save_fig, set_seed, write_json,
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

    # HVG 图：均值-离散度，标出被选中的
    #
    # **注意 scanpy 版本的 API 差异。** `sc.pl.highly_variable_genes` 在
    # 1.12 里**没有 `return_fig`**（传了会 TypeError），它画到当前 figure 上。
    # 所以这里用 plt.gcf() 取，而不是靠返回值。
    import matplotlib.pyplot as plt
    sc.pl.highly_variable_genes(adata, show=False)
    fig = plt.gcf()
    fig.suptitle(f"HVG ({hvg_flavor_used}, n={n_hvg})")
    save_fig(cfg, "02-02-01-unit1-hvg-selection", fig)

    # ---- 2. 子集到 HVG ------------------------------------------------------
    work = adata[:, adata.var["highly_variable"]].copy()

    # ---- 3. PCA -------------------------------------------------------------
    sc.tl.pca(work, svd_solver="arpack", random_state=cfg["analysis"]["seed"])
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

    var_ratio = work.uns["pca"]["variance_ratio"]
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
        sc.tl.pca(work, svd_solver="arpack", random_state=cfg["analysis"]["seed"])
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
                ax.legend(fontsize=6, markerscale=3, loc="outside upper right")
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
        "n_pcs": int(cfg["reduce"]["n_pcs"]),
        "pca_variance_ratio_top10": [round(float(v), 5) for v in var_ratio[:10]],
        "pca_cumvar_top10": [round(float(v), 5) for v in np.cumsum(var_ratio[:10])],
        "integration": integ_record,
        "use_rep": use_rep,
        "status": "ok",
    }
    write_json(res_dir / "integration_status.json", status)
    return status


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    t0 = time.time()
    try:
        run_02_integrate(cfg)
        record_step(cfg, "integrate", "ok", time.time() - t0)
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "integrate", "failed", time.time() - t0, message=str(e))
        raise
