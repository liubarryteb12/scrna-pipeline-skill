#!/usr/bin/env python3
"""
03_cluster_annotate.py — 邻居图、UMAP、Leiden 聚类、marker 基因、细胞类型打分

**注释是打分提示，不是结论。** 见 assets/celltype_markers.yml 开头的说明。
每个簇的 assignment 都带 `score_margin`（第一名与第二名的差）——
margin 小的 assignment 不该被当成结论，而这一点只有把 margin 写出来
才看得出来。只报一个类型名等于把不确定性藏起来。

**聚类分辨率是有后果的选择。** resolution 决定簇的粒度，而下游所有
marker / 注释 / 通讯分析都建立在这个粒度上。这里额外跑一个
resolution 扫描，把"多少个簇"随分辨率怎么变落盘 —— 让读者能判断
选 1.0 是不是合理，而不是只能接受它。
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

from common import (df_to_records, ensure_dirs, load_config, log_info,  # noqa: E402
                    log_warn, parse_args, record_step, save_fig, set_seed, write_json)


def load_signatures(cfg: dict) -> dict:
    p = Path(__file__).resolve().parent.parent / "assets" / "celltype_markers.yml"
    with open(p, encoding="utf-8") as fh:
        doc = yaml.safe_load(fh)
    name = (cfg.get("analysis") or {}).get("celltype_markers", "default")
    sigs = doc.get("signatures", {})
    if name not in sigs:
        raise KeyError(f"celltype_markers='{name}' 不在 {p.name} 里；"
                       f"可选: {', '.join(sorted(sigs))}")
    return sigs[name]


def score_celltypes(adata, sig) -> tuple:
    """
    对每个簇算每种细胞类型的平均打分。

    用 scanpy 的 `score_genes`（对每个细胞的签名基因做相对表达打分），
    再按簇取均值 —— 比"簇内平均表达"稳，因为它校正了基因的表达水平与
    技术噪声（用的是随机对照基因集）。

    **返回 margin。** 只报 argmax 会把"T 细胞 0.31 vs CD4 T 0.30"和
    "T 细胞 0.31 vs 上皮 0.02"报成同一个结论，而它们完全不同。
    """
    celltypes = sig["celltypes"]
    present, missing = {}, {}
    for ct, d in celltypes.items():
        genes = [g for g in d.get("markers", []) if g in adata.var_names]
        if genes:
            present[ct] = genes
        missing[ct] = [g for g in d.get("markers", []) if g not in adata.var_names]

    if not present:
        return None, None, {"reason": "签名里的基因一个都不在数据里"}

    for ct, genes in present.items():
        sc.tl.score_genes(adata, genes, score_name=f"score_{ct}",
                          random_state=42, use_raw=True)

    score_cols = [f"score_{ct}" for ct in present]
    per_cell = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
    per_cell.columns = list(present.keys())

    # 每个簇：第一名、第二名、margin
    rows = []
    for cl in per_cell.index:
        s = per_cell.loc[cl].sort_values(ascending=False)
        top = s.index[0]
        margin = float(s.iloc[0] - s.iloc[1]) if len(s) > 1 else float("nan")
        rows.append({
            "cluster": str(cl),
            "assigned": top,
            "top_score": round(float(s.iloc[0]), 4),
            "runner_up": s.index[1] if len(s) > 1 else None,
            "runner_up_score": round(float(s.iloc[1]), 4) if len(s) > 1 else None,
            "score_margin": round(margin, 4),
            # margin 小于 0.05 时第一名与第二名基本无差别
            "assignment_confident": bool(margin >= 0.05),
        })
    assign = pd.DataFrame(rows)
    diag = {
        "n_celltypes_scored": len(present),
        "celltypes": sorted(present.keys()),
        "missing_markers": {k: v for k, v in missing.items() if v},
        "n_clusters_low_margin": int((~assign["assignment_confident"]).sum()),
    }
    return per_cell, assign, diag


def run_03_cluster_annotate(cfg: dict) -> dict:
    ensure_dirs(cfg)
    set_seed(cfg)
    data_dir = Path(cfg["output"]["data_dir"])
    res_dir = Path(cfg["output"]["results_dir"])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    adata = sc.read_h5ad(data_dir / "integrated.h5ad")
    integ = {}
    try:
        import json
        with open(res_dir / "integration_status.json", encoding="utf-8") as fh:
            integ = json.load(fh)
    except Exception:  # noqa: BLE001
        pass
    use_rep = integ.get("use_rep", "X_pca")
    log_info(f"读入 {adata.n_obs} 细胞 x {adata.n_vars} HVG（use_rep={use_rep}）")

    rd = cfg["reduce"]
    n_pcs = int(rd["n_pcs"])
    n_pcs_use = rd.get("n_pcs_use") or n_pcs

    # ---- 1. 邻居 + UMAP -----------------------------------------------------
    sc.pp.neighbors(adata, n_neighbors=int(rd["n_neighbors"]),
                    n_pcs=int(n_pcs_use), use_rep=use_rep,
                    random_state=cfg["analysis"]["seed"])
    sc.tl.umap(adata, random_state=cfg["analysis"]["seed"])
    log_info(f"UMAP 完成（n_neighbors={rd['n_neighbors']}, n_pcs={n_pcs_use}）")

    # ---- 2. 分辨率扫描 ------------------------------------------------------
    # 让"选 1.0"变成有依据的决定，而不是默认值。
    scan_res = [0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0]
    scan_rows = []
    for r in scan_res:
        sc.tl.leiden(adata, resolution=r, key_added=f"leiden_r{r}",
                     flavor="igraph", n_iterations=2, directed=False,
                     random_state=cfg["analysis"]["seed"])
        scan_rows.append({"resolution": r,
                          "n_clusters": int(adata.obs[f"leiden_r{r}"].nunique())})
    scan = pd.DataFrame(scan_rows)
    scan.to_csv(res_dir / "cluster_resolution_scan.csv", index=False)
    log_info("分辨率扫描: " + ", ".join(f"{r.resolution}->{r.n_clusters}" for r in scan.itertuples()))

    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    ax.plot(scan["resolution"], scan["n_clusters"], "o-", color="#2C7FB8")
    ax.axvline(float(rd["resolution"]), color="#B2182B", ls="--", lw=1,
               label=f"used: {rd['resolution']}")
    ax.set_xlabel("Leiden resolution"); ax.set_ylabel("number of clusters")
    ax.set_title("Cluster count vs resolution", fontsize=10)
    ax.legend(fontsize=8)
    save_fig(cfg, "cluster_resolution_scan", fig)

    # ---- 3. 用配置的分辨率定稿 ----------------------------------------------
    res_used = float(rd["resolution"])
    sc.tl.leiden(adata, resolution=res_used, key_added="leiden",
                 flavor="igraph", n_iterations=2, directed=False,
                 random_state=cfg["analysis"]["seed"])
    n_clusters = int(adata.obs["leiden"].nunique())
    log_info(f"最终聚类: {n_clusters} 个簇（resolution={res_used}）")

    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    xy = adata.obsm["X_umap"]
    cats = adata.obs["leiden"].astype(str).values
    for c in sorted(set(cats), key=lambda x: int(x) if x.isdigit() else x):
        m = cats == c
        ax.scatter(xy[m, 0], xy[m, 1], s=4, alpha=0.75, label=c)
        cx, cy = xy[m, 0].mean(), xy[m, 1].mean()
        ax.text(cx, cy, c, fontsize=9, weight="bold",
                ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75))
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    ax.set_title(f"Leiden clusters (n={n_clusters}, resolution={res_used})", fontsize=10)
    ax.legend(fontsize=6, markerscale=2.5, loc="center left", bbox_to_anchor=(1.0, 0.5))
    save_fig(cfg, "umap_clusters", fig)

    # ---- 4. Marker 基因 -----------------------------------------------------
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon",
                            use_raw=True, random_state=cfg["analysis"]["seed"])
    top_n = int((cfg.get("analysis") or {}).get("top_markers", 25))
    frames = []
    for g in adata.obs["leiden"].cat.categories:
        df = sc.get.rank_genes_groups_df(adata, group=g).head(top_n)
        df.insert(0, "cluster", str(g))
        frames.append(df)
    markers = pd.concat(frames, ignore_index=True)
    markers.to_csv(res_dir / "markers_all.csv", index=False)
    log_info(f"marker 表: {len(markers)} 行（每簇前 {top_n}）")

    # marker 点图：每簇取前 3 个
    top3 = (markers.sort_values(["cluster", "scores"], ascending=[True, False])
            .groupby("cluster", observed=True).head(3)["names"].unique().tolist())
    top3 = [g for g in top3 if g in adata.raw.var_names][:40]
    if top3:
        sc.pl.dotplot(adata, top3, groupby="leiden", use_raw=True, show=False,
                      standard_scale="var")
        fig = plt.gcf()
        fig.suptitle("Top markers per cluster", fontsize=10)
        save_fig(cfg, "markers_dotplot", fig)

    # ---- 5. 细胞类型打分 ----------------------------------------------------
    sig = load_signatures(cfg)
    per_cell, assign, diag = score_celltypes(adata, sig)
    annot_record = {"signature": (cfg.get("analysis") or {}).get("celltype_markers", "default"),
                    "signature_description": sig.get("description"),
                    "signature_note": sig.get("note")}
    if assign is None:
        annot_record.update({"status": "not_possible", **diag})
        log_warn(f"细胞类型打分不可行: {diag.get('reason')}")
    else:
        per_cell.to_csv(res_dir / "celltype_scores.csv")
        assign.to_csv(res_dir / "celltype_annotation.csv", index=False)
        # 把 assignment 写回 obs（下游通讯/轨迹要用）
        m = dict(zip(assign["cluster"], assign["assigned"]))
        adata.obs["celltype"] = adata.obs["leiden"].astype(str).map(m).astype("category")
        annot_record.update({"status": "ok", **diag,
                             "assignments": df_to_records(assign)})
        log_info(f"细胞类型打分: {diag['n_celltypes_scored']} 种类型；"
                 f"{diag['n_clusters_low_margin']}/{n_clusters} 个簇 margin<0.05（assignment 不确定）")
        if diag["n_clusters_low_margin"] > 0:
            log_warn("margin<0.05 的簇其 assignment 不应被当成结论 —— 见 celltype_annotation.csv")

        # 打分热图
        fig, ax = plt.subplots(figsize=(max(6.0, 0.45 * len(per_cell.columns) + 3),
                                        max(3.0, 0.32 * len(per_cell) + 1.6)))
        im = ax.imshow(per_cell.values, aspect="auto", cmap="viridis")
        ax.set_xticks(range(len(per_cell.columns)))
        ax.set_xticklabels(per_cell.columns, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(per_cell)))
        ax.set_yticklabels([str(i) for i in per_cell.index], fontsize=7)
        ax.set_xlabel("cell type signature"); ax.set_ylabel("cluster")
        ax.set_title("Mean signature score per cluster", fontsize=10)
        fig.colorbar(im, ax=ax, label="score")
        save_fig(cfg, "celltype_scores_heatmap", fig)

    # ---- 6. 落盘 ------------------------------------------------------------
    out = data_dir / "clustered.h5ad"
    adata.write_h5ad(out)
    log_info(f"已写出 {out}")

    status = {
        "dataset_id": cfg["dataset_id"],
        "n_cells": int(adata.n_obs),
        "n_clusters": n_clusters,
        "resolution": res_used,
        "resolution_scan": df_to_records(scan),
        "n_neighbors": int(rd["n_neighbors"]),
        "n_pcs_used": int(n_pcs_use),
        "use_rep": use_rep,
        "n_markers_rows": int(len(markers)),
        "cluster_sizes": {str(k): int(v) for k, v in
                          adata.obs["leiden"].value_counts().sort_index().items()},
        "annotation": annot_record,
        "status": "ok",
    }
    write_json(res_dir / "cluster_status.json", status)
    return status


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    t0 = time.time()
    try:
        run_03_cluster_annotate(cfg)
        record_step(cfg, "cluster_annotate", "ok", time.time() - t0)
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "cluster_annotate", "failed", time.time() - t0, message=str(e))
        raise
