#!/usr/bin/env python3
"""
03_cluster_annotate.py — 邻居图、UMAP、Leiden 聚类、marker 基因、细胞类型打分

**注释是打分提示，不是结论。** 见 assets/celltype_markers.yml 开头的说明。
每个簇的 assignment 都带 `score_margin`（第一名与第二名的差）——
margin 小的 assignment 不该被当成结论，而这一点只有把 margin 写出来
才看得出来。只报一个类型名等于把不确定性藏起来。

**两条独立的注释路径**（文档 §2.4）：
  1. **marker 签名打分**（`score_genes` + 簇均值）→ `celltype_annotation.csv`
  2. **CellTypist 预训练模型**（文档指定的自动注释工具）→ `celltypist_labels.csv`

两者一致时结论更可信，不一致时那个簇值得人看 —— 所以这里**量化一致率**
而不是只报其中一条。CellTypist 的模型是运行期下载的，拿不到就记
`model_unavailable`，不退回、不假装。

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
                    log_warn, parse_args, record_step, result_status_of, save_fig, set_seed, write_json, W_DOUBLE, W_ONE_HALF, W_SINGLE, mm, plot_marker_dotplot, build_marker_dotplot_figure, PAL,)


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


# ---------------------------------------------------------------------------
# CellTypist（文档 §2.4 指定的自动注释工具）
# ---------------------------------------------------------------------------
# **与 marker 签名打分是两条独立的路。** marker 打分只看"签名基因在簇里
# 平均高不高"；CellTypist 用一个在几十万细胞上训过的逻辑回归模型，看
# **全部基因**的加权组合。两者一致时结论更可信，不一致时那个簇值得人看。
#
# **模型文件是运行期下载的**（~50MB，落在 ~/celltypist/）。CI 里可能拿不到，
# 拿不到就记 not_available，**不退回、不假装**。
#
# **输入必须是对数化的全基因集。** 本流水线在 02_integrate 里把 `adata.raw`
# 设成了全基因集，之后子集到了 HVG —— 直接把 HVG 子集喂给模型会丢掉
# 模型依赖的大部分基因，标签会退化成噪声。所以这里显式用 `.raw`。
CELLTYPIST_DEFAULT_MODEL = "Immune_All_Low.pkl"


def try_celltypist(adata, cfg: dict, log=log_info):
    """用 CellTypist 预训练模型给细胞打标签。

    返回 `(labels_series | None, info)`。**任何异常都吞掉并写进 info** ——
    拿不到模型不该让整步失败，但必须让人看见。
    """
    info = {"attempted": True, "status": None, "reason": "",
            "model": (cfg.get("analysis") or {}).get("celltypist_model",
                                                     CELLTYPIST_DEFAULT_MODEL)}
    try:
        import celltypist
        from celltypist import models as ctm
    except Exception as exc:  # noqa: BLE001
        info["status"] = "package_missing"
        info["reason"] = f"celltypist 未安装（{type(exc).__name__}: {exc}）"
        log(f"CellTypist 不可用：{info['reason']}")
        return None, info

    # **版本要拿得到，不能写 "unknown"。** `celltypist.__version__` 不一定存在
    # （它没在 `__init__.py` 里保证导出），而"用的是哪个版本的模型/代码"
    # 是复现的前提 —— 这正是模块零 §0.3 要记的东西。所以退回发行版元数据。
    info["version"] = getattr(celltypist, "__version__", None)
    if not info["version"]:
        try:
            from importlib import metadata as _md
            info["version"] = _md.version("celltypist")
        except Exception:  # noqa: BLE001
            info["version"] = "unknown"

    # ---- 模型文件：先看本地有没有，没有才下载 --------------------------------
    #
    # **`models_path` 是 `str`，不是 `Path`。** celltypist 1.7.1 的
    # `models.py:19` 是 `models_path = os.path.join(data_path, "models")`，
    # 对它做 `/` 直接抛
    # `TypeError: unsupported operand type(s) for /: 'str' and 'str'`。
    #
    # 实测（run 35486399043）这个 TypeError 被下面那个笼统的
    # `except Exception` 接住，于是记成了
    # **「模型拿不到 —— CI 可能无外网」** —— 一个纯本地代码/API 版本 bug
    # 被写成了网络问题。**归因错了比报错更糟**：下一个人会去查 runner 的
    # 出网策略、换镜像、加超时，而真正要改的是这一行。
    #
    # 所以这里把失败拆成三类，各自记自己的原因：
    #   1. 路径构造失败（本地代码 / celltypist API 变了）→ status="failed"
    #   2. 下载抛异常（网络）                            → "model_unavailable"
    #   3. 下载没抛异常但文件仍不在（名字不在清单里）     → "model_unavailable"
    #
    # 第 3 类必须单独查：`download_models` 内部把每个模型的下载异常
    # **吞掉只打日志**（models.py:512-517），所以"下载失败"并不总是抛出来。
    from pathlib import Path as _Path

    try:
        models_dir = _Path(ctm.models_path)      # str -> Path，见上面的说明
    except Exception as exc:  # noqa: BLE001
        info["status"] = "failed"
        info["reason"] = (f"celltypist.models_path 解析失败"
                          f"（{type(exc).__name__}: {exc}）—— 这是 celltypist "
                          f"API 与调用方不匹配，**不是网络问题**")
        log_warn(f"CellTypist 不可用：{info['reason']}")
        return None, info

    info["models_dir"] = str(models_dir)
    model_file = models_dir / info["model"]
    if model_file.exists():
        info["downloaded"] = False
    else:
        log(f"CellTypist 模型 {info['model']} 不在本地，尝试下载…")
        try:
            ctm.download_models(model=info["model"])
            info["downloaded"] = True
        except Exception as exc:  # noqa: BLE001
            info["status"] = "model_unavailable"
            info["reason"] = (f"模型 {info['model']} 下载失败"
                              f"（{type(exc).__name__}: {exc}）—— 模型服务器 "
                              f"celltypist.cog.sanger.ac.uk 不可达或超时")
            log_warn(f"CellTypist 不可用：{info['reason']}")
            return None, info
        if not model_file.exists():
            info["status"] = "model_unavailable"
            info["reason"] = (f"download_models 未抛异常，但 {model_file} 仍不存在 —— "
                              f"该名字可能不在 celltypist 的模型清单里"
                              f"（模型清单见 https://celltypist.cog.sanger.ac.uk/models/models.json）")
            log_warn(f"CellTypist 不可用：{info['reason']}")
            return None, info
    info["model_bytes"] = int(model_file.stat().st_size)

    try:
        # 用全基因集，不是 HVG 子集 —— 见上面的说明
        full = adata.raw.to_adata() if adata.raw is not None else adata.copy()
        if adata.raw is None:
            info["used_raw"] = False
            log_warn("adata.raw 为空 —— CellTypist 只能用当前（可能是 HVG）矩阵，"
                     "标签可靠性下降")
        else:
            info["used_raw"] = True
        info["n_genes_input"] = int(full.n_vars)

        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pred = celltypist.annotate(full, model=info["model"],
                                       majority_voting=True)
        labels = pred.predicted_labels
        col = ("majority_voting" if "majority_voting" in labels.columns
               else "predicted_labels")
        info["status"] = "ok"
        info["label_column"] = col
        info["n_labels"] = int(labels[col].nunique())
        log(f"CellTypist 完成：{info['n_labels']} 种标签"
            f"（模型 {info['model']} v{info['version']}，"
            f"{info['n_genes_input']} 基因输入）")
        return labels[col].astype(str), info
    except Exception as exc:  # noqa: BLE001
        info["status"] = "failed"
        info["reason"] = f"{type(exc).__name__}: {exc}"
        log_warn(f"CellTypist 跑失败：{info['reason']}")
        return None, info


def compare_annotations(assign: pd.DataFrame, ct_labels, adata, log=log_info) -> dict:
    """量化 marker 打分与 CellTypist 在**簇层面**的一致程度。

    簇层面的比较才有意义：marker 打分给的是每个簇一个标签，而
    CellTypist 给的是每个细胞一个标签。先按簇取众数再比。
    """
    out = {"compared": False}
    if assign is None or ct_labels is None:
        return out
    try:
        df = pd.DataFrame({
            "cluster": adata.obs["leiden"].astype(str).values,
            "celltypist": ct_labels.reindex(adata.obs_names).values,
        }).dropna()
        # 每簇的众数标签 + 该标签占比（占比低说明这个簇本身不纯）
        maj = (df.groupby("cluster")["celltypist"]
                 .agg(lambda s: s.value_counts().index[0]))
        purity = (df.groupby("cluster")["celltypist"]
                    .agg(lambda s: float(s.value_counts().iloc[0] / len(s))))
        own = dict(zip(assign["cluster"].astype(str), assign["assigned"]))
        common = sorted(set(own) & set(maj.index))
        agree = [c for c in common if str(own[c]).lower() == str(maj[c]).lower()]
        out.update({
            "compared": True,
            "n_clusters": len(common),
            "n_agree_exact": len(agree),
            "agreement_frac": round(len(agree) / len(common), 3) if common else None,
            "per_cluster": [
                {"cluster": c, "marker_signature": str(own[c]),
                 "celltypist_majority": str(maj[c]),
                 "celltypist_purity": round(float(purity[c]), 3),
                 "agree": str(own[c]).lower() == str(maj[c]).lower()}
                for c in common
            ],
        })
        log(f"两种注释在簇层面一致 {len(agree)}/{len(common)}"
            f"（{out['agreement_frac']}）—— 不一致的簇值得人工看")
    except Exception as exc:  # noqa: BLE001
        out["reason"] = f"对比失败：{type(exc).__name__}: {exc}"
        log_warn(out["reason"])
    return out


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

    fig, ax = plt.subplots(figsize=(W_SINGLE, mm(58)))
    ax.plot(scan["resolution"], scan["n_clusters"], "o-", color=PAL["primary"])
    ax.axvline(float(rd["resolution"]), color=PAL["highlight"], ls="--", lw=1,
               label=f"used: {rd['resolution']}")
    ax.set_xlabel("Leiden resolution"); ax.set_ylabel("number of clusters")
    ax.set_title("Cluster count vs resolution")
    fig.legend(fontsize=8, ncol=1, loc="outside right center")
    save_fig(cfg, "02-03-01-unit1-cluster-resolution-scan", fig)

    # ---- 3. 用配置的分辨率定稿 ----------------------------------------------
    res_used = float(rd["resolution"])
    sc.tl.leiden(adata, resolution=res_used, key_added="leiden",
                 flavor="igraph", n_iterations=2, directed=False,
                 random_state=cfg["analysis"]["seed"])
    n_clusters = int(adata.obs["leiden"].nunique())
    log_info(f"最终聚类: {n_clusters} 个簇（resolution={res_used}）")

    fig, ax = plt.subplots(figsize=(W_ONE_HALF, mm(84)))
    xy = adata.obsm["X_umap"]
    cats = adata.obs["leiden"].astype(str).values
    for c in sorted(set(cats), key=lambda x: int(x) if x.isdigit() else x):
        m = cats == c
        ax.scatter(xy[m, 0], xy[m, 1], s=4, alpha=0.75, label=c)
        cx, cy = xy[m, 0].mean(), xy[m, 1].mean()
        # 直接标注簇号 = 参考规范要的"颜色之外的第二条线索"
        ax.text(cx, cy, c, weight="bold",
                ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75))
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
    ax.set_title(f"Leiden clusters (n={n_clusters}, resolution={res_used})")
    fig.legend(fontsize=6, markerscale=2.5, ncol=1, loc="outside right center")
    save_fig(cfg, "02-03-02-unit1-umap-clusters", fig)

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
    # **上限按图宽算，不是随手取 40。** 双栏 183 mm 下每个基因约 7 mm，
    # 再多标签就挤成一片。原来取 40 会画出 323 mm 宽的图 —— 装不进任何
    # 期刊的一页。
    top3 = [g for g in top3 if g in adata.raw.var_names][:24]
    if top3:
        # **手工画 dotplot，不用 `sc.pl.dotplot`。** 用户 2026-09-24 第四轮
        # 指出四个问题，前两个是根因级的：
        # ① **标度矛盾** —— 副标题写 z-score，但 scanpy 的
        #    `standard_scale="var"` 实际是**逐基因 min-max 归一化到 0-1**
        #    （读 `_prepare_dot_data` 源码确认），图例 0-1 与"z-score"矛盾。
        #    z-score 均值应为 0、有正有负，不可能全为正。
        #    → 改为**真 z-score**（按基因跨簇标准化），RdBu_r 对称色标：
        #      0 = 簇间平均水平、红 = 高于均值、蓝 = 低于均值。
        # ② **基因名被裁成 C100A8** —— 数据里是 S100A8（markers_all.csv 可查），
        #    左缘裁切把 S 的左半吃掉了。→ 自绘布局，标签区显式留宽。
        # ③ 副标题截断、Y 轴没标 Cluster → 自绘布局 + 显式 ylabel。
        # ④ 图例圆点粘连 → 大小图例单独轴、间距显式给定。
        # 用 **log 化的全基因集**（adata.raw，AGENTS 规则 2）：marker 多为低表达
        raw = adata.raw.to_adata() if adata.raw is not None else adata
        sub = raw[:, top3]
        X = np.asarray(sub.X.todense()) if hasattr(sub.X, "todense") else np.asarray(sub.X)
        groups = adata.obs["leiden"].astype(str).values
        ug = sorted(set(groups), key=lambda v: int(v))
        # frac = 表达细胞比例（>0 计表达，与 scanpy 默认 expression_cutoff 一致）
        frac = np.zeros((len(ug), len(top3)))
        mean_expr = np.zeros((len(ug), len(top3)))
        for ri in range(len(ug)):
            blk = X[groups == ug[ri]]
            frac[ri] = (blk > 0).mean(axis=0)
            mean_expr[ri] = blk.mean(axis=0)
        # **按基因做 z-score**（均值 0、方差 1）；常数列（全同值）记 0
        sd = mean_expr.std(axis=0, ddof=0)
        sd[sd == 0] = 1.0
        zmat = (mean_expr - mean_expr.mean(axis=0)) / sd

        frac_df = pd.DataFrame(frac, index=ug, columns=top3)
        z_df = pd.DataFrame(zmat, index=ug, columns=top3)

        # **图例放主图下方的横带（Seurat do_DotPlot 范式，用户第七/八轮反馈）。**
        # 整图构建抽在 `common.build_marker_dotplot_figure` —— 前七轮把这段
        # 内联在脚本里、验证脚本又照抄一份，三份镜像不同步，导致"改了没效果"
        # 与七轮返工（详见该函数 docstring）。
        fig, size_handles = build_marker_dotplot_figure(
            frac_df, z_df,
            group_label="Leiden cluster",
            title="Top markers per cluster",
            subtitle="rows = Leiden clusters (identities: celltype_annotation.csv)\ndot size = fraction of cells expressing the gene")

        save_fig(cfg, "02-03-03-unit1-markers-dotplot", fig)
        # **灰底 marker UMAP 网格**（差距清单 #19，文献范式）：
        # dotplot 给"哪个簇高表达"，灰底 UMAP 给"高表达在哪块区域" ——
        # 空间/连续区域上的特异性一眼可见（灰底=不表达）。
        # 单图原则：每个基因一张（02-03-03 unit2..unit9）。
        DYNAMIC_FIG_BASES = {"03": 8}
        FIG_BASE = "-".join(["02", "03", "03", "unit"])
        xy = adata.obsm["X_umap"]
        for gi_, g in enumerate(top3[:8], start=2):
            try:
                expr = np.asarray(adata.raw[:, g].X.todense()).ravel()
            except Exception:  # noqa: BLE001
                expr = np.asarray(adata.raw[:, g].X).ravel()
            fig_m, ax_m = plt.subplots(figsize=(W_SINGLE, mm(58)))
            ax_m.scatter(xy[:, 0], xy[:, 1], c=np.ones(len(expr)),
                         s=3, cmap="Greys", vmin=0, vmax=2, alpha=0.25)
            sm = ax_m.scatter(xy[expr > 0, 0], xy[expr > 0, 1],
                              c=expr[expr > 0], s=4, cmap="viridis")
            fig_m.colorbar(sm, ax=ax_m, shrink=0.8, pad=0.02,
                           fraction=0.046, label="expression")
            ax_m.set_title(f"{g} on UMAP (grey = not detected)")
            ax_m.set_xlabel("UMAP1"); ax_m.set_ylabel("UMAP2")
            save_fig(cfg, FIG_BASE + str(gi_) + "-marker-" + g.lower(), fig_m)

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
        # 宽度夹在 [单栏半, 双栏]：类型少时不至于太空，类型多时也不会
        # 画出装不进一页的图（标签已旋转 45°）
        fig, ax = plt.subplots(figsize=(min(W_DOUBLE, max(W_ONE_HALF, 0.45 * len(per_cell.columns) + 3)),
                                        max(3.0, 0.32 * len(per_cell) + 1.6)))
        im = ax.imshow(per_cell.values, aspect="auto", cmap="viridis")
        ax.set_xticks(range(len(per_cell.columns)))
        ax.set_xticklabels(per_cell.columns, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(per_cell)))
        ax.set_yticklabels([str(i) for i in per_cell.index], fontsize=7)
        ax.set_xlabel("cell type signature"); ax.set_ylabel("cluster")
        ax.set_title("Mean signature score per cluster")
        fig.colorbar(im, ax=ax, label="score")
        save_fig(cfg, "02-03-04-unit1-celltype-scores-heatmap", fig)

    # ---- 6. CellTypist 自动注释（文档 §2.4，与 marker 打分并行）------------
    ct_labels, ct_info = try_celltypist(adata, cfg)
    ct_cmp = {"compared": False}
    if ct_labels is not None:
        pd.DataFrame({"cell": ct_labels.index, "celltypist": ct_labels.values}) \
            .to_csv(res_dir / "celltypist_labels.csv", index=False)
        ct_cmp = compare_annotations(assign, ct_labels, adata)
        if ct_cmp.get("compared"):
            ct_cmp["interpretation"] = (
                "marker 签名打分只看签名基因在簇里的平均水平；CellTypist 用"
                "预训练模型看全部基因的加权组合 —— **两者是独立证据**。"
                "一致的簇可信度更高；不一致的簇不应只报其中一个。")
    annot_record["celltypist"] = ct_info
    annot_record["celltypist_vs_marker"] = ct_cmp

    # ---- 7. 落盘 ------------------------------------------------------------
    out = data_dir / "clustered.h5ad"
    adata.write_h5ad(out)
    log_info(f"已写出 {out}")

    # ---- §0.2 跨部分交接的**产出侧**：这一份 h5ad 就是 Part 3 的参考 -------
    #
    # Part 3（空间转录组）的 `deconvolution.reference: h5ad` 要的就是
    # **这一份文件**：带细胞类型标签的单细胞 h5ad。契约有三条：
    #
    #   1. `layers['counts']` 必须是**原始计数** —— Part 3 用它的
    #      类型均值当参考谱，而 NNLS 解的是线性混合，log 值会破坏线性；
    #   2. `obs` 里要有细胞类型列（本仓库写的是 `celltype`），
    #      Part 3 通过 `celltype_key` 指名要哪一列；
    #   3. 基因集是**交集**：Part 3 只用两边共同基因，所以参考若是
    #      HVG 子集，参考谱覆盖的基因就跟着变少。
    #
    # **前两条不满足时 Part 3 会静默算错或报 KeyError，第三条会被静默
    # 接受。** 所以这里把契约状态**主动记下来**（而不是等 Part 3 去发现），
    # 验收也据此检查 —— 交接的两端都要能被核对。
    _ct_col = "celltype" if "celltype" in adata.obs.columns else None
    _has_counts = "counts" in adata.layers
    _n_hvg = int(adata.n_vars)
    part3_ref = {
        "path": str(out),
        "contract": {
            "layers['counts']": _has_counts,
            "celltype_column": _ct_col,
        },
        "n_cells": int(adata.n_obs),
        "n_genes": _n_hvg,
        "counts_layer_source": ("02_integrate 存的原始计数（HVG 子集）"
                                if _has_counts else None),
        "how_part3_uses_it": (
            "Part 3 的 `05_deconvolution.py` 在 `reference: h5ad` 下读"
            "`layers['counts']` + `celltype_key`，按类型求均值得到参考谱 S，"
            "再用 NNLS 解每个 spot 的组成；它会把这次交接记进"
            "`run_manifest.json` 的 `cross_language`（含丢失字段）"),
        "limitations": ([
            "**基因集是 HVG 子集** —— 本文件来自 `integrated.h5ad`"
            f"（{_n_hvg} 个高变基因），不是全基因集。Part 3 的参考谱"
            "因此只覆盖两边的共同基因，共同基因太少时 Part 3 会直接报错",
        ] if _n_hvg < 10000 else []) + ([
            "**没有 counts 层** —— Part 3 会退回用 `.X`，而 `.X` 是 log 后的值，"
            "解出的比例没有意义（Part 3 会 WARN 并记进 lost_fields）",
        ] if not _has_counts else []) + ([
            "**obs 里没有细胞类型列** —— Part 3 配 `celltype_key` 时会 KeyError",
        ] if _ct_col is None else []),
    }
    if not _has_counts or _ct_col is None:
        log_warn(f"Part 3 参考契约不完整：counts 层={_has_counts}，"
                 f"细胞类型列={_ct_col} —— 见 cluster_status.json 的 part3_reference")

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
        # §0.2 的产出侧契约（Part 2 → Part 3）
        "part3_reference": part3_ref,
        "status": "ok",
    }
    write_json(res_dir / "cluster_status.json", status)
    return status


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    t0 = time.time()
    try:
        res = run_03_cluster_annotate(cfg)
        # E-56：接住返回值。本步骤 `not_possible`（打分不可行）是早退路径。
        record_step(cfg, "cluster_annotate", "ok", time.time() - t0,
                    result_status=result_status_of(res))
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "cluster_annotate", "failed", time.time() - t0,
                    message=str(e), result_status="failed")
        raise
