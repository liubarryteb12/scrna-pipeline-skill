#!/usr/bin/env python3
"""
06_communication.py — 细胞间通讯（配体-受体）

**没有 LIANA/CellChat 时怎么办。** 那两个包是这类分析的标准工具，但它们
（以及底层的数据库）依赖较重、版本敏感。这里做的是**数据库驱动的
配体-受体共表达打分**，并**明确标注它不是 LIANA/CellChat**。

做法：
  1. 用内置的配体-受体对（assets/ligand_receptor.yml）
  2. 对每对 (发送簇, 接收簇)：配体在发送簇的平均表达 x 受体在接收簇的平均表达
  3. 用置换检验给出经验 p 值（打乱簇标签）

**必须说清楚的局限：**
  - 共表达不等于通讯。没有空间信息时，两个细胞类型"能通讯"只是说
    它们分别表达了配体和受体，不代表它们在组织里相邻。
  - 表达量是稳态丰度，不等于蛋白水平，也不等于分泌量。
  - 打分是启发式，不是 LIANA 的 consensus rank aggregate。
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
                    log_warn, parse_args, record_step, save_fig, set_seed, write_json, W_ONE_HALF,)

N_PERMUTATIONS = 200


def load_lr_pairs() -> list:
    p = Path(__file__).resolve().parent.parent / "assets" / "ligand_receptor.yml"
    with open(p, encoding="utf-8") as fh:
        doc = yaml.safe_load(fh)
    return doc.get("pairs", []) or []


def get_full_expression(adata):
    """
    取**全基因集**的表达矩阵，而不是 HVG 子集。

    **这是本步骤最容易踩的坑，实测踩过。** 上游为了降维把数据子集到了
    2000 个高变基因，而配体/受体基因大多是低表达的（细胞因子、趋化因子），
    几乎不会进 HVG。实测：用 HVG 子集时 38 对里只有 **3 对**可用；
    用全基因集时是 **27 对** —— 差 9 倍。

    后果是"没找到显著通讯"变成**假阴性**，而假阴性看起来和"真的没有通讯"
    一模一样。所以这里强制用 `.raw`（全基因集），并在缺失时明确报错
    而不是退回 HVG 悄悄继续。
    """
    if adata.raw is None:
        raise RuntimeError(
            "adata.raw 为空 —— 无法取全基因集表达。\n"
            "  为什么必须停: 配体/受体基因绝大多数不是高变基因，"
            "在 HVG 子集上做通讯分析会得到大量假阴性（实测 3/38 vs 27/38）。\n"
            "  修法: 确认 02_integrate.py 在子集到 HVG 之前设了 adata.raw。")
    src = adata.raw
    X = src.X
    X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    return X, list(src.var_names)


def mean_expression(genes: list, mask, X, var_names: list) -> float:
    """给定基因集与细胞掩码，返回平均表达（log 后）。"""
    lookup = {g: i for i, g in enumerate(var_names)}
    idx = [lookup[g] for g in genes if g in lookup]
    if not idx:
        return 0.0
    sub = X[np.ix_(mask, idx)]
    return float(sub.mean()) if sub.size else 0.0


def run_06_communication(cfg: dict) -> dict:
    ensure_dirs(cfg)
    set_seed(cfg)
    data_dir = Path(cfg["output"]["data_dir"])
    res_dir = Path(cfg["output"]["results_dir"])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    comm = cfg.get("communication") or {}
    if not comm.get("enabled", True):
        status = {"dataset_id": cfg["dataset_id"], "status": "disabled",
                  "reason": "配置 communication.enabled=false"}
        write_json(res_dir / "communication_status.json", status)
        log_info("细胞通讯分析已按配置关闭")
        return status

    adata = sc.read_h5ad(data_dir / "clustered.h5ad")
    group_key = "celltype" if "celltype" in adata.obs.columns else "leiden"
    groups = sorted(adata.obs[group_key].astype(str).unique())
    if comm.get("clusters"):
        groups = [g for g in groups if g in [str(c) for c in comm["clusters"]]]
    if len(groups) < 2:
        status = {"dataset_id": cfg["dataset_id"], "status": "not_applicable",
                  "reason": f"只有 {len(groups)} 个簇/细胞类型，通讯分析需要至少 2 个"}
        write_json(res_dir / "communication_status.json", status)
        log_warn(status["reason"])
        return status

    pairs = load_lr_pairs()
    # **用全基因集，不是 HVG 子集** —— 见 get_full_expression 的说明。
    X, var_names = get_full_expression(adata)
    lookup = set(var_names)
    log_info(f"表达矩阵（全基因集）: {X.shape[0]} 细胞 x {X.shape[1]} 基因")

    present = []
    for pr in pairs:
        lig = [g for g in pr.get("ligand", []) if g in lookup]
        rec = [g for g in pr.get("receptor", []) if g in lookup]
        if lig and rec:
            present.append({"name": pr["name"], "ligand": lig, "receptor": rec,
                            "ligand_missing": [g for g in pr["ligand"] if g not in lookup],
                            "receptor_missing": [g for g in pr["receptor"] if g not in lookup]})
    log_info(f"配体-受体对: {len(present)}/{len(pairs)} 的基因在数据里存在")

    if not present:
        status = {"dataset_id": cfg["dataset_id"], "status": "no_pairs_in_data",
                  "reason": f"数据库里 {len(pairs)} 对，没有一对的配体与受体基因同时出现在数据里"}
        write_json(res_dir / "communication_status.json", status)
        log_warn(status["reason"])
        return status

    masks = {g: (adata.obs[group_key].astype(str) == g).values for g in groups}
    rng = np.random.default_rng(cfg["analysis"]["seed"])
    n_cells = adata.n_obs

    # **只抽取涉及的配体/受体基因列，不要在 13714 基因的全矩阵上反复索引。**
    # 未优化时 200 次置换 x 673 组合要跑 4.5 分钟；抽成 ~100 列的小矩阵后
    # 降到几秒。置换检验的统计含义完全不变 —— 只是不再搬运用不到的列。
    needed = sorted({g for pr in present for g in pr["ligand"] + pr["receptor"]})
    col_of = {g: i for i, g in enumerate(needed)}
    col_idx = [var_names.index(g) for g in needed]
    small = np.ascontiguousarray(X[:, col_idx])
    log_info(f"置换用子矩阵: {small.shape[0]} 细胞 x {small.shape[1]} 个配体/受体基因")

    def score(lig_genes, rec_genes, s_mask, d_mask) -> float:
        li = [col_of[g] for g in lig_genes if g in col_of]
        ri = [col_of[g] for g in rec_genes if g in col_of]
        if not li or not ri:
            return 0.0
        a = small[np.ix_(s_mask, li)].mean()
        b = small[np.ix_(d_mask, ri)].mean()
        return float(a * b)

    rows = []
    for pr in present:
        for src in groups:
            for dst in groups:
                if src == dst:
                    continue
                obs_score = score(pr["ligand"], pr["receptor"], masks[src], masks[dst])
                if obs_score <= 0:
                    continue
                # 置换：打乱细胞标签，看这个分数有多容易随机出现
                n_s, n_d = int(masks[src].sum()), int(masks[dst].sum())
                null = np.empty(N_PERMUTATIONS)
                for i in range(N_PERMUTATIONS):
                    perm = rng.permutation(n_cells)
                    null[i] = score(pr["ligand"], pr["receptor"],
                                    perm[:n_s], perm[n_s:n_s + n_d])
                p = float((np.sum(null >= obs_score) + 1) / (N_PERMUTATIONS + 1))
                rows.append({
                    "pair": pr["name"], "sender": src, "receiver": dst,
                    "ligand": ",".join(pr["ligand"]), "receptor": ",".join(pr["receptor"]),
                    "score": round(obs_score, 5), "null_mean": round(float(null.mean()), 5),
                    "p_value": round(p, 5), "n_permutations": N_PERMUTATIONS,
                })

    if not rows:
        status = {"dataset_id": cfg["dataset_id"], "status": "no_signal",
                  "reason": "所有配体-受体对在所有簇对上的打分都为 0"}
        write_json(res_dir / "communication_status.json", status)
        log_warn(status["reason"])
        return status

    res = pd.DataFrame(rows).sort_values("score", ascending=False)
    res.to_csv(res_dir / "cell_communication.csv", index=False)

    # 多重检验：这里报 BH 校正后的值。**不校正的话几百个组合里
    # 一定有一堆 p<0.05，而那只是组合数多。**
    from statsmodels.stats.multitest import multipletests
    res["p_adj_bh"] = multipletests(res["p_value"], method="fdr_bh")[1].round(5)
    res.to_csv(res_dir / "cell_communication.csv", index=False)
    n_sig = int((res["p_adj_bh"] < 0.05).sum())
    log_info(f"通讯打分: {len(res)} 个 (配体受体, 发送, 接收) 组合，"
             f"BH 校正后 {n_sig} 个 p<0.05")

    # 热图：发送 x 接收 的总分
    top = res.head(min(20, len(res)))
    fig, ax = plt.subplots(figsize=(max(W_ONE_HALF, 0.55 * len(groups) + 2.2),
                                    max(3.4, 0.42 * len(top) + 1.8)))
    mat = top.pivot_table(index="pair", columns="receiver", values="score",
                          aggfunc="sum").fillna(0.0)
    im = ax.imshow(mat.values, aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(mat.columns)))
    ax.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(mat)))
    ax.set_yticklabels(mat.index, fontsize=7)
    ax.set_xlabel("receiver"); ax.set_ylabel("ligand-receptor pair")
    ax.set_title(f"Top {len(top)} LR pairs by score")
    fig.colorbar(im, ax=ax, label="score")
    save_fig(cfg, "communication_heatmap", fig)

    status = {
        "dataset_id": cfg["dataset_id"],
        "status": "ok",
        "group_key": group_key,
        "n_groups": len(groups),
        "n_pairs_in_database": len(pairs),
        "n_pairs_usable": len(present),
        "n_genes_used": int(X.shape[1]),
        "gene_set_note": ("用的是**全基因集**（adata.raw），不是 HVG 子集。"
                          "配体/受体基因大多是低表达的细胞因子/趋化因子，"
                          "几乎不进 HVG —— 实测在 HVG 子集上 38 对里只有 3 对可用，"
                          "全基因集上是 27 对，差 9 倍"),
        "n_combinations": int(len(res)),
        "n_significant_bh": n_sig,
        "n_permutations": N_PERMUTATIONS,
        "top_pairs": df_to_records(res.head(15)),
        "method": ("数据库驱动的配体-受体共表达打分 + 簇标签置换检验 + BH 校正。"
                   "**不是 LIANA/CellChat** —— 那两个包未安装"),
        "limitations": [
            "共表达不等于通讯：没有空间信息时，只能说两类细胞分别表达了配体和受体",
            "表达量是稳态丰度，不等于蛋白水平，也不等于分泌量",
            "打分是启发式，不是 LIANA 的 consensus rank aggregate",
            "置换检验打乱的是细胞标签，保留了每种细胞类型的细胞数",
            f"**内置库只有 {len(pairs)} 对**（免疫为主），远少于 CellChatDB 的数千对；"
            "覆盖不全时『没找到显著通讯』是假阴性，不是真的没有通讯",
        ],
    }
    write_json(res_dir / "communication_status.json", status)
    return status


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    t0 = time.time()
    try:
        run_06_communication(cfg)
        record_step(cfg, "communication", "ok", time.time() - t0)
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "communication", "failed", time.time() - t0, message=str(e))
        raise
