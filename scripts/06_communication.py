#!/usr/bin/env python3
"""
06_communication.py — 细胞间通讯（配体-受体）

**两条路都跑，都如实标注。**

文档 §2.7 把 LIANA 定为主工具、CellChat 为辅。本仓库的做法：

  1. **LIANA rank_aggregate**（装得上就跑）→ `liana_results.csv`
     LIANA 自带 consensus 资源（CellChatDB + CellPhoneDB 等）。
  2. **自建共表达打分**（始终跑）→ `cell_communication.csv`
     内置的配体-受体对 + 置换检验 + BH 校正。

**为什么并存而不是替换：** 自建打分的产物已被下游和验收引用，直接换成
LIANA 会让"现有结果"无从对比。两个都产出，并**量化两者一致性**
（Spearman + top-N 重叠）—— 不一致本身就是发现，不是谁错了。

**CellChat 是 R 包**，需要 rpy2 + R，与本仓库"云端不装 R"的设计冲突，
所以未使用，理由写进 status 的 limitations。

**必须说清楚的局限：**
  - 共表达不等于通讯。没有空间信息时，两个细胞类型"能通讯"只是说
    它们分别表达了配体和受体，不代表它们在组织里相邻。
  - 表达量是稳态丰度，不等于蛋白水平，也不等于分泌量。
  - 自建打分是启发式，不是 LIANA 的 consensus rank aggregate。
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
                    log_warn, parse_args, record_step, save_fig, set_seed,
                    write_json, W_DOUBLE, W_ONE_HALF,)

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


# ---------------------------------------------------------------------------
# LIANA（文档 §2.7 指定的主工具）
# ---------------------------------------------------------------------------
# **文档把 LIANA 定为主工具、CellChat 为辅。** 本仓库原先只有自建的
# 共表达打分，产物里必须写"不是 LIANA" —— 因为它没有 consensus rank
# aggregate，也没有多方法一致性这一层。
#
# 装得上就跑，装不上/跑不动就退回自建打分，**两条路都如实记录**。
# 不要因为"退回的路也能出图"就把方法写成 LIANA。
#
# **为什么并存而不是替换：** 自建打分的产物（cell_communication.csv）
# 已经被下游和验收引用，直接换成 LIANA 会让"现有结果"无从对比。
# 两个都产出，并**量化两者是否一致** —— 不一致本身就是发现。
LIANA_N_PERMS = 100


def try_liana(adata, group_key: str, cfg: dict, log=log_info):
    """尝试用 LIANA 跑 rank_aggregate。

    返回 `(res_df | None, info)`。**任何异常都吞掉并写进 info** ——
    LIANA 跑不动不该让整步失败，但必须让人看见它跑不动。
    """
    info = {"attempted": True, "status": None, "reason": "",
            "n_perms": LIANA_N_PERMS}
    try:
        import liana as li
    except Exception as exc:  # noqa: BLE001
        info["status"] = "package_missing"
        info["reason"] = f"liana 未安装（{type(exc).__name__}: {exc}）"
        log(f"LIANA 不可用：{info['reason']}")
        return None, info

    info["version"] = getattr(li, "__version__", "unknown")

    # LIANA 在 1.x 里把 `method` 改名成 `mt`。两个都试，别赌版本。
    fn = None
    for path in (("mt", "rank_aggregate"), ("method", "rank_aggregate")):
        obj = li
        try:
            for attr in path:
                obj = getattr(obj, attr)
            fn = obj
            info["api"] = ".".join(path)
            break
        except AttributeError:
            continue
    if fn is None:
        info["status"] = "api_not_found"
        info["reason"] = "liana 里找不到 mt.rank_aggregate 或 method.rank_aggregate"
        log_warn(info["reason"])
        return None, info

    import warnings
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fn(adata, groupby=group_key, use_raw=True,
               n_perms=LIANA_N_PERMS, n_jobs=1, verbose=False,
               seed=cfg["analysis"]["seed"])
        res = adata.uns.get("liana_res")
        if res is None:
            info["status"] = "no_result"
            info["reason"] = "liana 跑完但 uns['liana_res'] 不存在"
            log_warn(info["reason"])
            return None, info
        res = res.copy()
        info["status"] = "ok"
        info["n_rows"] = int(len(res))
        info["columns"] = sorted(map(str, res.columns))
        log(f"LIANA rank_aggregate 完成：{len(res)} 行，版本 {info['version']}")
        return res, info
    except Exception as exc:  # noqa: BLE001
        info["status"] = "failed"
        info["reason"] = f"{type(exc).__name__}: {exc}"
        log_warn(f"LIANA 跑失败，退回自建共表达打分：{info['reason']}")
        return None, info


def compare_with_liana(liana_res, own: pd.DataFrame, log=log_info) -> dict:
    """量化自建打分与 LIANA 的一致性。

    **两种方法给出不同排序是常态，不是 bug** —— LIANA 的 rank aggregate
    综合了 7 种方法的秩，而自建打分只是一个表达量乘积。所以这里报
    Spearman 相关与 top-N 重叠率，让"差多少"变成数字而不是印象。

    返回的字典直接进 status JSON。
    """
    out = {"compared": False}
    if liana_res is None or own is None or len(own) == 0:
        return out
    try:
        # LIANA 的列名在版本间变过，按优先级找"综合分"那一列
        score_col = None
        for c in ("aggregate_rank", "consensus_score", "magnitude_rank",
                  "expr_prod", "lr_means"):
            if c in liana_res.columns:
                score_col = c
                break
        if score_col is None:
            out["reason"] = f"liana_res 里没有可识别的分数列（有 {list(liana_res.columns)[:8]}）"
            return out

        lr = liana_res.copy()
        lr["_pair"] = (lr["ligand_complex"].astype(str) + "^" +
                       lr["receptor_complex"].astype(str))
        lr["_combo"] = (lr["_pair"] + "|" + lr["source"].astype(str) + "->" +
                        lr["target"].astype(str))
        lr_score = lr.groupby("_combo")[score_col].mean()

        ow = own.copy()
        ow["_pair"] = ow["ligand"].astype(str) + "^" + ow["receptor"].astype(str)
        ow["_combo"] = (ow["_pair"] + "|" + ow["sender"].astype(str) + "->" +
                        ow["receiver"].astype(str))
        # 自建打分越大越强；LIANA 的 *_rank 越小越强，要翻向
        ascending = score_col.endswith("_rank")
        ow_score = ow.groupby("_combo")["score"].mean()
        if ascending:
            lr_score = -lr_score

        common = ow_score.index.intersection(lr_score.index)
        out["compared"] = True
        out["score_column_used"] = score_col
        out["rank_direction"] = "越小越强（已翻向）" if ascending else "越大越强"
        out["n_common_combinations"] = int(len(common))
        out["n_own_only"] = int(len(ow_score.index.difference(lr_score.index)))
        out["n_liana_only"] = int(len(lr_score.index.difference(ow_score.index)))
        if len(common) >= 5:
            from scipy.stats import spearmanr
            rho, p = spearmanr(ow_score.loc[common], lr_score.loc[common])
            out["spearman_rho"] = round(float(rho), 4)
            out["spearman_p"] = float(p)
        # top-N 重叠
        for n in (10, 25, 50):
            if len(common) >= n:
                a = set(ow_score.loc[common].nlargest(n).index)
                b = set(lr_score.loc[common].nlargest(n).index)
                out[f"top{n}_overlap"] = len(a & b)
                out[f"top{n}_overlap_frac"] = round(len(a & b) / n, 3)
        log(f"自建 vs LIANA：共同组合 {out['n_common_combinations']}，"
            f"Spearman rho={out.get('spearman_rho')}，"
            f"top25 重叠 {out.get('top25_overlap')}/25")
    except Exception as exc:  # noqa: BLE001
        out["reason"] = f"对比失败：{type(exc).__name__}: {exc}"
        log_warn(out["reason"])
    return out


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

    # 热图：配体-受体对 x 接收细胞类型 的总分
    #
    # **按"配体-受体对"选前 N，不是按三元组选。**
    # `res` 的每一行是 (配体受体对, 发送, 接收) 三元组，所以
    # `res.head(20)` 拿到的是 20 个**组合**；再 pivot 到 pair x receiver，
    # 行数就塌成"这 20 个组合里出现过几个不同的 pair" —— 实测 pbmc3k 只有 **3** 个。
    # 于是图上出现 3 行、标题却写着 "Top 20 LR pairs"，而且画布约 8/8 是空的 0 值区。
    # **判据要匹配图在问的那件事**：这张图的行是 pair，就该按 pair 排序。
    # （同一批数据里 B2M_CD8A 的单个组合分就很高，把另外两个 pair 全挤出去了。）
    pair_score = res.groupby("pair")["score"].sum().sort_values(ascending=False)
    n_pairs = min(20, len(pair_score))
    top_pairs = list(pair_score.head(n_pairs).index)
    top = res[res["pair"].isin(top_pairs)]
    mat = top.pivot_table(index="pair", columns="receiver", values="score",
                          aggfunc="sum").fillna(0.0)
    # 行序按总分从高到低，与 `pair_score` 一致 —— 否则行序是 pivot 的字母序，
    # "按分数取的前 N"这句话在图上就看不出来。
    mat = mat.reindex([p for p in top_pairs if p in mat.index])

    # **高度按实际画出来的行数算，不按请求数算。** 原来用 `len(top)`（=20）算，
    # 而实际只有 3 行 -> 每行 3.4 英寸，格子被拉成巨大的纯色块。
    # 同时**必须夹到 W_DOUBLE**：实测原来算出 259.1 mm，装不进任何期刊的一栏
    # （规则 13）。
    fig_h = min(W_DOUBLE, max(2.6, 0.16 * len(mat) + 1.9))
    fig_w = min(W_DOUBLE, max(W_ONE_HALF, 0.55 * len(groups) + 2.2))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(mat.values, aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(mat.columns)))
    ax.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(mat)))
    ax.set_yticklabels(mat.index, fontsize=7)
    ax.set_xlabel("receiver"); ax.set_ylabel("ligand-receptor pair")
    # 标题写**实际画出来的**行数。原来写的是 `len(top)` = 组合数，
    # 标题对图上内容的描述是错的。
    ax.set_title(f"Top {len(mat)} ligand-receptor pairs by summed score")
    fig.colorbar(im, ax=ax, label="score")
    save_fig(cfg, "communication_heatmap", fig)
    log_info(f"通讯热图: {len(mat)} 个配体受体对 x {mat.shape[1]} 个接收类型"
             f"（{len(res)} 个组合 -> {len(pair_score)} 个不同 pair，取前 {n_pairs}）")

    # ---- LIANA（文档 §2.7 指定的主工具）------------------------------------
    # **自建打分照常产出**（cell_communication.csv 已被下游引用），
    # LIANA 另存一份并量化两者一致性 —— 直接替换会让"现有结果"无从对比。
    liana_res, liana_info = try_liana(adata, group_key, cfg)
    liana_cmp = {"compared": False}
    if liana_res is not None:
        liana_res.to_csv(res_dir / "liana_results.csv", index=False)
        liana_cmp = compare_with_liana(liana_res, res)
        if liana_cmp.get("compared"):
            liana_cmp["interpretation"] = (
                "LIANA 的 rank aggregate 综合了多种方法的秩，自建打分只是"
                "表达量乘积 —— **两者排序不同是预期内的，不是谁错了**。"
                "这里报出来是为了让『差多少』有数字，而不是让读者以为"
                "两个工具在回答同一个问题。")

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
        "liana": liana_info,
        "liana_vs_builtin": liana_cmp,
        "method": ("**两条路都跑了**：(1) 自建数据库驱动的配体-受体共表达打分"
                   "+ 簇标签置换检验 + BH 校正 → cell_communication.csv；"
                   + (f"(2) LIANA rank_aggregate v{liana_info.get('version')} "
                      f"→ liana_results.csv"
                      if liana_info.get("status") == "ok"
                      else f"(2) LIANA **未能运行**（{liana_info.get('status')}："
                           f"{str(liana_info.get('reason'))[:120]}）")),
        "primary_tool_per_spec": "LIANA（文档 §2.7）",
        "limitations": [
            "共表达不等于通讯：没有空间信息时，只能说两类细胞分别表达了配体和受体",
            "表达量是稳态丰度，不等于蛋白水平，也不等于分泌量",
            "自建打分是启发式，不是 LIANA 的 consensus rank aggregate",
            "置换检验打乱的是细胞标签，保留了每种细胞类型的细胞数",
            f"**内置库只有 {len(pairs)} 对**（免疫为主），远少于 CellChatDB 的数千对；"
            "覆盖不全时『没找到显著通讯』是假阴性，不是真的没有通讯",
            ("LIANA 用的是它自带的 consensus 资源（CellChatDB + CellPhoneDB 等），"
             "与内置库不是同一套配体-受体对 —— 两者的组合数不可直接比较"
             if liana_info.get("status") == "ok" else
             "LIANA 未运行，本轮只有自建打分这一条路"),
            "CellChat（文档 §2.7 的辅工具）是 R 包，需要 rpy2 + R —— "
            "本仓库不装 R，故未使用",
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
