#!/usr/bin/env python3
"""
04_pseudobulk_de.py — 拟bulk 差异表达（pydeseq2）

**为什么必须做拟bulk而不是直接在细胞层面跑检验。**
把每个细胞当一个独立样本，是把"细胞数"当成了"样本量"。
10000 个细胞来自 3 个供体，自由度不是 9997 而是 2 —— 而 Wilcoxon
会给出 p = 1e-300 这样的数字。这是单细胞差异分析里最严重的系统性错误，
**而且它给出的假阳性看起来极其显著**。

**样本量的单位是"生物学重复"（供体/病人），不是细胞。**

所以这里：
  1. 按 `design.sample_key` × `design.group_key` × 细胞类型 把计数加起来
  2. 每个 (样本, 细胞类型) 组合是一个拟bulk 样本
  3. 用 pydeseq2 做负二项检验
  4. **细胞数过少的组合要剔除并记录** —— 3 个细胞凑出来的拟bulk 是噪声

没有 `design.sample_key` 时记 `not_configured` 并写明原因，不静默跳过。
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
from scipy import sparse  # noqa: E402

from common import (df_to_records, ensure_dirs, load_config, log_info,  # noqa: E402
                    log_warn, parse_args, record_step, result_status_of,
                    save_fig, set_seed, write_json)

# 一个 (样本, 细胞类型) 组合至少要有这么多细胞才拿去做拟bulk
MIN_CELLS_PER_PSEUDOBULK = 10
# 一个细胞类型至少要在这么多个样本里出现，才能做组间比较
MIN_SAMPLES_PER_CELLTYPE = 2
# 每组至少这么多个生物学重复 —— 1 个重复算不出组内方差
MIN_REPLICATES_PER_GROUP = 2


def aggregate_pseudobulk(adata, sample_key: str, celltype_key: str):
    """
    按 (样本, 细胞类型) 加总原始计数。

    **必须用 counts layer，不能用 X。** X 已经 normalize + log1p 过，
    而 DESeq2 的负二项模型要求原始计数 —— 把 log 值喂进去不会报错，
    只会给出错的离散度估计。

    **没有 counts 层时不再静默退回 `adata.X`**（审计 S6）：原实现
    docstring 写着"必须用 counts"，下一行却 `else adata.X` 且不留任何
    记录 —— 于是「这一轮的离散度估计不可信」这件事在产物里完全看不见。
    现在返回一个 `counts_source` 说明，由调用方判成失败。
    """
    has_counts = "counts" in adata.layers
    if not has_counts:
        # 不在这里抛异常：调用方要把这件事**写进 status** 再决定怎么办，
        # 抛出去会被 `__main__` 的 except 记成步骤崩溃，反而看不出根因。
        return None, None, {
            "reason": ("`clustered.h5ad` 没有 `layers['counts']` —— "
                       "DESeq2 的负二项模型要求原始计数，喂 log 值不会报错"
                       "但离散度估计是错的。**拒绝用 `.X` 冒充计数**；"
                       "请确认上游 `01_qc` / `02_integrate` 保留了 counts 层"),
            "counts_source": "missing",
        }
    X = adata.layers["counts"]
    X = X.toarray() if sparse.issparse(X) else np.asarray(X)

    obs = adata.obs[[sample_key, celltype_key]].copy()
    obs.columns = ["sample", "celltype"]
    obs = obs.astype(str)
    obs["_row"] = np.arange(len(obs))

    mats, meta = {}, []
    for (s, ct), grp in obs.groupby(["sample", "celltype"], observed=True):
        idx = grp["_row"].values
        if len(idx) < MIN_CELLS_PER_PSEUDOBULK:
            continue
        mats[(s, ct)] = X[idx, :].sum(axis=0)
        meta.append({"sample": s, "celltype": ct, "n_cells": int(len(idx))})

    if not mats:
        return None, None, {"reason": f"没有任何 (样本, 细胞类型) 组合达到 "
                                      f"{MIN_CELLS_PER_PSEUDOBULK} 个细胞",
                            "counts_source": "layers['counts']"}
    keys = list(mats)
    mat = np.vstack([mats[k] for k in keys]).astype(np.float64)
    meta = pd.DataFrame(meta).set_index(pd.Index([f"{s}|{ct}" for s, ct in keys]))
    return mat, meta, {"n_pseudobulk_samples": len(keys),
                       "counts_source": "layers['counts']"}


def run_04_pseudobulk_de(cfg: dict) -> dict:
    ensure_dirs(cfg)
    set_seed(cfg)
    data_dir = Path(cfg["output"]["data_dir"])
    res_dir = Path(cfg["output"]["results_dir"])

    design = cfg.get("design") or {}
    sample_key = design.get("sample_key")
    group_key = design.get("group_key")
    ref_group = design.get("reference_group")

    status = {"dataset_id": cfg["dataset_id"], "status": "not_configured"}
    if not sample_key or not group_key:
        status["reason"] = (
            "design.sample_key / design.group_key 为空 —— 拟bulk DE 需要知道"
            "「哪些细胞来自同一个生物学重复」和「怎么分组」。"
            "单细胞数据的样本量单位是供体/病人，不是细胞；"
            "没有这两列就无法构造拟bulk，也不应退化成细胞层面的检验"
            "（那会把细胞数当成样本量，给出 p=1e-300 的假阳性）")
        log_warn(f"跳过拟bulk DE: {status['reason']}")
        write_json(res_dir / "pseudobulk_status.json", status)
        return status

    try:
        import pydeseq2  # noqa: F401
    except ImportError as e:
        status["status"] = "package_missing"
        status["reason"] = f"pydeseq2 未安装: {e}"
        log_warn(status["reason"])
        write_json(res_dir / "pseudobulk_status.json", status)
        return status

    adata = sc.read_h5ad(data_dir / "clustered.h5ad")
    for k, name in ((sample_key, "sample_key"), (group_key, "group_key")):
        if k not in adata.obs.columns:
            status["status"] = "column_missing"
            status["reason"] = (f"design.{name}='{k}' 不在 obs 里；"
                                f"实际列: {', '.join(sorted(adata.obs.columns)[:40])}")
            log_warn(status["reason"])
            write_json(res_dir / "pseudobulk_status.json", status)
            return status

    celltype_key = "celltype" if "celltype" in adata.obs.columns else "leiden"
    log_info(f"拟bulk: sample_key={sample_key}, group_key={group_key}, "
             f"celltype_key={celltype_key}")

    mat, meta, agg = aggregate_pseudobulk(adata, sample_key, celltype_key)
    if mat is None:
        # `counts_source == "missing"` 是**上游契约被破坏**，不是"这批数据
        # 不适合做拟bulk" —— 两者必须长得不一样（审计 S6）。
        if agg.get("counts_source") == "missing":
            status.update({"status": "missing_counts", **agg})
        else:
            status.update({"status": "no_usable_pseudobulk", **agg})
        log_warn(f"拟bulk 不可行: {agg['reason']}")
        write_json(res_dir / "pseudobulk_status.json", status)
        return status

    meta["group"] = [str(adata.obs.loc[adata.obs[sample_key].astype(str) == s, group_key]
                         .astype(str).iloc[0]) for s in meta["sample"]]
    log_info(f"拟bulk 矩阵: {mat.shape[0]} 个 (样本 x 细胞类型) x {mat.shape[1]} 基因")

    # ---- 每个细胞类型单独做 --------------------------------------------------
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    results, skipped = [], []
    contrasts_used = {}
    for ct, sub in meta.groupby("celltype", observed=True):
        if sub["group"].nunique() < 2:
            skipped.append({"celltype": ct, "reason": "只有一个分组取值"})
            continue
        counts = sub["group"].value_counts()
        thin = counts[counts < MIN_REPLICATES_PER_GROUP]
        if len(thin) > 0:
            skipped.append({"celltype": ct,
                            "reason": f"组内生物学重复不足 {MIN_REPLICATES_PER_GROUP}: "
                                      f"{thin.to_dict()}"})
            continue
        if len(sub) < MIN_SAMPLES_PER_CELLTYPE * 2:
            skipped.append({"celltype": ct, "reason": f"拟bulk 样本只有 {len(sub)} 个"})
            continue

        idx = [meta.index.get_loc(i) for i in sub.index]
        cts = pd.DataFrame(mat[idx, :].astype(int),
                           index=sub.index,
                           columns=adata.var_names)
        # DESeq2 不接受全零基因
        cts = cts.loc[:, cts.sum(axis=0) > 0]

        # ---- contrast 的分子/分母必须**显式定死**，不能靠行序（审计 S5）----
        #
        # 原实现 `str(sub["group"].unique()[0])` 取的是「第一次出现的取值」，
        # 而 `sub` 的顺序来自 `meta.groupby(...)` —— 也就是 obs 的行序。
        # **重排细胞顺序会让全部 log2FoldChange 变号，而 padj 一个都不变**
        # （对比方向翻转是符号对称的），产物看起来完全正常。
        #
        # 所以这里把分组水平按**排序后的字典序**定死（可复现、与行序无关），
        # 再把 numerator / reference 都记进 status 供事后核对。
        levels = sorted(sub["group"].astype(str).unique())
        if ref_group is not None and str(ref_group) not in levels:
            # 指定的参考组在这个细胞类型里不存在 —— DESeq2 会抛错，
            # 而原实现的 `except` 会把它吞成 skipped（审计 M20）。
            # 这里显式记成"跳过 + 原因"，与"算失败"区分开。
            skipped.append({
                "celltype": ct,
                "reason": (f"指定的 reference_group='{ref_group}' 在该细胞类型里"
                           f"不存在（实际水平: {levels}）")})
            log_warn(f"  {ct} 跳过：reference_group='{ref_group}' 不在 {levels}")
            continue
        numerator = levels[-1] if ref_group is None else [
            lv for lv in levels if lv != str(ref_group)][0]
        reference = levels[0] if ref_group is None else str(ref_group)
        coldata = pd.DataFrame({"group": sub["group"].astype(str).values},
                               index=sub.index)

        try:
            dds = DeseqDataSet(counts=cts, metadata=coldata, design="~group",
                               refit_cooks=True, quiet=True)
            dds.deseq2()
            stat = DeseqStats(dds, contrast=["group", numerator, reference],
                              quiet=True)
            stat.summary()
            df = stat.results_df.reset_index().rename(columns={"index": "gene"})
            df.insert(0, "celltype", ct)
            results.append(df)
            contrasts_used[ct] = {"numerator": numerator, "reference": reference,
                                  "levels": levels}
            log_info(f"  {ct}: {len(sub)} 个拟bulk 样本, "
                     f"contrast={numerator} vs {reference}, "
                     f"{(df['padj'] < 0.05).sum() if 'padj' in df else 0} 个 padj<0.05")
        except Exception as e:  # noqa: BLE001
            skipped.append({"celltype": ct, "reason": f"{type(e).__name__}: {e}"})
            log_warn(f"  {ct} 的 DESeq2 失败: {type(e).__name__}: {e}")

    if results:
        allres = pd.concat(results, ignore_index=True)
        allres.to_csv(res_dir / "pseudobulk_de.csv", index=False)
        meta.to_csv(res_dir / "pseudobulk_samples.csv")
        status["n_genes_tested"] = int(len(allres))
        status["n_celltypes_tested"] = len(results)
        status["status"] = "ok"
    else:
        status["status"] = "no_celltype_testable"

    status.update({
        "sample_key": sample_key, "group_key": group_key,
        "celltype_key": celltype_key, "reference_group": ref_group,
        # 实际用到的对比方向（审计 S5）：没有它就无法事后核对
        # `log2FoldChange` 的符号指的是哪个方向。
        "contrast_used": contrasts_used,
        "contrast_rule": ("分组水平按**字典序**定死；`reference_group` 为空时"
                          "numerator=最大的水平、reference=最小的水平。"
                          "**不看 obs 行序** —— 旧实现取 `unique()[0]`，"
                          "重排细胞会让全部 log2FoldChange 变号而 padj 不变"),
        "counts_source": agg.get("counts_source", "layers['counts']"),
        "min_cells_per_pseudobulk": MIN_CELLS_PER_PSEUDOBULK,
        "min_replicates_per_group": MIN_REPLICATES_PER_GROUP,
        "n_pseudobulk_samples": int(mat.shape[0]),
        "skipped_celltypes": skipped,
        "design_note": ("样本量单位是生物学重复（供体/病人），不是细胞 —— "
                        "把细胞当样本会给出 p=1e-300 的假阳性"),
    })
    write_json(res_dir / "pseudobulk_status.json", status)
    return status


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    t0 = time.time()
    try:
        res = run_04_pseudobulk_de(cfg)
        # E-56：接住返回值。本步骤有 4 条早退路径（`not_configured` /
        # `package_missing` / `column_missing` / `no_usable_pseudobulk` /
        # `no_celltype_testable`），全都不抛异常。
        record_step(cfg, "pseudobulk_de", "ok", time.time() - t0,
                    result_status=result_status_of(res))
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "pseudobulk_de", "failed", time.time() - t0,
                    message=str(e), result_status="failed")
        raise
