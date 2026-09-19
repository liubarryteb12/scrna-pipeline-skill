#!/usr/bin/env python3
"""
tools/selftest_pseudobulk.py — 拟bulk DE 代码路径自检

**这是代码测试，不是分析。** 它验证 DESeq2 那一段能跑通、能产出结果表，
仅此而已。它产出的 p 值**没有生物学含义** —— 分组是按细胞随机切的。

为什么需要它：CI 用的 pbmc3k 是**单样本**数据集，`04_pseudobulk_de.py`
会正确地记 `not_configured` 并跳过。那验证了"正确跳过"，但**没有验证
真正的 DESeq2 代码路径**。而"没跑到的代码"和"跑对了的代码"在 CI 日志里
长得一模一样 —— 这是本仓库反复踩到的同一类问题（Part 1 规则 24）。

所以：用一个明确标注的自检把那段代码跑到。

用法：
    python tools/selftest_pseudobulk.py --config assets/config.pbmc3k.yml
退出码：0 = 代码路径可用；1 = 失败
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts" / "lib"))
sys.path.insert(0, str(REPO / "scripts"))

from common import load_config, log_info, log_warn  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    import numpy as np
    import pandas as pd
    import scanpy as sc

    cfg = load_config(args.config)
    data_dir = Path(cfg["output"]["data_dir"])
    src = data_dir / "clustered.h5ad"
    if not src.exists():
        print(f"跳过：{src} 不存在（先跑 00-03）")
        return 0

    adata = sc.read_h5ad(src)
    log_info(f"自检读入 {adata.n_obs} 细胞 x {adata.n_vars} 基因")

    # 造 6 个"样本" x 2 个"分组"。**随机切分，没有生物学含义。**
    rng = np.random.default_rng(20260919)
    n = adata.n_obs
    sample = rng.choice([f"S{i}" for i in range(6)], size=n)
    group = np.where(np.isin(sample, ["S0", "S1", "S2"]), "A", "B")
    adata.obs["_selftest_sample"] = pd.Categorical(sample)
    adata.obs["_selftest_group"] = pd.Categorical(group)

    # 只留最大的 3 个簇，避免碎簇拉长运行时间
    top = adata.obs["leiden"].value_counts().head(3).index.tolist()
    sub = adata[adata.obs["leiden"].isin(top)].copy()
    log_info(f"自检子集: {sub.n_obs} 细胞，簇 {top}")

    # 直接复用生产代码的聚合函数 —— 自检要测的就是它
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "pb", REPO / "scripts" / "04_pseudobulk_de.py")
    pb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pb)

    mat, meta, agg = pb.aggregate_pseudobulk(sub, "_selftest_sample", "leiden")
    if mat is None:
        log_warn(f"聚合失败: {agg}")
        return 1
    log_info(f"聚合出 {mat.shape[0]} 个拟bulk 样本 x {mat.shape[1]} 基因")

    meta["group"] = [str(sub.obs.loc[sub.obs["_selftest_sample"].astype(str) == s,
                                     "_selftest_group"].astype(str).iloc[0])
                     for s in meta["sample"]]

    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    ct = meta["celltype"].iloc[0]
    sel = meta[meta["celltype"] == ct]
    idx = [meta.index.get_loc(i) for i in sel.index]
    cts = pd.DataFrame(mat[idx, :].astype(int), index=sel.index,
                       columns=sub.var_names)
    cts = cts.loc[:, cts.sum(axis=0) > 0]
    coldata = pd.DataFrame({"group": sel["group"].values}, index=sel.index)

    dds = DeseqDataSet(counts=cts, metadata=coldata, design="~group",
                       refit_cooks=True, quiet=True)
    dds.deseq2()
    stat = DeseqStats(dds, contrast=["group", "A", "B"], quiet=True)
    stat.summary()
    res = stat.results_df
    n_sig = int((res["padj"] < 0.05).sum()) if "padj" in res else 0
    log_info(f"DESeq2 自检通过: {len(res)} 个基因，{n_sig} 个 padj<0.05（无生物学含义）")

    outdir = Path(cfg["output"]["results_dir"]) / "selftest"
    outdir.mkdir(parents=True, exist_ok=True)
    res.reset_index().rename(columns={"index": "gene"}).to_csv(
        outdir / "pseudobulk_selftest_de.csv", index=False)

    import json
    (outdir / "README.txt").write_text(
        "本目录是**代码自检**产物，不是分析结果。\n"
        "分组是按细胞随机切的，p 值没有生物学含义。\n"
        "它的唯一用途是验证 04_pseudobulk_de.py 的 DESeq2 代码路径能跑通。\n",
        encoding="utf-8")
    log_info(f"自检产物写到 {outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
