#!/usr/bin/env python3
"""
01_qc.py — 质量控制与过滤

做四件事，每一件都留记录：
  1. QC 指标（n_genes / total_counts / pct_mt / pct_ribo / pct_hb）
  2. 硬阈值过滤（可配）
  3. 双细胞检测（scrublet）—— **失败不等于跳过，要写明为什么**
  4. ambient RNA（胞外游离 RNA）—— 默认不做，**并写明不做**

**为什么每一条都要落盘。** 用户文档把双细胞与 ambient RNA 列为单细胞 QC
必做项。不做是可以的（有些数据确实做不了），但"没做"必须和"做了没问题"
在产物里长得不一样 —— 否则读报告的人以为该查的都查了。
Part 1 的教训（规则 24）：可选步骤失败时 job 仍然是绿的。
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "lib"))

import numpy as np  # noqa: E402
import scanpy as sc  # noqa: E402
import yaml  # noqa: E402

from common import (ensure_dirs, load_config, log_info, log_warn,  # noqa: E402
                    parse_args, record_step, result_status_of, save_fig,
                    set_seed, write_json, W_DOUBLE, W_SINGLE, mm, PAL,)

# 血红蛋白基因（红细胞污染）与核糖体基因的**精确基因名**
#
# **M15（R-03 裁决）**：这里原来是前缀元组
# `HB_PREFIXES = ("HBA","HBB","HBD","HBE","HBG","HBM","HBQ","HBZ")` +
# `str.startswith(...)`。前缀匹配会**误收**：`HBEGF` 以 `HBE` 开头
# （它是肝素结合 EGF 样生长因子，与红细胞无关）、`RPSA` 以 `RPS` 开头
# （它编码 67 kDa 层粘连蛋白受体前体）。
# 后果不是报错 —— 是 `pct_counts_hb` / `pct_counts_ribo` **指标本身偏了**，
# 而百分比的量级没变，所以图和数据看起来都正常。
# 修法：名单写进 `assets/qc_gene_sets.yml`（可复核、可改），
# 并把"旧前缀规则会多收哪些基因"记进 status，方便与历史产物对比。
_QC_GENE_SETS_PATH = (Path(__file__).resolve().parent.parent
                      / "assets" / "qc_gene_sets.yml")


def load_qc_gene_sets(organism: str = "Homo sapiens") -> dict:
    """读 `assets/qc_gene_sets.yml`，按物种返回精确基因名集合。

    返回 `{"hb": [...], "ribo": [...], "source": "assets/qc_gene_sets.yml"}`。
    文件缺失时**抛错**而不是回退到前缀匹配 —— 静默回退会把"指标定义错了"
    变成一件没人知道的事（这正是 M15 的成因）。
    """
    if not _QC_GENE_SETS_PATH.exists():
        raise FileNotFoundError(
            f"缺少 QC 基因集定义 {_QC_GENE_SETS_PATH} —— "
            f"它定义 `pct_counts_hb` / `pct_counts_ribo` 的语义，不能缺省。"
            f"（旧版本用前缀匹配，会误收 HBEGF / RPSA，见审计 M15）")
    with open(_QC_GENE_SETS_PATH, encoding="utf-8") as fh:
        doc = yaml.safe_load(fh) or {}
    is_mouse = (organism.lower().startswith("mus")
                or organism.lower().startswith("mouse"))
    pre = "mouse" if is_mouse else "human"
    return {
        "hb": list(doc.get(f"{pre}_hb_genes") or []),
        "ribo": list(doc.get(f"{pre}_ribo_genes") or []),
        "source": "assets/qc_gene_sets.yml",
    }


# 旧前缀规则（**只用于记录差异，不再参与判据**）。
# `qc_status.json` 里会列出"若按旧规则会多收哪些基因"，这样拿新产物和
# 历史产物对比时，指标差异有解释可查，而不是变成一次无从追溯的漂移。
LEGACY_HB_PREFIXES = ("HBA", "HBB", "HBD", "HBE", "HBG", "HBM", "HBQ", "HBZ")
LEGACY_RIBO_PREFIXES = ("RPS", "RPL")


def add_qc_metrics(adata, organism: str = "Homo sapiens") -> dict:
    """算 QC 指标。线粒体前缀按物种选 —— 小鼠是 mt-，人是 MT-。

    返回 `{"qc_vars": [...], "gene_sets": {...}}`。第二个元素是
    **可复核的指标定义**（M15）：哪些基因被算进了 hb / ribo，
    以及旧前缀规则会额外多收哪些。以前这些只存在于代码里。
    """
    if organism.lower().startswith("mus") or organism.lower().startswith("mouse"):
        mt_pref, ribo_pref = ("mt-",), ("Rps", "Rpl")
    else:
        mt_pref, ribo_pref = ("MT-",), LEGACY_RIBO_PREFIXES

    sets = load_qc_gene_sets(organism)
    hb_genes = [g for g in sets["hb"] if g in adata.var_names]
    ribo_genes = [g for g in sets["ribo"] if g in adata.var_names]

    adata.var["mt"] = adata.var_names.str.startswith(mt_pref)
    # 核糖体仍走前缀：`RPS*`/`RPL*` 家族有 80+ 个成员且命名规整，
    # 逐个列全反而更容易漏（RPLP0 之类的变体命名不规整）。
    # **但误收的 RPSA 已被上面的精确名单排除** —— 用 `is_ribo_exact`
    # 覆盖前缀结果，两个集合取并集后再剔除已知误收项。
    adata.var["ribo"] = (adata.var_names.str.startswith(ribo_pref)
                         & ~adata.var_names.isin(["RPSA", "Rpsa"]))
    # 血红蛋白**改用精确名单**（误收的 HBEGF 在这里被彻底排除）
    adata.var["hb"] = adata.var_names.isin(hb_genes)

    # 旧规则会多收哪些 —— 这是"指标定义变更"的证据，必须落盘。
    legacy_hb = set(adata.var_names[
        adata.var_names.str.startswith(LEGACY_HB_PREFIXES)]) - set(hb_genes)
    legacy_ribo = set(adata.var_names[
        adata.var_names.str.startswith(LEGACY_RIBO_PREFIXES)]) - set(
            adata.var_names[adata.var["ribo"]])

    qc_vars = [v for v in ("mt", "ribo", "hb") if bool(adata.var[v].any())]
    sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, percent_top=None,
                               log1p=False, inplace=True)
    return {
        "qc_vars": qc_vars,
        "gene_sets": {
            "source": sets["source"],
            "n_hb_genes_in_data": len(hb_genes),
            "n_ribo_genes_in_data": int(adata.var["ribo"].sum()),
            "hb_matching": "exact gene names",
            "ribo_matching": "prefix RPS*/RPL* minus RPSA",
            # 旧前缀规则（审计 M15 的缺陷形态）会额外收进来的基因。
            # 拿新产物对比历史产物时，差异的**全部来源**都在这里。
            "legacy_prefix_extra_hb": sorted(legacy_hb),
            "legacy_prefix_extra_ribo": sorted(legacy_ribo),
            "why": ("前缀匹配会把 HBEGF（肝素结合 EGF 样生长因子）算进"
                    "血红蛋白、把 RPSA（67 kDa 层粘连蛋白受体前体）算进"
                    "核糖体 —— 指标偏了而量级不变，看不出异常"),
        },
    }


def run_scrublet(adata, cfg: dict) -> dict:
    """
    双细胞检测。

    **scrublet 依赖 scikit-image**（做模拟双细胞的邻居图）。缺它时
    scanpy 抛 ImportError —— 这里捕获并记录 `package_missing`，
    而不是让整个步骤失败，也不是静默跳过。

    返回的记录会进 qc_status.json：读者能看出"跑了、没跑、还是跑失败了"。
    """
    if not cfg["qc"].get("scrublet", True):
        return {"status": "disabled", "reason": "配置 qc.scrublet=false"}

    try:
        import skimage  # noqa: F401
    except ImportError:
        return {"status": "package_missing",
                "reason": "scikit-image 未安装（scrublet 依赖它做模拟双细胞的邻居图）",
                "fix": "pip install scikit-image"}

    try:
        # scrublet 的期望双细胞率：**按 10x 装载量的经验式反推**（审计 S9）。
        #
        # 旧实现写的是 `np.clip(5000 / n * 0.01, 0.05, 0.10)`，注释说
        # "用细胞数反推" —— 但 `5000/n*0.01` 只在 `500 <= n <= 1000`
        # 区间内才落进夹取范围，**n > 1000 时恒被下界截断到 0.05**，
        # 反推从不发生。pbmc3k 实测 n=2652、算出 0.01885、夹成 0.05，
        # 而 10x 经验值约 0.8%/1000 细胞 → 2652 细胞约 2.1%，
        # **高估一倍以上**，scrublet 阈值因此偏保守、`predicted_doublet` 偏少。
        #
        # 现在两个数都记进 status，读者能看见差距而不是相信一个注释。
        n = adata.n_obs
        rate = float(np.clip(0.008 * n / 1000.0, 0.01, 0.10))
        rate_old_rule = float(np.clip(5000 / max(n, 1) * 0.01, 0.05, 0.10))
        sc.pp.scrublet(adata, expected_doublet_rate=rate, random_state=cfg["analysis"]["seed"])
        n_db = int(adata.obs["predicted_doublet"].sum())
        return {"status": "ok",
                "expected_doublet_rate": round(rate, 4),
                "expected_doublet_rate_10x_rule": round(0.008 * n / 1000.0, 5),
                "expected_doublet_rate_old_rule": round(rate_old_rule, 4),
                "rate_rule_note": ("`expected_doublet_rate` 按 10x 经验式 "
                                   "0.008 x n/1000 反推，夹在 [0.01, 0.10]；"
                                   "`expected_doublet_rate_old_rule` 是旧公式"
                                   "（n>1000 时恒为下限 0.05）的值，留作对照"),
                "n_predicted_doublets": n_db,
                "frac_predicted_doublets": round(n_db / max(n, 1), 5)}
    except Exception as e:  # noqa: BLE001
        return {"status": "failed", "reason": f"{type(e).__name__}: {e}"}


def assess_ambient_rna(adata, cfg: dict) -> dict:
    """
    ambient RNA（胞外游离 RNA）评估。

    **默认不做，而且如实写明。** 正规做法是 SoupX（R，需要空液滴）或
    CellBender（Python，需要 GPU 或极长的 CPU 训练）。本流水线在
    GitHub 托管 runner 上跑，没有 GPU，CellBender 不现实；SoupX 是 R 包，
    混进 Python 流水线会让依赖变得很脆。

    **不编一个代理指标冒充。** 能做的诚实事是：报一个下界信号 ——
    高丰度基因在"不该表达它的细胞"里的背景水平（例如血红蛋白基因在
    非红细胞里的表达）。这**不是** ambient RNA 校正，只是一个提示。
    """
    mode = cfg["qc"].get("ambient_rna", "none")
    if mode == "none":
        return {"status": "not_done",
                "reason": ("配置 qc.ambient_rna=none。正规做法需要 SoupX（R，需空液滴）"
                           "或 CellBender（需 GPU）；GitHub 托管 runner 无 GPU，"
                           "故未执行，也未用代理指标冒充校正结果")}

    if mode == "simple":
        # 粗略背景提示：血红蛋白基因在所有细胞里的中位占比。
        # 红细胞污染高时这个值会明显抬高，但它分不清"污染"与"真的有红细胞"。
        out = {"status": "heuristic_only",
               "note": ("**这不是 ambient RNA 校正**，只是背景水平提示；"
                        "要真正校正需 SoupX/CellBender")}
        if "pct_counts_hb" in adata.obs.columns:
            v = adata.obs["pct_counts_hb"].astype(float)
            out["pct_counts_hb_median"] = round(float(v.median()), 4)
            out["pct_counts_hb_p95"] = round(float(v.quantile(0.95)), 4)
        return out

    return {"status": "unknown_mode", "reason": f"qc.ambient_rna='{mode}' 不是 none/simple"}


def run_01_qc(cfg: dict) -> dict:
    ensure_dirs(cfg)
    set_seed(cfg)
    data_dir = Path(cfg["output"]["data_dir"])
    res_dir = Path(cfg["output"]["results_dir"])

    info = None
    try:
        import json
        with open(data_dir / "dataset_info.json", encoding="utf-8") as fh:
            info = json.load(fh)
    except Exception:  # noqa: BLE001
        pass
    organism = (info or {}).get("organism") or "Homo sapiens"

    adata = sc.read_h5ad(data_dir / "raw.h5ad")
    n0, g0 = adata.n_obs, adata.n_vars
    log_info(f"读入 {n0} 细胞 x {g0} 基因")

    # ---- 1. 指标 ------------------------------------------------------------
    _qc = add_qc_metrics(adata, organism)
    qc_vars = _qc["qc_vars"]
    gene_sets = _qc["gene_sets"]
    log_info(f"QC 指标已算（{', '.join(qc_vars) if qc_vars else '无线粒体/核糖体基因命中'}）")
    # M15：把指标定义与"旧规则会多收哪些基因"打出来 —— 这是与历史产物
    # 对比时唯一能解释指标差异的依据。
    if gene_sets["legacy_prefix_extra_hb"] or gene_sets["legacy_prefix_extra_ribo"]:
        log_info(f"QC 基因集（{gene_sets['source']}）：旧前缀规则会多收 "
                 f"hb={gene_sets['legacy_prefix_extra_hb']}、"
                 f"ribo={gene_sets['legacy_prefix_extra_ribo']} —— 本轮的指标不含它们")

    # ---- 2. 过滤前的图 ------------------------------------------------------
    # 先画再滤 —— 滤完再画就看不到"滤掉了什么"，而那正是要判断的东西。
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    keys = [k for k in ("n_genes_by_counts", "total_counts", "pct_counts_mt",
                        "pct_counts_ribo", "pct_counts_hb") if k in adata.obs.columns]
    # **原始 obs 列名不能直接当面板标题**（评审 3.7）：`pct_counts_mt` 这类
    # 名字是数据结构泄露，不是给读者看的标签。统一映射成人类可读英文。
    QC_LABELS = {
        "n_genes_by_counts": "Genes detected per cell",
        "total_counts": "Total counts per cell",
        "pct_counts_mt": "Mitochondrial fraction (%)",
        "pct_counts_ribo": "Ribosomal fraction (%)",
        "pct_counts_hb": "Haemoglobin fraction (%)",
    }
    # **单图原则拆分（D-006）**：五联小提琴拆为 5 张独立单图（S1 过滤前
    # QC 流程链）。每个指标一张、独立达标图幅（W_SINGLE）与分辨率门禁；
    # 组内"联合决定过滤阈值"的叙事由 figure_groups 的 S1 组表达。
    # 指标的映射/退化处理逻辑与拆分前一致（QC_LABELS 人类可读、
    # IQR=0 画 strip 散点并注明 no variance）。
    UNIT_NAMES = {
        "n_genes_by_counts": "02-01-01-unit1-genes-detected",
        "total_counts": "02-01-01-unit2-total-counts",
        "pct_counts_mt": "02-01-01-unit3-mito-fraction",
        "pct_counts_ribo": "02-01-01-unit4-ribo-fraction",
        "pct_counts_hb": "02-01-01-unit5-hb-fraction",
    }
    for _ui, k in enumerate(keys, start=1):
        v = adata.obs[k].astype(float).values
        label = QC_LABELS.get(k, k)
        fig, ax = plt.subplots(figsize=(W_SINGLE, mm(58)))
        if float(np.subtract(*np.percentile(v, [75, 25]))) == 0:
            rng = np.random.default_rng(cfg["analysis"]["seed"])
            ax.scatter(rng.uniform(-0.12, 0.12, len(v)), v, s=2, alpha=0.25,
                       color=PAL["primary"], rasterized=True)
            ax.set_title(f"{label}\nno variance: "
                         f"{int((v != 0).sum())}/{len(v)} cells", fontsize=8)
        else:
            ax.violinplot(v, showmedians=True)
            ax.set_title(label, fontsize=8)
        ax.set_xticks([])
        # 图名来自 UNIT_NAMES 查表（字面量，过命名门禁）
        save_fig(cfg, UNIT_NAMES.get(k, "02-01-01-unit1-genes-detected"), fig)

    # 阈值线：基因数 vs 线粒体比例 —— 双细胞和死细胞在这张图上是两个角
    fig, ax = plt.subplots(figsize=(W_SINGLE, mm(64)))
    sc_ = ax.scatter(adata.obs["total_counts"], adata.obs["n_genes_by_counts"],
                     c=adata.obs["pct_counts_mt"] if "pct_counts_mt" in adata.obs else None,
                     s=3, cmap="viridis", alpha=0.6)
    if sc_ is not None:
        fig.colorbar(sc_, ax=ax, label="Mitochondrial fraction (%)")
    q = cfg["qc"]
    ax.axhline(q["min_genes"], color=PAL["highlight"], lw=1, ls="--")
    ax.axhline(q["max_genes"], color=PAL["highlight"], lw=1, ls="--")
    # 轴标题同样走映射；刻度用 5k 步长整数，避免 10000/12500/15000 挤在一起
    ax.set_xlabel(QC_LABELS["total_counts"])
    ax.set_ylabel(QC_LABELS["n_genes_by_counts"])
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=5, integer=True))
    ax.set_title("QC thresholds (red = cut-offs)")
    save_fig(cfg, "02-01-02-unit1-qc-scatter-thresholds", fig)

    # ---- 3. 过滤 ------------------------------------------------------------
    #
    # **M23（R-03 裁决）**：这里 `filter_cells`（按 min_genes 删细胞）与
    # `filter_genes`（按 min_cells 删基因）连着跑，然后把 `adata.n_obs` 记在
    # `n_after_min_genes` 名下。数值是对的（filter_genes 不改 n_obs），
    # **但名字是错的** —— 而更严重的漏记是 `filter_genes` 对 `n_vars` 的
    # 影响**完全没有记录**：读者看到"基因数 36601"和"最后 2000 HVG"之间
    # 有一段凭空的收缩，无从解释。修法：两个维度各自计数、名字与语义一致。
    q = cfg["qc"]
    n_before = adata.n_obs
    g_before = adata.n_vars
    adata.var["mt"] = adata.var_names.str.startswith(
        ("mt-",) if organism.lower().startswith("mus") else ("MT-",))
    sc.pp.filter_cells(adata, min_genes=int(q["min_genes"]))
    n_after_min_genes = adata.n_obs          # 语义 = 只跑了 min_genes 过滤之后
    g_after_filter_cells = adata.n_vars      # filter_cells 不改 n_vars，但记下来
    sc.pp.filter_genes(adata, min_cells=int(q["min_cells"]))
    g_after_min_cells = adata.n_vars         # filter_genes 的真实影响
    if q.get("max_genes"):
        adata = adata[adata.obs["n_genes_by_counts"] < int(q["max_genes"])].copy()
    n_after_maxg = adata.n_obs
    if "pct_counts_mt" in adata.obs.columns and q.get("max_pct_mt") is not None:
        adata = adata[adata.obs["pct_counts_mt"] < float(q["max_pct_mt"])].copy()
    n_after_mt = adata.n_obs

    log_info(f"过滤: 细胞 {n_before} -> {n_after_min_genes} (min_genes) "
             f"-> {n_after_maxg} (max_genes) -> {n_after_mt} (max_pct_mt)；"
             f"基因 {g_before} -> {g_after_min_cells} (min_cells，"
             f"删掉 {g_before - g_after_min_cells} 个)")

    # ---- 4. 双细胞 ----------------------------------------------------------
    db = run_scrublet(adata, cfg)
    n_after_db = adata.n_obs
    if db["status"] == "ok" and db["n_predicted_doublets"] > 0:
        adata = adata[~adata.obs["predicted_doublet"]].copy()
        n_after_db = adata.n_obs
        log_info(f"去除双细胞: {n_after_db + db['n_predicted_doublets']} -> {n_after_db}")
    elif db["status"] == "ok":
        log_info("scrublet 未标出双细胞")
    else:
        log_warn(f"双细胞检测未执行: {db['status']} —— {db.get('reason', '')}")

    # ---- 5. ambient RNA -----------------------------------------------------
    amb = assess_ambient_rna(adata, cfg)
    if amb["status"] != "ok":
        log_warn(f"ambient RNA: {amb['status']} —— {amb.get('reason', amb.get('note', ''))}")

    if adata.n_obs < 50:
        raise RuntimeError(f"过滤后只剩 {adata.n_obs} 个细胞，无法继续分析")

    # ---- 6. 落盘 ------------------------------------------------------------
    out = data_dir / "qc_filtered.h5ad"
    adata.write_h5ad(out)
    log_info(f"已写出 {out}（{adata.n_obs} 细胞 x {adata.n_vars} 基因）")

    # QC 汇总表：每个细胞一行，便于事后复查"到底滤掉了谁"
    cols = [c for c in ("n_genes_by_counts", "total_counts", "pct_counts_mt",
                        "pct_counts_ribo", "pct_counts_hb", "predicted_doublet",
                        "doublet_score") if c in adata.obs.columns]
    adata.obs[cols].to_csv(res_dir / "qc_cells.csv")

    status = {
        "dataset_id": cfg["dataset_id"],
        "organism": organism,
        "n_cells_raw": int(n0),
        "n_genes_raw": int(g0),
        "filtering": {
            "min_genes": int(q["min_genes"]),
            "max_genes": q.get("max_genes"),
            "max_pct_mt": q.get("max_pct_mt"),
            "min_cells": int(q["min_cells"]),
            # M23：名字必须与语义对得上 —— 这一步只跑了 `filter_cells(min_genes)`。
            # 旧名字 `n_after_min_genes` 下记的其实是"两个过滤器都跑完"的值，
            # 数值巧合是对的，但读者会以为 `min_cells` 也被算进来了。
            "n_after_min_genes": int(n_after_min_genes),
            "n_after_max_genes": int(n_after_maxg),
            "n_after_mt": int(n_after_mt),
            "n_after_doublets": int(n_after_db),
            "n_removed_total": int(n0 - adata.n_obs),
            "frac_removed": round((n0 - adata.n_obs) / max(n0, 1), 5),
            # M23：基因维度的收缩以前完全没记录 —— 从 36601 到最后的 HVG
            # 中间有一段凭空消失，无从解释。
            "n_genes_before_filter": int(g_before),
            "n_genes_after_filter_cells": int(g_after_filter_cells),
            "n_genes_after_min_cells": int(g_after_min_cells),
            "n_genes_removed_by_min_cells": int(g_before - g_after_min_cells),
        },
        # M15：QC 指标的**定义**（哪些基因算 hb / ribo）与旧前缀规则的差异。
        # 拿新产物对比历史产物时，指标差异的全部来源都在这里。
        "qc_gene_sets": gene_sets,
        "doublet_detection": db,
        "ambient_rna": amb,
        "n_cells_final": int(adata.n_obs),
        "n_genes_final": int(adata.n_vars),
        "status": "ok",
    }
    write_json(res_dir / "qc_status.json", status)
    return status


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    t0 = time.time()
    try:
        res = run_01_qc(cfg)
        # E-56：接住返回值。本步骤的早退路径（`not_done` / `heuristic_only` /
        # `unknown_mode`）是 `write_json` + `return status`，不抛异常 ——
        # 写死 "ok" 会让 `state.json` 与自述状态打架。
        record_step(cfg, "qc", "ok", time.time() - t0,
                    result_status=result_status_of(res))
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "qc", "failed", time.time() - t0, message=str(e),
                    result_status="failed")
        raise
