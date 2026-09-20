#!/usr/bin/env python3
"""
08_virtual_perturbation.py — 虚拟敲除 / 过表达（文档 §1.7 / §1.8）

**这是一个"保留框架"，不是已验证的因果预测。** 规范把 §1.7/§1.8 标为保留框架、
主语言 Python，并把候选靶基因的产出留给 Part 1（R）。所以这里做三件事：

  1. 记录三个点名工具的**实际可用性**（含不可用的确切原因）
  2. 用**现有调控网络**做一个一阶、单跳的 in-silico 扰动，给出可核对的数
  3. 把方法学限定写全 —— 这份表的每一列能被读成什么、不能被读成什么

## 三个点名工具在本流水线的真实状态

**这不是"没装"，是"装不了"，理由都是实测出来的：**

| 工具 | 状态 | 原因 |
|---|---|---|
| `scTenifoldKnk` | **不可用** | R/CRAN 包，**不在 PyPI**（`scTenifoldKnk` 与 `sctenifoldknk` 都查过） |
| `PerturbNet` | **不可用** | PyPI 上 0.0.2 / 0.0.3 钉 `requires_python='<3.8,>=3.7'`，0.0.3b0/b1 钉 `'<3.11,>=3.10'` —— **没有任何一版支持 CI 的 3.12** |
| `RegVelo` | **不可用** | 能装（`regvelo>=0.4.2`，`py>=3.10`），但它要 RNA velocity 的 spliced/unspliced 层，本数据没有；且拉入 `torch` + `scvi-tools` |

所以本步**不假装跑过它们**，而是用现有 GRN（`07_grn.py` 的
`tf_regulon_edges.csv`）做一个明确标为一阶近似的扰动。

## 算法（一阶、单跳）

对候选基因 G、细胞类型 C：

    Δz_G  = −z̄(G, C)              （敲除：G 从当前水平降到 0）
    Δz_G  = +overexpress_sd        （过表达：抬高 s 个标准差）
    Δz_Y  = w(G,Y) · Δz_G          对 G 的每个靶基因 Y，w 是相关系数

即：把 G 的变化量按边权线性传给它的直接靶基因。报三个量：

  - `perturbation_magnitude` = ‖Δz‖₂ —— 预测的影响**有多大**
  - `signature_alignment`    = Spearman(Δz, 疾病签名 logFC) —— 影响**朝哪个方向**
    （正 = 敲除把转录组推向疾病态，负 = 推离）
  - `n_targets`              = 有多少靶基因在数据里测到

## 这份表不能被读成什么

1. **不是因果。** 边来自共表达相关（`07_grn.py` 没有 motif 剪枝），
   相关不等于调控，更不等于"敲掉它下游就会这样变"。
2. **只有一阶、单跳。** 不走多跳传播 —— 在相关网络上做多跳会放大噪声，
   看起来像"网络效应"其实只是把相关系数乘了几遍。
3. **不是 scTenifoldKnk。** 它做流形对齐 + 张量分解；这里是线性一阶近似。
4. **不是 PerturbNet。** 它是深度生成模型；这里没有任何学习到的成分。
5. **表达层面，不是蛋白层面。** TF 的 mRNA 与其活性经常不相关。
6. `signature_alignment` 的可靠性**上界是 Part 1 签名的质量**。
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
from scipy.stats import spearmanr  # noqa: E402

from common import (df_to_records, ensure_dirs, load_config, log_info,  # noqa: E402
                    log_warn, parse_args, record_cross_language,
                    record_decision, record_step, save_fig, set_seed,
                    write_json, W_ONE_HALF, W_SINGLE, mm,)

# 候选基因没有外部靶基因表时，用调控子按簇特异性排序取前 N 个
DEFAULT_TOP_N = 20
# 一个基因至少要有几个靶基因在数据里才评估
MIN_TARGETS = 5
# 过表达抬高多少个标准差
DEFAULT_OVEREXPRESS_SD = 1.0

# 规范 §1.7/§1.8 点名的工具，以及实测的不可用原因。
# **键必须存在、值记 null 还是记原因，是两件事** —— 见 common.py 清单层约定。
TOOLS = [
    ("scTenifoldKnk", "R/CRAN 包，不在 PyPI（scTenifoldKnk / sctenifoldknk 都查过）"),
    ("PerturbNet", "PyPI 上所有版本都钉 Python 上限（<3.8 或 <3.11），CI 用 3.12 装不上"),
    ("RegVelo", "需 RNA velocity 的 spliced/unspliced 层，本数据没有；且拉入 torch + scvi-tools"),
]


def get_full_expression(adata):
    """取全基因集表达（同 06/07 —— 转录因子大多低表达，不在 HVG 里）。"""
    if adata.raw is None:
        raise RuntimeError("adata.raw 为空 —— 无法取全基因集表达做扰动模拟")
    src = adata.raw
    X = src.X
    X = X.toarray() if sparse.issparse(X) else np.asarray(X)
    return X, list(src.var_names)


def probe_tools() -> dict:
    """逐个探测工具能不能 import。**探测结果进状态文件，不靠猜。**"""
    import importlib.util
    out = {}
    for name, why in TOOLS:
        try:
            found = importlib.util.find_spec(name) is not None
        except (ImportError, ValueError):
            found = False
        out[name] = {"available": bool(found),
                     "reason": None if found else why}
    return out


def load_targets(cfg: dict, reg: pd.DataFrame) -> tuple:
    """候选靶基因：优先 Part 1 交接的 CSV（§0.2 只走 CSV），否则用内部调控子。

    返回 (DataFrame[gene, logfc], source_str)。
    """
    pert = cfg.get("perturbation") or {}
    spec = pert.get("targets_csv")
    if spec:
        p = Path(spec)
        if not p.is_absolute():
            # 相对路径按 data_dir 解析 —— 跨部分交接的文件落在这里
            p = Path(cfg["output"]["data_dir"]) / spec
        if p.exists():
            df = pd.read_csv(p)
            # 列名宽松一点：gene / target / symbol 都接受
            gcol = next((c for c in ("gene", "target", "symbol")
                         if c in df.columns), None)
            if gcol is None:
                log_warn(f"Part 1 交接表 {p.name} 里没有 gene/target/symbol 列 —— 回退到内部调控子")
            else:
                # logFC 列名大小写不定，统一找一遍
                fcol = next((c for c in df.columns
                             if c.lower() in ("logfc", "log2fc", "log_fc")), None)
                out = pd.DataFrame({"gene": df[gcol].astype(str)})
                out["logfc"] = (pd.to_numeric(df[fcol], errors="coerce")
                                if fcol else np.nan)
                out = out.dropna(subset=["gene"]).drop_duplicates("gene")
                log_info(f"候选靶基因来自 Part 1 交接表 {p.name}：{len(out)} 个"
                         + (f"，含 logFC（{fcol}）" if fcol else "，**无 logFC 列**"))
                # §0.2 跨语言转换记录。**这是本仓库唯一一处真正的跨语言交接** ——
                # 上游是 R（Part 1 的 `09_export_targets.R`），下游是 Python。
                #
                # 记的是"丢了什么"，不是"传了什么"。
                #
                # **`lost` 里不能出现 `gene` 或 logFC 那一列** —— 它们是
                # **留下的**。第一版写成"除 gene 外的所有列"，于是把
                # `logfc` 也报成了"丢失" —— 而它明明是唯一跨部分传过来的
                # 数值列。**把留下的说成丢掉的，方向和事实正好相反**，
                # 比不记更糟。CI 里的 `tools/selftest_part1_handoff.py`
                # 现在会断言这一点。
                #
                # 列名大小写不定（`logFC` / `log2FC` / `log_fc`），所以
                # **改名**本身也是信息：记进 `note`，否则下游会以为
                # 上游那一列本来就叫 `logfc`。
                _kept = {gcol} | ({fcol} if fcol else set())
                _lost = [c for c in df.columns if c not in _kept]
                _rename = (f"`{fcol}` → `logfc`（列名别名归一）"
                           if fcol and fcol != "logfc" else "")
                if fcol is None:
                    _lost.append("logFC 列（上游有，但列名不在这批别名里 → "
                                 "signature_alignment 会是 NaN）")
                record_cross_language(
                    cfg,
                    src=f"Part 1 (R) {p.name}", dst="08_virtual_perturbation 候选表",
                    fmt="csv",
                    before={"n_rows": int(len(df)), "n_cols": int(df.shape[1]),
                            "columns": [str(c) for c in df.columns]},
                    after={"n_genes": int(len(out)),
                           "n_with_logfc": int(out["logfc"].notna().sum()),
                           "kept_columns": ["gene", "logfc"]},
                    lost=_lost,
                    note=("**只走 CSV，不做对象级转换**（§0.2 禁止 sceasy，"
                          "本仓库也不用 zellkonverter / anndata2ri）。"
                          "`logfc` 是**唯一**跨部分传过来的数值列；"
                          "Part 1 的模块归属、来源计数等没有交接，"
                          "所以这边无法按来源分层看扰动结果。"
                          + (f" 列名归一：{_rename}。" if _rename else "")))
                return out, f"part1_handoff:{p.name}"
        else:
            log_warn(f"配置指定的 Part 1 交接表不存在: {p} —— 回退到内部调控子")

    n = int(pert.get("top_n", DEFAULT_TOP_N))
    head = reg.head(n)
    out = pd.DataFrame({"gene": head["tf"].astype(str)})
    out["logfc"] = np.nan
    log_info(f"候选靶基因来自内部调控子（按簇特异性取前 {len(out)} 个）—— "
             f"**没有 Part 1 签名，signature_alignment 会是 NaN**")
    # §0.2：**这一轮没有发生跨语言交接，必须说出来。**
    #
    # 清单的 `cross_language` 这时是空数组，而**空数组和"忘了记"长得一模一样**
    # （和 geo 的 `record_decision()` 定义了却没调用点是同一个病）。
    # 默认配置 `targets_csv: null` 走的正是这个分支，所以这条尤其重要。
    #
    # 注意这里记的是**决策**（说明为什么没有），不是往 `cross_language` 里
    # 塞一条假的转换记录 —— 后者会让"本轮交接了"变成假话。
    record_decision(
        cfg, "cross_language",
        "§0.2 跨语言转换：这一轮有没有 Part 1 → Part 2 的 CSV 交接？",
        ("**没有。** 配置里 `perturbation.targets_csv` 为空或文件不存在，"
         "本轮走的是回退分支（内部调控子），`cross_language` 为空是"
         "**本轮配置的结果**，不是遗漏"),
        evidence=(f"`targets_csv={spec!r}`"
                  + ("（已配置但文件不存在）" if spec else "（未配置）")
                  + "。**交接代码路径本身是好的** —— 由 "
                  "`tools/selftest_part1_handoff.py` 每轮 CI 单独跑到并断言"
                  "（造一份假交接表，验证 `record_cross_language()` 真的"
                  "写进了清单、且上游多出来的列进了 `lost_fields`）。"
                  "所以这条不是「没实现」，是「本轮没配」"))
    return out, "internal_top_regulons"


def perturb_one(gene: str, w: np.ndarray, tgt_idx: np.ndarray,
                delta_z: float, logfc: np.ndarray, shared: np.ndarray) -> dict:
    """一阶扰动。w 与 tgt_idx 是 G 的靶基因权重与列号，一一对应。"""
    dz = w * delta_z                       # 每个靶基因的预测变化
    mag = float(np.sqrt(np.sum(dz ** 2)))  # L2 范数
    out = {"perturbation_magnitude": round(mag, 4),
           "mean_abs_delta_z": round(float(np.mean(np.abs(dz))), 4),
           "max_abs_delta_z": round(float(np.max(np.abs(dz))), 4),
           "n_targets": int(len(dz))}
    # 方向：预测变化 vs 疾病签名 logFC。**共享基因少于 3 个不报相关** ——
    # 两个点永远能连成一条线，那不是相关性。
    if shared.size >= 3 and not np.all(np.isnan(logfc[shared])):
        ok = ~np.isnan(logfc[shared])
        if ok.sum() >= 3:
            rho, p = spearmanr(dz[shared][ok], logfc[shared][ok])
            out["signature_alignment"] = round(float(rho), 4)
            out["signature_alignment_p"] = round(float(p), 4)
            out["n_signature_shared"] = int(ok.sum())
        else:
            out["signature_alignment"] = None
            out["signature_alignment_p"] = None
            out["n_signature_shared"] = int(ok.sum())
    else:
        out["signature_alignment"] = None
        out["signature_alignment_p"] = None
        out["n_signature_shared"] = int(shared.size)
    return out


def run_08_virtual_perturbation(cfg: dict) -> dict:
    ensure_dirs(cfg)
    set_seed(cfg)
    data_dir = Path(cfg["output"]["data_dir"])
    res_dir = Path(cfg["output"]["results_dir"])
    status_path = res_dir / "virtual_perturbation_status.json"

    # **开跑前先删掉自己的状态文件**（AGENTS.md 规则 14）：
    # 崩溃时旧文件会冒充本轮结果。
    if status_path.exists():
        status_path.unlink()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pert = cfg.get("perturbation") or {}
    tools = probe_tools()

    if not pert.get("enabled", True):
        status = {"dataset_id": cfg["dataset_id"], "status": "disabled",
                  "reason": "配置 perturbation.enabled=false", "tools": tools}
        write_json(status_path, status)
        log_info("虚拟扰动分析已按配置关闭")
        return status

    # ---- 0. 工具可用性 ------------------------------------------------------
    for name, info in tools.items():
        if info["available"]:
            log_info(f"  {name}: 可用")
        else:
            log_warn(f"  {name}: 不可用 —— {info['reason']}")

    # ---- 1. 读输入 ----------------------------------------------------------
    edges_path = res_dir / "tf_regulon_edges.csv"
    if not edges_path.exists():
        status = {"dataset_id": cfg["dataset_id"], "status": "not_configured",
                  "reason": "缺少 tf_regulon_edges.csv —— 先跑 07_grn（它可能被关闭或失败）",
                  "tools": tools}
        write_json(status_path, status)
        log_warn(status["reason"])
        return status

    reg_path = res_dir / "tf_regulons.csv"
    reg = pd.read_csv(reg_path) if reg_path.exists() else pd.DataFrame({"tf": []})

    adata = sc.read_h5ad(data_dir / "clustered.h5ad")
    # **列名实测是 `celltype`（无下划线）**，03_cluster_annotate.py 就是这么写的。
    # 三种写法都认，是为了万一将来改名不至于静默退回 leiden ——
    # 退回 leiden 不报错，只是把"细胞类型"换成了"簇编号"，
    # 而下游解读完全变了。
    ckey = next((c for c in ("celltype", "cell_type", "celltype_label",
                             "leiden") if c in adata.obs.columns), None)
    if ckey is None:
        status = {"dataset_id": cfg["dataset_id"], "status": "not_configured",
                  "reason": "clustered.h5ad 里既没有 celltype / cell_type 也没有 leiden 列",
                  "tools": tools}
        write_json(status_path, status)
        log_warn(status["reason"])
        return status
    labels = adata.obs[ckey].astype(str).values
    groups = sorted(set(labels))
    if ckey == "leiden":
        log_warn("分组键退回到 leiden（没有细胞类型注解）—— "
                 "结果按簇编号报，解读时要知道这一点")

    X, var_names = get_full_expression(adata)
    lookup = {g: i for i, g in enumerate(var_names)}
    log_info(f"表达矩阵（全基因集）: {X.shape[0]} 细胞 x {X.shape[1]} 基因；"
             f"分组键 {ckey}（{len(groups)} 组）")

    # 逐基因标准化（z 分数）—— 让 Δz 与边权（相关系数）同尺度
    mu = X.mean(axis=0, keepdims=True)
    sd = X.std(axis=0, keepdims=True)
    sd[sd == 0] = 1.0
    Z = (X - mu) / sd

    edges = pd.read_csv(edges_path)
    log_info(f"调控子边表: {len(edges)} 条边，{edges['tf'].nunique()} 个 TF")

    targets, target_source = load_targets(cfg, reg)

    # 疾病签名（可选）
    sig = {}
    if "logfc" in targets.columns and targets["logfc"].notna().any():
        sig = dict(zip(targets["gene"], targets["logfc"]))
        sig = {g: float(v) for g, v in sig.items()
               if isinstance(v, float) and not np.isnan(v) and g in lookup}
    logfc_arr = np.full(len(var_names), np.nan)
    for g, v in sig.items():
        logfc_arr[lookup[g]] = v

    # ---- 2. 逐个候选基因 × 细胞类型做扰动 ------------------------------------
    overexpress_sd = float(pert.get("overexpress_sd", DEFAULT_OVEREXPRESS_SD))
    rows = []
    skipped = {}
    for gene in targets["gene"]:
        if gene not in lookup:
            skipped[gene] = "不在表达矩阵里"
            continue
        sub = edges[edges["tf"] == gene]
        if sub.empty:
            skipped[gene] = "在调控子边表里没有作为 TF 的出边"
            continue
        # 只保留在数据里测到的靶基因
        keep = [(t, float(w)) for t, w in zip(sub["target"], sub["corr"])
                if t in lookup]
        if len(keep) < MIN_TARGETS:
            skipped[gene] = f"可用靶基因只有 {len(keep)} 个（< {MIN_TARGETS}）"
            continue
        tgt_names = [t for t, _ in keep]
        w = np.array([x for _, x in keep], dtype=float)
        tgt_idx = np.array([lookup[t] for t in tgt_names])
        shared = tgt_idx[~np.isnan(logfc_arr[tgt_idx])]

        z_gene = Z[:, lookup[gene]]
        # **网络内在敏感度：与表达无关。** ko_magnitude = |Δz_G| × 这个值，
        # 所以只报 ko_magnitude 会把两件事混在一起 —— 排在前面的往往只是
        # "该基因在这个细胞类型里表达高"，而不是"它在网络里被连得紧"。
        # 分开报，读者才能自己判断是哪种。
        net_sens = float(np.sqrt(np.sum(w ** 2)))
        for c in groups:
            m = labels == c
            if m.sum() < 3:
                continue
            zbar = float(z_gene[m].mean())
            ko = perturb_one(gene, w, tgt_idx, -zbar, logfc_arr, shared)
            oe = perturb_one(gene, w, tgt_idx, +overexpress_sd, logfc_arr, shared)
            rows.append({
                "gene": gene, "cell_type": c, "n_cells": int(m.sum()),
                "mean_z_in_celltype": round(zbar, 4),
                "network_sensitivity": round(net_sens, 4),
                "ko_magnitude": ko["perturbation_magnitude"],
                "ko_mean_abs_delta_z": ko["mean_abs_delta_z"],
                "ko_max_abs_delta_z": ko["max_abs_delta_z"],
                "ko_signature_alignment": ko["signature_alignment"],
                "ko_signature_alignment_p": ko["signature_alignment_p"],
                "oe_magnitude": oe["perturbation_magnitude"],
                "oe_signature_alignment": oe["signature_alignment"],
                "n_targets": ko["n_targets"],
                "n_signature_shared": ko["n_signature_shared"],
            })

    if not rows:
        status = {"dataset_id": cfg["dataset_id"], "status": "no_candidates",
                  "reason": f"{len(targets)} 个候选基因没有一个满足 >= {MIN_TARGETS} 个可用靶基因",
                  "skipped": skipped, "tools": tools, "target_source": target_source}
        write_json(status_path, status)
        log_warn(status["reason"])
        return status

    df = pd.DataFrame(rows)
    df = df.sort_values("ko_magnitude", ascending=False)
    df.to_csv(res_dir / "virtual_perturbation.csv", index=False)
    log_info(f"虚拟扰动: {df['gene'].nunique()} 个候选基因 x {df['cell_type'].nunique()} 个细胞类型"
             f" = {len(df)} 行")

    # 每个候选基因影响最大的细胞类型（这是这张表最该看的一列）
    best = (df.loc[df.groupby("gene")["ko_magnitude"].idxmax()]
              [["gene", "cell_type", "ko_magnitude", "network_sensitivity",
                "ko_signature_alignment"]]
              .rename(columns={"cell_type": "most_affected_celltype"})
              .sort_values("ko_magnitude", ascending=False))
    best.to_csv(res_dir / "virtual_perturbation_top.csv", index=False)
    log_info(f"影响最大的候选基因: {best['gene'].iloc[0]} "
             f"（最敏感细胞类型 {best['most_affected_celltype'].iloc[0]}，"
             f"‖Δz‖={best['ko_magnitude'].iloc[0]:.3f}）")

    # ---- 3. 出图 ------------------------------------------------------------
    top = best.head(12)
    fig, ax = plt.subplots(figsize=(min(W_ONE_HALF, max(W_SINGLE, 0.30 * len(top) + 3.0)),
                                    mm(70)))
    ax.barh(range(len(top))[::-1], top["ko_magnitude"].values,
            color="#2C7FB8", alpha=0.85)
    ax.set_yticks(range(len(top))[::-1])
    ax.set_yticklabels([f"{g}  ({c})" for g, c in
                        zip(top["gene"], top["most_affected_celltype"])], fontsize=7)
    ax.set_xlabel(r"first-order knockout effect  $\|\Delta z\|_2$")
    ax.set_title("Virtual knockout: predicted effect size\n"
                 "(one-hop linear propagation on a co-expression GRN)",
                 fontsize=9)
    save_fig(cfg, "virtual_perturbation_effect", fig)

    n_align = int(df["ko_signature_alignment"].notna().sum())
    status = {
        "dataset_id": cfg["dataset_id"],
        "status": "ok",
        "n_candidates": int(df["gene"].nunique()),
        "n_cell_types": int(df["cell_type"].nunique()),
        "n_rows": int(len(df)),
        "n_with_signature_alignment": n_align,
        "target_source": target_source,
        "group_key": ckey,
        "top_by_effect": df_to_records(best.head(10)),
        "skipped": skipped,
        "tools": tools,
        "tools_all_unavailable": all(not v["available"] for v in tools.values()),
        "method": ("一阶单跳线性传播：把候选基因的 z 分数变化按共表达边权"
                   "传给它的直接靶基因。**不是 scTenifoldKnk，也不是 PerturbNet**"),
        "limitations": [
            "**不是因果预测**：边来自共表达相关（07_grn 没有 motif 剪枝），"
            "相关不等于调控",
            "**只有一阶、单跳**：不走多跳传播 —— 在相关网络上做多跳会放大噪声，"
            "看起来像网络效应，其实只是把相关系数乘了几遍",
            "不是 scTenifoldKnk（流形对齐 + 张量分解），不是 PerturbNet（深度生成模型）",
            "表达层面而非蛋白层面：TF 的 mRNA 与其活性经常不相关",
            "线性传播假设网络就是有效因果模型，本数据无法验证这个假设",
            ("**`ko_magnitude` 里混了两个因素**：它是 `|Δz_G| × network_sensitivity`，"
             "所以排在前面的基因往往只是**在该细胞类型里表达高**，"
             "而不是在网络里被连得紧。要看后者请用 `network_sensitivity` 列"
             "（它与表达无关）"),
            (f"signature_alignment 的可靠性上界是 Part 1 签名的质量；"
             f"本轮候选来自 {target_source}"
             + ("，**没有疾病签名，该列为空**" if n_align == 0 else "")),
        ],
    }
    write_json(status_path, status)
    return status


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    t0 = time.time()
    try:
        run_08_virtual_perturbation(cfg)
        record_step(cfg, "virtual_perturbation", "ok", time.time() - t0)
    except Exception as e:  # noqa: BLE001
        record_step(cfg, "virtual_perturbation", "failed", time.time() - t0,
                    message=str(e))
        raise
