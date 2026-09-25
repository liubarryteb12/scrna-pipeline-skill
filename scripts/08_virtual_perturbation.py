#!/usr/bin/env python3
"""
08_virtual_perturbation.py — 虚拟敲除 / 过表达（文档 §1.7 / §1.8）

**这是一个"保留框架"，不是已验证的因果预测。** 规范把 §1.7/§1.8 标为保留框架、
主语言 Python，并把候选靶基因的产出留给 Part 1（R）。所以这里做三件事：

  1. 记录三个点名工具的**实际可用性**（含不可用的确切原因）
  2. 用**两套引擎**做 in-silico 扰动，并给出两法的一致性
  3. 把方法学限定写全 —— 这份表的每一列能被读成什么、不能被读成什么

## 两套引擎，并列落盘

| 引擎 | 是什么 | 产物 |
|---|---|---|
| `first_order` | 在共表达 GRN 上做一阶单跳线性传播（本仓库自己实现） | `virtual_perturbation.csv` / `virtual_perturbation_top.csv` |
| `tenifold` | **真 · scTenifoldKnk**（R/CRAN），流形对齐 + 张量分解 | `virtual_perturbation_tenifold.csv` / `virtual_perturbation_tenifold_top.csv` |

**为什么两套都留着，而不是用 tenifold 替换掉一阶近似**：两者回答的不是同一个
问题，而且代价差两个数量级。

- 一阶近似是 **O(边数)**、**逐细胞类型**的（每个基因 × 每种细胞类型一个数）；
- scTenifoldKnk 要建 10 个网络 + 张量分解 + 流形对齐，是 **分钟级**的，
  而且**给不出细胞类型分辨率** —— 它整份数据只建一个网络，输出是
  「扰动基因 × 网络基因」的距离矩阵，没有细胞类型这一维。

所以 tenifold 的结果是**每基因一行**的全局量，与一阶近似的「基因 × 细胞类型」
不在同一个索引空间里。硬把它们拼成一张表只会让人以为 tenifold 也有细胞类型
分辨率。**两法各自落盘，另加一张一致性表**（按基因排名的 Spearman）。

## 三个点名工具在本流水线的真实状态

**PerturbNet 与 RegVelo 是"装不了"，scTenifoldKnk 曾经是"不在 PyPI"——
现在改走 R 通道，所以它可用了。** 这三条都是实测出来的：

| 工具 | 状态 | 原因 |
|---|---|---|
| `scTenifoldKnk` | **可用（走 R）** | 不在 PyPI，但在 **CRAN**（v1.1，`NeedsCompilation: no`）。所以不经 PyPI，改由 `lib/tenifold_knk.R` 调 `Rscript` 跑；CI 装 R + 该包 |
| `PerturbNet` | **不可用** | PyPI 上 0.0.2 / 0.0.3 钉 `requires_python='<3.8,>=3.7'`，0.0.3b0/b1 钉 `'<3.11,>=3.10'` —— **没有任何一版支持 CI 的 3.12** |
| `RegVelo` | **不可用** | 能装（`regvelo>=0.4.2`，`py>=3.10`），但它要 RNA velocity 的 spliced/unspliced 层，本数据没有；且拉入 `torch` + `scvi-tools` |

**"不在 PyPI" ≠ "不可用"。** 这条曾经被当成不可用的理由写进状态文件 ——
理由是**真的**（PyPI 上确实没有），但结论**越界了**：一个 R 包能不能用，
取决于 CI 里有没有 R，而不是取决于它有没有 Python wheel。
K-01b 的教训记在这里：**把"某个通道里没有"写成"不存在"，会让一个
本来可解的问题看起来无解。**

## 引擎一：一阶、单跳（`first_order`）

对候选基因 G、细胞类型 C：

    Δz_G  = −z̄(G, C)              （敲除：G 从当前水平降到 0）
    Δz_G  = +overexpress_sd        （过表达：抬高 s 个标准差）
    Δz_Y  = w(G,Y) · Δz_G          对 G 的每个靶基因 Y，w 是相关系数

即：把 G 的变化量按边权线性传给它的直接靶基因。报三个量：

  - `perturbation_magnitude` = ‖Δz‖₂ —— 预测的影响**有多大**
  - `signature_alignment`    = Spearman(Δz, 疾病签名 logFC) —— 影响**朝哪个方向**
    （正 = 敲除把转录组推向疾病态，负 = 推离）
  - `n_targets`              = 有多少靶基因在数据里测到

## 引擎二：scTenifoldKnk（`tenifold`）

R 脚本 `lib/tenifold_knk.R` 调 `scTenifoldKnk::scTenifoldKnk()`。本步负责
**按它要求的形状**准备输入并解读输出：

  - 矩阵方向必须是 **genes × cells**（基因行、细胞列）—— 源码 `R/scTenifoldKnk.R` L10
    逐字写着。**写反了不会报错**，只会算出一个垃圾网络而所有状态字段看着正常。
    所以 meta 里记了每个基因的总计数，R 侧用 `rowSums` 复核（转置后行和会变成
    文库大小，数量级都不一样）。
  - `transcriptomeWide = TRUE`：建一次 WT 网络、批量扰动全部候选基因。
    默认的单基因模式**每扰动一个基因就重建一次网络**，慢 N 倍。
  - `qc = FALSE`：`01_qc.py` 已经做过 QC，这里再滤一遍会**改变细胞集**，
    那么两套引擎用的就不是同一批细胞，数不可比。

## 这份表不能被读成什么

1. **不是因果。** 边来自共表达相关（`07_grn.py` 没有 motif 剪枝），
   相关不等于调控，更不等于"敲掉它下游就会这样变"。这一条对两套引擎
   都成立 —— scTenifoldKnk 也是从共表达网络出发的。
2. **一阶近似只有一阶、单跳。** 不走多跳传播 —— 在相关网络上做多跳会放大噪声，
   看起来像"网络效应"其实只是把相关系数乘了几遍。
3. **tenifold 没有细胞类型分辨率。** 它整份数据只建一个网络。本步**不**把它的
   距离按细胞类型重新加权 —— 那会造出一个既不是 scTenifoldKnk、也不是本仓库
   一阶近似的新方法，然后借 scTenifoldKnk 的名字发出去。
4. **不是 PerturbNet。** 它是深度生成模型；两套引擎里都没有学习到的成分。
5. **表达层面，不是蛋白层面。** TF 的 mRNA 与其活性经常不相关。
6. `signature_alignment` 的可靠性**上界是 Part 1 签名的质量**。
7. **两法一致性只说明排序像不像**，不说明哪个对。它们共享同一个上游
   （`07_grn.py` 的共表达边），一致性高可能只是**共享了同一个偏差**。
"""

from __future__ import annotations

import json
import shutil
import subprocess
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
                    log_warn, PAL, parse_args, probe_r_packages,
                    record_cross_language, record_decision, record_step,
                    save_fig, set_seed, write_json, W_ONE_HALF, W_SINGLE, mm,)

# 候选基因没有外部靶基因表时，用调控子按簇特异性排序取前 N 个
DEFAULT_TOP_N = 20
# 一个基因至少要有几个靶基因在数据里才评估
MIN_TARGETS = 5
# 过表达抬高多少个标准差
DEFAULT_OVEREXPRESS_SD = 1.0

# R 引擎脚本（真 scTenifoldKnk）。与 08 同属一个模块，放在 lib/ 下。
TENIFOLD_R = Path(__file__).resolve().parent / "lib" / "tenifold_knk.R"

# 引擎选择。`both` 是默认 —— 两法并列落盘才有得比。
DEFAULT_ENGINE = "both"
ENGINES = ("tenifold", "first_order", "both")

# ---- scTenifoldKnk 的参数默认值 -------------------------------------------
# **网络基因数是唯一真正要紧的那个旋钮。** scTenifoldNet 的 README 给了一张
# 实测运行时间表（cells × genes）：
#
#     300 × 1,000 =  3.5 min / 0.4 GB
#   1,000 × 1,000 =  4.3 min
#   1,000 × 5,000 =  2 h 52 min / 9.2 GB
#   7,500 × 7,500 = 10 h 16 min / 22.6 GB
#
# 本数据集是 2652 细胞 × 13714 基因 —— **全基因集绝对跑不动**，
# 而且 CI 的总上限是 45 分钟。所以必须限定基因子集：
# 取候选基因（TFs 及其靶基因）里表达方差最高的 `max_genes` 个，
# 并**保证全部候选敲除基因都在里面**（否则包会 stop）。
DEFAULT_TENIFOLD = {
    "max_genes": 600,     # 网络基因数上限（见上面的实测运行时间表）
    "n_net": 10,          # 包默认
    "n_cells": 500,       # 每个网络抽多少细胞（有放回）；包默认
    "n_comp": 3,          # PC 回归的主成分数；包默认
    "td_k": 3,            # 张量分解秩；包默认
    "ma_n_dim": 2,        # 流形对齐维度；包默认
    "n_cores": 1,         # CI 钉单线程（与 workflow 的 env 一致）
    "timeout_sec": 1500,  # 单个 R 进程上限；CI job 上限 45 min，留足余量
}

# 规范 §1.7/§1.8 点名的工具，以及实测的不可用原因。
# **键必须存在、值记 null 还是记原因，是两件事** —— 见 common.py 清单层约定。
#
# `scTenifoldKnk` 的 reason 只在**探测失败**时用（R 没装 / 包没装），
# 措辞要能区分"CI 没配 R"和"这个包真的不存在"—— 这两种情况的处置完全不同。
TOOLS = [
    ("scTenifoldKnk", "R 通道不可用（Rscript 不在 PATH 上，或 CRAN 包 scTenifoldKnk 未安装）"),
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


def _probe_r_package(pkg: str) -> dict:
    """探测一个 R 包能不能用。**不能走 `importlib.util.find_spec`** ——
    那是 Python 的模块查找器，对 R 包永远返回 False。第一版就是这么把
    scTenifoldKnk 判成不可用的：理由写的是"不在 PyPI"（真的），
    结论却是"不可用"（越界了）。R 包能不能用，取决于 CI 里有没有 R。

    实现转调 `common.probe_r_packages()` —— **同一件事只留一份代码**。
    这里原本自己拼了一次 `Rscript -e requireNamespace`，而 `capture_versions`
    那条通道查的是 `importlib.metadata`，于是同一轮里
    `virtual_perturbation_status.json` 说"装了 1.1"、
    `run_manifest.json` 说"没装"（两份产物互相矛盾，而后者是复现依据）。
    两份实现的差别不在写法，在**问错了对象** —— 修法是把问的对象统一成 R。
    """
    info = probe_r_packages([pkg])
    version = info["packages"].get(pkg)
    if version is not None:
        return {"available": True, "reason": None,
                "version": version, "r_version": info["r_version"]}
    if info["reason"]:
        reason = info["reason"]
    else:
        reason = (f"R {info['r_version']} 装上了，但 CRAN 包 {pkg} 没装上 —— "
                  f"检查 CI 的 R 依赖步骤")
    return {"available": False, "reason": reason,
            "version": None, "r_version": info["r_version"]}


def probe_tools() -> dict:
    """逐个探测工具能不能用。**探测结果进状态文件，不靠猜。**

    Python 包走 `importlib.util.find_spec`；R 包走 `Rscript -e requireNamespace`。
    两条通道分开探 —— 混在一起会让"R 没装"看起来像"包不存在"。
    """
    import importlib.util
    out = {}
    for name, why in TOOLS:
        if name == "scTenifoldKnk":
            info = _probe_r_package("scTenifoldKnk")
            if not info["available"]:
                info["reason"] = f"{why}。具体: {info['reason']}"
            out[name] = info
            continue
        try:
            found = importlib.util.find_spec(name) is not None
        except (ImportError, ValueError):
            found = False
        out[name] = {"available": bool(found),
                     "reason": None if found else why,
                     "version": None, "r_version": None}
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
    # **先按 TF 去重，再取前 N 行。** `tf_regulons.csv` 是**行级**排名：同一个
    # TF 在不同簇上各有一行（217 行 / 201 个唯一 TF，15 个 TF 出现多次）。
    # `reg.head(n)` 是行级截断，于是 20 个候选里只有 17 个唯一 —— TBX21 /
    # IRF1 / EZH2 各出现两次。后果不只是候选表多几行：
    #   - `virtual_perturbation.csv` 出现 21 行重复的 (gene, cell_type)；
    #   - status 的 `n_rows` 报 140（真值 119），虚报 18%；
    #   - 传给 R 的 target 列表是 20 条，而实际只有 17 个基因被敲。
    # 去重后**不丢基因、还会多出 3 个**（E2F3 / ELK4 / TCF3）—— 因为名额
    # 原本被重复行占了。这正是 AGENTS 规则 23.4 记的那件事，当时只修了图，
    # 消费侧这条漏了。
    #
    # `drop_duplicates` 保首行，而 `reg` 已按 `cluster_specificity` 降序 ——
    # 所以留下的是该 TF 特异性最高的那一行，不是任意一行。
    reg_u = reg.drop_duplicates(subset=["tf"])
    head = reg_u.head(n)
    out = pd.DataFrame({"gene": head["tf"].astype(str)})
    out["logfc"] = np.nan
    if len(reg_u) < len(reg):
        log_info(f"候选表按 TF 去重：{len(reg)} 行 -> {len(reg_u)} 个唯一 TF"
                 f"（取前 {len(out)} 个）")
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


def open_full_counts(cfg: dict, adata) -> object:
    """取**全基因集**的原始计数矩阵（稀疏），细胞顺序对齐 `adata`。

    **不能用 `clustered.h5ad` 的 `layers["counts"]`。** 那一层确实是
    `02_integrate.py` L50 在子集化**之前**存的原始计数，但 L104 的
    `work = adata[:, adata.var["highly_variable"]].copy()` 把整个对象缩到了
    2000 个 HVG，层也跟着缩。而候选敲除基因来自 `07_grn.py` 的**全基因集**
    推断（它读 `adata.raw`，13714 个基因）—— 拿 2000 基因的计数矩阵去喂
    scTenifoldKnk，绝大多数候选基因根本不在矩阵里，R 侧会直接 `stop()`。

    全基因集的原始计数只在 `qc_filtered.h5ad` 里（`01_qc.py` L249 落盘，
    尚未子集化）。`01_qc.py` 全程没有 normalize/log1p，所以 `.X` 就是计数。
    """
    qc_path = Path(cfg["output"]["data_dir"]) / "qc_filtered.h5ad"
    if not qc_path.exists():
        raise RuntimeError(
            f"{qc_path} 不存在 —— 全基因集原始计数只在这份文件里。"
            f"clustered.h5ad 的 layers['counts'] 是 2000 个 HVG 的子集，"
            f"拿它喂 scTenifoldKnk 会让绝大多数候选基因不在网络里")
    qc = sc.read_h5ad(qc_path)
    if list(qc.obs_names) != list(adata.obs_names):
        # 细胞顺序不一致会让"哪个细胞属于哪个类型"整体错位。**错位不报错**,
        # 只会算出看着正常的垃圾 —— 所以这里显式对齐，并把这件事说出来。
        qc = qc[adata.obs_names].copy()
        log_warn("qc_filtered.h5ad 与 clustered.h5ad 的细胞顺序不同 —— 已按后者重排")
    X = qc.X if "counts" not in qc.layers else qc.layers["counts"]
    if not sparse.issparse(X):
        X = sparse.csr_matrix(np.asarray(X))
    return X, list(qc.var_names)


def _gene_variance(X_sparse) -> np.ndarray:
    """稀疏矩阵上按列算方差（不做 densify）。

    `X.var(axis=0)` 对 scipy 稀疏矩阵不存在；`np.var(X.toarray())` 会把
    2652 x 13714 展成 291 MB 的 float64。用 E[X²] − E[X]²。
    """
    n = X_sparse.shape[0]
    mean = np.asarray(X_sparse.mean(axis=0)).ravel()
    sq = np.asarray(X_sparse.multiply(X_sparse).mean(axis=0)).ravel()
    return np.maximum(sq - mean ** 2, 0.0)


def select_tenifold_genes(var_names, cand_genes, max_genes: int,
                          variance=None) -> list:
    """挑送进 scTenifoldKnk 的基因子集。

    **这是本引擎唯一真正要紧的选择**，因为它直接决定跑不跑得完
    （见 `DEFAULT_TENIFOLD` 上方的实测运行时间表）。规则：

      1. **全部候选敲除基因必须在里面** —— 包对不在网络里的 `gKO` 是
         `stop()` 而不是跳过，少一个就整体失败。
      2. 其余名额给**表达方差最高**的基因。传了 `variance`（与 var_names
         同序）就按它排；没传则退化为 var_names 原顺序（此时"方差最高"
         这句话不成立，所以会记一条 warn）。
      3. 顺序按**输入 var_names 的原顺序**输出 —— R 侧断言行名逐位相符，
         顺序在这里定下来，后面不能重排。

    返回基因名列表（长度 >= 候选数，通常 <= max_genes）。
    """
    names = list(var_names)
    cand = [g for g in cand_genes if g in set(names)]
    if not cand:
        return []
    keep = set(cand)
    if len(cand) >= max_genes:
        # 候选本身就超了上限：**不能砍候选**（砍了就少敲一个基因，
        # 而下游会以为那个基因"跑了但没结果"）。宁可超上限，也要把
        # 候选敲全 —— 上限是时间预算，候选完整性是正确性。
        log_warn(f"候选基因 {len(cand)} 个已超过网络基因上限 {max_genes} —— "
                 f"不砍候选（砍了会让下游以为该基因跑了却没结果），"
                 f"按 {len(cand)} 个建网络")
        return [g for g in names if g in keep]
    if variance is None:
        log_warn("select_tenifold_genes 没有收到方差 —— 补足名额时用的是 "
                 "var_names 原顺序，不是表达方差")
        pool = names
    else:
        order = np.argsort(np.asarray(variance))[::-1]
        pool = [names[int(i)] for i in order]
    for g in pool:
        if len(keep) >= max_genes:
            break
        keep.add(g)
    return [g for g in names if g in keep]


def _sig6(v: float) -> float:
    """保留 6 位**有效数字**（不是 6 位小数）。

    **这是 K-01b 第一版的一个真 bug 的修法。** 当时写的是 `round(v, 6)`，
    而 scTenifoldKnk 的流形距离量级是 **1e-9 ~ 1e-3** —— 6 位小数把
    `1.48e-08` 归成 `0.0`、把 `1.09e-06` 压成 `1e-06`。后果：

      - 15 个基因里 8 个 mean 变成 `0.0`、6 个并列 `1e-06`，**排序被摧毁**
        —— 而排序正是这张表存在的全部意义；
      - 出图用的就是这列（`ax.barh(top["tenifold_mean_distance"])`），
        于是 5 根柱子长度为零、6 根等长，图上看不出谁强谁弱；
      - 一致性系数被改动：原始 mean 算 Spearman rho=0.575（p=0.0249），
        round6 后 0.5634（p=0.0287）—— 落盘的是后者。

    `round()` 的位数是**绝对**刻度，对跨数量级的量必然出事；有效数字是
    **相对**刻度，这正是"只想压缩显示、不想改变排序"时该用的东西。
    一阶近似那边用 `round(x, 4)` 没事，因为它的量级是 O(1)~O(10)。

    实现用 `:.6g`（6 位有效数字），转回 float 让 pandas 写出来是数值而不是
    字符串。非有限值原样返回。
    """
    f = float(v)
    if not np.isfinite(f):
        return f
    return float(f"{f:.6g}")


def compress_tenifold_distances(dist: pd.DataFrame, rmeta: dict) -> pd.DataFrame:
    """把「扰动基因 × 网络基因」的距离矩阵压成每基因一行的可解读量。

    scTenifoldKnk 的 `perturbationDistances[g, y]` 是「敲掉 g 之后，基因 y
    在流形里移动了多远」。**没有细胞类型这一维** —— 它整份数据只建一个网络。
    所以这里报的是全局量，列名里明确带 `tenifold_` 前缀，防止与一阶近似的
    细胞类型分辨率混为一谈（两边都有 `gene` 列，但分辨率含义不同）。

    **空敲除**：R 侧把出度为 0 的候选基因整行置了 NA（见 `tenifold_knk.R`
    的 `zero_outdegree_targets`）。这里用 meta 里那份**结构证据**
    （`empty_knockout_genes`）打标记，而不是反过来从「整行都是 NA」去推断 ——
    推断会把「R 侧因为别的原因给了 NA」也误标成空敲除。

    抽成独立函数是为了能**不装 R 就测**（见 `tools/selftest_tenifold.py`）：
    空敲除这条保护一旦失效，失败方式是"输出一排看起来正常的数"，
    没有别的机会发现。
    """
    empty_ko = set(rmeta.get("empty_knockout_genes") or [])
    outdeg = rmeta.get("target_outdegree") or {}
    dmat = dist.to_numpy(dtype=float)
    rows = []
    for i, g in enumerate(dist.index.astype(str)):
        v = dmat[i]
        ok = ~np.isnan(v)
        is_empty = g in empty_ko
        # 空敲除的行必须报 None，不能报 0.0 —— 0.0 读起来是"敲除没有影响"，
        # 而真相是这个网络表达不了这个扰动。
        if ok.sum() == 0 or is_empty:
            rows.append({"gene": g, "tenifold_n_genes_scored": 0,
                         "tenifold_mean_distance": None,
                         "tenifold_max_distance": None,
                         "tenifold_top_gene": None,
                         "tenifold_top_distance": None,
                         "tenifold_target_outdegree": outdeg.get(g),
                         "tenifold_empty_knockout": bool(is_empty)})
            continue
        vv = np.where(ok, v, -np.inf)
        j = int(np.argmax(vv))
        rows.append({
            "gene": g,
            "tenifold_n_genes_scored": int(ok.sum()),
            "tenifold_mean_distance": _sig6(float(np.nanmean(v))),
            "tenifold_max_distance": _sig6(float(np.nanmax(v))),
            "tenifold_top_gene": str(dist.columns[j]),
            "tenifold_top_distance": _sig6(float(v[j])),
            "tenifold_target_outdegree": outdeg.get(g),
            "tenifold_empty_knockout": False,
        })
    tdf = pd.DataFrame(rows)
    # 按平均距离降序 —— 与一阶近似的 `ko_magnitude` 降序方向一致，
    # 这样两边的"top N"是可比的排名。NaN 排在最后（pandas 默认
    # na_position='last'），正好让空敲除的行沉底。
    return tdf.sort_values("tenifold_mean_distance", ascending=False)


def run_tenifold_engine(cfg: dict, X_counts, var_names, targets: pd.DataFrame,
                        out_dir: Path, tparams: dict) -> dict:
    """跑真 scTenifoldKnk，落盘距离矩阵 + meta。

    `X_counts` 必须是**全基因集原始计数**（稀疏，列与 `var_names` 一一对应），
    由 `open_full_counts()` 提供。
    返回一个 dict（无论成败），由调用点决定怎么记进状态文件。
    """
    probe = _probe_r_package("scTenifoldKnk")
    if not probe["available"]:
        return {"status": "unavailable", "reason": probe["reason"]}

    cand = [g for g in targets["gene"].astype(str) if g in set(var_names)]
    if not cand:
        return {"status": "no_candidates",
                "reason": "候选基因一个都不在计数矩阵的基因集里"}

    genes = select_tenifold_genes(var_names, cand, int(tparams["max_genes"]),
                                  variance=_gene_variance(X_counts))
    missing = [g for g in cand if g not in set(genes)]
    if missing:
        # select_tenifold_genes 已保证不砍候选，走到这里说明有 bug
        return {"status": "failed",
                "reason": f"内部错误：候选基因 {missing} 没进网络基因集"}

    pos = {g: i for i, g in enumerate(var_names)}
    gi = [pos[g] for g in genes]
    sub = X_counts[:, gi]
    sub = sub.toarray() if sparse.issparse(sub) else np.asarray(sub)
    sub = sub.T                                     # genes x cells —— 关键
    log_info(f"scTenifoldKnk: 网络 {sub.shape[0]} 基因 x {sub.shape[1]} 细胞，"
             f"敲除 {len(cand)} 个候选")

    tmp = out_dir / "_tenifold_tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    mat_path = tmp / "counts_matrix.txt"
    meta_path = tmp / "counts_meta.json"
    targets_path = tmp / "targets.txt"

    # 裸文本矩阵：第 1 行基因名，其后每行一个基因的细胞计数。
    # **不能写 CSV** —— 10^6 量级的整数，`read.csv` 要几十秒，而它唯一的
    # 好处（自带列名）我们不需要，因为基因顺序记在 meta 里。更要紧的是
    # `read.csv` 会把细胞条码里的 `-`/`.` 改写成 `X` 前缀。
    with open(mat_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(genes) + "\n")
        for row in sub:
            fh.write("\t".join(str(int(v)) for v in row) + "\n")
    write_json(meta_path, {
        "n_genes": int(sub.shape[0]), "n_cells": int(sub.shape[1]),
        "genes": genes,
        # **方向见证**：每个基因的总计数。R 侧用 rowSums 复核 ——
        # 转置后行和会变成文库大小，数量级都不一样。只查维度拦不住方阵转置。
        "gene_sums": [float(v) for v in sub.sum(axis=1)],
    })
    targets_path.write_text("\n".join(cand) + "\n", encoding="utf-8")

    cmd = [shutil.which("Rscript"), str(TENIFOLD_R),
           "--input", str(mat_path), "--meta", str(meta_path),
           "--targets", str(targets_path), "--out_dir", str(tmp),
           "--n_net", str(tparams["n_net"]), "--n_cells", str(tparams["n_cells"]),
           "--n_comp", str(tparams["n_comp"]), "--td_k", str(tparams["td_k"]),
           "--ma_n_dim", str(tparams["ma_n_dim"]), "--n_cores", str(tparams["n_cores"])]
    log_info(f"调用 Rscript: {' '.join(cmd[:2])} … "
             f"(nNet={tparams['n_net']} nCells={tparams['n_cells']} "
             f"nComp={tparams['n_comp']} tdK={tparams['td_k']})")
    t0 = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=int(tparams["timeout_sec"]))
    except subprocess.TimeoutExpired:
        return {"status": "timeout",
                "reason": (f"scTenifoldKnk 超过 {tparams['timeout_sec']} s 未结束"
                           f"（网络 {sub.shape[0]} 基因）。调小 "
                           f"perturbation.tenifold.max_genes")}
    elapsed = time.time() - t0
    # R 的 stderr 直接进 CI 日志 —— 这是不走 rpy2 的主要好处之一
    for line in (p.stdout or "").splitlines():
        log_info(f"  [R] {line}")
    if p.returncode != 0:
        tail = "\n".join((p.stderr or "").strip().splitlines()[-8:])
        return {"status": "failed",
                "reason": f"Rscript 退出码 {p.returncode}",
                "stderr_tail": tail, "elapsed_sec": round(elapsed, 1)}

    dist_csv = tmp / "tenifold_perturbation_distances.csv"
    rmeta_path = tmp / "tenifold_meta.json"
    if not dist_csv.exists():
        return {"status": "failed",
                "reason": "R 脚本退出码 0 但没有写出 tenifold_perturbation_distances.csv",
                "elapsed_sec": round(elapsed, 1)}
    dist = pd.read_csv(dist_csv, index_col=0)
    rmeta = json.loads(rmeta_path.read_text(encoding="utf-8")) \
        if rmeta_path.exists() else {}

    # 输入矩阵已经读进内存了，删掉它 —— `_tenifold_tmp` 在 results/ 下，
    # 而 results/ 是整目录上传的，不删就会把一个几 MB 的可重建文本矩阵
    # 一起塞进 artifact。**失败时不删**（早返回路径走不到这里）：
    # 那种时候这份矩阵正是"R 到底吃到了什么"的唯一证据。
    shutil.rmtree(tmp, ignore_errors=True)

    # ---- 把「扰动基因 × 网络基因」的距离压成每基因一行的可解读量 ----------
    # 逻辑在 `compress_tenifold_distances()` 里（抽出来是为了能不装 R 就测）。
    tdf = compress_tenifold_distances(dist, rmeta)
    n_empty = int(tdf["tenifold_empty_knockout"].sum())
    if n_empty:
        log_warn(f"scTenifoldKnk: {n_empty}/{len(tdf)} 个候选基因在去噪网络中"
                 f"出度为 0（{'、'.join(sorted(set(rmeta.get('empty_knockout_genes') or [])))}）"
                 f"—— 它们的距离行已置空，"
                 f"**不是「效应为 0」，而是这个网络表达不了该扰动**")
    tdf.to_csv(out_dir / "virtual_perturbation_tenifold.csv", index=False)
    tdf.head(10).to_csv(out_dir / "virtual_perturbation_tenifold_top.csv",
                        index=False)
    # 完整距离矩阵也留一份：它是唯一的原始产物，压成每基因一行会丢掉
    # "这个基因影响了哪些下游"的信息
    dist.to_csv(out_dir / "tenifold_perturbation_distances.csv")

    return {"status": "ok", "elapsed_sec": round(elapsed, 1),
            "n_genes_network": int(dist.shape[1]),
            "n_targets": int(dist.shape[0]),
            "n_empty_knockout": n_empty,
            # 从 meta 取，不从 `compress_tenifold_distances` 的局部变量取 ——
            # 那个变量随函数抽出去一起走了（抽出时这里漏改过一次，
            # `py_compile` 查不出来，只有 tenifold 真跑通那一刻才会
            # NameError）。
            "empty_knockout_genes": sorted(rmeta.get("empty_knockout_genes") or []),
            "genes": genes, "targets": cand,
            "meta": rmeta,
            "table": tdf}


def compare_engines(fo: pd.DataFrame, td: pd.DataFrame) -> dict:
    """两套引擎的一致性。**只比排名，不比数值** —— 两个量的单位根本不同
    （一阶是 z 分数的 L2 范数，tenifold 是流形上的欧氏距离），
    比数值等于比苹果和橘子。

    一阶近似是「基因 × 细胞类型」，所以每个基因先按**最大** ko_magnitude
    取一行（与 `virtual_perturbation_top.csv` 的取法一致），再与 tenifold
    的每基因一行对齐。
    """
    if fo is None or td is None or fo.empty or td.empty:
        return {"comparable": False, "reason": "两套引擎没有同时产出结果"}
    # 空敲除的基因（距离为空）不能进相关性 —— 它们的"距离"是浮点噪声，
    # 算进去等于往相关系数里掺随机数。这里显式剔掉并记账，而不是靠
    # dropna() 悄悄少几个点。
    td_valid = td[td["tenifold_mean_distance"].notna()]
    excluded = sorted(set(td["gene"]) - set(td_valid["gene"]))
    fo_best = (fo.loc[fo.groupby("gene")["ko_magnitude"].idxmax()]
                 [["gene", "ko_magnitude"]])
    m = fo_best.merge(td_valid[["gene", "tenifold_mean_distance"]], on="gene",
                      how="inner").dropna()
    out = {"comparable": True, "n_shared_genes": int(len(m)),
           "engines": ["first_order", "scTenifoldKnk"]}
    if excluded:
        out["excluded_empty_knockout"] = excluded
        out["excluded_reason"] = ("这些基因在去噪网络里出度为 0，scTenifoldKnk "
                                  "对它们的输出是浮点噪声，不参与相关性")
    if len(m) < 3:
        # 两个点永远能连成一条线 —— 少于 3 个共享基因不报相关系数
        out["reason"] = f"只有 {len(m)} 个共享基因，少于 3 个不报相关系数"
        return out
    rho, p = spearmanr(m["ko_magnitude"], m["tenifold_mean_distance"])
    out["spearman_rho"] = round(float(rho), 4)
    out["spearman_p"] = round(float(p), 4)
    out["note"] = ("**只说明排序像不像，不说明哪个对。** 两法共享同一个上游"
                   "（07_grn 的共表达边），一致性高可能只是共享了同一个偏差")
    # 两个排名的 top 5 重合度：比单一相关系数更好读
    a = list(fo_best.sort_values("ko_magnitude", ascending=False)["gene"][:5])
    b = list(td_valid.sort_values("tenifold_mean_distance",
                                  ascending=False)["gene"][:5])
    out["first_order_top5"] = a
    out["tenifold_top5"] = b
    out["top5_overlap"] = sorted(set(a) & set(b))
    return out


def run_first_order_engine(targets, edges, lookup, Z, logfc_arr, labels, groups,
                           overexpress_sd: float) -> tuple:
    """引擎一：一阶单跳线性传播。返回 `(rows, skipped)`。

    对每个候选基因 G、每种细胞类型 C：
        Δz_G = −z̄(G,C)（敲除）/ +overexpress_sd（过表达）
        Δz_Y = w(G,Y) · Δz_G
    """
    rows, skipped = [], {}
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
    return rows, skipped


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

    # ---- 引擎选择 -----------------------------------------------------------
    engine = str(pert.get("engine", DEFAULT_ENGINE)).strip().lower()
    if engine not in ENGINES:
        # **拼错引擎名不能静默退回默认** —— 那会让"我配了 tenifold"变成
        # 一句假话，而状态文件里看不出差别。
        raise ValueError(f"perturbation.engine={engine!r} 不认识，可选 {list(ENGINES)}")
    tparams = dict(DEFAULT_TENIFOLD)
    tparams.update(pert.get("tenifold") or {})
    log_info(f"扰动引擎: {engine}"
             + (f"（tenifold 网络基因上限 {tparams['max_genes']}）"
                if engine in ("tenifold", "both") else ""))

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

    # ---- 2a. 引擎一：一阶单跳（基因 × 细胞类型）-----------------------------
    overexpress_sd = float(pert.get("overexpress_sd", DEFAULT_OVEREXPRESS_SD))
    rows, skipped = [], {}
    if engine in ("first_order", "both"):
        rows, skipped = run_first_order_engine(
            targets, edges, lookup, Z, logfc_arr, labels, groups, overexpress_sd)

    # ---- 2b. 引擎二：真 scTenifoldKnk（R，整份数据一个网络）-----------------
    tenifold = None
    if engine in ("tenifold", "both"):
        Xc, cvar = open_full_counts(cfg, adata)
        tenifold = run_tenifold_engine(cfg, Xc, cvar, targets, res_dir, tparams)
        if tenifold["status"] == "ok":
            log_info(f"scTenifoldKnk 完成：{tenifold['elapsed_sec']} s，"
                     f"{tenifold['n_targets']} 个候选 x "
                     f"{tenifold['n_genes_network']} 个网络基因")
        else:
            log_warn(f"scTenifoldKnk 未产出结果（{tenifold['status']}）—— "
                     f"{tenifold.get('reason', '')}")

    # 两个引擎都没东西可报才算"没有候选"
    if not rows and (tenifold is None or tenifold["status"] != "ok"):
        status = {"dataset_id": cfg["dataset_id"], "status": "no_candidates",
                  "reason": f"{len(targets)} 个候选基因没有一个满足 >= {MIN_TARGETS} 个可用靶基因",
                  "engine": engine,
                  "skipped": skipped, "tools": tools, "target_source": target_source}
        if tenifold is not None:
            status["tenifold"] = {k: v for k, v in tenifold.items() if k != "table"}
        write_json(status_path, status)
        log_warn(status["reason"])
        return status

    if rows:
        df = pd.DataFrame(rows)
        df = df.sort_values("ko_magnitude", ascending=False)
        df.to_csv(res_dir / "virtual_perturbation.csv", index=False)
        log_info(f"虚拟扰动（一阶）: {df['gene'].nunique()} 个候选基因 x "
                 f"{df['cell_type'].nunique()} 个细胞类型 = {len(df)} 行")

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
    else:
        df, best = pd.DataFrame(), pd.DataFrame()
        log_info("一阶引擎本轮没有产出（引擎=%s）" % engine)

    # ---- 2c. 两法一致性 ------------------------------------------------------
    # `td_ok` = tenifold **真的跑出了结果**（决定 engines_used / 方法串 /
    # 一致性对比）；`use_td` = **图里画的是 tenifold**（决定画哪张图）。
    # 两者必须分开：全部候选都是空敲除时 tenifold 确实跑了（该记进
    # engines_used），但图得退回一阶 —— 把"跑了"和"画了"当成一件事，
    # 会让状态文件里 tenifold 消失、方法串也不再提它。
    td_ok = tenifold is not None and tenifold["status"] == "ok"
    cmp = None
    if rows and td_ok:
        cmp = compare_engines(df, tenifold["table"])
        log_info(f"两法排名一致性: Spearman rho="
                 f"{cmp.get('spearman_rho')}（{cmp['n_shared_genes']} 个共享基因）")

    # ---- 3. 出图 ------------------------------------------------------------
    #
    # 一张图只能画一个引擎 —— 两法的量纲不同（一阶是 z 分数的 L2 范数，
    # tenifold 是流形欧氏距离），画在同一根 x 轴上会让人以为可以直接比大小。
    # 默认画 tenifold（真方法），它没跑出来时退回一阶，并在标题里写明是哪个。
    #
    # 空敲除的行距离是空的（见 `tenifold_knk.R` 的 `zero_outdegree_targets`），
    # 它们**不能进图** —— 柱长会变成 0，读起来正是"敲除没有影响"这个我们要
    # 避免的误读。若全部候选都是空敲除，图会变成一张白板，而**白板能通过
    # 所有既有图门禁**（图非空/图名/图幅/dpi/配色都只查"有没有产出"），
    # 所以这里显式退回一阶的图。
    td_plot = (tenifold["table"][tenifold["table"]["tenifold_mean_distance"].notna()]
               if td_ok else pd.DataFrame())
    use_td = td_ok and not td_plot.empty
    if td_ok and td_plot.empty:
        log_warn("scTenifoldKnk 跑通了，但**所有候选基因都是空敲除**"
                 "（在去噪网络里出度为 0），没有可画的效应 —— 图退回一阶引擎；"
                 "这不是「效应都很小」，是网络装不下这批候选")
    if use_td:
        top = td_plot.head(12)
        vals = top["tenifold_mean_distance"].to_numpy(dtype=float)
        ys = np.arange(len(top))[::-1]
        fig, ax = plt.subplots(
            figsize=(min(W_ONE_HALF, max(W_SINGLE, 0.30 * len(top) + 3.0)), mm(70)))
        # **点，不是柱子。** 十个候选的距离跨 3~4 个数量级（1.1e-09 ~ 4.5e-06），
        # 柱子的长度是"从基线量起"的，零基线在这里没有意义 —— 画线性柱时
        # 弱的一半会变成零长（看起来正是"敲除没有影响"这个我们要避免的
        # 误读），画对数柱时"柱长"又会被读成倍数。位置编码允许非零起点，
        # 所以改成对数轴上的点：排名看得见，且轴标签写明是对数。
        ax.plot(vals, ys, "o", linestyle="none", color=PAL["primary"],
                ms=5.0, markeredgecolor="white", markeredgewidth=0.6)
        ax.set_xscale("log")
        ax.set_yticks(ys)
        ax.set_yticklabels([str(g) for g in top["gene"]], fontsize=7)
        ax.set_ylim(-0.6, len(top) - 0.4)
        # 左右各留出标注文字的净空（数值标在点的右侧）。**不用 `margins()`** ——
        # 对数轴上它的伸缩是按倍率算的，量级跨度大时会留出一大段空白。
        ax.set_xlim(vals.min() * 0.45, vals.max() * 4.5)
        # 每个点右边标出真实数值 —— 对数轴上看不出绝对量级，
        # 而这列数字是这张图唯一的定量产出。
        for y, v in zip(ys, vals):
            ax.annotate(f"{v:.2g}", (v, y), xytext=(4, 0),
                        textcoords="offset points", va="center",
                        fontsize=6.5, color=PAL["muted"])
        ax.set_xlabel("scTenifoldKnk perturbation distance (manifold, log scale)")
        ax.set_title("Virtual knockout: predicted effect size\n"
                     "(scTenifoldKnk: tensor-decomposed network + manifold alignment;"
                     " mean distance over the 600-gene network)",
                     fontsize=9)
    else:
        top = best.head(12)
        fig, ax = plt.subplots(
            figsize=(min(W_ONE_HALF, max(W_SINGLE, 0.30 * len(top) + 3.0)), mm(70)))
        ax.barh(range(len(top))[::-1], top["ko_magnitude"].values,
                color=PAL["primary"], alpha=0.85)
        ax.set_yticks(range(len(top))[::-1])
        ax.set_yticklabels([f"{g}  ({c})" for g, c in
                            zip(top["gene"], top["most_affected_celltype"])], fontsize=7)
        ax.set_xlabel(r"first-order knockout effect  $\|\Delta z\|_2$")
        ax.set_title("Virtual knockout: predicted effect size\n"
                     "(one-hop linear propagation on a co-expression GRN)",
                     fontsize=9)
    save_fig(cfg, "02-08-01-unit1-virtual-perturbation-effect", fig)

    # ---- 4. 状态 ------------------------------------------------------------
    # 这里用的是 `td_ok`（**跑了**），不是 `use_td`（**画了**）——
    # 全部候选都是空敲除时图退回一阶，但 tenifold 确实跑了、也确实产出了
    # 距离表，它必须出现在 engines_used 和方法串里，否则状态文件会谎称
    # "本轮没跑 scTenifoldKnk"。
    n_align = int(df["ko_signature_alignment"].notna().sum()) if rows else 0
    engines_used = []
    if rows:
        engines_used.append("first_order")
    if td_ok:
        engines_used.append("tenifold")

    method_bits = []
    if rows:
        method_bits.append(
            "**一阶单跳线性传播**：把候选基因的 z 分数变化按共表达边权传给它的"
            "直接靶基因；有细胞类型分辨率（基因 x 细胞类型）")
    if td_ok:
        method_bits.append(
            f"**scTenifoldKnk {tenifold['meta'].get('engine_version', '')}**（R/CRAN）："
            f"对 {tenifold['n_genes_network']} 个基因的计数矩阵做"
            f"主成分回归建网络 -> 张量分解去噪 -> 逐个敲除候选基因 -> "
            f"流形对齐量扰动距离；**没有细胞类型分辨率**（整份数据只建一个网络）")
    if not td_ok and engine in ("tenifold", "both"):
        method_bits.append(
            f"**scTenifoldKnk 本轮没有产出结果**（{tenifold['status'] if tenifold else '未执行'}："
            f"{(tenifold or {}).get('reason', '')}）—— 状态里的 `tenifold` 字段记了原因")
    method_bits.append("**不是 PerturbNet**（深度生成模型），也不是因果预测")
    method = "；".join(method_bits)

    limitations = []
    if rows:
        limitations += [
            "**一阶引擎不是因果预测**：边来自共表达相关（07_grn 没有 motif 剪枝），"
            "相关不等于调控",
            "**一阶引擎只有一阶、单跳**：不走多跳传播 —— 在相关网络上做多跳会放大噪声，"
            "看起来像网络效应，其实只是把相关系数乘了几遍",
            ("**一阶引擎的 `ko_magnitude` 里混了两个因素**：它是 "
             "`|Δz_G| x network_sensitivity`，所以排在前面的基因往往只是"
             "**在该细胞类型里表达高**，而不是在网络里被连得紧。"
             "要看后者请用 `network_sensitivity` 列（它与表达无关）"),
        ]
    if td_ok:
        limitations += [
            ("**scTenifoldKnk 没有细胞类型分辨率**：它整份数据只建一个网络，"
             "输出是「扰动基因 x 网络基因」的全局距离。本仓库**没有**把距离按"
             "细胞类型重新加权 —— 那会造出一个既不是 scTenifoldKnk、"
             "也不是本仓库一阶近似的新方法，然后借它的名字发出去"),
            (f"**scTenifoldKnk 的网络基因集是子集**：{tenifold['n_genes_network']} 个"
             f"（上限 {tparams['max_genes']}），不是全基因集 —— 距离只在网络内可比。"
             f"包内固定 `set.seed(1)`，与本仓库 analysis.seed 无关"),
            ("**scTenifoldKnk 的输入是原始计数**（`qc_filtered.h5ad`，全基因集），"
             "且 `qc=FALSE` —— 本仓库 01_qc 已做过 QC，再滤会改变细胞集，"
             "两法就不可比了"),
        ]
        if tenifold.get("n_empty_knockout"):
            limitations.append(
                f"**scTenifoldKnk 有 {tenifold['n_empty_knockout']} 个候选基因是"
                f"「空敲除」**（{'、'.join(tenifold.get('empty_knockout_genes') or [])}）："
                f"它们在去噪后的网络里**出度为 0**，而包的敲除方式是把网络里该基因的"
                f"那一行清零 —— 清一行本来就全 0 的行等于什么都没敲，`KO` 与 `WT` "
                f"逐位相同，返回的\"距离\"只剩浮点噪声（实测 ~1e-16）。"
                f"本仓库把这些行置空并标了 `tenifold_empty_knockout` 列，"
                f"**它们不表示「敲除没有影响」，只表示这个网络表达不了该扰动**；"
                f"要拿到真实效应，该基因需要先进入网络（调大 "
                f"`perturbation.tenifold.max_genes`）或换网络推断方式")
    limitations += [
        "表达层面而非蛋白层面：TF 的 mRNA 与其活性经常不相关",
        ("两法共享同一个上游（07_grn 的共表达边），一致性高可能只是"
         "**共享了同一个偏差**；一致性只说明排序像不像，不说明哪个对"),
        (f"signature_alignment 的可靠性上界是 Part 1 签名的质量；"
         f"本轮候选来自 {target_source}"
         + ("，**没有疾病签名，该列为空**" if n_align == 0 else "")),
    ]

    status = {
        "dataset_id": cfg["dataset_id"],
        "status": "ok",
        "engine": engine,
        "engines_used": engines_used,
        "n_candidates": int(targets["gene"].nunique()),
        "n_cell_types": int(df["cell_type"].nunique()) if rows else 0,
        "n_rows": int(len(df)),
        "n_with_signature_alignment": n_align,
        "target_source": target_source,
        "group_key": ckey,
        "top_by_effect": df_to_records(best.head(10)) if rows else [],
        "skipped": skipped,
        "tools": tools,
        "tools_all_unavailable": all(not v["available"] for v in tools.values()),
        "method": method,
        "limitations": limitations,
    }
    if tenifold is not None:
        status["tenifold"] = {k: v for k, v in tenifold.items() if k != "table"}
    if cmp is not None:
        status["engine_comparison"] = cmp
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
