#!/usr/bin/env python3
"""
main_analysis.py — 编排全部步骤 + 验收清单

**验收清单是这一步存在的理由。** 单细胞流水线的失败模式大多是
"job 绿了但结果不对"：某个可选步骤静默跳过、某张图是空白的、
某个中间产物没写出来。这些在 CI 日志里和成功长得一模一样。

所以最后强制检查：
  - 必需步骤是否都 ok
  - 必需产物文件是否都在（且非空）
  - 图是否都真的画出了东西（交给 tools/check_figures.mjs 做像素级检查）

**任何必需项失败 → 退出码 1 → CI 红。** 可选步骤失败只在报告里标出。
"""

from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts" / "lib"))

from common import (capture_versions, init_manifest, load_config,  # noqa: E402
                    log_error, log_info, log_warn, manifest_path,
                    manifest_summary, named_tools_note, parse_args,
                    probe_named_tools, read_json, read_manifest, read_state,
                    record_decision, record_human_review, record_input,
                    record_params, record_step, set_orchestrated, write_json)

# (步骤 id, 模块文件, 函数名, 是否必需, 中文名)
STEPS = [
    ("fetch",            "00_fetch.py",             "run_00_fetch",            True,  "取数与校验"),
    ("qc",               "01_qc.py",                "run_01_qc",               True,  "质量控制"),
    ("integrate",        "02_integrate.py",         "run_02_integrate",        True,  "标准化与降维"),
    ("cluster_annotate", "03_cluster_annotate.py",  "run_03_cluster_annotate", True,  "聚类与注释"),
    ("pseudobulk_de",    "04_pseudobulk_de.py",     "run_04_pseudobulk_de",    False, "拟bulk 差异表达"),
    ("trajectory",       "05_trajectory.py",        "run_05_trajectory",       False, "轨迹推断"),
    ("communication",    "06_communication.py",     "run_06_communication",    False, "细胞通讯"),
    ("grn",              "07_grn.py",               "run_07_grn",              False, "转录因子调控"),
    ("virtual_perturbation", "08_virtual_perturbation.py",
     "run_08_virtual_perturbation", False, "虚拟敲除/过表达（§1.7/§1.8 保留框架）"),
]

# 每步会写的状态文件。**跑之前先删掉** —— 否则步骤崩溃时旧文件还在，
# 下游"产物存在"检查读的是**上一轮的**结果，会给出虚假的通过。
# 实测踩过：07_grn 因漏 import 崩了，而 grn_status.json 是上一轮的，
# 验收照样 40 项全绿。
STEP_STATUS_FILES = {
    "qc": "qc_status.json",
    "integrate": "integration_status.json",
    "cluster_annotate": "cluster_status.json",
    "pseudobulk_de": "pseudobulk_status.json",
    "trajectory": "trajectory_status.json",
    "communication": "communication_status.json",
    "grn": "grn_status.json",
    "virtual_perturbation": "virtual_perturbation_status.json",
}

# 文档 §2「本部分人工复核节点」。**默认 pending，不是 confirmed** ——
# 自动化流水线不能替人签字，把未确认的节点记成已确认，等于把复核节点
# 变成摆设。验收里作为**可见但不阻断**的项列出（required 只标"这节点
# 是否适用本数据集"）。
HUMAN_REVIEW_NODES = [
    ("celltype_labels",      "细胞类型注释最终标签",              True),
    ("cluster_resolution",   "聚类分辨率选择依据",                True),
    ("pseudobulk_design",    "拟bulk 差异分析设计",               False),
    ("trajectory_direction", "拟时序轨迹方向确认（marker 验证）",  True),
    ("trajectory_branches",  "拟时序分支点的生物学解释",           True),
    # §1.7/§1.8 的人工复核节点：**虚拟扰动的靶基因在生物学上是否讲得通**。
    # 这是整个虚拟扰动框架里唯一能挡住"算出来一个数就当真"的机制 ——
    # 一阶网络模型必然会给出一份排名，排名本身不含任何合理性判据。
    ("virtual_perturbation_targets", "虚拟扰动靶基因的生物学合理性", False),
]

# 需要登记哈希的输入（相对 data_dir）。(文件名, 中文说明, 是否必需)
INPUT_FILES = [
    ("raw.h5ad",          "原始计数矩阵", True),
    ("dataset_info.json", "数据集元信息",  True),
    ("qc_filtered.h5ad",  "QC 后矩阵",    True),
    ("integrated.h5ad",   "整合后矩阵",    True),
    ("clustered.h5ad",    "聚类后矩阵",    True),
]

# 必需产物（相对 results_dir）。(文件名, 中文说明, 是否必需)
REQUIRED_FILES = [
    ("state.json",                 "步骤状态",           True),
    ("qc_status.json",             "QC 结论",            True),
    ("integration_status.json",    "降维与整合结论",      True),
    ("cluster_status.json",        "聚类与注释结论",      True),
    ("markers_all.csv",            "marker 基因表",       True),
    ("cluster_resolution_scan.csv", "分辨率扫描",         True),
    ("celltype_annotation.csv",    "细胞类型注释",        True),
    ("celltype_scores.csv",        "细胞类型打分矩阵",     True),
    ("qc_cells.csv",               "每细胞 QC 指标",      True),
    ("pseudobulk_status.json",     "拟bulk 状态",         False),
    ("trajectory_status.json",     "轨迹状态",            False),
    ("communication_status.json",  "通讯状态",            False),
    ("grn_status.json",            "GRN 状态",            False),
    ("virtual_perturbation_status.json", "虚拟扰动状态",    False),
]

# 必需的图（相对 figures_dir，不含扩展名）
REQUIRED_FIGURES = [
    ("02-01-01-unit1-qc-violin-before",          "过滤前 QC 分布"),
    ("02-01-02-unit1-qc-scatter-thresholds",     "QC 阈值散点"),
    ("02-02-01-unit1-hvg-selection",             "高变基因选择"),
    ("02-02-02-unit1-pca-variance-ratio",        "PCA 方差解释"),
    ("02-03-01-unit1-cluster-resolution-scan",   "分辨率扫描曲线"),
    ("02-03-02-unit1-umap-clusters",             "UMAP 聚类图"),
    ("02-03-03-unit1-markers-dotplot",           "marker 点图"),
    ("02-03-04-unit1-celltype-scores-heatmap",   "细胞类型打分热图"),
]


def load_step_fn(module_file: str, fn_name: str):
    spec = importlib.util.spec_from_file_location(
        module_file.replace(".py", ""), REPO / "scripts" / module_file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return getattr(mod, fn_name)


def has_file(p: Path) -> bool:
    return p.exists() and p.stat().st_size > 0


def run_all(cfg: dict, only: list = None) -> int:
    set_orchestrated(True)
    res_dir = Path(cfg["output"]["results_dir"])
    res_dir.mkdir(parents=True, exist_ok=True)

    log_info("=" * 68)
    log_info(f"单细胞流水线 | 数据集 {cfg['dataset_id']}")
    log_info(f"结果目录 {res_dir}")
    log_info("=" * 68)

    # ---- 模块零：建立本轮运行清单（§0.3 / §0.4）-----------------------------
    # **必须在任何步骤之前建，且先清掉上一轮** —— 清单描述的是本轮。
    # 清单和 state.json 分开：state 记"跑没跑成"（每步重写），
    # 清单记"在什么条件下跑出来的"（证据，写入后不该再变）。
    init_manifest(cfg)
    capture_versions(cfg)
    record_params(cfg, {
        "seed": cfg.get("analysis", {}).get("seed"),
        "qc": cfg.get("qc", {}),
        "cluster": cfg.get("cluster", {}),
        "integration": cfg.get("integration", {}),
        "trajectory": cfg.get("trajectory", {}),
        # **这三段原先漏了。** §0.3 要求"全部参数含 seed"入清单 ——
        # 漏掉的可选步骤参数意味着那几步的结果无法被复现，
        # 而清单看起来是完整的（缺的是键，不是值）。
        "communication": cfg.get("communication", {}),
        "grn": cfg.get("grn", {}),
        "perturbation": cfg.get("perturbation", {}),
    })
    for node, label, req in HUMAN_REVIEW_NODES:
        record_human_review(cfg, node, required=req, status="pending",
                            note=f"{label} —— 需人工确认，本轮自动化未确认")

    # ---- §2 点名工具的缺口登记 ------------------------------------------
    #
    # **为什么放在清单里，而不是各步骤的 status 文件里：**
    # "R 包一个都跑不了"是本仓库**整轮运行**的属性（CI 没有 R + rpy2），
    # 不是某一步的属性。放进清单只写一次，不会出现"6 个状态文件里
    # 有 5 个写了、1 个漏了"这种半对半错的状态。
    #
    # 它必须存在，因为产物**看不出来**：02 有归一化、04 有差异表、
    # 05 有四条轨迹、07 有调控子 —— 每一步都有东西，
    # 所以"点名的方法一个都没用上"这件事得自己说出来。
    named = probe_named_tools(log=log_warn)
    record_decision(
        cfg, "named_tools",
        "文档 §2.1–§2.8 点名的工具，哪些真的用上了？",
        (f"{sum(1 for i in named.values() if i['available'])}/{len(named)} 个可用；"
         "本仓库实际使用的是内置/替代实现"),
        evidence=named_tools_note(),
    )
    record_params(cfg, {"named_tools": named})
    log_info("")
    log_info("§2 点名工具的使用情况（**多数没用上，这是缺口不是已覆盖**）：")
    for tool, info in named.items():
        mark = "可用" if info["available"] else "未使用"
        log_info(f"  [{mark}] {info['section']} {tool}（{info['kind']}）"
                 f"—— {info['reason'][:66]}")
    log_info(f"运行清单：{manifest_path(cfg)}")

    failed_required = []
    for sid, mfile, fn, required, label in STEPS:
        if only and sid not in only:
            log_info(f"--- 跳过 {label} ({sid})：不在 --steps 里")
            continue
        log_info("")
        log_info(f"--- {label} ({sid})" + ("" if required else "  [可选]"))
        # 先删本步的状态文件：崩溃时不留旧文件冒充本轮结果
        stale = res_dir / STEP_STATUS_FILES.get(sid, "")
        if sid in STEP_STATUS_FILES and stale.exists():
            stale.unlink()
        t0 = time.time()
        try:
            fn = load_step_fn(mfile, fn)
            fn(cfg)
            record_step(cfg, sid, "ok", time.time() - t0, required=required)
        except Exception as e:  # noqa: BLE001
            msg = f"{type(e).__name__}: {e}"
            record_step(cfg, sid, "failed", time.time() - t0,
                        message=msg, required=required)
            if required:
                log_error(f"{label} 失败（必需）: {msg}")
                failed_required.append((sid, label, msg))
            else:
                # **可选步骤失败不让 job 变红，但必须显眼。**
                # Part 1 的教训：可选步骤失败时 job 仍然是绿的。
                log_warn(f"{label} 失败（可选，不影响 job 结论）: {msg}")

    # ---- 模块零：登记输入哈希（§0.4）----------------------------------------
    # 放在所有步骤之后 —— 可选步骤的产物这轮有没有，跑完才知道。
    data_dir = Path(cfg["output"]["data_dir"])
    for fname, desc, req in INPUT_FILES:
        e = record_input(cfg, data_dir / fname, label=f"{desc} ({fname})",
                         required=req)
        if e["status"] == "missing" and req:
            log_warn(f"输入缺失：{desc} ({fname})")
    msum = manifest_summary(cfg)
    log_info(f"输入登记 {msum['n_inputs']} 项"
             + (f"，其中缺失 {len(msum['inputs_missing'])} 项"
                if msum["inputs_missing"] else "，全部就位"))

    # ---- 验收清单 ----------------------------------------------------------
    log_info("")
    log_info("=" * 68)
    log_info("验收清单")
    log_info("=" * 68)

    fig_dir = Path(cfg["output"]["figures_dir"])
    checks = []

    for sid, mfile, fn, required, label in STEPS:
        if only and sid not in only:
            continue
        st = {s["id"]: s for s in read_state(cfg).get("steps", [])}.get(sid, {})
        status = st.get("status", "not_run")
        ok = (status == "ok")
        # **`failed` 不等于"正确地跳过"。**
        # `not_configured`/`not_applicable`/`disabled`/`not_run` 是**设计如此**
        # （如没有分组信息就不做拟bulk DE），可选步骤这样算通过。
        # 但 `failed` 是**崩了** —— 上一版把 failed 也归进 ok=True，
        # 结果是 07_grn 因漏 import 崩溃而验收全绿、产物还是旧的。
        if not required and status in ("not_configured", "not_applicable",
                                       "disabled", "not_run"):
            ok = True
        checks.append({"item": f"步骤 {label}", "ok": ok,
                       "required": required, "detail": status})

    for fname, desc, required in REQUIRED_FILES:
        ok = has_file(res_dir / fname)
        checks.append({"item": f"产物 {desc} ({fname})", "ok": ok,
                       "required": required,
                       "detail": "存在" if ok else "**缺失或为空**"})

    for fname, desc in REQUIRED_FIGURES:
        ok = has_file(fig_dir / f"{fname}.png")
        checks.append({"item": f"图 {desc} ({fname}.png)", "ok": ok,
                       "required": True,
                       "detail": "存在" if ok else "**缺失**"})

    # ---- 模块零：运行清单（§0.3 / §0.4）-------------------------------------
    # 清单缺项不是"分析错了"，而是"这轮跑出来的东西没法追溯"。
    # **人工复核未确认不算失败** —— 默认就是 pending，那是设计如此；
    # 把它判成 FAIL 会让每个 job 都红，反而没人看。但必须可见。
    checks.append({
        "item": "运行清单存在（run_manifest.json）",
        "ok": msum.get("present", False), "required": True,
        "detail": (f"{msum.get('n_versions', 0)} 个包版本、"
                   f"{msum.get('n_inputs', 0)} 项输入、"
                   f"{msum.get('n_decisions', 0)} 条决策"
                   if msum.get("present") else "**缺失**"),
    })
    if msum.get("present"):
        checks.append({
            "item": "版本记录非空（pip freeze 全量）",
            "ok": msum["n_versions"] >= 20, "required": True,
            "detail": f"{msum['n_versions']} 个已安装包",
        })
        checks.append({
            "item": "输入哈希已登记且必需项无缺失",
            # **只看 required 的缺失。** 可选项（如 clinical.csv）本来就可以
            # 不存在，算进来会让没有该文件的数据集全部误判失败。
            "ok": msum["n_inputs"] >= len(INPUT_FILES)
                  and not msum["inputs_missing_required"],
            "required": True,
            "detail": (f"{msum['n_inputs']} 项"
                       + (f"，必需缺失 {','.join(msum['inputs_missing_required'])}"
                          if msum["inputs_missing_required"] else "，必需项齐全")
                       + (f"（可选缺失 {','.join(msum['inputs_missing'])}）"
                          if msum["inputs_missing"] else "")),
        })
        pend = msum["human_review_pending"]
        checks.append({
            "item": f"人工复核节点待确认（{len(pend)} 个，不阻断 job）",
            "ok": True, "required": False,
            "detail": (", ".join(pend) if pend else "全部已确认"),
        })
        # ---- §0.2 跨语言转换：**空数组必须被解释** --------------------------
        #
        # 本仓库唯一一处真正的跨语言交接是 08 读 Part 1 的 `part2_targets.csv`。
        # 默认配置 `targets_csv: null` 走回退分支，所以 `cross_language` 是空的
        # —— 而**空数组和"忘了记"长得一模一样**（geo 的 `record_decision()`
        # 定义了却没调用点，就是这么藏了一整轮）。
        #
        # 判据不是"必须有记录"，而是"**空的话必须有解释**"：
        # 要么 `cross_language` 非空，要么决策链里有一条说明为什么。
        # 将来真加了跨语言步骤却忘了记，这条会立刻变红。
        #
        # **读全量清单而不是 `msum`** —— `manifest_summary()` 只给计数。
        #
        # 用 `read_manifest(cfg)` 而不是 `read_json(res_dir /
        # "run_manifest.json")`：两者等价，但**文件名只该有一处** ——
        # `MANIFEST_NAME` 改了而这里写死字符串的话，这条检查会静默地
        # 读一个不存在的文件、`_dec_nodes` 变空，然后报"没解释"。
        _full = read_manifest(cfg)
        _cl = msum.get("n_cross_language", 0)
        _dec_nodes = {d.get("node") for d in (_full.get("decisions") or [])}
        checks.append({
            "item": "§0.2 跨语言交接已登记或有解释",
            "ok": _cl > 0 or "cross_language" in _dec_nodes,
            "required": False,
            "detail": (f"{_cl} 条跨语言转换记录（交接表 + 丢失字段）"
                       if _cl else
                       ("0 条，**但决策链里已说明本轮没有交接**"
                        if "cross_language" in _dec_nodes else
                        "**0 条且没有任何解释** —— 读者无法区分"
                        "『本轮没配』和『忘了记』")),
        })

    # 可选步骤的"没做"要在报告里可见 —— 不能只是绿
    for sid, fname in (("pseudobulk_de", "pseudobulk_status.json"),
                       ("trajectory", "trajectory_status.json"),
                       ("communication", "communication_status.json"),
                       ("grn", "grn_status.json"),
                       ("virtual_perturbation", "virtual_perturbation_status.json")):
        d = read_json(res_dir / fname)
        if not d:
            continue
        st = d.get("status")
        note = d.get("reason") or d.get("method") or ""
        checks.append({"item": f"可选步骤状态 {sid} = {st}", "ok": True,
                       "required": False, "detail": str(note)[:150]})

    # ---- 内容级检查：状态文件在 ≠ 结果是对的 --------------------------------
    #
    # Part 1 的教训：TF 第一次跑时 tf_status.json 正常产出、验收全绿，
    # 而 figure_written 其实是 false（图一张没出），同时 121 个样本因列名
    # 写错全部匹配失败、组间比较静默为空。
    #
    # 所以这里**读 status 里的真实字段**，而不是只看文件在不在。
    tj = read_json(res_dir / "trajectory_status.json") or {}
    if tj.get("status") == "ok":
        n_m = int(tj.get("n_methods") or 0)
        checks.append({
            "item": "轨迹用了 >=2 种方法交叉验证",
            "ok": n_m >= 2, "required": True,
            "detail": (f"{n_m} 种：{','.join(tj.get('methods_ok', {}).keys())}"
                       if n_m >= 2 else
                       f"**只有 {n_m} 种** —— 单一方法的拟时序是某个算法的一次"
                       f"输出，不是数据里的结构"),
        })
        cv = tj.get("cross_validated_methods") or []
        mean_rho = tj.get("method_correlation_mean_offdiag")
        checks.append({
            "item": "轨迹方法间一致性已量化（且排除方向参考）",
            "ok": bool(cv) and mean_rho is not None,
            "required": True,
            "detail": (f"交叉验证 {len(cv)} 种，平均 rho={mean_rho:+.4f}，"
                       f"参考方法 {tj.get('direction_reference_method')}"
                       if bool(cv) and mean_rho is not None else
                       "**缺 cross_validated_methods 或一致性数值**"),
        })
        checks.append({
            "item": "轨迹方向来源已写明",
            "ok": bool(tj.get("direction_source")),
            "required": True,
            "detail": str(tj.get("direction_source") or "**缺失**"),
        })
        checks.append({
            "item": "scVelo 的不可得已如实记录",
            "ok": bool((tj.get("scvelo") or {}).get("status")),
            "required": False,
            "detail": str((tj.get("scvelo") or {}).get("status", "**缺失**")),
        })
        checks.append({
            "item": "轨迹方法学限定已写明（>=4 条）",
            "ok": len(tj.get("limitations") or []) >= 4,
            "required": True,
            "detail": f"{len(tj.get('limitations') or [])} 条",
        })
        for fn, desc in (("02-05-02-unit1-trajectory-method-correlation", "方法一致性矩阵"),
                         ("02-05-03-unit1-trajectory-modules", "沿轨迹基因模块")):
            ok = has_file(fig_dir / f"{fn}.png")
            checks.append({"item": f"图 {desc} ({fn}.png)", "ok": ok,
                           "required": True,
                           "detail": "存在" if ok else "**缺失**"})

    # ---- 文档指定工具的落地情况（§2.4 CellTypist / §2.7 LIANA）-------------
    #
    # **不阻断 job，但没跑成就要红字显示。** 这两个都是 pip 能装的，
    # 装了就该跑；没跑成（模型下载失败 / 版本不兼容）是环境问题不是分析
    # 错了，所以 required=False。但绝不能混在绿字里 ——
    # 规则 4：可选步骤的"没做"必须和"做了没问题"长得不一样。
    cj = read_json(res_dir / "cluster_status.json") or {}
    ann = cj.get("annotation") or {}
    ct = ann.get("celltypist") or {}
    ct_cmp = ann.get("celltypist_vs_marker") or {}

    # ---- §0.2 交接的**产出侧**：Part 3 拿 `clustered.h5ad` 当参考 ----------
    #
    # Part 3 的 `deconvolution.reference: h5ad` 要的就是这一份文件。
    # **交接的两端都要能被核对** —— 只检查消费侧（Part 3 记没记）的话，
    # 产出侧悄悄丢掉 counts 层不会被任何人发现，而 Part 3 会静默退回
    # 用 `.X`（log 值）去求参考谱，解出的比例没有意义。
    #
    # 所以这里检查产出侧的契约：**必须有 counts 层**（NNLS 要计数），
    # **必须有细胞类型列**（Part 3 的 `celltype_key` 要指名一列）。
    # 两条都不阻断 job（Part 3 可以不用这份文件），但缺了要红字显示。
    _p3 = cj.get("part3_reference") or {}
    _p3c = _p3.get("contract") or {}
    _p3_missing = [k for k, v in (("layers['counts']", _p3c.get("layers['counts']")),
                                  ("obs 细胞类型列", _p3c.get("celltype_column")))
                   if not v]
    checks.append({
        "item": "§0.2 Part 3 参考导出的契约（counts 层 + 细胞类型列）",
        "ok": bool(_p3) and not _p3_missing, "required": False,
        "detail": (f"契约完整：counts 层 + `{_p3c.get('celltype_column')}` 列，"
                   f"{_p3.get('n_cells')} 细胞 x {_p3.get('n_genes')} 基因"
                   f"（Part 3 配 `reference: h5ad` 即可直接用）"
                   if _p3 and not _p3_missing else
                   (f"**契约缺**：{_p3_missing} —— Part 3 拿这份文件当参考时"
                    f"会静默算错（退回 log 值）或 KeyError"
                    if _p3 else
                    "**没有 part3_reference 记录** —— cluster_status.json 是"
                    "上一轮的旧文件，或这一步没跑完")),
    })
    checks.append({
        "item": "CellTypist 自动注释（§2.4）",
        "ok": ct.get("status") == "ok", "required": False,
        "detail": (f"ok：{ct.get('n_labels')} 种标签，模型 {ct.get('model')}，"
                   f"{ct.get('n_genes_input')} 基因输入"
                   if ct.get("status") == "ok" else
                   f"**未跑成**（{ct.get('status')}：{str(ct.get('reason'))[:120]}）"
                   " —— marker 打分仍在，但少了第二条独立证据"),
    })
    if ct.get("status") == "ok":
        checks.append({
            "item": "CellTypist 与 marker 注释的一致性已量化（§2.4）",
            "ok": bool(ct_cmp.get("compared")), "required": False,
            "detail": (f"簇层面一致 {ct_cmp.get('n_agree_exact')}/"
                       f"{ct_cmp.get('n_clusters')}"
                       f"（{ct_cmp.get('agreement_frac')}）"
                       if ct_cmp.get("compared") else
                       f"**未对比**：{ct_cmp.get('reason')}"),
        })

    cm = read_json(res_dir / "communication_status.json") or {}
    li = cm.get("liana") or {}
    li_cmp = cm.get("liana_vs_builtin") or {}
    checks.append({
        "item": "LIANA rank_aggregate（§2.7 指定的主工具）",
        "ok": li.get("status") == "ok", "required": False,
        "detail": (f"ok：{li.get('n_rows')} 行，v{li.get('version')}，"
                   f"api {li.get('api')}"
                   if li.get("status") == "ok" else
                   f"**未跑成**（{li.get('status')}：{str(li.get('reason'))[:120]}）"
                   " —— 自建共表达打分仍在，但它不是 consensus rank aggregate"),
    })
    if li.get("status") == "ok":
        checks.append({
            "item": "LIANA 与自建打分的差异已量化（§2.7）",
            "ok": bool(li_cmp.get("compared")), "required": False,
            "detail": (f"共同组合 {li_cmp.get('n_common_combinations')}，"
                       f"Spearman rho={li_cmp.get('spearman_rho')}，"
                       f"top25 重叠 {li_cmp.get('top25_overlap')}/25"
                       if li_cmp.get("compared") else
                       f"**未对比**：{li_cmp.get('reason')}"),
        })

    # ---- §1.7 / §1.8 虚拟扰动（保留框架）------------------------------------
    #
    # 这一节的**重点是"没做什么"**。规范点名的三个工具在本环境全都装不了，
    # 而"装不了"和"没装"是两件事 —— 前者有确切原因，后者是疏忽。
    # 所以逐个列出原因，并且把"用的是自建一阶近似"这件事写在最显眼处。
    vp = read_json(res_dir / "virtual_perturbation_status.json") or {}
    vp_tools = vp.get("tools") or {}
    if vp_tools:
        unavail = {k: v.get("reason") for k, v in vp_tools.items()
                   if not v.get("available")}
        checks.append({
            "item": f"§1.7/§1.8 点名工具不可用的原因已逐个记录（{len(unavail)}/{len(vp_tools)} 个不可用）",
            # 判据是"每个不可用的都有原因"，不是"全都不可用" ——
            # 将来某个工具能装了，这条应该自动变成 PASS 而不是 FAIL。
            "ok": all(bool(r) for r in unavail.values()) if unavail else True,
            "required": False,
            "detail": ("；".join(f"{k}：{str(r)[:70]}" for k, r in unavail.items())
                       if unavail else "三个工具都可用"),
        })
    if vp.get("status") == "ok":
        checks.append({
            "item": "虚拟扰动的候选来源已写明（Part 1 交接 vs 内部回退）",
            "ok": bool(vp.get("target_source")), "required": False,
            "detail": (f"{vp.get('target_source')}，"
                       f"{vp.get('n_candidates')} 个候选 x "
                       f"{vp.get('n_cell_types')} 个细胞类型"
                       + ("（**内部回退：没有 Part 1 签名，"
                          "signature_alignment 为空是预期的**）"
                          if str(vp.get("target_source", "")).startswith("internal")
                          else "")),
        })
        checks.append({
            "item": "虚拟扰动的方法学限定已写明（>=5 条）",
            "ok": len(vp.get("limitations") or []) >= 5, "required": False,
            "detail": f"{len(vp.get('limitations') or [])} 条",
        })
        checks.append({
            "item": "虚拟扰动明确声明不是 scTenifoldKnk / PerturbNet",
            # 这一条是防"报了个数就被当成因果预测"的。方法串里必须出现
            # 这两个否定，否则读者会以为跑的是规范点名的工具。
            "ok": all(t in str(vp.get("method", ""))
                      for t in ("scTenifoldKnk", "PerturbNet")),
            "required": False,
            "detail": str(vp.get("method", ""))[:120],
        })
        _vp_fig = "02-08-01-unit1-virtual-perturbation-effect.png"
        ok_fig = has_file(fig_dir / _vp_fig)
        checks.append({"item": f"图 虚拟敲除效应 ({_vp_fig})",
                       "ok": ok_fig, "required": False,
                       "detail": "存在" if ok_fig else "**缺失**"})

    # ---- §2 点名工具的缺口登记 -------------------------------------------
    #
    # **判据是"理由写了没有"，不是"工具跑了没有"。**
    # 将来某个工具能装了（比如 CI 换成带 R 的镜像），这条应该依然 PASS，
    # 而不是因为 `available=False` 就变红 —— 那会把"如实记录"惩罚成失败。
    #
    # 必须存在的原因是产物**看不出来**：每一步都有东西产出，
    # 而 §2 点名的 R 包（SCTransform / scran / DESeq2 / Monocle3 /
    # Slingshot / CellChat / SoupX）一个都跑不了。
    # 用 `read_manifest(cfg)` 而不是拼 `res_dir / "run_manifest.json"` ——
    # 文件名只该有一处（`MANIFEST_NAME`）。`read_manifest` 在文件不存在时
    # 已经返回 `{}`，所以这里不需要 `has_file` 判断。
    msum = read_manifest(cfg)
    named = (msum.get("params") or {}).get("named_tools")
    if not isinstance(named, dict) or not named:
        checks.append({
            "item": "§2 点名工具的缺口已登记（清单 named_tools）",
            "ok": False, "required": False,
            "detail": "清单里没有 named_tools —— 读者会以为 §2 点名的方法都用上了",
        })
    else:
        no_reason = [t for t, i in named.items() if not (i or {}).get("reason")]
        no_kind = [t for t, i in named.items() if not (i or {}).get("kind")]
        checks.append({
            "item": "§2 点名工具的缺口已登记（清单 named_tools）",
            "ok": not no_reason and not no_kind,
            "required": False,
            "detail": (f"{len(named)} 个点名工具已登记，"
                       f"{sum(1 for i in named.values() if i.get('available'))} 个当前可用"
                       if not no_reason and not no_kind else
                       f"无理由 {no_reason}；无分类 {no_kind}"),
        })
        # 顶名的包必须被标成 name_taken —— 这是最危险的一类
        squat = [t for t, i in named.items() if (i or {}).get("kind") == "name_taken"]
        checks.append({
            "item": "PyPI 同名无关包已标为 name_taken（edgeR / Slingshot）",
            "ok": set(squat) >= {"edgeR", "Slingshot"},
            "required": False,
            "detail": f"已标注: {sorted(squat)}" if squat else "**一个都没标**",
        })
        # 决策链也要留痕（§0.4）
        dec = msum.get("decisions") or []
        checks.append({
            "item": "点名工具的使用情况进了决策链（§0.4）",
            "ok": any((d or {}).get("node") == "named_tools" for d in dec),
            "required": False,
            "detail": f"decisions 里 {len(dec)} 条，named_tools "
                      f"{'在' if any((d or {}).get('node') == 'named_tools' for d in dec) else '**不在**'}",
        })

    n_fail = 0
    for c in checks:
        mark = "PASS" if c["ok"] else "FAIL"
        if not c["ok"] and c["required"]:
            n_fail += 1
        if c["ok"] and c["required"]:
            log_info(f"  [{mark}] {c['item']}")
        elif not c["ok"]:
            log_error(f"  [{mark}] {c['item']}  {c['detail']}")
        else:
            log_info(f"  [{mark}] {c['item']}  {c['detail']}")

    summary = {
        "dataset_id": cfg["dataset_id"],
        "n_checks": len(checks),
        "n_failed_required": n_fail,
        "checks": checks,
        "steps_failed_required": [{"id": i, "label": l, "error": m}
                                  for i, l, m in failed_required],
        "manifest": msum,
        "verdict": "ok" if n_fail == 0 else "failed",
    }
    write_json(res_dir / "acceptance.json", summary)

    log_info("")
    if n_fail == 0:
        log_info(f"验收通过：{len(checks)} 项检查，0 项必需失败")
    else:
        log_error(f"验收失败：{n_fail} 项必需检查未通过")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    args = parse_args()
    cfg = load_config(args.config)
    only = [s.strip() for s in args.steps.split(",")] if args.steps else None
    sys.exit(run_all(cfg, only))
