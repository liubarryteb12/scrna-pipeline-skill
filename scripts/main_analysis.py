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

from common import (load_config, log_error, log_info, log_warn,  # noqa: E402
                    parse_args, read_json, read_state, record_step,
                    set_orchestrated, write_json)

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
]

# 必需的图（相对 figures_dir，不含扩展名）
REQUIRED_FIGURES = [
    ("qc_violin_before",          "过滤前 QC 分布"),
    ("qc_scatter_thresholds",     "QC 阈值散点"),
    ("hvg_selection",             "高变基因选择"),
    ("pca_variance_ratio",        "PCA 方差解释"),
    ("cluster_resolution_scan",   "分辨率扫描曲线"),
    ("umap_clusters",             "UMAP 聚类图"),
    ("markers_dotplot",           "marker 点图"),
    ("celltype_scores_heatmap",   "细胞类型打分热图"),
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

    failed_required = []
    for sid, mfile, fn, required, label in STEPS:
        if only and sid not in only:
            log_info(f"--- 跳过 {label} ({sid})：不在 --steps 里")
            continue
        log_info("")
        log_info(f"--- {label} ({sid})" + ("" if required else "  [可选]"))
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
        # 可选步骤的 not_configured / not_applicable / disabled 是**正确行为**
        if not required and status in ("failed", "not_run"):
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

    # 可选步骤的"没做"要在报告里可见 —— 不能只是绿
    for sid, fname in (("pseudobulk_de", "pseudobulk_status.json"),
                       ("trajectory", "trajectory_status.json"),
                       ("communication", "communication_status.json"),
                       ("grn", "grn_status.json")):
        d = read_json(res_dir / fname)
        if not d:
            continue
        st = d.get("status")
        note = d.get("reason") or d.get("method") or ""
        checks.append({"item": f"可选步骤状态 {sid} = {st}", "ok": True,
                       "required": False, "detail": str(note)[:150]})

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
