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
import re
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts" / "lib"))

from common import (capture_versions, classify_step_result, init_manifest,  # noqa: E402
                    load_config, log_error, log_info, log_warn, manifest_path,
                    manifest_summary, named_tools_note, parse_args,
                    probe_named_tools, read_json, read_manifest, read_state,
                    record_decision, record_human_review, record_input,
                    record_params, record_step, result_status_of,
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

# **M10（R-03 裁决）**：这张表的键必须与 `STEPS` 的 id 集合一致。
#
# 旧代码 `stale = res_dir / STEP_STATUS_FILES.get(sid, "")` 用的是**默认空串**
# —— `res_dir / ""` 就是 `res_dir` 本身，一旦上面那个 `if sid in STEP_STATUS_FILES`
# 守卫被改动（比如将来有人加了新步骤却忘了登记状态文件），`unlink()` 就会
# 作用在**结果目录**上，`IsADirectoryError` 被 `except` 接住 ⇒ 步骤莫名失败。
# 守卫现在挡着，所以是**潜伏**缺陷，但"靠守卫恰好成立"不是设计。
#
# 这里在导入时就把两个集合对齐 —— 将来加步骤忘了登记，**启动即报错**。
#
# `fetch` 是唯一的例外，且是**设计如此**：它不写 `res_dir/*_status.json`，
# 写的是 `data_dir/dataset_info.json`（见 `00_fetch.py:338`）—— 取数步骤的
# 产物是数据本身，"状态"没有独立载体。所以它进豁免集，而不是给它编一个
# 空文件名（编了就会在跑前删掉一个不该删的东西）。
_STEPS_WITHOUT_STATUS_FILE = {"fetch"}

_MISSING_STATUS_FILES = sorted(
    {s[0] for s in STEPS} - set(STEP_STATUS_FILES) - _STEPS_WITHOUT_STATUS_FILE)
if _MISSING_STATUS_FILES:
    raise RuntimeError(
        f"STEPS 里的步骤既没有登记状态文件、也不在豁免集里: "
        f"{_MISSING_STATUS_FILES} —— 没有登记的话跑前删旧状态文件的逻辑"
        f"会落到结果目录上（审计 M10）。请在 STEP_STATUS_FILES 里补上，"
        f"或（确实不写状态文件时）加进 _STEPS_WITHOUT_STATUS_FILE 并写明理由。")

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
#
# **M11（R-03 裁决）：这张表不再单独承担判据。** 它是一张**手抄**的表，
# 而手抄的表会漂移 —— 实测漏了 `02-05-01-unit1-paga-graph`、
# `02-05-04-unit1/unit2/unit3`、`02-05-05-unit1` 与全部 `02-06-*`，
# 其中 `02-06-01-unit1-communication-heatmap` 在 `06_communication.py` 里
# 写出却**没有任何检查看得见**。漂移的方向恰好是"新加的图不在表里"，
# 也就是把 E-48 那个盲区原样再造一遍。
#
# 所以真正的判据换成 `declared_figures()`（从源码扫声明，见下），
# 这张表退化成**说明文字**：给每张图一个人话标题。表里缺条目不再是盲区
# （扫出来的图名不依赖它），但会让报告少一句解释 —— 所以下面有一条
# 检查专门盯"扫出来的图有没有说明"。
REQUIRED_FIGURES = [
    ("02-01-01-unit1-genes-detected",     "过滤前 QC: genes detected"),
    ("02-01-01-unit2-total-counts",       "过滤前 QC: total counts"),
    ("02-01-01-unit3-mito-fraction",      "过滤前 QC: mito fraction"),
    ("02-01-01-unit4-ribo-fraction",      "过滤前 QC: ribo fraction"),
    ("02-01-01-unit5-hb-fraction",        "过滤前 QC: hb fraction"),
    ("02-01-02-unit1-qc-scatter-thresholds",     "QC 阈值散点"),
    ("02-02-01-unit1-hvg-selection",             "高变基因选择"),
    ("02-02-02-unit1-pca-variance-ratio",        "PCA 方差解释"),
    ("02-02-03-unit1-batch-mixing",              "批次混合前后对比"),
    ("02-03-01-unit1-cluster-resolution-scan",   "分辨率扫描曲线"),
    ("02-03-02-unit1-umap-clusters",             "UMAP 聚类图"),
    ("02-03-03-unit1-markers-dotplot",           "marker 点图"),
    ("02-03-04-unit1-celltype-scores-heatmap",   "细胞类型打分热图"),
    ("02-05-01-unit1-paga-graph",                "PAGA 连接图（M11 补）"),
    ("02-05-02-unit1-trajectory-method-correlation", "轨迹方法一致性"),
    ("02-05-03-unit1-trajectory-modules-heatmap",    "轨迹模块热图"),
    ("02-05-03-unit2-trajectory-module-profiles",    "轨迹模块轮廓"),
    ("02-05-04-unit1-pseudotime-consensus",      "共识拟时序 UMAP"),
    ("02-05-04-unit2-pseudotime-dpt",            "DPT 拟时序 UMAP"),
    ("02-05-04-unit3-celltype-on-umap",          "细胞类型 UMAP（M6，条件产出）"),
    ("02-05-04-unit4-pseudotime-principal-path", "拟时序主路径+root"),
    ("02-05-05-unit1-pseudotime-by-cluster",     "各簇拟时序箱线图"),
    ("02-05-05-unit2-pseudotime-ridgeline",      "拟时序山脊图"),
    ("02-06-01-unit1-communication-heatmap",     "细胞通讯热图（M11 补）"),
    # Q-26 / E-48 补：这四条原来**不在这张表里**，所以"该有的图没有"这一整类
    # 问题没有任何检查看得见。02-07-01 的 5 张图（本图 + unit2..5 单 TF 面板）
    # 因 `figsize` 三元素元组从未产出过，而验收 70 项全绿。
    ("02-07-01-unit1-tf-activity-vs-pseudotime", "TF 活性沿拟时序"),
    ("02-07-02-unit1-tf-activity-heatmap",       "TF 活性热图"),
    ("02-07-03-unit1-tf-specificity-scatter",    "TF 特异性散点"),
    ("02-08-01-unit1-virtual-perturbation-effect", "虚拟扰动效应"),
]
_FIG_DESC = {nm: desc for nm, desc in REQUIRED_FIGURES}

PART = "02"


def _strip_comments(src: str) -> str:
    """逐行剥注释，**引号内不剥** —— 与 `tools/check_fig_names.mjs` 的
    `stripComments` 同义。

    口径必须一致：门禁层用 JS 那份扫"声明了哪些图"，验收层用这份扫，
    两边算法不同就会对同一份源码给出不同的图名集合，而**没有任何东西
    能发现它们不一致**（门禁绿、验收也绿）。

    与 E-62 同源：纯文本扫描器必须先把非代码区域抹掉，否则一个写注释里的
    图名会让靠括号配平的扫描器一路吞到文件尾。
    """
    out = []
    for line in src.split("\n"):
        q, cut = None, None
        for i, c in enumerate(line):
            if q:
                if c == q:
                    q = None
            elif c in "\"'":
                q = c
            elif c == "#":
                cut = i
                break
        out.append(line if cut is None else line[:cut])
    return "\n".join(out)


def _script_sources():
    for f in sorted((REPO / "scripts").glob("[0-9][0-9]_*.py")):
        # **必须排除 main_analysis.py**：本文件的验收层会引用图名，
        # 那是"检查对象"不是"出图声明"。不排除的话验收会要求自己
        # 引用过的每张图都存在，把口径搞反。
        if f.name == "main_analysis.py":
            continue
        yield f, _strip_comments(f.read_text(encoding="utf-8"))


_FIG_NAME_RE = re.compile(
    rf"^{PART}-\d{{2}}-\d{{2}}-unit\d+-[a-z0-9]+(?:-[a-z0-9]+)*$")


def declared_figures() -> list:
    """扫出**声明要出**的静态图名（字符串字面量）。M11。

    声明从源码扫出来而不是手抄一张表 —— 手抄的表会漂移，而漂移的方向
    恰好是"新加的图不在表里"，也就是把 E-48 那个盲区原样再造一遍。

    含 `{ }` 的模板串跳过（运行时拼名，由 `DYNAMIC_FIG_BASES` 声明豁免）。
    """
    out = []
    pat = re.compile(rf'"{PART}-[^"]*"')
    for _f, src in _script_sources():
        for m in pat.finditer(src):
            nm = re.sub(r"\.(pdf|png)$", "", m.group(0)[1:-1])
            if "{" in nm or "}" in nm:
                continue
            if _FIG_NAME_RE.match(nm):
                out.append(nm)
    return sorted(set(out))


def dynamic_fig_bases() -> dict:
    """扫出 `DYNAMIC_FIG_BASES = {"<图号>": <张数>}` 声明。

    键补全成完整前缀 `02-<模块号>-<图号>`。这是**槽位上限**，不是精确值：
    少出合法（只有 top 8 个 marker 基因时 8 个槽位里出不满），
    但**一张都没有说明那段循环整段没跑**。
    """
    out = {}
    pat = re.compile(r"DYNAMIC_FIG_BASES\s*=\s*\{([^}]*)\}")
    for f, src in _script_sources():
        m = pat.search(src)
        if not m:
            continue
        for pair in m.group(1).split(","):
            if ":" not in pair:
                continue
            k, v = pair.split(":", 1)
            key = k.strip().strip("\"'")
            try:
                out[f"{PART}-{f.name[:2]}-{key}"] = int(v.strip())
            except ValueError:
                continue
    return out


# 条件产出的图：**图名 -> (状态文件, 判据路径, 能力就位时的取值, 说明)**。
#
# **M6（R-03 裁决）**：`05_trajectory.py` 的 unit3 在 `celltype` 列缺失时
# 回退成按 `leiden` 着色，然后**跳过保存**（图名账目不含 leiden 回退）。
# 旧实现把这件事只写进 `log_warn` —— 日志是过程性的，事后没人读得到。
# 这里把它变成可判定的：`celltype` 列存在时这张图**变成必需**，
# 不存在时允许缺失但**必须可见**（报告里列出豁免原因）。
#
# 语义与 spatial 的 `CONDITIONAL_FIGURES` 一致：**不是"已知缺陷白名单"**，
# 它写的是"为什么可以没有"，且是**自愈**的 —— `celltype` 列哪天有了，
# 这张图立刻自动变成必需。
CONDITIONAL_FIGURES = {
    "02-05-04-unit3-celltype-on-umap": (
        "cluster_status.json", ("annotation", "status"), "ok",
        "03 步骤的细胞类型注释没成功时没有 `celltype` 列，"
        "此时 unit3 会退回 leiden 且按约定不落盘"),
    "02-02-03-unit1-batch-mixing": (
        "integration_status.json", ("has_batch_key",), True,
        "单样本数据没有批次，没有「前后对比」可画 —— 这张图本就不该存在"),
}


def _dig(obj, path):
    """按路径取值；任一层缺失返回 `None`（**不抛异常**）。"""
    cur = obj
    for k in path:
        if isinstance(cur, dict):
            cur = cur.get(k)
        elif isinstance(cur, list) and isinstance(k, int) and 0 <= k < len(cur):
            cur = cur[k]
        else:
            return None
    return cur


def chk(cid: str, kind: str, ok: bool, detail: str,
        severity: str = "required") -> dict:
    """构造一条验收项（Q-27 范式，与 spatial 侧同签名）。

    **第二个位置参数是 `kind` 不是 `severity`** —— 两者同名同型、位置相邻，
    传错不会报错，只会把 severity 值写进 kind 槽而检查仍然"看起来正常"
    （E-53 自查抓到的正是这个）。所以这里两个参数都写成关键字更安全的
    形式：`kind` 在前、`severity` 有默认值且在末位。

    `severity` 三档：
      - `"required"`：失败即验收红（默认）；
      - `"content"`：内容正确性，失败即红但计数分开；
      - `"info"`：只可见，永不判红（如"这张图缺中文说明"）。
    """
    return {"id": cid, "kind": kind, "item": f"[{kind}] {cid}",
            "ok": bool(ok), "required": severity != "info",
            "severity": severity, "detail": detail}


def load_step_fn(module_file: str, fn_name: str):
    spec = importlib.util.spec_from_file_location(
        module_file.replace(".py", ""), REPO / "scripts" / module_file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return getattr(mod, fn_name)


def has_file(p: Path) -> bool:
    return p.exists() and p.stat().st_size > 0


# 嵌套 status 里，哪些取值算"崩了"。其余（`not_applied` / `not_done` /
# `not_configured` / `skipped` …）都是**设计如此地没做**，只可见、不阻断。
#
# **注意这个集合只管嵌套字段，不要和 `common.STEP_ABORT_VALUES` 合并**：
# 顶层返回值里的 `bad_root` / `insufficient_methods` 是"这一步没做成"，
# 而嵌套字段里同名或近义的值往往只是"某个可选能力没成"——内置引擎照样出了
# 结果（`08_virtual_perturbation.py` 的 `tenifold.status` 可以是 `timeout`
# 或 `no_candidates`，顶层仍是 `ok`）。判红会让每个 job 都红，
# 反而没人看（台账元规则 ④：假阳性比假阴性更危险）。
NESTED_FAILED_VALUES = ("failed", "error", "fail")


def _iter_nested_status(obj, path=()):
    """递归产出所有名为 `status` 的字段：`(路径, 值, 同级 reason)`。

    Q-26 / E-48：`grn_status.json` 的**顶层** `status` 是 `"ok"`，而里面
    `regulon_vs_pseudotime.status` 是 `"failed"` —— 五张图从未产出，而
    `acceptance.json` 70 项全绿。步骤验收的判据是 `ok = (status == "ok")`，
    只看顶层，所以嵌套字段**没有任何一条检查看得见**。

    结构上这和"`_content_overflow()` 打了 WARN 没人读"是同一个错误：
    `except` 把异常降级成了一个**没人读的字段**。所以这里把它读出来。
    """
    if isinstance(obj, dict):
        if "status" in obj:
            reason = obj.get("reason") or obj.get("message") or ""
            yield path, obj["status"], reason
        for k, v in obj.items():
            if isinstance(v, (dict, list)):
                yield from _iter_nested_status(v, path + (str(k),))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            if isinstance(v, (dict, list)):
                yield from _iter_nested_status(v, path + (str(i),))


def _vp_engine_consistent(vp: dict) -> bool:
    """状态里的引擎清单与方法串必须自洽。

    原来这条查的是"方法串里必须出现 'scTenifoldKnk' 和 'PerturbNet' 两个词"，
    用来证明**没跑**这两个工具。K-01b 把真 scTenifoldKnk 接进来之后，那条
    判据会**把正确的结果判红** —— 方法串里现在必须出现 scTenifoldKnk 来
    *说明跑了它*，而不是来否认它。所以判据换成自洽性：

      * 记了 `engines_used` 就必须逐个在方法串里被点名；
      * 没跑 tenifold 时，方法串必须写明**为什么没跑**（而不是假装跑了）；
      * PerturbNet 始终没跑，方法串必须仍然否掉它。
    """
    if not vp:
        return False
    used = vp.get("engines_used") or []
    method = str(vp.get("method", ""))
    names = {"first_order": "一阶", "tenifold": "scTenifoldKnk"}
    if not all(names.get(e, e) in method for e in used):
        return False
    if "tenifold" not in used:
        td = vp.get("tenifold") or {}
        # 没跑 tenifold：要么配置里就没选它，要么它失败并留了原因
        if vp.get("engine") in ("tenifold", "both") and not td.get("reason"):
            return False
    return "PerturbNet" in method


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
    # 溢出记录同样要**跑前清空**（与状态文件同理）：否则这一轮修好了，
    # 上一轮留下的记录还在，验收会把已经修好的图判红。
    _ovf_stale = res_dir / "figure_overflow.json"
    if _ovf_stale.exists():
        _ovf_stale.unlink()
    for sid, mfile, fn, required, label in STEPS:
        if only and sid not in only:
            log_info(f"--- 跳过 {label} ({sid})：不在 --steps 里")
            continue
        log_info("")
        log_info(f"--- {label} ({sid})" + ("" if required else "  [可选]"))
        # 先删本步的状态文件：崩溃时不留旧文件冒充本轮结果。
        # M10：**显式取键再判空**，不用 `.get(sid, "")` —— 空串会让
        # `res_dir / ""` 等于结果目录本身。上面导入时的断言已保证键存在，
        # 这里仍然写成显式判空，两层各管一件事（断言管"配置一致"，
        # 判空管"这一行不会删到目录"）。
        fname = STEP_STATUS_FILES.get(sid)
        if fname:
            stale = res_dir / fname
            if stale.exists() and stale.is_file():
                stale.unlink()
        t0 = time.time()
        try:
            fn = load_step_fn(mfile, fn)
            # **返回值必须接住**（E-56）。步骤函数有一条"跑完了、但结果是
            # 『没做成』"的返回路径（`bad_root` / `insufficient_methods` /
            # `no_candidates` / `no_tfs_in_data` / `no_usable_pseudobulk` …），
            # 它们是 `write_json` + `log_warn` + `return status`，**不抛异常**。
            # 旧代码写死 `record_step(..., "ok", ...)` 并丢弃返回值 ——
            # 于是 `state.json` 记 `ok`、验收记绿，而这一步什么也没产出。
            # 这是 E-48 事故**未修完的另一半**。
            res = fn(cfg)
            record_step(cfg, sid, "ok", time.time() - t0, required=required,
                        result_status=result_status_of(res))
        except Exception as e:  # noqa: BLE001
            msg = f"{type(e).__name__}: {e}"
            record_step(cfg, sid, "failed", time.time() - t0,
                        message=msg, required=required,
                        result_status="failed")
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

        # ---- 步骤函数自己返回的 status 也要看（E-56）-----------------------
        #
        # `run_all` 现在把 `fn(cfg)` 的返回值记进了 `result_status`。**光记不读
        # 等于没记** —— 这正是 E-48 的形态（`except` 把异常降级成一个没人读的
        # 字段）。`bad_root` / `insufficient_methods` / `no_candidates` /
        # `no_tfs_in_data` / `no_usable_pseudobulk` 都是"步骤跑完了、但结果是
        # 『没做成』"：不抛异常、`state.json` 记 `ok`、验收记绿，而这一步
        # 什么也没产出 —— 最容易读成成功的一种。
        rs = st.get("result_status")
        if rs is None or str(rs).lower() == "ok":
            continue
        verdict = classify_step_result(rs)
        abort = (verdict == "abort")
        checks.append({
            "item": f"步骤 {label} 自述状态 = {rs}",
            "ok": not abort,
            "required": abort,
            "detail": (f"**步骤跑完了但结果是失败**：{sid} 返回 status={rs!r}"
                       if abort else
                       f"设计如此地没做（{sid} 返回 status={rs!r}），可见不阻断"),
        })

    for fname, desc, required in REQUIRED_FILES:
        ok = has_file(res_dir / fname)
        checks.append({"item": f"产物 {desc} ({fname})", "ok": ok,
                       "required": required,
                       "detail": "存在" if ok else "**缺失或为空**"})

    # 图属于哪一步，由图名里的模块号决定（`02-07-01-…` 的 `07` → `07_grn.py`
    # → `grn`）。**图不能无条件要求产出** —— 可选步骤没配时它本就不该有图
    # （如 `trajectory.enabled: false` 时 `02-05-*` 一张都不该有）。
    # 上一版把 `02-05-04/02-05-05` 无条件写成 required=True，等于一旦关掉
    # 轨迹，验收必红 —— 那是判据错了，不是产物错了。
    #
    # **M11（R-03 裁决）**：判据从"手抄表里的图在不在"换成
    # "**源码声明过的图在不在**"。手抄表会漂移，实测漏了 `02-05-01` /
    # `02-05-04-unit1..3` / `02-05-05-unit1` / `02-06-01` —— 其中
    # `02-06-01-unit1-communication-heatmap` 在 `06_communication.py` 里
    # 正常写出，却**没有任何检查看得见**。漂移方向恰好是"新加的图不在表里"，
    # 就是把 E-48 那个盲区原样再造一遍。
    _mod2sid = {mfile[:2]: sid for sid, mfile, _f, _r, _l in STEPS}
    _step_status = {s["id"]: s.get("status", "not_run")
                    for s in read_state(cfg).get("steps", [])}
    _step_required = {sid: req for sid, _m, _f, req, _l in STEPS}

    fig_names = sorted(p.stem for p in fig_dir.glob("*.png"))
    declared = declared_figures()
    dyn = dynamic_fig_bases()
    declared_set = set(declared)

    def _owner_of(fname: str):
        mod = fname.split("-")[1] if fname.count("-") >= 1 else ""
        owner = _mod2sid.get(mod)
        return owner, (_step_status.get(owner, "not_run") if owner else "not_run")

    # ---- 静态图：声明过的每一张都要在（可选步骤没跑则豁免但可见）------------
    missing, waived, skipped_optional = [], [], []
    for fname in declared:
        owner, owner_status = _owner_of(fname)
        if has_file(fig_dir / f"{fname}.png"):
            continue
        cond = CONDITIONAL_FIGURES.get(fname)
        if cond:
            st_file, path, ready_val, why = cond
            p = res_dir / st_file
            got = None
            if p.exists():
                try:
                    got = _dig(json.loads(p.read_text(encoding="utf-8")), path)
                except Exception:  # noqa: BLE001
                    got = None
            if got != ready_val:
                waived.append(f"{fname}（{why}；{st_file} "
                              f"{'.'.join(str(x) for x in path)}={got!r}）")
                continue
        if not _step_required.get(owner, False) and owner_status != "ok":
            # 可选步骤没跑（或配置关闭）→ 不要求这张图，但**必须可见**
            skipped_optional.append(f"{fname}（步骤 {owner} = {owner_status}）")
            continue
        missing.append(fname)

    _notes = []
    if waived:
        _notes.append(f"{len(waived)} 张条件图本轮不适用：{waived}")
    if skipped_optional:
        _notes.append(f"{len(skipped_optional)} 张属于未运行的可选步骤："
                      f"{skipped_optional}")
    _tail = ("；".join(_notes)) if _notes else ""
    checks.append(chk("figures:declared", "required", not missing,
                      (f"源码声明 {len(declared)} 张静态图，全部产出"
                       + (f"；{_tail}" if _tail else "")
                       if not missing else
                       f"**声明了但没产出** {missing}"
                       + (f"（{_tail}）" if _tail else "")
                       + f" —— 实际产出 {len(fig_names)} 张")))

    # ---- 动态图名：每组前缀至少 1 张 ----------------------------------------
    dyn_missing = []
    for base, n_slots in dyn.items():
        got = [nm for nm in fig_names
               if nm.startswith(base + "-") and nm not in declared_set]
        if not got:
            dyn_missing.append(f"{base}（声明 {n_slots} 个槽位，实际 0 张）")
    checks.append(chk("figures:dynamic", "required", not dyn_missing,
                      (f"{len(dyn)} 组动态图名共 {sum(dyn.values())} 个槽位，"
                       f"各自至少产出 1 张"
                       if not dyn_missing else
                       f"**动态图名整组没产出** {dyn_missing} —— "
                       f"槽位是上限不是精确值，少出合法，"
                       f"但**一张都没有说明那段循环整段没跑**")))

    # 保留计数作为**下限兜底**：声明扫描本身失效时（如源码结构大改导致
    # 一条字面量都扫不到）这条还能拦住"一张图都没有"。
    checks.append(chk("figures:count", "required", len(fig_names) >= 8,
                      f"{len(fig_names)} 张图（要求 >=8）"))

    # 扫出来的图名必须有说明 —— 否则报告里会出现一个只有文件名、没人知道
    # 它想表达什么的条目。**这条不判红**（`severity="info"`）：
    # `REQUIRED_FIGURES` 只是说明文字表，缺一条说明不代表产物有问题，
    # 判红会让"加了新图"这件好事变成一次 CI 失败。
    _no_desc = [nm for nm in declared if nm not in _FIG_DESC]
    checks.append(chk("figures:documented", "required", True,
                      (f"{len(declared)} 张声明图都有中文说明"
                       if not _no_desc else
                       f"{len(_no_desc)} 张声明图缺中文说明（不影响正确性，"
                       f"但报告里只有文件名）：{_no_desc}"),
                      severity="info"))

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
        # **旧判据是 `msum["n_versions"] >= 20`（L12），一个没有任何依据的
        # 魔数。** 它有两个毛病：
        #   1. 20 从哪来说不清 —— 实测本仓装 96 个包，20 只是"看起来够多"；
        #   2. 它把「枚举成功」和「枚举到多少个」混成一个量。枚举器半路抛
        #      异常时 `n_versions` 仍可能凑够 21 个，验收照样绿，而
        #      `versions` 里缺的正是后来要用来复现的那几个包。
        # 现在判「枚举这个动作成功没有」——那是一个事实，不是阈值。
        #
        # **三态，不是二态。** `versions_enumeration` 有三个取值：
        #   `"ok"`     —— 枚举成功；
        #   `"failed"` —— 枚举器抛了异常（这是本轮唯一该判红的）；
        #   `"unknown"`—— 清单里**没有这个字段**，即清单由旧版本代码写出。
        # 把 `unknown` 判红会让"读一份历史 artifact"变成失败，而那不是
        # 任何人的缺陷；把 `unknown` 判绿又会让"枚举从没跑过"混进通过里 ——
        # 与 E-64 同一条教训：**"没跑"和"跑了没问题"必须长得不一样**。
        # 所以 `unknown` 只可见（`required=False`），红只留给 `failed`。
        _enum = msum.get("versions_enumeration", "unknown")
        _enum_ok = (_enum == "ok")
        checks.append({
            "item": "版本枚举成功（importlib.metadata 未抛异常）",
            "ok": _enum != "failed", "required": _enum == "failed",
            "detail": (f"{msum.get('n_versions', 0)} 个已安装包"
                       if _enum_ok else
                       (f"**枚举失败**：{msum.get('versions_enumeration_error')}"
                        f" —— versions 只有 {msum.get('n_versions', 0)} 项，"
                        f"不足以复现本轮"
                        if _enum == "failed" else
                        f"**无法判断**：清单里没有 `versions_enumeration` "
                        f"字段（该清单由旧版本代码写出，不记枚举成败）；"
                        f"versions 有 {msum.get('n_versions', 0)} 项。"
                        f"本轮代码写出的清单会带这个字段。")),
        })
        # 关键工具解析率**可见但不阻断**：本仓 CI 只装 §2 的一个子集，
        # 实测 22 个关键工具里 8 个解析出来（含 scTenifoldKnk 1.1），
        # 其余是"查过了，没装"。把它判红会让每轮都红，但完全不报又会
        # 让"关键工具一个都没记上"从验收里消失 —— 所以要求**至少有一个
        # 解析出来**，并把未解析的名单完整列出。
        _nkr, _nkt = msum.get("n_key_resolved", 0), msum.get("n_key_total", 0)
        _unres = msum.get("key_unresolved") or []
        checks.append({
            "item": "关键工具版本有记录（§0.3）",
            "ok": _nkt > 0 and _nkr > 0, "required": True,
            "detail": (f"{_nkr}/{_nkt} 个关键工具解析出版本"
                       + (f"；未解析（查过了，未安装）: {', '.join(_unres)}"
                          if _unres else "")),
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
    #
    # **但"崩了"和"没做"必须分开** —— 见下面独立的内嵌 status 扫描（Q-26）。
    #
    # **M12（R-03 裁决）**：这里原来写死 `"ok": True` —— 无论 `st` 是什么
    # 都恒为真。后果：`ok=True` 的检查在报告里是 `[PASS]`，于是
    # "`trajectory` 崩了"和"`trajectory` 正常跑完"打印得一模一样。
    # 恒为真的量比没有这个量更糟 —— 它看起来像一条检查。
    #
    # 修法：`ok` 反映真实取值，`severity="info"` 保证**永不判红**
    # （可选步骤没做是设计如此，判红会让每个 job 都红，台账元规则 ④）。
    # 这样"崩了"在报告里显示成 `[INFO] ... ok=False`，与 `ok=True` 可区分，
    # 而不改变退出码。
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
        # `classify_step_result` 把自述状态归成 ok / abort / skip ——
        # 只有 `abort` 算"崩了"，其余（`not_applied` / `package_missing` …）
        # 是设计如此地没做。这里只借它算 `ok` 的真假，`severity="info"`
        # 保证退出码不受影响。
        _ok = classify_step_result(st) != "abort"
        checks.append(chk(f"step_status:{sid}", "info", _ok,
                          f"可选步骤状态 {sid} = {st!r}；{str(note)[:150]}",
                          severity="info"))

    # ---- 内嵌 status 扫描（Q-26 / E-48，本轮新增）----------------------------
    #
    # 上一版没有任何检查看得见**嵌套**的失败：`grn_status.json` 的顶层
    # `status` 是 `"ok"`，而里面 `regulon_vs_pseudotime.status` 是 `"failed"`
    # （`figsize` 三元素元组 → 5 张图从未产出），`acceptance.json` 70 项全绿。
    # 步骤验收的判据是 `ok = (status == "ok")`，只看顶层。
    #
    # 判据依据（真 artifact `scrna-results-55` 实测的嵌套取值分布）：
    #   {'ok': 5, 'failed': 1, 'not_applied': 1, 'not_done': 2}
    # 其中 `not_applied`（单样本无批次，不做整合）与 `not_done`（输入只有一套
    # 计数，不做 RNA 速率）是**设计如此地没做**，必须放行；只有 `failed` 是崩了。
    #
    # 结构上这和"`_content_overflow()` 打了 WARN 没人读"是同一个错误：
    # `except` 把异常降级成了一个**没人读的字段**。这里把它读出来。
    #
    # **E-56：这里以前写着 `if p`，把顶层路径滤掉了。** 顶层调用
    # `_iter_nested_status(d)` 时 `path` 是空元组 —— **falsy** —— 所以
    # `p and ...` 恒为假，**顶层 `status` 一个都扫不到**。而
    # `pseudobulk_status.json` / `trajectory_status.json` 这些文件的失败
    # 恰恰写在**顶层**（`05_trajectory.py` 的 `bad_root` /
    # `insufficient_methods` 整份文件就只有顶层一个 `status`）。姊妹仓库
    # spatial 的同构判据用 `key = ".".join(path) + ".status" if path else "status"`
    # 把顶层写成 `"status"` 而不是丢掉 —— 这里照抄那个写法。
    for stf in sorted(res_dir.glob("*status.json")):
        d = read_json(stf)
        if not isinstance(d, dict):
            continue
        bad = []
        for p, v, r in _iter_nested_status(d):
            if str(v).lower() not in NESTED_FAILED_VALUES:
                continue
            key = ".".join(p) + ".status" if p else "status"
            bad.append((f"{stf.name} → {key}", str(r)[:200]))
        for where, why in bad:
            checks.append({
                "item": f"{where} = failed",
                "ok": False, "required": True,
                "detail": f"**内嵌失败**：{why or '(无 reason)'}",
            })

    # ---- 内容超出画布（Q-26 / E-49，本轮新增）--------------------------------
    #
    # `save_fig` 里的 `_content_overflow()` **早就检测到了**
    # `{'width_overflow_frac': 0.3523}`、也**打了 WARN**，而 `02-08-01` 的标题
    # 两侧仍被静默裁掉 35%（左端只剩 `dKnk:`、右端断在 `the 60`）—— 因为
    # **告警没有消费者**。`savefig.bbox=standard` 下超出的部分直接被裁，
    # 文件照样生成：图名合规、图幅合规（1606 px 对应 136 mm 没超）、墨迹正常
    # （被裁的图照样有墨），四条现有门禁全绿。
    #
    # 所以判据强度从**警告级**（`log_warn`，会被忽略）升到**产物级**
    # （写进 `figure_overflow.json`，这里读它并判红）。
    ovf = read_json(res_dir / "figure_overflow.json") or {}
    ovf_figs = ovf.get("figures") or {}
    if ovf_figs:
        for fname_o, info in sorted(ovf_figs.items()):
            checks.append({
                "item": f"图 {fname_o} 内容超出画布",
                "ok": False, "required": True,
                "detail": f"**会被静默裁掉**：{info}",
            })
    else:
        checks.append({
            "item": "没有图的内容超出画布", "ok": True, "required": False,
            "detail": "figure_overflow.json 无记录（或无图被检测出溢出）",
        })

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
        # **02-05-03 已按单图原则拆成两张**（原双面板违反 D-006）：
        # unit1 = 热图、unit2 = 模块曲线。两张都要在，缺一张即判红。
        # 改名时**必须同步这里** —— 否则验收会查一个不存在的文件而静默变 false。
        for fn, desc in (("02-05-02-unit1-trajectory-method-correlation", "方法一致性矩阵"),
                         ("02-05-03-unit1-trajectory-modules-heatmap", "沿轨迹基因模块热图"),
                         ("02-05-03-unit2-trajectory-module-profiles", "各模块拟时序曲线")):
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
        # **"恒为 0 的一致率"必须显示成红的**（审计 S3 / 台账 E-58）。
        # 两个词表不相交时字符串全等恒为 0，而它看起来像一个结论
        # （"两条路完全不一致"），实际只是词表没映射。
        # 所以判据分三层：没对比 → 红；有对比但**一个簇都没映射上** → 红
        # （映射表缺了或过时了，正是旧缺陷的形态）；有映射 → 报真的一致率。
        _cmp_ok = bool(ct_cmp.get("compared"))
        _n_map = ct_cmp.get("n_mapped")
        _no_map = _cmp_ok and (_n_map == 0) and (ct_cmp.get("n_clusters") or 0) > 0
        checks.append({
            "item": "CellTypist 与 marker 注释的一致性已量化（§2.4）",
            "ok": _cmp_ok and not _no_map, "required": False,
            "detail": (
                f"簇层面一致 {ct_cmp.get('n_agree_mapped')}/{_n_map}"
                f"（{ct_cmp.get('agreement_frac_mapped')}）"
                f"；{ct_cmp.get('n_unmapped')} 个簇的词表无映射，不计入分母"
                if _cmp_ok and not _no_map else
                (f"**{ct_cmp.get('n_clusters')} 个簇一个都没映射上** —— "
                 f"`{ct_cmp.get('mapping_source')}` 缺条目或词表已变；"
                 "此时字符串全等恒为 0，**不能读成『两条路不一致』**"
                 if _no_map else
                 f"**未对比**：{ct_cmp.get('reason')}")),
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
    # 这一节的**重点是"没做什么"**。规范点名的三个工具里，PerturbNet 与
    # RegVelo 在本环境装不了，而"装不了"和"没装"是两件事 —— 前者有确切
    # 原因，后者是疏忽。所以逐个列出原因，并且把"用的是自建一阶近似"这件事
    # 写在最显眼处。**scTenifoldKnk 从 K-01b 起已经真的跑了**（R 引擎，
    # 见 `scripts/lib/tenifold_knk.R`），它的可用性由探针实测，不再是
    # 一条写死的"不可用"。
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
            "item": "虚拟扰动声明的方法与状态里的引擎一致",
            # 这一条是防"报了个数就被当成因果预测"的。原来它只查方法串里
            # 有没有"不是 scTenifoldKnk"两个词 —— 那个判据在真跑通 R 引擎
            # 之后会**把正确的结果判红**（K-01b 换掉了它）。现在的判据是
            # **自洽性**：状态里记了哪些引擎，方法串就必须点名哪些引擎，
            # 且必须仍然否掉没跑的 PerturbNet。
            "ok": _vp_engine_consistent(vp), "required": False,
            "detail": (f"engines_used={vp.get('engines_used')}；"
                       + str(vp.get("method", ""))[:100]),
        })
        _vp_td = vp.get("tenifold") or {}
        if _vp_td:
            _td_ok = _vp_td.get("status") == "ok"
            checks.append({
                "item": "scTenifoldKnk（R 引擎）的产出或失败原因已记录",
                # 失败也必须留下原因 —— 否则"没跑"和"跑了没结果"长得一样。
                "ok": bool(_vp_td.get("reason")) or _td_ok,
                "required": False,
                "detail": (f"网络 {_vp_td.get('n_genes_network')} 基因，"
                           f"{_vp_td.get('elapsed_sec')} s，"
                           f"版本 {(_vp_td.get('meta') or {}).get('engine_version')}"
                           if _td_ok else
                           f"{_vp_td.get('status')}：{str(_vp_td.get('reason'))[:80]}"),
            })
        # 空敲除（出度为 0 的候选基因）必须被点名。
        #
        # 包的敲除方式是"把网络里该基因那一行清零"；该基因出度为 0 时，
        # 清一行全 0 的行等于没敲，`KO` 与 `WT` 逐位相同，返回的"距离"
        # 只剩浮点噪声（实测 ~1e-16）。**它不报错、不给 NA**，读表的人
        # 只会得出"敲除这个基因没有影响"—— 而这正是错的。
        #
        # 所以判据是：只要 R 侧数出了空敲除，状态文件里就必须有它的
        # 数量和基因名，供下游（和读报告的人）区分"效应为 0"与
        # "这个网络表达不了该扰动"。
        if _vp_td.get("status") == "ok":
            _n_empty = int(_vp_td.get("n_empty_knockout") or 0)
            _empty_genes = _vp_td.get("empty_knockout_genes") or []
            checks.append({
                "item": "scTenifoldKnk 的空敲除（网络里出度为 0 的候选）已点名",
                "ok": (_n_empty == 0) or (len(_empty_genes) == _n_empty),
                "required": False,
                "detail": ("没有空敲除" if _n_empty == 0 else
                           f"{_n_empty} 个：{'、'.join(_empty_genes)}"
                           f"（**不是效应为 0，是该网络表达不了该扰动**）"),
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
    #
    # **变量名是 `mfull` 不是 `msum`（L13）。** 旧版这一段复用同一个 `msum`
    # 先装 `manifest_summary(cfg)`（计数摘要，键是 `n_versions` /
    # `inputs_missing` 这类）再装 `read_manifest(cfg)`（全量清单，键是
    # `versions` / `params` / `decisions`）。两个结构**键完全不同**，
    # 而 Python 对"读一个不存在的键"的默认行为是抛 `KeyError` ——
    # 或者更糟：`.get()` 静默给 `None`，检查项就变成一个恒真的空判据。
    # 同名的代价在下一次改动时才兑现：有人加一行 `msum["versions"]`
    # 会得到 `KeyError`，而加一行 `msum.get("versions")` 会**静默通过**。
    # 所以两个结构必须有两个名字。
    mfull = read_manifest(cfg)
    named = (mfull.get("params") or {}).get("named_tools")
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
        dec = mfull.get("decisions") or []
        checks.append({
            "item": "点名工具的使用情况进了决策链（§0.4）",
            "ok": any((d or {}).get("node") == "named_tools" for d in dec),
            "required": False,
            "detail": f"decisions 里 {len(dec)} 条，named_tools "
                      f"{'在' if any((d or {}).get('node') == 'named_tools' for d in dec) else '**不在**'}",
        })

    n_fail = 0
    for c in checks:
        mark = "PASS" if c["ok"] else ("FAIL" if c["required"] else "INFO")
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
        "n_info": sum(1 for c in checks if c.get("severity") == "info"),
        "checks": checks,
        "steps_failed_required": [{"id": i, "label": l, "error": m}
                                  for i, l, m in failed_required],
        "manifest": mfull,
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
