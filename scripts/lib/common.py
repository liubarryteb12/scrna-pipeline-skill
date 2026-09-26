"""
lib/common.py — 单细胞流水线公共库

约定与 geo-normal-pipeline-skill 保持一致（同一套心智模型）：
  - `parse_args()` **故意没有默认配置** —— 多数据集下静默默认到其中某一个，
    正是"跑错数据集"的来源。
  - 日志走 log_info/log_warn/log_error，不用裸 print。
  - 路径一律从 cfg 派生（dataset_id -> results/<id>、data/<id>），不硬编码。
  - 步骤失败必须抛异常让编排器记录，不要 try/except 后静默继续。
  - 可选步骤的失败要在自己的状态 JSON 里留 reason。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

# ============================================================================
# 确定性 —— **必须在 import numpy / scanpy / numba 之前执行**
# ============================================================================
# numba 在**首次 JIT 时**读这些环境变量决定线程数与线程层，晚设无效；
# OpenBLAS 在**运行期**按 CPU 型号选 SIMD 内核，也在 import 时定型。
# Part 1 的教训：浮点末位分叉有两个独立来源 —— 线程调度（多线程归约的
# 求和顺序）和 OpenBLAS 运行期按 CPU 型号选 SIMD 内核。这里都钉住。
#
# **M22（R-03 裁决）**：这一段原来写在 `import numpy as np` **之后**，
# 而注释自己写着"必须在 import scanpy/numba 之前设" —— **注释陈述的行为
# 与代码实际行为相反**（同 E-58 防复发②：S6 的"必须用 counts layer"、
# S9 的"用细胞数反推"、S10 的"判断整个 X"三处都是同一形态）。
# 本机直跑时这一层是**唯一**一层防护（CI 里 workflow 另有 job 级 env），
# 所以它失效时 CI 看不出来。同时补上原来**缺失**的 `NUMBA_NUM_THREADS`
# 与 `NUMBA_THREADING_LAYER`（AGENTS 规则 20.3 点名过；CI workflow 里
# 设了、common.py 里没设 ⇒ 本地直跑与 CI 不等价）。
#
# 顺序上必须放在所有第三方 import 之前 —— 下面 `import numpy as np` 等
# 一律排在这一段之后。`import os` 是标准库、不触发这些库的初始化，可以
# 留在原位（文件头）。
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
           "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("NUMBA_THREADING_LAYER", "workqueue")
os.environ.setdefault("OPENBLAS_CORETYPE", "Haswell")
os.environ.setdefault("PYTHONHASHSEED", "0")

import numpy as np  # noqa: E402  —— 必须排在上面那段 env 之后
import yaml  # noqa: E402


# ============================================================================
# 日志
# ============================================================================
_ORCHESTRATED = False


def set_orchestrated(flag: bool = True) -> None:
    """被 main_analysis.py source 时置 True，抑制重复的步骤横幅。"""
    global _ORCHESTRATED
    _ORCHESTRATED = flag


def orchestrated() -> bool:
    return _ORCHESTRATED


def _log(level: str, msg: str) -> None:
    ts = time.strftime("%H:%M:%S")
    print(f"[{ts}] {level:<7} {msg}", flush=True)


def log_info(msg: str) -> None:
    _log("INFO", msg)


def log_warn(msg: str) -> None:
    _log("WARN", msg)


def log_error(msg: str) -> None:
    _log("ERROR", msg)


# ============================================================================
# 参数与配置
# ============================================================================
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="单细胞转录组流水线",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    # **没有 default。** 见模块 docstring。
    p.add_argument("--config", required=True, help="配置文件路径（assets/config.<id>.yml）")
    p.add_argument("--steps", default=None,
                   help="只跑指定步骤，逗号分隔（调试用）。默认跑全部")
    return p.parse_args(argv)


def load_config(path: str) -> dict:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"配置文件不存在: {path}")
    with open(p, encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    if not isinstance(cfg, dict):
        raise ValueError(f"配置文件顶层必须是映射: {path}")
    if not cfg.get("dataset_id"):
        raise ValueError("配置缺少 dataset_id")

    did = str(cfg["dataset_id"])
    out = cfg.setdefault("output", {})
    # 目录由 dataset_id 派生 —— 一个数据集一个产物目录，互不覆盖。
    out.setdefault("results_dir", f"results/{did}")
    out.setdefault("data_dir", f"data/{did}")
    out.setdefault("figures_dir", f"results/{did}/figures")

    ana = cfg.setdefault("analysis", {})
    ana.setdefault("seed", 20260919)
    # **300 dpi 是投稿图的底线，不是"够用就行"。**
    # 原来默认 150：183 mm 宽的图在 150 dpi 下只有 1080 px，放大或印刷后
    # 字形和细线都发虚 —— 而"发虚"从图注上完全看不出来，文件大小也正常。
    #
    # **注意这只是兜底，改它并不够。**
    # `assets/config.pbmc3k.yml` 里显式写了 `figure_dpi`，而 `setdefault`
    # 只在键缺失时才生效 —— 所以那份配置根本用不到这个默认值。
    # **这是实测踩过的**：只改了这里，CI 跑绿，但从 artifact 里读 PNG 头
    # 发现宽度仍是 907 px（= 6.05 in × 150），也就是**根本没生效**。
    # 日志里看不出任何异常 —— 只有把 artifact 拿下来量像素才能发现。
    # 现在两处都是 300。
    # （姊妹项目 geo-normal-pipeline-skill 的 config 恰好没有这个键，
    # 所以那边只改 `save_pdf()` 的函数默认值就够了 —— 那是巧合，不是通例。）
    ana.setdefault("figure_dpi", 300)
    return cfg


def ensure_dirs(cfg: dict) -> None:
    for k in ("results_dir", "data_dir", "figures_dir"):
        Path(cfg["output"][k]).mkdir(parents=True, exist_ok=True)


def set_seed(cfg: dict) -> int:
    """把全局随机种子钉死。**紧挨着随机调用设置**，与上游消耗了多少随机数无关。"""
    seed = int(cfg["analysis"]["seed"])
    random.seed(seed)
    np.random.seed(seed)
    # **样式必须在这里应用，不能等到 save_fig。** rcParams 只在 figure
    # 创建时被读取 —— 等图建好了再 plt.style.use，那张图仍然是旧样式。
    # 每个步骤脚本第一句都是 set_seed(cfg)，所以这里是唯一的正确位置。
    apply_style(cfg)
    return seed


def _scrub_nonfinite(obj):
    """把 `NaN` / `Infinity` 递归换成 `None`，并返回 `(新对象, 命中数)`。

    **M9（R-03 裁决）**：`json.dump` 默认 `allow_nan=True`，会把 `NaN` /
    `Infinity` **原样写进文件** —— 那不是合法 JSON，`node`/`jq`/R 侧读者
    一律解析失败。而 Python 自己 `json.load` 读得回来，所以**本地看不出来**
    （M17 的 `score_margin: NaN` 就是这么落盘的）。

    换成 `None`（JSON `null`）而不是删键：**"这个量算不出来"本身是信息**，
    删掉它读者只会以为没这个字段。同时把命中数报出来，让它可见。
    """
    if isinstance(obj, dict):
        out, hits = {}, 0
        for k, v in obj.items():
            nv, h = _scrub_nonfinite(v)
            out[k], hits = nv, hits + h
        return out, hits
    if isinstance(obj, (list, tuple)):
        out, hits = [], 0
        for v in obj:
            nv, h = _scrub_nonfinite(v)
            out.append(nv)
            hits += h
        return (out if isinstance(obj, list) else tuple(out)), hits
    if isinstance(obj, float) and not math.isfinite(obj):
        return None, 1
    return obj, 0


def write_json(path, obj) -> None:
    """写 JSON。**保证写出的是合法 JSON**（M9）—— 非有限浮点先换成 `null`。

    用了 `allow_nan=False` 作为**兜底断言**：即使 `_scrub_nonfinite` 漏了
    某个非有限值（如藏在自定义对象里），`json.dump` 会抛 `ValueError`
    而不是**悄悄写出一个非法文件**。宁可这一轮失败，也不要下游读者解析失败。
    """
    Path(path).parent.mkdir(parents=True, exist_ok=True)

    def default(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            f = float(o)
            return None if not math.isfinite(f) else f
        if isinstance(o, (np.bool_,)):
            return bool(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, Path):
            return str(o)
        return str(o)

    clean, n_nonfinite = _scrub_nonfinite(obj)
    if n_nonfinite:
        log_warn(f"write_json: {path} 有 {n_nonfinite} 个 NaN/Infinity，"
                 f"已写成 null（非法 JSON 会让 node/jq/R 读者解析失败）")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(clean, fh, ensure_ascii=False, indent=2,
                  default=default, allow_nan=False)


def read_json(path):
    """读 JSON。文件不存在 → `None`；**内容坏了 → 抛异常**（M8）。

    原实现 `except Exception: return None` 把两种完全不同的情况混成一个
    `None`：① 这轮没这个文件（正常，可选步骤没跑）；② **文件在、但内容是
    坏的**（写坏了 / 截断了 / 上次运行留下的半截文件）。

    ② 的后果不是"少一条数据"，而是**静默丢掉一条检查** ——
    消费点 `main_analysis.py` 写的是 `d = read_json(...)` / `if not d: continue`，
    所以一个损坏的 `*_status.json` 会让那条验收项**直接消失**，报告里
    看不出少了什么。这比"报错"糟得多（E-62/E-63/E-64 同一形态：
    假阴性把自己藏了起来）。

    所以：**不存在 → None；存在但解析失败 → 抛**。需要"坏也当没有"的
    调用点显式用 `read_json_or_none()`，让意图在调用处可见。
    """
    p = Path(path)
    if not p.exists():
        return None
    with open(p, encoding="utf-8") as fh:
        return json.load(fh)


def read_json_or_none(path):
    """`read_json` 的宽容版：**任何**失败都返回 `None`，但会 `log_warn`。

    只给"这一步的可选产物、坏了大不了当没跑"的调用点用。判据产物
    （`*_status.json` 等）一律用 `read_json` —— 坏掉必须炸出来。
    """
    try:
        return read_json(path)
    except Exception as e:  # noqa: BLE001
        log_warn(f"read_json_or_none: {path} 存在但解析失败"
                 f"（{type(e).__name__}: {e}）—— 当作不存在处理")
        return None


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ============================================================================
# 步骤状态（state.json）
# ============================================================================
def state_path(cfg: dict) -> Path:
    return Path(cfg["output"]["results_dir"]) / "state.json"


def read_state(cfg: dict) -> dict:
    return read_json(state_path(cfg)) or {"steps": []}


# ---------------------------------------------------------------------------
# 步骤自述状态的翻译（E-56）
# ---------------------------------------------------------------------------
# 步骤函数有一条"跑完了、但结果是『没做成』"的返回路径：它们
# `write_json(状态文件, status)` + `log_warn(...)` + `return status`，
# **不抛异常**。旧编排器写死 `record_step(cfg, sid, "ok", ...)` 并丢弃返回值 ——
# 于是 `state.json` 记 `ok`、验收记绿，而这一步什么也没产出。
# 这是 E-48 事故**未修完的另一半**。
STEP_ABORT_VALUES = (
    # 通用崩溃
    "failed", "error", "fail",
    # 05 轨迹：配置的根簇不在簇列表里 / 上游没产出 leiden / 成功方法不足 2 种
    "bad_root", "missing_clusters", "insufficient_methods",
    # 08 虚拟敲除：候选基因没有一个满足最小靶基因数
    "no_candidates",
    # 07 GRN：TF 一个都不在数据里 / 没有 TF 找到足够多的正相关靶
    "no_tfs_in_data", "no_regulons",
    # 06 通讯：数据库里没有一对配体受体同时在数据里 / 打分全 0
    "no_pairs_in_data", "no_signal",
    # 04 拟bulk：聚合不出可用拟bulk / 没有一个细胞类型可检验 / 配置的列不存在
    #   `missing_counts`（审计 S6）：上游没保留 counts 层 —— 这是**契约被破坏**，
    #   不是"这批数据不适合做拟bulk"。旧实现静默退回 log 后的 `.X`，
    #   离散度估计全错而产物看起来正常。
    "no_usable_pseudobulk", "no_celltype_testable", "column_missing",
    "missing_counts",
    # 04 拟bulk（审计 M1）：样本与分组不是一对一 —— 旧实现对每个样本取
    # `iloc[0]` 当分组，混了分组的样本被静默归到第一个细胞的组，
    # DESeq2 照样给出完整差异表而分组标签是错的。**契约被破坏要判红。**
    "mixed_group",
)

# **设计如此地没做**（配置关掉了 / 输入不支持 / 环境缺包）—— 只可见、不阻断。
STEP_SKIP_VALUES = (
    "disabled", "not_configured", "not_applicable", "not_run",
    "not_applied", "not_done", "not_available", "not_possible",
    "heuristic_only", "package_missing", "unavailable", "needs_reference",
    # 审计 M19/M20：**部分成功**。子集被跳过（细胞类型重复不足 / 某个
    # 可选能力没成）但主体产出了结果。它必须与 `ok` **长得不一样** ——
    # 否则"7 个里成了 2 个"和"7 个全成"在状态文件里无法区分。
    # 登记为 skip 而不是 abort：跳过本身是设计如此，判红会让每个 job 都红。
    "partial",
)

# 上面两个集合**只管步骤函数顶层返回值**。嵌套字段里的同一批词不能照搬：
# `08_virtual_perturbation.py` 的 `tenifold.status` 可能是 `timeout` /
# `no_candidates`，而内置引擎照样出了结果（顶层 `ok`）；`annotation.celltypist`
# 的 `model_unavailable`、`liana` 的 `api_not_found` 也只是"这个可选能力没成"。
# 把它们判红会让每个 job 都红 —— 假阳性比假阴性更危险（台账元规则 ④）。


def result_status_of(res) -> str:
    """步骤函数**返回值**里自述的 `status`（E-56）。

    为什么需要它：步骤函数有一条"跑完了、但结果是『没做成』"的返回路径 ——
    `05_trajectory.py` 的 `bad_root` / `insufficient_methods`、
    `08_virtual_perturbation.py` 的 `no_candidates`、`07_grn.py` 的
    `no_tfs_in_data` 都是 `write_json(...)` + `log_warn(...)` + `return status`，
    **不抛异常**。编排器只要不接住这个返回值，`record_step` 就会记 `ok`，
    验收记绿 —— 而这一步什么也没产出。

    返回 `"ok"` 的两种情形：① 返回的 dict 里没有 `status` 键（如
    `00_fetch.py` 返回的 `info`）；② 返回的不是 dict。两者都表示
    "函数正常跑完、没有自述失败"，与调用方的 `except` 分支互补。
    """
    if isinstance(res, dict):
        st = res.get("status")
        return str(st) if st is not None else "ok"
    return "ok"


def classify_step_result(result_status: str) -> str:
    """把一个步骤自述状态归成三类：`"ok"` / `"abort"` / `"skip"`（E-56）。

    - `"ok"`    —— 成功（含 `None` / 空）。
    - `"abort"` —— **语义上是「这一步没做成」**，判红（`severity="required"`）。
    - `"skip"`  —— **设计如此地没做**（配置关了 / 环境缺包），可见不阻断。
    - 未登记的取值也归 `"skip"` —— 宁可只可见也不判红。判据的原则是
      **只把确知是失败的判红**：一个没见过的词有可能是新加的"设计如此"，
      判红会让 job 变红、让人开始忽略告警（台账元规则 ④ 假阳性更危险）。
    """
    v = str(result_status or "").strip().lower()
    if not v or v == "ok":
        return "ok"
    if v in STEP_ABORT_VALUES:
        return "abort"
    return "skip"


def record_step(cfg: dict, step_id: str, status: str, seconds: float = None,
                message: str = "", required: bool = True,
                result_status: str = None) -> None:
    st = read_state(cfg)
    st.setdefault("steps", [])
    st["steps"] = [s for s in st["steps"] if s.get("id") != step_id]
    entry = {"id": step_id, "status": status, "required": bool(required),
             "message": message}
    if result_status is not None:
        # 步骤函数自己返回的 status —— 与 `status`（编排器记的"有没有崩"）
        # 是**两件事**：`status="ok"` + `result_status="bad_root"` 完全可能，
        # 那正是 E-56 要让它显形的那种组合。
        entry["result_status"] = str(result_status)
    if seconds is not None:
        entry["seconds"] = round(float(seconds), 1)
    st["steps"].append(entry)
    st["dataset_id"] = cfg["dataset_id"]
    st["updated_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    write_json(state_path(cfg), st)


# ============================================================================
# 运行清单（模块零规范层）
# ============================================================================
# 规范来源：用户整合文档「模块零：语言与运行时规范」。
#   §0.2 跨语言接口 —— 只走 CSV；每次转换记录维度/metadata/丢失字段
#   §0.3 版本记录   —— pip freeze 全量 + 文档点名的关键工具单独记版本
#   §0.3 随机种子   —— 所有随机过程固定种子并记录
#   §0.4 运行日志   —— 输入数据哈希、软件版本、关键参数、决策链、
#                      人工干预记录、跨语言转换记录
#
# 产物：results/<dataset_id>/run_manifest.json
#
# **为什么不塞进 state.json：** state.json 记的是"这一步跑没跑成"，每步重写；
# manifest 记的是"本轮是在什么条件下跑出来的"，是证据，写入后不该再变。
# 混在一起会让后者被前者覆盖。
#
# **诚实性要求（AGENTS.md 规则 4）：** 没做的分析、没装的工具、没确认的
# 复核节点，都要在 manifest 里留下痕迹，不能因为"不影响结论"就不写。
# `decisions` 记的是**实际的选择**，不是"应该怎么做"。

MANIFEST_NAME = "run_manifest.json"

# 文档 §2 点名的工具。**没装的记 None，不省略键** —— 键消失和"版本是 None"
# 看起来完全不同，后者才说明"这个工具本该有但没装"。
#
# **但这张表里的 R 包不能靠 `importlib.metadata` 判"没装"。** 那句话的前提是
# "本仓库 CI 里没有 R"，K-01b 之后不成立了：CI 装了 R，`scTenifoldKnk` 真跑了
# 6 分钟，而 `key_versions` 仍记 None —— 于是 `run_manifest.json`（复现依据）
# 与 `virtual_perturbation_status.json` 对**同一个工具**给出相反结论。
# 所以 R 包从这张表里挪出去，单独走 `R_KEY_PACKAGES` + `probe_r_packages()`
# （起一次 `Rscript` 问，而不是问 Python 的包数据库）。
KEY_PACKAGES = [
    # §2 核心 Python 包
    "scanpy", "anndata", "scvi-tools", "cellbender", "harmonypy", "scvelo",
    "celltypist", "pyscenic", "liana", "doubletdetection", "scrublet",
    # §2.6 拟时序
    "palantir", "scfates", "cytotrace",
    # §1.7/§1.8 虚拟扰动。PerturbNet / RegVelo 是 **Python** 包（PyPI 上的版本
    # 都钉死了 Python 上限），所以走上面那条 `importlib.metadata` 通道；
    # 记 None 是"查过了，装不上"，省略键是"没查"。
    # 确切原因写在 virtual_perturbation_status.json 的 tools 字段里。
    "PerturbNet", "RegVelo",
]

# R 包（CRAN / Bioconductor）。**版本只能问 R。** `importlib.metadata` 枚举的是
# Python 发行版，对 R 包原理上永远返回 None —— 那个 None 会被读成
# "查过了，装不上"，而它其实只是"问错了地方"。
#
# `scTenifoldKnk` 是 §1.8 点名的虚拟敲除工具，K-01b 起 CI 真的装了它并每轮跑
# （见 workflow 的 setup-r 步骤与 `scripts/lib/tenifold_knk.R`）。
R_KEY_PACKAGES = [
    "monocle3", "slingshot", "cellchat", "soupx", "scdblfinder",
    "scTenifoldKnk",
]


# ---- 文档点名、但本仓库用不了的工具（§2.1–§2.8）-----------------------------
#
# 与空间侧（`spatial-pipeline-skill/scripts/lib/common.py`）**同一套结构**，
# 四类 kind 同名同义。判据全部是实测的 PyPI metadata，不是推测。
#
# ## 为什么单细胞侧更需要这张表
#
# §2 点名的工具里 R 包特别多（SCTransform / scran / DESeq2 / edgeR /
# Monocle3 / Slingshot / CellChat / SoupX），而本仓库**没有 rpy2 路径** ——
# 规范说"Python (+R via rpy2)"，但 CI 里没装 R，所以这些一个都跑不了。
#
# 产物里却看不出来：`02_integrate.py` 有归一化、`04_pseudobulk_de.py` 有
# 差异表、`05_trajectory.py` 有四条轨迹 —— **每一步都有东西，
# 所以"点名的方法一个都没用上"这件事必须自己说出来。**
#
# ## `name_taken` —— 单细胞侧也有，而且更荒谬
#
#   pip install edgeR     -> "Redirect Microsoft Edge to your preferred browser"
#   pip install slingshot -> ElasticSearch 索引迁移
#
# `edgeR` 是 §2.5 点名的差异分析工具，`slingshot` 是 §2.6 点名的轨迹工具。
# 这两个名字**绝对不能进 requirements.txt**：装上了不会报错，
# 只会让 `import edgeR` 拿到一个浏览器重定向库。
#
# 判断依据是 summary / author / project_urls，不是"名字存不存在" ——
# `SingleR` 就是反例（见下）。
NAMED_TOOLS = {
    # ---- §2.1 原始矩阵校正 --------------------------------------------------
    "SoupX": dict(
        kind="r_package", section="§2.1",
        reason=("R 包（CRAN/GitHub），PyPI 上无同名包；本仓库 CI 不装 R + rpy2。"
                "它还需要空液滴（empty droplet）信息"),
    ),
    "CellBender": dict(
        kind="deps", section="§2.1",
        reason=("PyPI 有真包（0.4.0），但依赖 `torch` + `pyro-ppl>=1.8.4` —— "
                "GPU 导向的深度生成模型；GitHub 托管 runner 无 GPU，"
                "CPU 训练时间不现实"),
    ),
    # ---- §2.2 归一化 --------------------------------------------------------
    "SCTransform": dict(
        kind="r_package", section="§2.2",
        reason="R 包（Seurat 生态），PyPI 上无同名包；本仓库用 normalize_total + log1p",
    ),
    "scran": dict(
        kind="r_package", section="§2.2",
        reason="Bioconductor R 包，PyPI 上无同名包",
    ),
    # ---- §2.4 细胞类型注释 --------------------------------------------------
    "SingleR": dict(
        kind="needs_reference", section="§2.4",
        reason=("**这个不是顶名的**：PyPI 上的 `SingleR` 0.5.0 是 BiocPy/singler，"
                "R 那个算法的官方 Python 绑定（作者 Aaron Lun）。"
                "但它需要**带标签的参考数据集**（R 侧的 celldex / ImmGen / HPCA），"
                "本流水线没有 —— 而 CellTypist 自带可下载的预训练模型。"
                "所以第二条证据用 CellTypist，不是 SingleR"),
    ),
    # ---- §2.5 拟bulk 差异分析 -----------------------------------------------
    "DESeq2": dict(
        kind="r_package", section="§2.5",
        reason="Bioconductor R 包，PyPI 上无同名包；规范给的 rpy2 路径本仓库没有",
    ),
    "edgeR": dict(
        kind="name_taken", section="§2.5",
        reason=("**PyPI 上的 `edgeR` 是「Redirect Microsoft Edge to your preferred "
                "browser」（作者 Daniel Agans，github.com/phwelo/edger）—— "
                "和 Bioconductor 的 edgeR 毫无关系。** 真 edgeR 是 R 包。"
                "这个名字进了 requirements.txt 就会装进来一个浏览器重定向工具"),
    ),
    # ---- §2.6 轨迹 ----------------------------------------------------------
    "Monocle3": dict(
        kind="r_package", section="§2.6",
        reason="R 包（GitHub cole-trapnell-lab/monocle3），PyPI 上无同名包",
    ),
    "Slingshot": dict(
        kind="name_taken", section="§2.6",
        reason=("**PyPI 上的 `slingshot` 是「Index Migration for ElasticSearch」"
                "（作者 Thierry Jossermoz）—— 与 Bioconductor 的 Slingshot 无关。** "
                "真 Slingshot 是 R 包；本仓库用它的 Python 移植 scFates 代替"),
    ),
    "CytoTRACE2": dict(
        kind="not_on_pypi", section="§2.6",
        reason=("PyPI 上 `CytoTRACE` / `cytotrace` / `CytoTRACE2` 都查不到。"
                "本仓库**自行实现**了 CytoTRACE 的核心统计量（GCS，"
                "基因计数特征），产物里写明「这是自行实现，不是 CytoTRACE 包」"),
    ),
    "scVelo": dict(
        kind="needs_layers", section="§2.6",
        reason=("PyPI 有真包（0.3.4），依赖也不重 —— 但它需要 RNA velocity 的 "
                "**spliced / unspliced 层**，而 10x 的 filtered_feature_bc_matrix "
                "只有 counts，没有这两层。硬跑会把未剪接信息当 0"),
    ),
    # ---- §2.7 细胞通讯 ------------------------------------------------------
    "CellChat": dict(
        kind="r_package", section="§2.7",
        reason="R 包（GitHub JinmiaoChenLab/CellChat），PyPI 上无同名包；主用 LIANA",
    ),
    # ---- §2.8 基因调控网络 --------------------------------------------------
    "pySCENIC": dict(
        kind="needs_resources", section="§2.8",
        reason=("PyPI 有真包（pyscenic 0.12.1，依赖是 numba/dask/ctxcore 这些，"
                "装得上）。**卡住的不是 pip，是资源**：SCENIC 需要 cisTarget 的 "
                "motif 注释数据库（feather，GB 级）+ TF 列表，CI 上下载不现实。"
                "所以本仓库跑的是共表达推断，产物里显式写明「不是 SCENIC」"),
    ),
}


def probe_named_tools(log=None, only=None) -> dict:
    """把 NAMED_TOOLS 整理成可写进状态 JSON 的登记表。

    与空间侧同名同义。**不尝试 import** —— 这一节的结论是"没装/装不了"。
    用 `importlib.util.find_spec` 复核一次，发现"登记说过不了、环境里却能
    import"时 WARN（说明登记过期了）。
    """
    import importlib.util

    import_name = {
        "CellBender": "cellbender", "scran": "scran", "SingleR": "singler",
        "DESeq2": "DESeq2", "edgeR": "edgeR", "scVelo": "scvelo",
        "pySCENIC": "pyscenic", "SoupX": "SoupX", "SCTransform": "sctransform",
        "Monocle3": "monocle3", "Slingshot": "slingshot",
        "CytoTRACE2": "cytotrace", "CellChat": "cellchat",
    }
    out = {}
    for tool, meta in NAMED_TOOLS.items():
        if only is not None and tool not in only:
            continue
        mod = import_name.get(tool)
        avail = False
        if mod:
            try:
                avail = importlib.util.find_spec(mod) is not None
            except (ImportError, ValueError):
                avail = False
        if avail and log:
            log_warn(f"工具 {tool} 登记为不可用，但环境里能 import —— 登记需要更新")
        out[tool] = {
            "available": bool(avail),
            "kind": meta["kind"],
            "section": meta["section"],
            "reason": meta["reason"],
        }
    return out


def named_tools_note() -> str:
    """一句话说明本仓库为什么 §2 点名的方法多数没用上。"""
    return ("文档 §2.1–§2.8 点名的工具里，R 包（SoupX / SCTransform / scran / "
            "DESeq2 / edgeR / Monocle3 / Slingshot / CellChat）一个都用不了 —— "
            "CI 从 K-01b 起装了 R，但只为 §1.8 的 `scTenifoldKnk` 一个包装的，"
            "**没有 rpy2 通道**，所以这些 R 包仍然调不到；"
            "CellBender / scVelo / pySCENIC 有 PyPI 真包，"
            "但分别卡在 GPU、缺 spliced/unspliced 层、需要 GB 级 motif 数据库；"
            "**`edgeR` 与 `slingshot` 在 PyPI 上是同名无关包**（浏览器重定向 / "
            "ElasticSearch 迁移）。逐条理由见各状态文件的 named_tools 字段。"
            "**这是缺口，不是「已覆盖」。**")



def manifest_path(cfg: dict) -> Path:
    return Path(cfg["output"]["results_dir"]) / MANIFEST_NAME


def read_manifest(cfg: dict) -> dict:
    return read_json(manifest_path(cfg)) or {}


def init_manifest(cfg: dict, language: str = "python") -> dict:
    """建立本轮清单骨架。**会清掉上一轮的内容** —— 清单描述的是本轮。"""
    m = {
        "dataset_id": cfg.get("dataset_id"),
        "language": language,
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "seed": (cfg.get("analysis") or {}).get("seed"),
        "versions": {},
        "key_versions": {},
        "inputs": [],
        "params": {},
        "decisions": [],
        "human_review": [],
        "cross_language": [],
    }
    write_json(manifest_path(cfg), m)
    return m


def _manifest_append(cfg: dict, key: str, entry) -> None:
    m = read_manifest(cfg)
    m.setdefault(key, [])
    m[key].append(entry)
    m["dataset_id"] = cfg.get("dataset_id")
    m["updated_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    write_json(manifest_path(cfg), m)


def _norm_pkg(name: str) -> str:
    """PEP 503 归一化：包名大小写与 -/_/. 不敏感，且**连续分隔符折叠成一个**。

    **L19（R-03 裁决）**：原实现是
        `lower().replace("_", "-").replace(".", "-")`
    —— 只做了"换字符"，没做"折叠"。PEP 503 的规范式是
    `re.sub(r"[-_.]+", "-", name).lower()`，即 `foo__bar` 与 `foo-bar`
    是**同一个发行版**，而原实现给出 `foo--bar` 与 `foo-bar` 两个不同的串。
    后果是"同一个包被数成两个"：`key_versions` 里出现两条键，
    去重统计偏大，而**没有任何地方会因此报错** —— 只会让
    `n_versions >= 20` 这条验收在一个依赖很少的环境里**假通过**
    （少一个真包、多一条重复键，净效果抵消）。
    """
    return re.sub(r"[-_.]+", "-", str(name).strip().lower())


# R 包名允许的字符（CRAN 规范：字母开头，只含字母数字点）。
# 用来挡住把任意字符串拼进 `Rscript -e` 表达式 —— 那不是解析，是拼串。
_R_PKG_OK = re.compile(r"^[A-Za-z][A-Za-z0-9.]*$")


def probe_r_packages(pkgs) -> dict:
    """问 **R 自己** 这几个包装了没有。一次 `Rscript` 问完全部。

    **为什么不能用 `importlib.metadata` 或 `find_spec`。** 那两条查的都是
    Python 的包数据库 / 模块查找器，对 R 包原理上永远返回"没有"。K-01b 之前
    本仓库 CI 里确实没有 R，所以那个 `None` 恰好是对的；K-01b 起 CI 装了 R、
    `scTenifoldKnk` 真跑了 6 分钟，而 `key_versions` 仍记 `None` ——
    于是 `run_manifest.json`（**复现依据**）与
    `virtual_perturbation_status.json` 对同一个工具给出相反结论。
    错的不是那个 `None` 的值，是**问错了地方**。

    返回 `{"rscript", "r_version", "packages": {pkg: version|None}, "reason"}`。
    Rscript 不在 PATH 时 `reason` 写明"CI 没有装 R"，各包版本记 `None` ——
    这个 `None` 与"问了 R，R 说没装"含义不同，所以两个原因分开写。

    这是本仓库**唯一**一处"问 R 包版本"的实现：`08_virtual_perturbation.py`
    的 `_probe_r_package()` 转调这里，免得同一件事有两份代码（抄两份时
    验证的往往只是副本 —— 见 AGENTS 规则 16 的教训）。
    """
    names = [str(p) for p in pkgs]
    bad = [p for p in names if not _R_PKG_OK.match(p)]
    if bad:
        return {"rscript": None, "r_version": None, "packages": {p: None for p in names},
                "reason": f"包名不是合法的 R 标识符，拒绝拼进表达式: {bad}"}

    rscript = shutil.which("Rscript")
    if rscript is None:
        return {"rscript": None, "r_version": None,
                "packages": {p: None for p in names},
                "reason": ("Rscript 不在 PATH 上 —— 本机/本 CI 没有装 R"
                           "（见 workflow 的 setup-r 步骤）")}

    # 一次问完：每个包一行 `OK<TAB>名字<TAB>版本` 或 `MISSING<TAB>名字<TAB>`。
    quoted = ", ".join('"' + p + '"' for p in names)
    expr = (
        f"pk <- c({quoted})\n"
        "for (p in pk) {\n"
        "  if (requireNamespace(p, quietly = TRUE)) {\n"
        '    cat("OK\\t", p, "\\t", as.character(utils::packageVersion(p)), "\\n", sep = "")\n'
        "  } else {\n"
        '    cat("MISSING\\t", p, "\\t\\n", sep = "")\n'
        "  }\n"
        "}\n"
        'cat("RVERSION\\t", as.character(getRversion()), "\\n", sep = "")\n'
    )
    try:
        # **L10（R-03 裁决）**：这里用 `capture_output=True`（管道），而
        # `capture_versions` 的 docstring 说"管道捕获输出在受限沙箱里会
        # EPERM"。两句话曾经并列在同一份文件里、互相矛盾。
        #
        # 事实是：**那条约束只对 Python 包成立**（`pip freeze` 有等价的
        # 非子进程替代，所以没必要冒 EPERM 的险）；而问 R 包**只能**起
        # `Rscript`，它的答案只从 stdout 出来。所以这里必须捕获，
        # 并且在捕获失败时把失败**如实写进 reason**，而不是返回一个
        # 会被读成"R 说没装"的 None。两条通道的取舍不同，不是自相矛盾。
        p = subprocess.run([rscript, "-e", expr], capture_output=True,
                           text=True, timeout=180)
    except (OSError, subprocess.SubprocessError) as e:  # noqa: BLE001
        return {"rscript": rscript, "r_version": None,
                "packages": {p2: None for p2 in names},
                "reason": f"调用 Rscript 失败: {type(e).__name__}: {e}"}
    if p.returncode != 0:
        return {"rscript": rscript, "r_version": None,
                "packages": {p2: None for p2 in names},
                "reason": (f"Rscript 退出码 {p.returncode}: "
                           f"{(p.stderr or '').strip()[:200]}")}

    versions = {p2: None for p2 in names}
    r_version = None
    for line in (p.stdout or "").splitlines():
        parts = line.split("\t")
        if parts[0] == "RVERSION" and len(parts) > 1:
            r_version = parts[1].strip() or None
        elif parts[0] in ("OK", "MISSING") and len(parts) > 1:
            name = parts[1].strip()
            if name in versions and parts[0] == "OK" and len(parts) > 2:
                versions[name] = parts[2].strip() or None
    return {"rscript": rscript, "r_version": r_version, "packages": versions,
            "reason": None}


def capture_versions(cfg: dict, key_packages=None, extra: dict = None) -> dict:
    """§0.3 版本记录。

    Python 包全量走 `importlib.metadata` 枚举已安装发行版，**不起子进程** ——
    管道捕获输出在受限沙箱里会 EPERM，而这里拿到的信息与 `pip freeze` 等价。

    **R 包走另一条路**（`probe_r_packages`，起一次 `Rscript`）：`importlib`
    对 R 包永远返回 None，那个 None 会被读成"查过了，装不上"。两条通道
    分开走，且分开写原因 —— 混在一起会让"没装 R"看起来像"包不存在"。

    文档点名的关键工具单独放进 `key_versions`：全量 freeze 有几百行，
    关键工具淹没在里面。
    """
    from importlib import metadata as _md

    # 枚举失败**必须落盘**，不能只打一条 WARN（R-03 裁决 L12）。
    #
    # 旧写法 `except Exception as exc: log_warn(...)` 是 E-48 的形态：异常被
    # 降级成一行日志，而日志没有人读。后果是 `versions` 变成空字典（或只枚举
    # 到一半），验收项「版本记录非空」却仍然可能通过 —— 因为它的判据是
    # `n_versions >= 20`，一个**魔数**：既没人知道 20 从哪来，也不区分
    # "枚举了 96 个包"与"枚举器半路抛异常前凑够了 21 个"。
    # 现在把"枚举成功没有"本身记成结构化字段，验收直接读它。
    full = {}
    enum_error = None
    try:
        for d in _md.distributions():
            try:
                n = d.metadata["Name"]
            except Exception:
                continue
            if n:
                full[_norm_pkg(n)] = d.version
    except Exception as exc:
        enum_error = f"{type(exc).__name__}: {exc}"
        log_warn(f"枚举已安装包失败（{enum_error}）—— versions 会不完整")

    key = {}
    for p in (key_packages if key_packages is not None else KEY_PACKAGES):
        key[p] = full.get(_norm_pkg(p))

    # R 包：只有走默认清单时才探（调用方显式传 key_packages 时，
    # 说明它只要那几个，别往里塞东西）。
    r_probe = None
    if key_packages is None and R_KEY_PACKAGES:
        r_probe = probe_r_packages(R_KEY_PACKAGES)
        key.update(r_probe["packages"])

    for k, v in (extra or {}).items():
        key[k] = v

    m = read_manifest(cfg)
    m["versions"] = dict(sorted(full.items()))
    m["key_versions"] = key
    m["n_packages"] = len(full)
    # `versions_enumeration`: "ok" / "failed"。**与 `n_packages` 分开** ——
    # 计数是个下界，"枚举成功"是个事实，二者不能互相替代。
    m["versions_enumeration"] = "failed" if enum_error else "ok"
    if enum_error:
        m["versions_enumeration_error"] = enum_error
    # 关键工具解析率：**这才是"版本记录有没有用"的判据**（L12）。
    # `n_packages` 只看"装了 96 个包"，而 `key_versions` 里 14 个是 null
    # 完全不影响它 —— 实测基线 artifact 就是 n_versions=96 而 key 里
    # 只有 8/22 解析出来。真正的风险是"关键工具一个都没记上，而验收全绿"。
    m["n_key_resolved"] = sum(1 for v in key.values() if v)
    m["n_key_total"] = len(key)
    m["python"] = sys.version.split()[0]
    if r_probe is not None:
        # R 环境单独记一段：`key_versions` 是平铺的 name→version，
        # 读不出"哪些是 R 包、R 是什么版本、R 到底有没有"。
        m["r"] = {"rscript": r_probe["rscript"],
                  "r_version": r_probe["r_version"],
                  "packages": r_probe["packages"],
                  "reason": r_probe["reason"]}
    try:
        import platform
        m["platform"] = platform.platform()
    except Exception:
        pass
    m["dataset_id"] = cfg.get("dataset_id")
    m["updated_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    write_json(manifest_path(cfg), m)

    # 警告分两句：Python 的"没装"和 R 的"没装"原因完全不同，
    # 混成一句 `关键工具未安装（15/22）` 会把"问错了地方"也列进去。
    py_missing = sorted(p for p in key if p not in set(R_KEY_PACKAGES)
                        and key[p] is None)
    if py_missing:
        log_warn(f"关键 Python 工具未安装（{len(py_missing)}）：{', '.join(py_missing)}")
    if r_probe is not None:
        r_missing = sorted(p for p, v in r_probe["packages"].items() if v is None)
        if r_probe["reason"]:
            log_warn(f"R 包版本未取到：{r_probe['reason']}")
        elif r_missing:
            log_warn(f"关键 R 工具未安装（{len(r_missing)}，R "
                     f"{r_probe['r_version']}）：{', '.join(r_missing)}")
        else:
            log_info(f"关键 R 工具全部就位（{len(r_probe['packages'])} 个，"
                     f"R {r_probe['r_version']}）")
    if not py_missing and (r_probe is None or not r_probe["reason"]):
        log_info(f"关键工具版本已记录，共 {len(full)} 个已安装 Python 包")
    return key


def record_input(cfg: dict, path, label: str = "", required: bool = True) -> dict:
    """§0.4 输入数据哈希。

    文件不存在时**记 missing 而不是抛异常** —— 调用点未必知道某个输入
    这轮会不会产生（可选步骤的产物就是）。`required=True` 时 missing
    会在验收里被看见。
    """
    p = Path(path)
    entry = {"label": label or p.name, "path": str(p), "required": bool(required)}
    if p.exists() and p.is_file():
        entry["sha256"] = sha256_file(p)
        entry["bytes"] = p.stat().st_size
        entry["status"] = "present"
    else:
        entry["status"] = "missing"
    _manifest_append(cfg, "inputs", entry)
    return entry


def record_params(cfg: dict, params: dict) -> None:
    """§0.4 关键参数完整记录（含随机种子）。"""
    m = read_manifest(cfg)
    m.setdefault("params", {})
    m["params"].update(params or {})
    m["seed"] = (cfg.get("analysis") or {}).get("seed")
    m["dataset_id"] = cfg.get("dataset_id")
    m["updated_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    write_json(manifest_path(cfg), m)


def record_decision(cfg: dict, node: str, question: str, answer, evidence: str = "") -> None:
    """§0.4 Agent 决策链：从原始问题到最终结论的每一步推理。

    `evidence` 要写**支持这个选择的实际数字**，不是"因为这是通行做法"。
    没有量化依据的决策也要记，但 evidence 就写"没有量化依据"。
    """
    _manifest_append(cfg, "decisions", {
        "node": node, "question": question, "answer": answer,
        "evidence": evidence, "at": time.strftime("%Y-%m-%dT%H:%M:%S"),
    })


def record_human_review(cfg: dict, node: str, required: bool = True,
                        status: str = "pending", note: str = "") -> None:
    """§0.4 人工干预记录。

    status：pending（需确认，未确认）/ confirmed / overridden（人推翻了
    自动结果，note 写改成什么）/ not_needed（本数据集不涉及）。

    **默认 pending 而不是 confirmed。** 自动化流水线不能替人签字 ——
    把未确认的节点默认记成已确认，等于把复核节点变成摆设。
    """
    _manifest_append(cfg, "human_review", {
        "node": node, "required": bool(required), "status": status,
        "note": note, "at": time.strftime("%Y-%m-%dT%H:%M:%S"),
    })


def record_cross_language(cfg: dict, src: str, dst: str, fmt: str,
                          before: dict = None, after: dict = None,
                          lost=None, tool: str = "", note: str = "") -> None:
    """§0.2 跨语言转换记录。

    文档要求记录转换前后维度、metadata 字段数、丢失字段清单。桥接工具限定
    zellkonverter / anndata2ri，**禁止 sceasy**（维护状态差、metadata 丢失
    风险高）。

    本仓库与姊妹仓库之间只走 CSV，所以正常路径下 `before`/`after` 是
    行列数与列名集合；真正发生对象级转换时才填 `tool`。
    """
    _manifest_append(cfg, "cross_language", {
        "src": src, "dst": dst, "format": fmt, "tool": tool,
        "before": before or {}, "after": after or {},
        "lost_fields": sorted(lost or []),
        "note": note, "at": time.strftime("%Y-%m-%dT%H:%M:%S"),
    })


def manifest_summary(cfg: dict) -> dict:
    """给验收用的一行摘要。

    **必需项缺失和可选项缺失要分开报。** 可选项（如 `clinical.csv`）本来就
    允许不存在，把它算进"缺失"会让每个没有该文件的数据集都判失败 ——
    那是把"设计如此"当成"出错了"。两者都必须**可见**，但只有必需项判失败。
    """
    m = read_manifest(cfg)
    if not m:
        return {"present": False}
    ins = m.get("inputs") or []
    miss = [i for i in ins if i.get("status") == "missing"]
    key = m.get("key_versions") or {}
    return {
        "present": True,
        "n_versions": len(m.get("versions") or {}),
        # **"枚举成功"与"枚举到多少个"是两件事**（L12）。前者是事实，
        # 后者是下界 —— 枚举器半路抛异常时 `n_versions` 仍可能非零。
        "versions_enumeration": m.get("versions_enumeration", "unknown"),
        "versions_enumeration_error": m.get("versions_enumeration_error"),
        # 关键工具的解析率。`n_versions = 96` 但 `key_versions` 里 22 个只有
        # 8 个解析出来，是基线 artifact 的真实状态 —— 只看 `n_versions`
        # 完全看不见这件事。
        "n_key_resolved": m.get("n_key_resolved",
                                sum(1 for v in key.values() if v)),
        "n_key_total": m.get("n_key_total", len(key)),
        "key_unresolved": sorted(k for k, v in key.items() if not v),
        "n_inputs": len(ins),
        "inputs_missing": sorted(i.get("label") or "?" for i in miss),
        "inputs_missing_required": sorted(
            i.get("label") or "?" for i in miss if i.get("required")
        ),
        "n_decisions": len(m.get("decisions") or []),
        "human_review_pending": sorted(
            h["node"] for h in (m.get("human_review") or [])
            if h.get("status") == "pending"
        ),
        "n_cross_language": len(m.get("cross_language") or []),
    }


# ============================================================================
# 出图样式与调色板
# ============================================================================
# 规范来源：scientific-agent-skills/skills/scientific-visualization
#   assets/publication.mplstyle 与 assets/color_palettes.py（K-Dense，MIT）。
# 本仓库把样式文件放在 assets/publication.mplstyle，两处偏离写在该文件头部。
#
# **颜色不能是唯一线索。** Okabe-Ito 只是把颜色本身做成色盲友好；
# 分类图上仍要加 marker / 线型 / 直接标注，否则灰度打印就全糊了。
PAL = {
    # Okabe-Ito（Wong, Nature Methods 8:441, 2011）里**能做实心标记的**几个。
    #
    # **原来的八个标准色里有两个用不了**，理由是实测的，不是偏好：
    #   * `yellow` #F0E442 —— 白底对比度只有 **1.32:1**，实心圆点在白底上
    #     几乎看不见（`tools/check_palette.mjs` 会拦下来）；
    #   * `sky_blue` #56B4E9 —— 与 `blue` 色相只差约 5°，远小于 15° 判据，
    #     同一张图上会被读成同一个颜色。
    #
    # 所以这两个**从 PAL 里删掉了**：留着等于暗示"可以用"，
    # 而下一个想用它们的人只会得到一张看不清的图。理由记在这里。
    "orange": "#E69F00",
    "green": "#009E73",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "black": "#000000",
    # 语义别名 —— 脚本里用语义名，换配色时只改这里
    "primary": "#0072B2",     # 主序列（原 #2C7FB8）
    "highlight": "#D55E00",   # 阈值线 / 强调（原 #B2182B）
    "muted": "#999999",       # 次要参照（如随机基线）
    # 循环色扩到 10 色时补的（挑法见 PAL_CYCLE 的说明）
    "cyan": "#17BECF",
    "maroon": "#A50F15",
    "indigo": "#5B4FCF",
    "grey": "#A8A8A8",
}
# **分类循环色必须 >= 类别数，否则会撞色 —— 而撞色不报错。**
#
# 原来只有 5 色，而实测 `resolution=1.0` 给 **10 个簇**，于是
# 0≡5 / 1≡6 / 2≡7 / 3≡8 / 4≡9：图上簇 5 与簇 0 是**上下相邻的两个同色团块**，
# 读者看不出分界在哪；`trajectory_modules` 的 M0 与 M5 同理。
#
# 前 5 个**保持原样**（blue / vermillion / green / purple / black）——
# 这样 <=5 个类别的图（占多数）外观不变，改动的影响面最小。
#
# 后 5 个是**按 `tools/check_palette.mjs` 的三条判据挑出来的，不是凭眼睛选的**：
# 两两色相 >=15°、三种色盲（protanopia / deuteranopia / tritanopia）下
# OKLab 距离 >=0.05、白底对比度 >=2.0。
#
# 挑的过程值得记下来，因为**它证明了色空间是饱和的**：试了 13 个候选色，
# 12 个都被拦下，而且每个只差一项 ——
#
# | 候选 | 被什么拦下 |
# |---|---|
# | `#3D3D00` 暗橄榄 | protanopia 下与 maroon 撞（Δ=0.018）|
# | `#7F7F7F` 中灰 | deuteranopia 下与 green 撞（Δ=0.033）|
# | `#1F5C3A` 暗绿 | 色相与 green 只差 9.9° |
# | `#4A2A00` 暗棕 | 色相与 orange 只差 10.1° |
# | `#008080` 青绿 | tritanopia 下与 blue 撞，且色相与 cyan 只差 11.7° |
# | `#9467BD` 紫罗兰 | protanopia 下与 blue 撞（Δ=0.019）|
# | `#8C564B` 棕 | 色相与 maroon 只差 5.2° |
# | `#D3D3D3` 浅灰 | 对比度 1.50:1 |
# | `#B4B4B4` / `#B0B0B0` / `#ADADAD` | protanopia 下与 cyan 撞 |
# | `#A0A0A0` | deuteranopia 下与 purple 撞 |
#
# **规律：红绿色盲把 20°–110° 的暖色区压成一条轴，只剩亮度能区分。**
# 暖色区已有三个亮度级（maroon 0.460 / vermillion 0.621 / orange 0.753），
# 第 4 个暖色无论放哪个亮度都会撞上其中之一。
# 所以最后补的是**无彩色**（grey）：它不受色相判据约束，
# 且亮度 0.72 与所有彩色都拉得开。**grey 排最后** ——
# 它最不显眼，让它承担第 10 个簇而不是第 8 个。
#
# **10 色不是"够用"，是这三条判据下的上限。** 类别再多就要靠
# marker / 直接标注（这张 `umap_clusters` 已经在质心标了簇号）——
# 见上面那句"颜色不能是唯一线索"。
PAL_CYCLE = [PAL["blue"], PAL["vermillion"], PAL["green"], PAL["purple"],
             PAL["black"], PAL["orange"], PAL["cyan"], PAL["maroon"],
             PAL["indigo"], PAL["grey"]]

_STYLE_APPLIED = False


def apply_style(cfg: dict = None) -> None:
    """
    应用出版级样式。**幂等**，重复调用无副作用。

    样式文件在 assets/publication.mplstyle。找不到时退回手工设几个
    关键 rcParam 并警告 —— 静默用 matplotlib 默认样式会让图看起来
    "能出"但不符合任何投稿规范。
    """
    global _STYLE_APPLIED
    if _STYLE_APPLIED:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    style = Path(__file__).resolve().parent.parent.parent / "assets" / "publication.mplstyle"
    if style.exists():
        plt.style.use(str(style))
    else:
        log_warn(f"找不到样式文件 {style}，退回手工设置（图不符合投稿规范）")
        matplotlib.rcParams.update({
            "figure.constrained_layout.use": True,
            "savefig.bbox": "standard",
            "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
            "axes.spines.top": False, "axes.spines.right": False,
            "pdf.fonttype": 42, "ps.fonttype": 42,
        })
    # 分类循环色显式设一遍，不依赖样式文件被正确解析
    matplotlib.rcParams["axes.prop_cycle"] = matplotlib.cycler(color=PAL_CYCLE)
    # 顺序色标：viridis 感知均匀且对色盲友好
    matplotlib.rcParams["image.cmap"] = "viridis"
    _STYLE_APPLIED = True


# ============================================================================
# 出图
# ============================================================================
def mm(*vals: float):
    """
    毫米 → 英寸。**投稿图的尺寸单位是毫米，不是英寸。**

    常用宽度（Nature 的规范，见参考 skill 的 journal_requirements.md）：
    89 mm 单栏，183 mm 双栏，120-136 mm 单栏半。
    直接用英寸写 figsize 会让"单栏图"到底多宽变成一个没人检查的猜测。
    """
    if len(vals) == 1:
        return vals[0] / 25.4
    return tuple(v / 25.4 for v in vals)


# 三种标准宽度，单位英寸。绘图脚本一律用它们，不写裸英寸。
W_SINGLE = mm(89)      # 单栏
W_ONE_HALF = mm(136)   # 单栏半（Nature 允许 120-136 mm）
W_DOUBLE = mm(183)     # 双栏（= 满版宽）


def grid_figsize(n_cols: int, n_rows: int, *,
                 col_mm: float = 11.4, row_mm: float = 8.6,
                 min_w_mm: float = 136.0, max_w_mm: float = 183.0,
                 min_h_mm: float = 60.0, max_h_mm: float = 200.0):
    """按"列数 x 行数"算热图尺寸，**全程用毫米**，且**高度也封顶**。

    **L3 / L4（R-03 裁决）**：原来的写法是

        figsize=(min(W_DOUBLE, max(W_ONE_HALF, 0.45 * len(cols) + 3)),
                 max(3.0, 0.32 * len(rows) + 1.6))

    宽度是 `mm()` 出来的**英寸**值（被 `W_DOUBLE` 封顶），高度是**裸英寸、
    且不封顶**。两个单位混在同一个元组里，而 matplotlib 两个都收 ——
    **单位混用不会报错**，它只是把"装不装得进一页"变成一个没人检查的
    猜测。实测：24 个簇时高度 `0.32*24+1.6 = 9.28` 英寸 = **236 mm**，
    超过 A4 的可用高度，投稿时被排版工具缩放或裁掉。

    行/列的上限之外**不再线性增长**：超过一页就只能靠缩字号或拆图，
    线性拉长只会把图拉出页面。
    """
    w_mm = min(max_w_mm, max(min_w_mm, col_mm * max(1, n_cols) + 20.0))
    h_mm = min(max_h_mm, max(min_h_mm, row_mm * max(1, n_rows) + 14.0))
    return (mm(w_mm), mm(h_mm))


def _content_overflow(fig) -> dict:
    """
    检查内容有没有超出画布（= 被裁掉）。

    为什么需要这个检查：`savefig.bbox` 从 "tight" 改成 "standard" 之后，
    装不下的标签会被**直接裁掉**，而图文件照样生成、`check_figures.mjs`
    照样报"有墨迹" —— 只有打开图才看得出来。

    constrained layout 正常情况下会把内容塞进画布；这个检查兜住
    "某个图用了 add_axes / 手工 GridSpec，constrained layout 管不到"的情况。
    """
    try:
        fig.canvas.draw()
        tb = fig.get_tightbbox(fig.canvas.get_renderer())
        if tb is None:
            return {}
        w, h = fig.get_size_inches()
        # 留 2% 容差：constrained layout 会把 pad 也算进去，少量溢出是正常的
        ow = (tb.width - w) / w
        oh = (tb.height - h) / h
        bad = {}
        if ow > 0.02:
            bad["width_overflow_frac"] = round(float(ow), 4)
        if oh > 0.02:
            bad["height_overflow_frac"] = round(float(oh), 4)
        return bad
    except Exception:  # noqa: BLE001
        # 检查本身失败不该让出图失败
        return {}



def plot_marker_dotplot(ax, frac, zmat, *, cmap="RdBu_r",
                        size_max=170, size_min=10, dot_edge=0.3):
    """
    手工画 marker dotplot：颜色 = 按基因做 z-score（跨簇可比偏离方向），
    点大小 = 表达该基因的细胞比例。

    **为什么不再用 `sc.pl.dotplot(..., standard_scale="var")`。**
    读 scanpy 源码（`_prepare_dot_data`）确认：`standard_scale="var"` 是
    **逐基因 min-max 归一化到 0–1**（减最小值、除最大值），**不是 z-score**。
    而副标题写"z-scored" —— 两者矛盾（用户 2026-09-24 指出的核心问题）。
    两个口径的科学含义不同：

    ================  ===============================================
    口径              含义
    ================  ===============================================
    min-max（旧）     0 = 该基因在所有簇里的最低表达，1 = 最高
                      —— 只有"相对排名"，没有偏离方向
    z-score（本版）   0 = 平均水平，正 = 高于均值，负 = 低于均值
                      —— 跨簇可比"哪个簇偏离更大"，且 0 有锚点含义
    ================  ===============================================

    自带的三个额外好处（都是用户逐条点名的）：
      * 标签不再被裁 —— 布局由本函数控制，不与 scanpy 的 grid 抢空间；
      * 色标对称 —— RdBu_r 以 0 为中点，正负偏离等权可见；
      * 图例间距 —— 点大小图例单独画、间距显式给定，不再粘连。

    :param ax: 目标 axes
    :param frac: DataFrame，index=分组（行），columns=基因（列），值 0..1
    :param zmat: DataFrame，同形状，值 = 按基因 z-score 后的值
    :param cmap: 色标（默认 RdBu_r：红=高表达，蓝=低表达，白=0）
    :param size_max: 100% 表达时点的面积（pt^2）
    :param size_min: 0% 表达时点的面积（仍画一个极小点，表示"测到但几乎不表达"）
    :param dot_edge: 点描边宽度
    :returns: (ScalarMappable, size_handles) —— 供调用方画共享色标与大小图例
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    genes = list(zmat.columns)
    groups = list(zmat.index)
    ng, ngrp = len(genes), len(groups)

    # **z 轴对称归一**：正负偏离等权，0 永远是白色 —— 这就是"色标含负值"
    # 且"白色 = 平均水平（不是 0 表达）"的来源。上限取全部 |z| 的 95 分位
    # 再放宽到 >=1.5，避免个别极端基因把色标撑得两极分化。
    vmax = max(1.5, float(np.nanpercentile(np.abs(zmat.to_numpy()), 95)))
    norm = Normalize(vmin=-vmax, vmax=vmax)
    cmap_obj = plt.get_cmap(cmap)

    xs = np.arange(ng)
    ys = np.arange(ngrp)
    for gi in range(ng):
        for ri in range(ngrp):
            f = float(frac.iloc[ri, gi])
            z = float(zmat.iloc[ri, gi])
            s = size_min + f * (size_max - size_min)
            ax.scatter(xs[gi], ys[ri], s=s,
                       color=cmap_obj(norm(z)),
                       edgecolor="black", linewidth=dot_edge, zorder=3)

    ax.set_xticks(xs)
    ax.set_xticklabels(genes, rotation=90, fontsize=7)
    ax.set_yticks(ys)
    ax.set_yticklabels(groups, fontsize=8)
    ax.set_xlim(-0.7, ng - 0.3)
    ax.set_ylim(ngrp - 0.5, -0.5)
    ax.tick_params(length=2, pad=2)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap_obj)
    # 大小图例句柄：25/50/75/100% 四档（用户第七轮：五档太挤）
    size_handles = [(f, size_min + f * (size_max - size_min)) for f in
                    (0.25, 0.5, 0.75, 1.0)]
    return sm, size_handles


def build_marker_dotplot_figure(frac_df, z_df, *, group_label, title,
                                subtitle, fig_width=None, fig_height_mm=126):
    """
    构建**整张** marker dotplot（主图 + 底部横带图例），返回 `(fig, size_handles)`。

    **为什么把它抽成函数（这是本图第八轮才做对的关键）。** 前七轮图例布局的
    代码是**内联在调用脚本里**的，而验证脚本又**照抄了一份** —— 同一段布局
    存在三份镜像（scrna 脚本 / spatial 脚本 / 本地验证脚本）。于是：

    * 改一处、另两处不同步 → 验证脚本测的不是生产代码，**"全 PASS"是假的**；
    * 每次微调都要改三处，而 matplotlib 的位置参数**看不到效果**，
      只能靠一轮轮烧 CI 去猜 —— 实测调了七轮、每轮都要你重新看图。

    抽成一个函数后，**脚本与验证脚本调用的是同一份代码**，改一次全体生效。

    **图例为什么放底部横带。** 用户提供的参考代码
    （`可参考代码/99.单细胞：自动注释`）用 `scCustomize::do_DotPlot(dot.scale=12)`
    —— Seurat 生态的标准范式：图例在主图**下方一条横带**，左半点大小、
    右半横向色标。前七轮把两块图例竖着塞进右侧 1/6 宽的窄列，两块图例加
    两个标题挤在 0.4 图高的竖条里，**没有不受挤的排法**（实测调四轮仍相撞）。
    改横带后主图横向延展到全宽，两个标题各在自己图例正上方，与点列物理分离。

    **必须 `set_layout_engine("none")`。** 本仓库全局开着
    `figure.constrained_layout.use: True`（AGENTS 规则 13），而它会**在每次
    draw 时重排手动设的坐标** —— 这正是前几轮"改了没效果、底部留白越调越多"
    的原因：`fig.add_axes([...])` 设的位置被 constrained layout 覆盖了。
    这张图的几何由本函数全权控制，所以显式关掉它。

    :param frac_df: DataFrame，index=分组，columns=基因，值 = 表达细胞比例
    :param z_df: DataFrame，同形状，值 = 按基因 z-score 后的表达
    :param group_label: y 轴标签（"Leiden cluster" / "Spatial domain"）
    :param title: suptitle（如 "Top markers per cluster"）
    :param subtitle: 主图标题第二行（轴语义说明）
    :param fig_width: 图宽（默认 W_DOUBLE）
    :param fig_height_mm: 图高（毫米）
    :returns: `(fig, size_handles)`
    """
    import matplotlib.pyplot as plt

    w = W_DOUBLE if fig_width is None else fig_width
    fig = plt.figure(figsize=(w, mm(fig_height_mm)))
    # **关掉 constrained layout** —— 见 docstring：它会覆盖下面所有手动坐标
    fig.set_layout_engine("none")

    # 主图占满上方（左右各留一点给 y 轴标签与右缘）
    ax = fig.add_axes([0.075, 0.315, 0.905, 0.535])
    sm, size_handles = plot_marker_dotplot(ax, frac_df, z_df)
    ax.set_xlabel("gene", fontsize=8)
    ax.set_ylabel(group_label, fontsize=8)
    fig.suptitle(title, fontsize=11, y=0.965)
    ax.set_title(subtitle, fontsize=8, pad=6, loc="left")

    # ---- 底部横带 · 左半：Percent Expressed (%)（四点横排）---------------
    lax = fig.add_axes([0.075, 0.045, 0.42, 0.13])
    lax.set_xlim(0, 1)
    lax.set_ylim(0, 1)
    lax.axis("off")
    lax.text(0.0, 0.88, "Percent Expressed (%)", ha="left", va="center",
             fontsize=7.5)
    for k, (f_, s_) in enumerate(size_handles):
        xx = 0.06 + k * 0.20
        lax.scatter([xx], [0.42], s=s_, color="gray",
                    edgecolor="black", linewidth=0.3)
        lax.text(xx, 0.02, f"{int(f_ * 100)}", ha="center", va="center",
                 fontsize=7.5)

    # ---- 底部横带 · 右半：Mean Expression（横向色标）--------------------
    cax = fig.add_axes([0.60, 0.085, 0.28, 0.035])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cax.text(0.0, 1.9, "Mean Expression", transform=cax.transAxes,
             ha="left", va="bottom", fontsize=7.5)
    # **L11（R-03 裁决）**：刻度原来写死 `[-1, 0, 1]`，而色标上限是
    # `vmax = max(1.5, p95|z|)`（见 `plot_marker_dotplot`）。p95|z| > 1.5
    # 时（数据里有少数强 marker 时很常见），真实的红端是 ±2.5 甚至 ±3，
    # 而刻度只标到 ±1 —— **色标看起来只有中间一小段有数据、两端"空着"**，
    # 读者会以为红/蓝端对应的是"极端值"，实际上 1 就已经是饱和色。
    # 按 `sm.norm.vmax` 生成刻度，让刻度与色标实际范围对得上。
    _vmax = float(sm.norm.vmax)
    if not math.isfinite(_vmax) or _vmax <= 0:
        _vmax = 1.0
    _vmax = round(_vmax, 1)
    cb.set_ticks([-_vmax, 0.0, _vmax])
    cb.set_ticklabels(["Low", "Mid", "High"])
    cb.ax.tick_params(labelsize=7, top=False, bottom=True,
                      labeltop=False, labelbottom=True)
    return fig, size_handles


def place_labels(ax, xs, ys, texts, fontsize=7, pad_px=2.0,
                 max_iters=400, seed=0):
    """
    在散点图上放**互不重叠**的文字标签（纯几何，不引第三方依赖）。

    **为什么需要它。** `adjustText` 不在本仓库依赖里（规则：不要临时引入），
    而"固定偏移 + 上下交替"这类纯参数法在**点挤成一条竖列**时必然失效 ——
    实测 `02-07-03-unit1-tf-specificity-scatter`：8 个 TF 的 x 都在 0–0.05，
    上下交替只把标签分成两层，同层内仍然互压（用户反馈"基因标签有重叠"）。

    做法：把每个标签候选位置按"离锚点由近及远"排序（右、左、上、下、
    四个对角共 8 个方向 × 多圈），**贪心**选第一个与已放标签不重叠的位置；
    都放不下就继续外扩。重叠判据用 matplotlib 的文本包围盒（渲染后实测，
    不是估算字宽）。

    **保证**：只要画布还有空间，返回的标签两两不重叠；实在放不下的会被
    推远（有引导线时仍可读）。返回实际使用的位置列表，便于日志记录。

    :param ax: 目标 axes（需已完成 scatter 且坐标范围已定）
    :param xs: 锚点 x（数据坐标）
    :param ys: 锚点 y（数据坐标）
    :param texts: 标签文字，与 xs/ys 等长
    :param fontsize: 标签字号（pt）
    :param pad_px: 标签之间要求的最小间隙（像素）
    :param max_iters: 每个标签最多尝试的候选位置数
    :param seed: 保留参数（本函数确定性，不用随机；留给调用方对齐口径）
    :returns: `[(x, y), ...]` 实际放置的**数据坐标**
    """
    import itertools

    fig = ax.figure
    # **必须先把 figure 的 dpi 对齐到"实际保存用的 dpi"再量包围盒。**
    #
    # 实测（2026-09-24）：本函数在创建 figure 时的默认 dpi（100）下量尺寸并
    # 排布，而 `save_fig()` 用 `dpi=300` 重新渲染 —— 文本的**字号是点、与 dpi
    # 无关，但包围盒是像素**，于是 300 dpi 下每个标签都放大 3 倍、而按点给的
    # 偏移没跟着放大，布局全散。实测 dpi=100 时 0 重叠、dpi=300 时
    # **28 对重叠 + 8 个越界** —— 正是 CI 图上"仍重叠"的原因。
    #
    # 所以：先设 dpi（与 `save_fig` 一致），再 draw、再量、再排。
    target_dpi = 300.0
    fig.set_dpi(target_dpi)
    fig.canvas.draw()  # 需要 renderer 才能量文本包围盒
    renderer = fig.canvas.get_renderer()

    # 候选方向：8 个方位，由近及远多圈外扩。
    # **优先向右/左上** —— 这类散点的点在左侧挤成竖列（x≈0），右侧是空的。
    dirs = [(1, 0), (1, 1), (1, -1), (0, 1), (0, -1),
            (-1, 1), (-1, -1), (-1, 0)]
    placed_boxes = []
    used = []
    # **标签必须留在坐标轴内。** 只判"标签之间不重叠"是不够的 ——
    # 实测（run 35974071689）：点在左上角挤成一团时，贪心把 8 个标签
    # 一路往右上推，**全部推出画布顶边**、还压住标题。
    # 所以加一条包含判据：标签包围盒必须完全落在 axes 内。
    ax_bb = ax.get_window_extent(renderer=renderer)

    def _mk(text, x, y, off_pt, dx, dy, arrow=False):
        """建一个标签。`off_pt` 是**以点为单位的偏移** —— 不是像素。"""
        kw = {}
        if arrow:
            kw["arrowprops"] = dict(arrowstyle="-", lw=0.4,
                                    color=PAL.get("muted", "#999999"))
        return ax.annotate(
            text, (x, y), xytext=off_pt, textcoords="offset points",
            fontsize=fontsize,
            ha="left" if dx > 0 else ("right" if dx < 0 else "center"),
            va="bottom" if dy > 0 else ("top" if dy < 0 else "center"),
            color=PAL.get("ink", "#1A1A1A"), **kw)

    for x, y, text in zip(xs, ys, texts):
        chosen = None
        best_outside = None   # 实在放不进时的次优（溢出最少的一个）
        for ring in itertools.count(1):
            if ring * len(dirs) > max_iters:
                break
            step_pt = 4.0 + (ring - 1) * 5.0   # 点（1/72 英寸），与 DPI 无关
            for dx, dy in dirs:
                # **`textcoords="offset points"` 要的是偏移量，不是绝对坐标。**
                # 实测踩过：这里原先传的是 `transData.transform()` 出来的
                # **像素绝对坐标**，于是标签被推到画布外几万像素处 ——
                # 包含判据全部拒绝，最后统统落进兜底分支挤在一起。
                off = (dx * step_pt, dy * step_pt)
                t = _mk(text, x, y, off, dx, dy)
                raw = t.get_window_extent(renderer=renderer)
                bb = raw.expanded(1 + pad_px / max(raw.width, 1),
                                  1 + pad_px / max(raw.height, 1))
                if any(bb.overlaps(o) for o in placed_boxes):
                    t.remove()
                    continue
                inside = (bb.x0 >= ax_bb.x0 and bb.x1 <= ax_bb.x1 and
                          bb.y0 >= ax_bb.y0 and bb.y1 <= ax_bb.y1)
                if not inside:
                    t.remove()
                    if best_outside is None:
                        over = (max(0, ax_bb.x0 - bb.x0) + max(0, bb.x1 - ax_bb.x1) +
                                max(0, ax_bb.y0 - bb.y0) + max(0, bb.y1 - ax_bb.y1))
                        best_outside = (over, off, dx, dy)
                    continue
                placed_boxes.append(bb)
                chosen = off
                used.append(off)
                break
            if chosen is not None:
                break
        if chosen is None:
            # **放不进就带引导线放最近处**，不静默丢弃、也不推到画布外。
            if best_outside is not None:
                _, off, dx, dy = best_outside
            else:
                off, dx, dy = (18.0, 0.0), 1, 0
            _mk(text, x, y, off, dx, dy, arrow=True)
            used.append(off)
    return used


def _record_figure_overflow(cfg: dict, name: str, bad: dict) -> None:
    """把"内容超出画布"落盘成产物，而不只是打一行 WARN。

    Q-26 / E-49 的教训：`_content_overflow()` **早就检测到了**
    `{'width_overflow_frac': 0.3523}`（与本地复现逐位吻合）、也**打了 WARN**，
    而 `02-08-01` 的标题两侧仍被静默裁掉 35% —— 因为**告警没有消费者**。
    "人读 CI 日志"已被证明会漏（E-30/E-31 同族，但这次连要读的那行字都是
    程序打出来的）。

    所以判据强度要从**警告级**升到**产物级**：写进 `figure_overflow.json`，
    由 `main_analysis.py` 的验收项读它并判红。文件写在 `figures_dir` 的上一级
    （即 `results/<dataset>/`），与其余 status/产物同级。
    """
    p = Path(cfg["output"]["figures_dir"]).parent / "figure_overflow.json"
    try:
        cur = json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}
    except Exception:
        cur = {}
    figs = cur.get("figures") or {}
    figs[name] = bad
    cur["figures"] = figs
    cur["n_overflow"] = len(figs)
    cur["threshold_note"] = (
        "超 2% 才记（_content_overflow 的容差）。width/height 键分别代表"
        "横向/纵向超出画布的比例；savefig.bbox=standard 下超出的部分会被"
        "**静默裁掉**，文件照样生成、图名图幅墨迹四条门禁全绿。")
    try:
        write_json(p, cur)
    except Exception as e:  # 落盘失败不该让出图失败，但必须吼一声
        log_warn(f"图 {name} 的溢出记录写不进 {p}：{type(e).__name__}: {e}")


def save_fig(cfg: dict, name: str, fig=None, tight: bool = False) -> list:
    """
    保存一张图为 PNG + PDF，并返回写出的路径。

    **PDF 与 PNG 都要出。** PDF 是矢量、可再编辑、字体以 Type 42 内嵌；
    PNG 用于快速像素级非空白检查（check_figures 解 PNG）。

    **`bbox_inches` 默认不再是 "tight"。** 参考规范
    （scientific-visualization）明确写着 tight 会改变输出的物理尺寸 ——
    投稿要求"单栏 89 mm"时，tight 出来的就不是 89 mm。装不下由
    constrained layout 解决，另有 `_content_overflow()` 兜底并把溢出
    **落盘**成 `figure_overflow.json`（供验收判红，见 `_record_figure_overflow`）。

    **不要在绘图代码里调用 plt.savefig 后不管 plt.close** —— 不关的话
    同一进程里后续的图会叠在旧 figure 上，产出"看着正常但内容错"的图。
    """
    import matplotlib.pyplot as plt

    apply_style(cfg)
    figdir = Path(cfg["output"]["figures_dir"])
    figdir.mkdir(parents=True, exist_ok=True)
    f = fig if fig is not None else plt.gcf()

    bad = _content_overflow(f)
    if bad:
        log_warn(f"图 {name} 的内容超出画布（{bad}）—— 标签可能被裁掉。"
                 f"调大 figsize 或改用 layout='constrained'")
        _record_figure_overflow(cfg, name, bad)

    written = []
    for ext in ("png", "pdf"):
        p = figdir / f"{name}.{ext}"
        f.savefig(p, dpi=cfg["analysis"]["figure_dpi"],
                  bbox_inches="tight" if tight else None)
        written.append(str(p))
    plt.close(f)
    return written


# ============================================================================
# 小工具
# ============================================================================
#
# **2026-09-26 删除三个零调用点的函数（R-03 裁决 L8 / L9）。**
#
# 删掉的是 `fix_dotplot_legends`（~180 行）、`annotate_commit`、`fmt_bytes`、
# `require_pkg` —— 四个名字在全仓库的**唯一**出现处就是它们的 `def` 行
# （逐个 grep 确认，见台账 E-66）。
#
# 为什么删而不是留着：
#
# * `fix_dotplot_legends` 是给 `sc.pl.dotplot` 的后处理链用的。E-36 当时写
#   「保留给仍在用 scanpy 路径的图」，但那个前提**已经不存在了** ——
#   `03_cluster_annotate.py` 改成手绘 `build_marker_dotplot_figure` 之后，
#   全仓 `sc.pl.*` 只剩 `05_trajectory.py` 的 `sc.pl.paga`，没有一处会产出
#   scanpy 的 DotPlot。**保留一个"等以后有人用"的函数，等于留一段没人跑过
#   的代码**：它下次被调用时，scanpy 的内部结构（`_plot_size_legend` /
#   `_plot_colorbar` 的 axes 布局）早已变过好几版，而没有任何门禁能发现。
#   真要再走 scanpy 路径，应当照当时的版本重写并当场标定，而不是复活这段。
# * `annotate_commit` 想做的是「把 commit 短哈希写进图名后缀」，但 geo 仓
#   那件事是**由 workflow 的 shell 做的**（`geo_analysis.yml` 的
#   "Write versioned copies of the figures" 步骤），scrna/spatial 两仓的
#   workflow 里根本没有对应步骤 —— 所以它不是"还没接线"，是**接线方式选错了
#   层**。留一个永远返回 `"local"`（CI 外）或 `sha[:7]`（CI 内）却没人调的
#   函数，只会让下一个人以为图名里已经有版本号。
# * `fmt_bytes` / `require_pkg` 同理：前者没有任何日志/状态文件用到字节数，
#   后者在本仓被 `importlib.metadata` 与 `probe_r_packages` 取代（geo 侧的
#   R 版 `require_pkg` 仍在用，那是**另一份代码**，不受此处影响）。
#
# 同一条教训的另一面写在 AGENTS 规则 15.1：判据的输入域要和它要防的缺陷
# 匹配。这里反过来 —— **没有判据覆盖的代码，删除比保留安全**。
#
# 注意：`spatial-pipeline-skill/scripts/lib/common.py` 里 `fix_dotplot_legends`
# 与 `require_pkg` **同样零调用点**，但那是 R-04 的范围，本次不动。

def df_to_records(df) -> list:
    """DataFrame -> JSON 安全的 records（NaN 变 None）。"""
    import pandas as pd
    if df is None or len(df) == 0:
        return []
    return json.loads(df.where(pd.notna(df), None).to_json(orient="records"))
