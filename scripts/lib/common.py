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
import os
import random
import sys
import time
from pathlib import Path

import numpy as np
import yaml

# ============================================================================
# 确定性
# ============================================================================
# **必须在 import scanpy/numba 之前设。** numba 在首次 JIT 时读这些环境变量
# 决定线程数，晚设无效。Part 1 的教训：浮点末位分叉有两个独立来源 ——
# 线程调度（多线程归约的求和顺序）和 OpenBLAS 运行期按 CPU 型号选 SIMD 内核。
# 这里两个都钉住。
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("OPENBLAS_CORETYPE", "Haswell")
os.environ.setdefault("PYTHONHASHSEED", "0")


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


def write_json(path, obj) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)

    def default(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.bool_,)):
            return bool(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, Path):
            return str(o)
        return str(o)

    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, ensure_ascii=False, indent=2, default=default)


def read_json(path):
    p = Path(path)
    if not p.exists():
        return None
    try:
        with open(p, encoding="utf-8") as fh:
            return json.load(fh)
    except Exception:
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


def record_step(cfg: dict, step_id: str, status: str, seconds: float = None,
                message: str = "", required: bool = True) -> None:
    st = read_state(cfg)
    st.setdefault("steps", [])
    st["steps"] = [s for s in st["steps"] if s.get("id") != step_id]
    entry = {"id": step_id, "status": status, "required": bool(required),
             "message": message}
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
KEY_PACKAGES = [
    # §2 核心 Python 包
    "scanpy", "anndata", "scvi-tools", "cellbender", "harmonypy", "scvelo",
    "celltypist", "pyscenic", "liana", "doubletdetection", "scrublet",
    # §2.6 拟时序
    "palantir", "scfates", "cytotrace",
    # §2 点名的 R 包（本仓库无 rpy2 路径，正常就是 None）
    "monocle3", "slingshot", "cellchat", "soupx", "scdblfinder",
    # §1.7/§1.8 虚拟扰动。**三个都装不上，但键必须留着** ——
    # 记 None 是"查过了，装不上"，省略键是"没查"。
    # 确切原因写在 virtual_perturbation_status.json 的 tools 字段里。
    "scTenifoldKnk", "PerturbNet", "RegVelo",
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
            "本仓库 CI 没有 R + rpy2；CellBender / scVelo / pySCENIC 有 PyPI 真包，"
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
    """PEP 503 归一化：包名大小写与 -/_/. 不敏感。"""
    return str(name).strip().lower().replace("_", "-").replace(".", "-")


def capture_versions(cfg: dict, key_packages=None, extra: dict = None) -> dict:
    """§0.3 版本记录。

    全量走 `importlib.metadata` 枚举已安装发行版，**不起子进程** ——
    管道捕获输出在受限沙箱里会 EPERM，而这里拿到的信息与 `pip freeze` 等价。

    文档点名的关键工具单独放进 `key_versions`：全量 freeze 有几百行，
    关键工具淹没在里面。
    """
    from importlib import metadata as _md

    full = {}
    try:
        for d in _md.distributions():
            try:
                n = d.metadata["Name"]
            except Exception:
                continue
            if n:
                full[_norm_pkg(n)] = d.version
    except Exception as exc:
        log_warn(f"枚举已安装包失败（{exc}）—— versions 会不完整")

    key = {}
    for p in (key_packages if key_packages is not None else KEY_PACKAGES):
        key[p] = full.get(_norm_pkg(p))
    for k, v in (extra or {}).items():
        key[k] = v

    m = read_manifest(cfg)
    m["versions"] = dict(sorted(full.items()))
    m["key_versions"] = key
    m["n_packages"] = len(full)
    m["python"] = sys.version.split()[0]
    try:
        import platform
        m["platform"] = platform.platform()
    except Exception:
        pass
    m["dataset_id"] = cfg.get("dataset_id")
    m["updated_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    write_json(manifest_path(cfg), m)

    missing = sorted(k for k, v in key.items() if v is None)
    if missing:
        log_warn(f"关键工具未安装（{len(missing)}/{len(key)}）：{', '.join(missing)}")
    else:
        log_info(f"关键工具全部就位（{len(key)} 个），共记录 {len(full)} 个已安装包")
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
    return {
        "present": True,
        "n_versions": len(m.get("versions") or {}),
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


def verticalize_dotplot_size_legend(fig, title=None):
    """
    把 `sc.pl.dotplot` 的**点大小图例从横排转成纵排**（约定 v2）。

    **为什么必须后处理。** scanpy 的 `DotPlot` 没有控制图例方向的参数
    （`legend()` 只收 `width` / `show_size_legend` / `colorbar_title`），
    而 `_plot_size_legend()` 内部把示例点画在 **x 轴**上
    （`scatter(arange(len(size))+0.5, repeat(0,...))`）—— 就是横排。
    实测两个仓库的 dotplot 都是"老问题图例横着排布、示例横向"
    （用户 2026-09-24 对 02-03-03 / 03-03-02 的反馈）。

    做法：找到那个大小图例 axes（特征是**只有 x 刻度标签、没有 y 刻度**，
    且含一个 PathCollection），把它清空后**按原尺寸竖着重画一遍**，
    刻度标签保留原值。这样示例点的相对大小不变，只是排布方向变了。

    **不伪造数据**：示例点的面积直接从原 scatter 的 `get_sizes()` 读回，
    不重新推算（重算要复制 scanpy 的 step 规则，容易与上游分叉）。

    :param fig: dotplot 所在的 figure
    :param title: 可选的新标题（原 `colorbar_title` 语义）
    :returns: 是否成功转换（False = 没找到符合特征的 axes，图保持原样）
    """
    import numpy as np

    for ax in fig.axes:
        # 大小图例的判据：有 x 刻度标签、**没有** y 刻度标签、含散点
        if ax.get_yticklabels() and any(t.get_text() for t in ax.get_yticklabels()):
            continue
        colls = [c for c in ax.collections if hasattr(c, "get_sizes")]
        if not colls:
            continue
        labels = [t.get_text() for t in ax.get_xticklabels()]
        if not labels or not all(labels):
            continue
        # 只认"看起来像百分比数字"的刻度，避免误伤其它图
        try:
            [float(s) for s in labels]
        except ValueError:
            continue

        sizes = colls[0].get_sizes()
        if len(sizes) != len(labels):
            continue

        # **先把 legend axes 撑高**：scanpy 给这块区域的高度是按"横排一行点"
        # 分配的，竖排 5 个点 + 5 个刻度标签塞不进原高度 —— 实测标签被裁成
        # `10/20/30/40` 的残影（run 35974091761）。
        # 做法：把这个 axes 的物理高度按点数放大，并把它的下边界上移，
        # 让它占住 legend 区里原本空着的上方空间（那里是 top_spacer）。
        pos = ax.get_position()
        n = len(sizes)
        # 每个点至少需要约 0.022 图高（含标签），上限不超过 0.55（别盖住主图）
        need_h = min(max(pos.height, 0.022 * n + 0.05), 0.55)
        ax.set_position([pos.x0, pos.y1 - need_h, pos.width, need_h])

        ax.clear()
        # 竖排：从下往上依次变大（与横排的左小右大语义一致）
        ys = np.arange(n)
        # **点贴近轴的左缘**：`set_yticklabels` 把标签画在轴**左侧**，
        # 原写法把点画在轴中央（x=0），于是"标签 ……… 点"中间留一大片空白。
        # 点移到左缘附近后，标签紧贴其左，一行读作"标签 — 点"。
        # （不能把标签放右侧：图例列本来就在最右，右侧没有空间，会被裁。）
        ax.scatter(np.full(n, 0.12), ys, s=sizes,
                   color="gray", edgecolor="black", linewidth=0.5, zorder=100)
        # 给上下留半个间距，否则首尾的点和标签贴边被裁
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(-0.6, n - 0.4)
        ax.set_yticks(ys)
        ax.set_yticklabels(labels, fontsize="small")
        ax.set_xticks([])
        ax.tick_params(axis="x", bottom=False, labelbottom=False)
        ax.tick_params(axis="y", left=False, labelleft=True, pad=1)
        for sp in ax.spines.values():
            sp.set_visible(False)
        if title:
            ax.set_title(title, fontsize="small", pad=3, loc="left")
        return True
    return False


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
    import matplotlib.transforms as mtransforms

    fig = ax.figure
    fig.canvas.draw()  # 需要 renderer 才能量文本包围盒
    renderer = fig.canvas.get_renderer()
    inv = ax.transData.inverted()

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

    for x, y, text in zip(xs, ys, texts):
        px, py = ax.transData.transform((x, y))
        chosen = None
        best_outside = None   # 实在放不进时的次优（离轴内最近的一个）
        for ring in itertools.count(1):
            if ring * len(dirs) > max_iters:
                break
            step = 6.0 + (ring - 1) * 7.0  # 像素：每圈外扩
            for dx, dy in dirs:
                cx, cy = px + dx * step, py + dy * step
                t = ax.annotate(
                    text, (x, y), xytext=(cx, cy), textcoords="offset pixels",
                    fontsize=fontsize,
                    ha="left" if dx > 0 else ("right" if dx < 0 else "center"),
                    va="bottom" if dy > 0 else ("top" if dy < 0 else "center"),
                    color=PAL.get("ink", "#1A1A1A"))
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
                        # 记录"溢出最少"的那个（后面兜底用）
                        over = (max(0, ax_bb.x0 - bb.x0) + max(0, bb.x1 - ax_bb.x1) +
                                max(0, ax_bb.y0 - bb.y0) + max(0, bb.y1 - ax_bb.y1))
                        best_outside = (over, cx, cy, dx, dy)
                    continue
                placed_boxes.append(bb)
                chosen = inv.transform((cx, cy))
                used.append(chosen)
                break
            if chosen is not None:
                break
        if chosen is None:
            # **放不进就带引导线放最近处**，不静默丢弃、也不推到画布外。
            if best_outside is not None:
                _, cx, cy, dx, dy = best_outside
                ha = "left" if dx > 0 else ("right" if dx < 0 else "center")
                va = "bottom" if dy > 0 else ("top" if dy < 0 else "center")
            else:
                cx, cy, ha, va = px + 30, py, "left", "center"
            ax.annotate(text, (x, y), xytext=(cx, cy), textcoords="offset pixels",
                        fontsize=fontsize, ha=ha, va=va,
                        color=PAL.get("ink", "#1A1A1A"),
                        arrowprops=dict(arrowstyle="-", lw=0.4,
                                        color=PAL.get("muted", "#999999")))
            used.append(inv.transform((cx, cy)))
    return used


def save_fig(cfg: dict, name: str, fig=None, tight: bool = False) -> list:
    """
    保存一张图为 PNG + PDF，并返回写出的路径。

    **PDF 与 PNG 都要出。** PDF 是矢量、可再编辑、字体以 Type 42 内嵌；
    PNG 用于快速像素级非空白检查（check_figures 解 PNG）。

    **`bbox_inches` 默认不再是 "tight"。** 参考规范
    （scientific-visualization）明确写着 tight 会改变输出的物理尺寸 ——
    投稿要求"单栏 89 mm"时，tight 出来的就不是 89 mm。装不下由
    constrained layout 解决，另有 `_content_overflow()` 兜底并告警。

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

    written = []
    for ext in ("png", "pdf"):
        p = figdir / f"{name}.{ext}"
        f.savefig(p, dpi=cfg["analysis"]["figure_dpi"],
                  bbox_inches="tight" if tight else None)
        written.append(str(p))
    plt.close(f)
    return written


def annotate_commit(cfg: dict) -> str:
    """把 commit 短哈希写进文件名后缀，便于把图与代码版本对上。"""
    sha = os.environ.get("GITHUB_SHA", "")
    return sha[:7] if sha else "local"


# ============================================================================
# 小工具
# ============================================================================
def fmt_bytes(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} GB"


def require_pkg(name: str, hint: str = "") -> None:
    try:
        __import__(name)
    except ImportError as e:
        extra = f"（{hint}）" if hint else ""
        raise RuntimeError(f"缺少依赖 {name}{extra}: {e}") from e


def df_to_records(df) -> list:
    """DataFrame -> JSON 安全的 records（NaN 变 None）。"""
    import pandas as pd
    if df is None or len(df) == 0:
        return []
    return json.loads(df.where(pd.notna(df), None).to_json(orient="records"))
