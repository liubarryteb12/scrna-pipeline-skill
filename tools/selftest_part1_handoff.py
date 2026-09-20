#!/usr/bin/env python3
"""
tools/selftest_part1_handoff.py — §0.2 跨语言交接代码路径自检

**这是代码测试，不是分析。** 它验证 `08_virtual_perturbation.py` 的
Part 1 → Part 2 CSV 交接那一段能跑通、能把 `cross_language` 记进清单，
仅此而已。它造的 CSV **不是真的 Part 1 输出** —— 基因名是从调控子里抄的，
logFC 是编的。

为什么需要它：`assets/config.pbmc3k.yml` 的 `perturbation.targets_csv`
是 `null`，所以 CI 每轮走的都是**回退分支**（`internal_top_regulons`）。
那验证了"正确回退"，但**没有验证真正的交接代码路径** ——
包括 `record_cross_language()` 到底有没有被调用。

而"没跑到的代码"和"跑对了的代码"在 CI 日志里长得一模一样。
本仓库已经踩过同一个坑两次（`07_grn` 漏 import 却验收全绿、
`record_decision()` 定义了三个仓库但两个从没被调用），所以这里主动把它跑到。

**这个自检抓到过一次真问题**：`record_cross_language()` 在 `common.py` 里
定义了很久，但**一次都没被调用** —— 清单的 `cross_language` 永远是空数组，
而空数组和"这一轮没有跨语言转换"长得一模一样。

用法：
    python tools/selftest_part1_handoff.py --config assets/config.pbmc3k.yml
退出码：0 = 代码路径可用；1 = 失败
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts" / "lib"))
sys.path.insert(0, str(REPO / "scripts"))

from common import load_config, log_info, read_manifest  # noqa: E402


def _load_08():
    spec = importlib.util.spec_from_file_location(
        "vp", REPO / "scripts" / "08_virtual_perturbation.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    import pandas as pd
    import tempfile

    cfg = load_config(args.config)

    # ---- 0. 全程在临时目录里跑，**不碰 data/ 和 results/** --------------------
    #
    # 交接表要是留在 `data/` 下，`check_artifact_paths.mjs` 会把它当成
    # "脚本写出了清单外的产物"；而 `cross_language` 条目要是写进**真的**清单，
    # 读者会以为这一轮真的发生了 Part 1 交接 —— **那正是这条自检要防的
    # "看起来一样"**。所以 data_dir 和 results_dir 都指向临时目录。
    with tempfile.TemporaryDirectory(prefix="selftest_handoff_") as tmp:
        tmpdir = Path(tmp)
        data_dir = tmpdir / "data"
        res_dir = tmpdir / "results"
        data_dir.mkdir(parents=True, exist_ok=True)
        res_dir.mkdir(parents=True, exist_ok=True)

        cfg_selftest = json.loads(json.dumps(cfg))      # 深拷一份，不改真配置
        cfg_selftest.setdefault("output", {})
        cfg_selftest["output"]["data_dir"] = str(data_dir)
        cfg_selftest["output"]["results_dir"] = str(res_dir)
        cfg_selftest.setdefault("perturbation", {})

        # ---- 1. 造一份"Part 1 交接表"----------------------------------------
        #
        # 列名故意和 Part 1 的 `09_export_targets.R` 对齐（`gene` + `logfc`），
        # 并**多带一列 `module_label`** —— 那是 Part 1 有、但按 §0.2 不交接的
        # 字段。自检要验证的正是"多出来的列被记进 `lost_fields` 了没有"。
        fixture = data_dir / "part2_targets_selftest.csv"
        genes = ["GATA2", "TBX21", "IRF7", "STAT1", "MS4A1", "CD79A",
                 "NKG7", "GZMB", "FCGR3A", "LYZ", "CCR7", "IL7R"]
        pd.DataFrame({
            "gene": genes,
            "logfc": [1.8, -0.9, 0.4, 0.0, -1.2, 0.7,
                      1.1, -0.3, 0.6, -0.8, 0.2, -0.5],
            "module_label": ["turquoise"] * 6 + ["blue"] * 6,
        }).to_csv(fixture, index=False)
        log_info(f"自检造出交接表 {fixture.name}：{len(genes)} 个基因，3 列")

        # ---- 2. 走真正的 load_targets() 交接分支 -----------------------------
        #
        # **直接复用生产代码** —— 自检要测的就是它，不是另写一遍。
        vp = _load_08()
        cfg_selftest["perturbation"]["targets_csv"] = fixture.name
        # 交接表已存在，`load_targets` 应该走 `part1_handoff` 分支而不是回退

        reg = pd.DataFrame({"tf": genes})     # 回退分支用不到，给个占位
        out, source = vp.load_targets(cfg_selftest, reg)
        log_info(f"load_targets -> source={source}，{len(out)} 个基因，"
                 f"含 logFC 的 {int(out['logfc'].notna().sum())} 个")

        ok = True

        # ---- 3. 断言：走了交接分支，而不是回退 ------------------------------
        if not source.startswith("part1_handoff:"):
            log_info(f"**FAIL**：应该走 part1_handoff 分支，实际是 {source} —— "
                     "说明交接表没被识别（列名匹配或路径解析坏了）")
            ok = False
        if len(out) != len(genes):
            log_info(f"**FAIL**：应该有 {len(genes)} 个基因，实际 {len(out)}")
            ok = False
        if int(out["logfc"].notna().sum()) != len(genes):
            log_info("**FAIL**：logFC 没有全部解析出来 —— "
                     "列名别名匹配坏了，下游 signature_alignment 会是 NaN")
            ok = False

        # ---- 4. 断言：cross_language 真的进了清单 ----------------------------
        #
        # **这一条才是自检的主要目的。** `record_cross_language()` 定义了但
        # 不调用时，清单里 `cross_language` 是空数组 —— 而空数组和
        # "这一轮没有跨语言转换"长得一模一样。
        m = read_manifest(cfg_selftest)
        cl = (m or {}).get("cross_language") or []
        if not cl:
            log_info("**FAIL**：清单里 cross_language 是空的 —— "
                     "record_cross_language() 没有被调用，或者没写进清单")
            ok = False
        else:
            e = cl[-1]
            log_info(f"cross_language 已登记：{e.get('src')} -> {e.get('dst')}"
                     f"（format={e.get('format')}）")
            for k in ("src", "dst", "format", "before", "after", "lost_fields"):
                if k not in e:
                    log_info(f"**FAIL**：cross_language 条目缺字段 {k}")
                    ok = False
            lost = e.get("lost_fields") or []
            log_info(f"  丢失字段（{len(lost)} 个）：{lost}")
            # **`module_label` 必须出现在丢失清单里。** 它就是"Part 1 有、
            # 但按 §0.2 不交接"的那一类 —— 不记下来，将来要用时没人知道
            # 该回 Part 1 的哪张表取。
            if not any("module_label" in str(x) for x in lost):
                log_info("**FAIL**：`module_label` 没有被记进 lost_fields —— "
                         "上游多出来的列被静默丢掉了")
                ok = False
            # **反过来也要断言：留下的不能出现在丢失清单里。**
            # 第一版把"除 gene 外的所有列"都算成丢失，于是 `logfc`
            # 被报成"丢了" —— 而它明明是唯一跨部分传过来的数值列。
            # **把留下的说成丢掉的，方向和事实正好相反**，比不记更糟。
            # 这条断言就是被那次云端日志逼出来的（CI 打出
            # `丢失字段（2 个）：['logfc', 'module_label']`）。
            _wrong = [x for x in lost if "logfc" in str(x).lower()
                      and "列名" not in str(x)]
            if _wrong:
                log_info(f"**FAIL**：{_wrong} 被报成丢失，但它是**留下的** —— "
                         "`after.kept_columns` 里明明写着 logfc")
                ok = False
            # before/after 的维度必须是数，不是占位符
            b, a = e.get("before") or {}, e.get("after") or {}
            if b.get("n_cols") != 3 or a.get("n_genes") != len(genes):
                log_info(f"**FAIL**：before/after 维度不对 —— "
                         f"before={b}, after={a}")
                ok = False
            if a.get("n_with_logfc") != len(genes):
                log_info(f"**FAIL**：after.n_with_logfc={a.get('n_with_logfc')}，"
                         f"应该是 {len(genes)}")
                ok = False

        # ---- 5. 第二例：列名别名 + 列序打乱 ---------------------------------
        #
        # 第一例的列名是"理想情况"。Part 1 真交接过来的表列名大小写不定
        # （`logFC` / `log2FC` / `log_fc`），列序也不保证。**别名匹配坏了
        # 不报错，只会让 `signature_alignment` 静默变成 NaN** ——
        # 而 NaN 和"真的没有相关性"长得一样。
        #
        # 这一例还断言**改名本身被记进 note**：上游那列叫 `log2FC`，
        # 下游叫 `logfc`，不说清楚下一个人会以为上游本来就叫 `logfc`。
        _fx2 = data_dir / "part2_targets_selftest_alias.csv"
        pd.DataFrame({
            "symbol": genes,                       # gene 的别名列名
            "n_sources": [3] * len(genes),         # 该丢的
            "log2FC": [1.8, -0.9, 0.4, 0.0, -1.2, 0.7,
                       1.1, -0.3, 0.6, -0.8, 0.2, -0.5],
        }).to_csv(_fx2, index=False)
        cfg_selftest["perturbation"]["targets_csv"] = _fx2.name
        out2, src2 = vp.load_targets(cfg_selftest, reg)
        log_info(f"别名例：source={src2}，{len(out2)} 个基因，"
                 f"含 logFC 的 {int(out2['logfc'].notna().sum())} 个")
        if not src2.startswith("part1_handoff:"):
            log_info(f"**FAIL**：`symbol` 列没被认成基因列（source={src2}）—— "
                     "Part 1 换个列名就静默回退了")
            ok = False
        if int(out2["logfc"].notna().sum()) != len(genes):
            log_info("**FAIL**：`log2FC` 没被认成 logFC 列 —— "
                     "下游 signature_alignment 会是 NaN，而且不报错")
            ok = False
        _m2 = read_manifest(cfg_selftest)
        _e2 = (_m2.get("cross_language") or [{}])[-1]
        _lost2 = _e2.get("lost_fields") or []
        log_info(f"  别名例丢失字段：{_lost2}")
        if not any("n_sources" in str(x) for x in _lost2):
            log_info("**FAIL**：`n_sources` 没被记进丢失字段")
            ok = False
        if any("log2fc" in str(x).lower() and "列名" not in str(x)
               for x in _lost2):
            log_info("**FAIL**：`log2FC` 被报成丢失，但它是留下的那一列")
            ok = False
        if "log2FC" not in str(_e2.get("note") or ""):
            log_info("**FAIL**：列名归一（`log2FC` → `logfc`）没有记进 note —— "
                     "下游会以为上游本来就叫 `logfc`")
            ok = False

    log_info("自检" + ("通过" if ok else "**失败**"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
