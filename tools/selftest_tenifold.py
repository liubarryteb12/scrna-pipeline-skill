#!/usr/bin/env python3
"""
tools/selftest_tenifold.py — scTenifoldKnk 结果处理自检（不需要装 R）

**这是代码测试，不是分析。** 它拿几张小假表喂给
`08_virtual_perturbation.py` 的 `compress_tenifold_distances()` 与
`compare_engines()`，验证「空敲除」那条保护真的生效。

## 为什么需要它

`scTenifoldKnk` 的 `transcriptomeWide` 模式敲除一个基因的方式是把 WT 网络的
**那一行清零**（`scTenifoldKnk.R` L206–L207：`KO <- WT; KO[g, ] <- 0`）。
如果该基因在网络里**出度本来就是 0**，清一行全 0 的行等于什么都没敲，
`KO` 与 `WT` 逐位相同，`manifoldAlignment` 对两个完全一样的网络对齐，
`dRegulation` 返回的"距离"只剩浮点噪声（实测 ~1e-16）。

**这个失败不报错、不给 NA、不给零** —— 它给出一排看起来完全正常的数。
读表的人只会得出"敲除 FOXM1 没有影响"，而真相是这个网络表达不了该扰动。
所有既有门禁（图非空 / 图名 / 图幅 / dpi / 配色 / 图例）都会通过，
因为它们检查的是"有没有产出"而不是"产出对不对"。

这条保护一旦失效，**没有别的机会发现它**。而 CI 上真跑一次 tenifold
要 ~6.5 分钟且依赖 R + P3M 二进制；本自检把这段逻辑单独拿出来，
用假数据在几毫秒内跑完，**在没有 R 的机器上也能跑**。

## 与 R 侧自检的分工

- `scripts/lib/tenifold_knk.R --selftest`：验证 **R 侧**（矩阵方向断言、
  `zero_outdegree_targets` 的出度计算、读写契约）。
- 本文件：验证 **Python 侧**（meta 里的结构证据有没有被正确消费成
  None / 标记列，以及空敲除有没有被排除在相关性之外）。

两边都必要：R 侧算出 `empty_knockout_genes` 而 Python 侧把它读丢，
症状同样是"输出一排正常的数"。

用法：
    python tools/selftest_tenifold.py
退出码：0 = 代码路径可用；1 = 失败
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts" / "lib"))
sys.path.insert(0, str(REPO / "scripts"))


def _load_08():
    spec = importlib.util.spec_from_file_location(
        "vp", REPO / "scripts" / "08_virtual_perturbation.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    import numpy as np
    import pandas as pd
    from scipy import sparse

    vp = _load_08()
    ok = True

    def _say(msg):
        print(f"[selftest_tenifold] {msg}")

    def _fail(msg):
        nonlocal ok
        _say(f"**FAIL**：{msg}")
        ok = False

    # ---- 1. 正常基因：距离压成均值 / 最大值 / 最强下游 ---------------------
    #
    # 真数据里的量级参考（pbmc3k，17×600 距离矩阵）：min 1.7e-18、
    # max 1.2e-03、mean 7.2e-07。这里用同量级的数，免得自检通过
    # 而真数据因为量级差几个数量级走进别的分支。
    genes = ["TBX21", "GATA2", "FOXM1"]
    net_genes = ["SRGN", "GZMB", "CD3D", "IL7R", "NKG7"]
    dist = pd.DataFrame(
        [[1.2e-03, 4.0e-04, 1.0e-06, 0.0, 3.0e-05],
         [9.0e-04, 2.0e-04, 5.0e-07, 0.0, 1.0e-05],
         [3.865e-16, 2.1e-16, 9.0e-17, 0.0, 1.5e-16]],
        index=genes, columns=net_genes)
    rmeta = {
        "empty_knockout_genes": ["FOXM1"],
        "target_outdegree": {"TBX21": 7, "GATA2": 4, "FOXM1": 0},
        "n_genes": 600,
    }

    tdf = vp.compress_tenifold_distances(dist, rmeta)

    if list(tdf["gene"]) != ["TBX21", "GATA2", "FOXM1"]:
        _fail(f"排序不对 —— 期望按平均距离降序 [TBX21, GATA2, FOXM1]，"
              f"实得 {list(tdf['gene'])}")
    # 空敲除行必须沉底（NaN 排在最后）。这一条如果坏了，FOXM1 会
    # 因为"距离是 0"排到最前面 —— 正好与事实相反。
    if tdf.iloc[-1]["gene"] != "FOXM1":
        _fail(f"空敲除的 FOXM1 没有沉底，排在第 "
              f"{list(tdf['gene']).index('FOXM1') + 1} 位 —— "
              "NaN 排序坏掉会让'什么都没敲'看起来像'影响最小'")

    r0 = tdf[tdf["gene"] == "TBX21"].iloc[0]
    if r0["tenifold_n_genes_scored"] != 5:
        _fail(f"TBX21 的 n_genes_scored={r0['tenifold_n_genes_scored']}，"
              f"应该是 5")
    if r0["tenifold_top_gene"] != "SRGN":
        _fail(f"TBX21 的最强下游是 {r0['tenifold_top_gene']!r}，应该是 SRGN")
    if abs(r0["tenifold_max_distance"] - 1.2e-03) > 1e-9:
        _fail(f"TBX21 的 max_distance={r0['tenifold_max_distance']}，"
              f"应该是 0.0012")
    if r0["tenifold_empty_knockout"]:
        _fail("TBX21 被误标成空敲除 —— 它在网络里有 7 条出边")
    if r0["tenifold_target_outdegree"] != 7:
        _fail(f"TBX21 的 target_outdegree={r0['tenifold_target_outdegree']}，"
              f"应该是 7（meta 里的结构证据没有被读出来）")

    # ---- 2. 空敲除：距离字段必须是 None，**不能是 0.0** -------------------
    #
    # 这是本自检存在的全部理由。`0.0` 和 `None` 在 CSV 里长得完全一样
    # （都是空或 0），但在读表的人眼里含义相反：
    #   0.0  → "敲除 FOXM1 没有影响"（一个**结论**）
    #   None → "这个网络表达不了该扰动"（一个**拒绝回答**）
    r2 = tdf[tdf["gene"] == "FOXM1"].iloc[0]
    for col in ("tenifold_mean_distance", "tenifold_max_distance",
                "tenifold_top_gene", "tenifold_top_distance"):
        v = r2[col]
        if v is not None and not (isinstance(v, float) and np.isnan(v)):
            _fail(f"空敲除的 FOXM1 报了 {col}={v!r} —— 必须是 None。"
                  f"报了 0.0 就等于说'敲除没有影响'，而真相是这个网络"
                  f"表达不了该扰动")
    if not r2["tenifold_empty_knockout"]:
        _fail("FOXM1 在 meta 的 empty_knockout_genes 里，"
              "但没被标成空敲除")
    if r2["tenifold_target_outdegree"] != 0:
        _fail(f"FOXM1 的 target_outdegree={r2['tenifold_target_outdegree']}，"
              f"应该是 0（出度为 0 正是它被判为空敲除的原因）")
    if r2["tenifold_n_genes_scored"] != 0:
        _fail(f"空敲除的 n_genes_scored={r2['tenifold_n_genes_scored']}，"
              f"应该是 0")

    # ---- 3. 判据必须是结构证据，不是「整行都是 NaN」------------------------
    #
    # 若把判据写成"这一行全是 NA 就是空敲除"，那么 R 侧因为**别的原因**
    # 给了 NA（比如某个下游基因没进网络、包返回 NaN）也会被误标成
    # 空敲除，然后它的真实原因被掩盖成"网络表达不了该扰动"。
    # 这里造一个"全 NA 但不在 empty_knockout_genes 里"的基因：
    # 距离字段照样该是 None（没数可报），但 **tenifold_empty_knockout
    # 必须是 False** —— 两种 NA 的原因不同，不能混成一个标记。
    dist_na = pd.DataFrame(
        [[1.0e-03, 2.0e-04, 3.0e-05],
         [np.nan, np.nan, np.nan]],
        index=["TBX21", "MYSTERY"], columns=["SRGN", "GZMB", "CD3D"])
    tdf_na = vp.compress_tenifold_distances(
        dist_na, {"empty_knockout_genes": [], "target_outdegree": {}})
    rna = tdf_na[tdf_na["gene"] == "MYSTERY"].iloc[0]
    if rna["tenifold_empty_knockout"]:
        _fail("MYSTERY 整行是 NA 但**不在** empty_knockout_genes 里，"
              "却被标成空敲除 —— 判据退化成了「整行都是 NA」，"
              "会把 R 侧别的 NA 原因误标成'网络表达不了该扰动'")
    if rna["tenifold_n_genes_scored"] != 0:
        _fail(f"全 NA 行的 n_genes_scored={rna['tenifold_n_genes_scored']}，"
              f"应该是 0")

    # ---- 4. 空敲除不能进相关性 -------------------------------------------
    #
    # 空敲除的"距离"是浮点噪声（~1e-16）。把它算进 Spearman 等于
    # 往相关系数里掺随机数 —— 而且它排在最后，对秩相关的影响不小。
    # 这里断言它被**显式剔除并记账**，而不是靠 dropna() 悄悄少几个点。
    fo = pd.DataFrame({
        "gene": ["TBX21", "GATA2", "FOXM1"],
        "ko_magnitude": [2.52, 5.41, 2.50],
    })
    cmp = vp.compare_engines(fo, tdf)
    if cmp.get("excluded_empty_knockout") != ["FOXM1"]:
        _fail(f"空敲除没有被排除在相关性之外："
              f"excluded_empty_knockout={cmp.get('excluded_empty_knockout')}，"
              f"应该是 ['FOXM1']")
    # 这条同样要挂在"确实排除了东西"上 —— 若一个都没排除，
    # 不写 excluded_reason 才是对的。
    if cmp.get("excluded_empty_knockout") and "excluded_reason" not in cmp:
        _fail("排除了空敲除却没写原因 —— 下一个人会以为只是少了个数据点")
    if cmp.get("n_shared_genes") != 2:
        _fail(f"n_shared_genes={cmp.get('n_shared_genes')}，"
              f"应该是 2（3 个基因里排掉 1 个空敲除）")
    # 只剩 2 个共享基因，少于 3 个不该报相关系数 ——
    # 两个点永远能连成一条线。
    #
    # **这条断言必须挂在 `n_shared_genes == 2` 上**：如果空敲除没被排除
    # （上面那条已经判红了），共享基因会变成 3 个，此时报相关系数是
    # **正确的**。不挂条件的话，自检会在一个正确行为上打出一条
    # "只有 2 个共享基因却报了相关系数" 的假消息 —— 而这句假话
    # 会把排查方向带偏（去看相关性计算，而真问题是空敲除没排除）。
    if cmp.get("n_shared_genes") == 2 and "spearman_rho" in cmp:
        _fail("只有 2 个共享基因却报了相关系数 —— 两个点永远能连成一条线")

    # ---- 5. 对照：没有空敲除时不出现 excluded 键 --------------------------
    #
    # 反过来也要断言。若 `excluded_empty_knockout` 无条件出现（哪怕是
    # 空数组），读表的人分不清"这轮没有空敲除"和"这轮忘了算"。
    tdf_clean = vp.compress_tenifold_distances(
        dist.iloc[:2], {"empty_knockout_genes": [], "target_outdegree": {}})
    cmp_clean = vp.compare_engines(fo[fo["gene"] != "FOXM1"], tdf_clean)
    if "excluded_empty_knockout" in cmp_clean:
        _fail("这轮没有空敲除，却写了 excluded_empty_knockout —— "
              "空数组和'忘了算'长得一样")

    # ---- 6. 端到端：跑一遍 run_tenifold_engine（用假的 Rscript）-------------
    #
    # **这一段是有来历的。** 把压行逻辑抽成 `compress_tenifold_distances()`
    # 之后，`run_tenifold_engine` 的 `return` 里还留着一句
    # `"empty_knockout_genes": sorted(empty_ko)` —— 而 `empty_ko` 已经
    # 跟着函数一起搬走了。`python -m py_compile` **查不出来**（它是合法
    # 语法，只是名字不存在），静态门禁也查不出来，本地没有 R 所以
    # 整条路径跑不到。它只会在 **CI 上 tenifold 真跑通那一刻** 抛
    # NameError，而那时 R 已经跑了 6.5 分钟。
    #
    # 前面 5 段测的是抽出来的那个纯函数，**碰不到 `return` 语句** ——
    # 所以它们全绿，缺陷照样在。这一段把 `subprocess.run` 换成假的，
    # 让 `run_tenifold_engine` 完整走到底。
    import json
    import shutil as _shutil
    import tempfile
    import types

    td_genes = ["TBX21", "GATA2", "FOXM1"]
    td_net = ["SRGN", "GZMB", "CD3D", "IL7R", "NKG7"]

    def _fake_run(cmd, **kw):
        """假 Rscript：把该写的两个产物写进 --out_dir，返回退出码 0。"""
        out = Path(cmd[cmd.index("--out_dir") + 1])
        rows = ["TBX21,1.2e-03,4.0e-04,1.0e-06,0.0,3.0e-05",
                "GATA2,9.0e-04,2.0e-04,5.0e-07,0.0,1.0e-05",
                # 空敲除：R 侧置的 NA
                "FOXM1,NA,NA,NA,NA,NA"]
        (out / "tenifold_perturbation_distances.csv").write_text(
            '"","' + '","'.join(td_net) + '"\n'
            + "\n".join(f'"{r.split(",")[0]}",'
                        + r.split(",", 1)[1] for r in rows) + "\n",
            encoding="utf-8")
        (out / "tenifold_meta.json").write_text(json.dumps({
            "engine": "scTenifoldKnk", "engine_version": "1.1",
            "n_genes": len(td_net), "n_cells": 20,
            "empty_knockout_genes": ["FOXM1"],
            "target_outdegree": {"TBX21": 7, "GATA2": 4, "FOXM1": 0},
        }), encoding="utf-8")
        return types.SimpleNamespace(returncode=0, stdout="TENIFOLD_OK\n",
                                     stderr="")

    with tempfile.TemporaryDirectory() as _td:
        _out = Path(_td)
        _orig_run = vp.subprocess.run
        _orig_which = vp.shutil.which
        _orig_probe = vp._probe_r_package
        try:
            vp.subprocess.run = _fake_run
            vp.shutil.which = lambda name: "Rscript"   # 本地没有 R，假一个
            vp._probe_r_package = lambda pkg: {
                "available": True, "reason": "", "version": "1.1",
                "r_version": "4.6.1"}
            # 计数矩阵要同时含**候选基因**（被敲的）与网络基因（被量的）——
            # `run_tenifold_engine` 先按 var_names 筛候选，再从中选网络基因集。
            # 只给网络基因的话候选一个都筛不出来，会早退成 no_candidates。
            _all = td_genes + td_net
            _X = sparse.csr_matrix(
                np.abs(np.random.default_rng(0).integers(
                    0, 50, size=(20, len(_all)))).astype(float))
            # 用 try 包住：`run_tenifold_engine` 里的 `return` 语句如果引用了
            # 已经不存在的名字（抽出函数时最容易漏的那种），这里会抛
            # NameError。裸抛也能让 CI 变红，但报出来是一长串 traceback，
            # 看不出"是 return 语句坏了"还是"R 侧真的失败了"。
            try:
                _res = vp.run_tenifold_engine(
                    {}, _X, _all,
                    pd.DataFrame({"gene": td_genes}), _out,
                    {"max_genes": 600, "n_net": 2, "n_cells": 10, "n_comp": 3,
                     "td_k": 3, "ma_n_dim": 2, "n_cores": 1,
                     "timeout_sec": 60})
            except Exception as _exc:  # noqa: BLE001
                _fail(f"run_tenifold_engine 抛异常："
                      f"{type(_exc).__name__}: {_exc} —— "
                      f"这一段端到端路径在本地没有 R 时跑不到，"
                      f"`py_compile` 也查不出这类名字错误")
                _res = None
        finally:
            vp.subprocess.run = _orig_run
            vp.shutil.which = _orig_which
            vp._probe_r_package = _orig_probe

        # 已经报过异常就不再补一条"没走到 ok" —— 那是同一件事的第二种说法，
        # 两条消息叠在一起反而让人以为是两个问题。
        if _res is None:
            pass                       # 异常已经报过，不再重复
        elif _res.get("status") != "ok":
            _fail(f"run_tenifold_engine 没走到 ok（status="
                  f"{_res.get('status')!r}, reason={_res.get('reason')!r}）—— "
                  f"整条端到端路径没被跑到，前面几段全绿也说明不了它没问题")
        else:
            if _res.get("n_empty_knockout") != 1:
                _fail(f"端到端 n_empty_knockout="
                      f"{_res.get('n_empty_knockout')}，应该是 1")
            if _res.get("empty_knockout_genes") != ["FOXM1"]:
                _fail(f"端到端 empty_knockout_genes="
                      f"{_res.get('empty_knockout_genes')!r}，应该是 ['FOXM1']")
            if _res.get("n_genes_network") != len(td_net):
                _fail(f"n_genes_network={_res.get('n_genes_network')}，"
                      f"应该是 {len(td_net)}")
            for _f in ("virtual_perturbation_tenifold.csv",
                       "virtual_perturbation_tenifold_top.csv",
                       "tenifold_perturbation_distances.csv"):
                if not (_out / _f).exists():
                    _fail(f"端到端跑完但没写出 {_f}")
            # 临时矩阵必须被清掉 —— `_tenifold_tmp` 在 results/ 下，
            # 而 results/ 是整目录上传的
            if (_out / "_tenifold_tmp").exists():
                _fail("端到端跑完后 `_tenifold_tmp` 还在 —— "
                      "它会被整个塞进 artifact")

    # ---- 7. 验收判据 `_vp_engine_consistent` ------------------------------
    #
    # 这条判据原来是"方法串里必须出现 scTenifoldKnk 和 PerturbNet 两个词"，
    # 用来证明**没跑**这两个工具。K-01b 之后它会**把正确的结果判红**，
    # 所以换成了自洽性检查 —— 而自洽性检查有个共同的风险：
    # 它可能松到"什么都通过"。下面正反两面都测。
    import importlib.util as _ilu
    _ma_spec = _ilu.spec_from_file_location(
        "ma", REPO / "scripts" / "main_analysis.py")
    ma = _ilu.module_from_spec(_ma_spec)
    _ma_spec.loader.exec_module(ma)

    _ok_method = "一阶近似 + scTenifoldKnk（R 引擎）；**不是 PerturbNet**"
    _base = {"engine": "both", "engines_used": ["first_order", "tenifold"],
             "method": _ok_method, "tenifold": {"status": "ok"}}

    if not ma._vp_engine_consistent(_base):
        _fail("两套引擎都跑了、方法串也点名了，判据却判红 —— "
              "这条会把正确的结果判红，比不判更糟")

    # 反例 1：跑了 tenifold 却不在方法串里点名
    _a = dict(_base, method="一阶近似；**不是 PerturbNet**")
    if ma._vp_engine_consistent(_a):
        _fail("跑了 tenifold 却没在方法串里点名，判据仍然通过 —— "
              "读状态的人会以为只有一阶近似")

    # 反例 2：tenifold 失败但没写原因
    _c = dict(_base, engines_used=["first_order"],
              method="一阶近似；**不是 PerturbNet**",
              tenifold={"status": "failed"})
    if ma._vp_engine_consistent(_c):
        _fail("tenifold 失败了却没写原因，判据仍然通过 —— "
              "'没跑'和'跑了没结果'会长得一样")

    # 反例 3：没否掉 PerturbNet
    _d = dict(_base, method="一阶近似 + scTenifoldKnk（R 引擎）")
    if ma._vp_engine_consistent(_d):
        _fail("方法串没否掉 PerturbNet，判据仍然通过")

    # 正例 2：tenifold 失败但写明了原因 —— 这是合法状态，必须通过
    _b = dict(_base, engines_used=["first_order"],
              method="一阶近似；tenifold 没跑（Rscript 退出码 1）；"
                     "**不是 PerturbNet**",
              tenifold={"status": "failed", "reason": "Rscript 退出码 1"})
    if not ma._vp_engine_consistent(_b):
        _fail("tenifold 失败且写明了原因，判据却判红 —— "
              "环境问题不该被当成分析错误")

    # ---- 8. `_sig6` 必须保住排序（K-01b 第二版修的那个真 bug）-------------
    #
    # 第一版写的是 `round(v, 6)`，而 tenifold 距离的量级是 1e-9 ~ 1e-3 ——
    # 6 位**小数**把 1.48e-08 归成 0.0、把 1.09e-06 压成 1e-06。
    # 后果不是"数字不好看"：**排序被摧毁**，而排序是这张表存在的全部意义
    # （图直接画这列）。下面这段的判据就是"压缩前后顺序必须一致"。
    _raw = [4.5093e-06, 1.4774e-06, 1.4413e-06, 1.1208e-06, 1.0944e-06,
            8.1719e-07, 5.5998e-07, 4.5300e-07, 3.8342e-07, 2.0413e-07,
            7.8427e-08, 1.4847e-08, 1.0909e-08, 8.2891e-09, 1.1387e-09]
    _sig = [vp._sig6(v) for v in _raw]
    if any(s == 0.0 for s in _sig):
        _fail(f"_sig6 把 {sum(1 for s in _sig if s == 0.0)} 个非零距离压成了 "
              f"0.0 —— 这正是 `round(v, 6)` 干的事，会让弱的一半在图上消失")
    if len(set(_sig)) != len(set(_raw)):
        _fail(f"_sig6 把 {len(set(_raw))} 个不同的距离压成了 "
              f"{len(set(_sig))} 个 —— 并列之后排序变成任意的")
    _ord = sorted(range(len(_raw)), key=lambda i: -_raw[i])
    _ord2 = sorted(range(len(_raw)), key=lambda i: -_sig[i])
    if _ord != _ord2:
        _fail(f"_sig6 改变了排序：{_ord} -> {_ord2}")
    # 反向断言：旧的写法**必须**被这段判据抓住。若哪天有人把 `_sig6`
    # 改回 `round(v, 6)` 而这里恰好也松了，这条就是最后一道网。
    if len(set(round(v, 6) for v in _raw)) == len(set(_raw)):
        _fail("这段判据本身失效了 —— `round(v, 6)` 竟然没压平这批数，"
              "说明测试数据换了量级，该断言已经证明不了任何事")
    # 非有限值原样返回（nan 不能变成字符串 'nan' 写进 CSV）
    if not np.isnan(vp._sig6(float("nan"))):
        _fail("_sig6(nan) 没原样返回 nan")
    if vp._sig6(0.0) != 0.0:
        _fail(f"_sig6(0.0)={vp._sig6(0.0)!r}，应该是 0.0")

    # ---- 9. `load_targets` 内部回退分支必须按 TF 去重 ----------------------
    #
    # `tf_regulons.csv` 是**行级**排名：同一个 TF 在不同簇上各有一行
    # （真数据 217 行 / 201 个唯一 TF）。`reg.head(n)` 是行级截断，于是
    # 20 个候选里只有 17 个唯一（TBX21 / IRF1 / EZH2 各两次），下游
    # `virtual_perturbation.csv` 出现 21 行重复的 (gene, cell_type)。
    #
    # 这条路径**不能靠图门禁发现** —— 图照样画得出来，只是 x 轴或
    # 候选表里多了几行。它只能在代码层测。
    _reg = pd.DataFrame({
        "tf": ["TBX21", "TBX21", "IRF1", "IRF1", "EZH2", "EZH2", "GATA2"],
        # 已按 cluster_specificity 降序（`drop_duplicates` 保首行，
        # 所以留下的必须是特异性最高的那一行）
        "cluster_specificity": [9.0, 3.0, 8.0, 2.0, 7.0, 1.0, 6.0],
    })
    # `load_targets` 会经 `record_decision()` 往清单里写东西，所以需要一个
    # 真目录。**单独开一个 tempdir**，不复用第 6 段那个 —— 那个 `with` 块
    # 已经退出、目录已被删掉，复用它会让 `write_json` 把目录重新建出来，
    # 而那个目录从此再也没人清理。
    with tempfile.TemporaryDirectory() as _td2:
        _cfg = {"output": {"results_dir": str(Path(_td2) / "res"),
                           "data_dir": str(Path(_td2) / "data")},
                "dataset_id": "selftest", "perturbation": {"top_n": 4}}
        _tgt, _src = vp.load_targets(_cfg, _reg)
    if _src != "internal_top_regulons":
        _fail(f"没配 targets_csv 时应该走内部回退，实得 source={_src!r}")
    if list(_tgt["gene"]) != ["TBX21", "IRF1", "EZH2", "GATA2"]:
        _fail(f"load_targets 没按 TF 去重：候选 = {list(_tgt['gene'])}，"
              f"应该是 ['TBX21', 'IRF1', 'EZH2', 'GATA2'] —— "
              f"重复行会占掉名额，去重后反而能多出候选")
    if len(set(_tgt["gene"])) != len(_tgt):
        _fail(f"候选表里仍有重复基因：{list(_tgt['gene'])}")

    _say("自检" + ("通过" if ok else "**失败**"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
