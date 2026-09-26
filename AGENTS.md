# 单细胞流水线的工程规则

> **治理层（2026-09-20 起）**：本仓库是 `scientific_agent_skill` 工作区三仓库之一，
> 受工作区治理层约束：一切产物只写工作区内；任务先登记在
> `governance/02_TASKLIST.md`；推送前跑 `node governance/hooks/pre-push.mjs`；
> checkpoint 台账见 `governance/04_CHECKPOINT_PLAN.md`；行为准则
> `governance/01_SPEC_v1.0.md`。冲突按规范 §7.4 报告裁决。

姊妹项目 `geo-normal-pipeline-skill/AGENTS.md` 的规则在这里同样适用。
本文只写**单细胞特有的**部分。

---

## 1. 样本量的单位是生物学重复，不是细胞

**永远不要在细胞层面跑差异检验。** 把每个细胞当一个独立样本，
是把细胞数当成了样本量。10000 个细胞来自 3 个供体，自由度是 2
不是 9997 —— 而 Wilcoxon 会给出 `p = 1e-300`。

这类假阳性**看起来极其显著**，是最容易骗过自己的错误。

正确做法：按 (样本 × 分组 × 细胞类型) 加总原始计数做拟bulk。
没有分组信息时记 `not_configured`，**不退化成细胞层面的检验**。

## 2. 需要特定基因的分析必须用 `adata.raw`（全基因集）

配体、受体、转录因子大多是低表达的，几乎不进 HVG。
实测：HVG 子集上 38 对配体-受体只有 3 对可用，全集上 27 对。

**假阴性看起来和"真的没有"一模一样**，所以宁可报错也不要悄悄退回 HVG。
`06_communication.py` 和 `07_grn.py` 都有 `get_full_expression()`，
`adata.raw` 为空时直接抛异常。

## 3. seurat_v3 口味的 HVG 必须在 normalize 之前算

`seurat_v3` 内部对**原始计数**做负二项拟合，所以要在
`normalize_total` + `log1p` 之前跑；`seurat` / `cell_ranger` 口味
相反，要在之后。

**顺序反了不报错，只是给出错的 HVG** —— 而 HVG 错了下游全部跟着错，
图上完全看不出来。

缺 `scikit-misc` 时降级到 `seurat` 口味，但**必须记录口味变了**
（`hvg_flavor_requested` vs `hvg_flavor_used`）—— 两种口味选出的
基因集合不同，下游 PCA/聚类/marker 全部跟着变。

## 4. 可选步骤的"没做"必须和"做了没问题"长得不一样

这是从 Part 1 继承的最重要一条。

`qc_status.json` 里 `doublet_detection.status` 有这些取值：
`ok` / `disabled` / `package_missing` / `failed`。
`ambient_rna.status` 有 `not_done` / `heuristic_only`。

**不要用一个代理指标冒充校正结果。** 做不了就写 `not_done` 加原因。
`ambient_rna` 的正规做法需要 SoupX（R，需空液滴）或 CellBender（需 GPU）；
GitHub 托管 runner 无 GPU，所以默认不做，并如实记录。

**"失败原因"必须是原因码，不能是布尔量（E-66 第二版 M2，2026-09-26 亲读
artifact 才抓到）。** 一个 `p_adj_available: false` 只能表达"没做校正"，而
"没做"至少有三种互不相干的处境。实测 `scrna-results-68` 的
`trajectory_status.json` 写着 `note: "statsmodels 不可用"`，而同一份 artifact
的 `run_manifest.json` 里 **`versions.statsmodels = 0.15.0`（装着呢）** ——
真因是配置里 `group_key: null`。**一个真实但错误的原因，比"没有原因"更糟**：
下一个人会去查依赖装没装、换镜像、加超时，而真正要改的只有那一行。

所以：处境有几种就写几个**原因码**（`not_configured` / `single_group` /
`no_pairs` / `import_failed` / `ok`），文案做成 `{码: 文案}[码]` 的映射表且
两两不同，布尔量**由码推导**（`p_adj_available = (ks_reason == "ok")`）而不是
与它并列 —— 两个独立量一定会漂开。**"数据里没有分组"与"配置漏了"是两件事。**

**标出来的状态必须真的有人消费（L6 第二半，自查 2026-09-26）。**
`05_trajectory.py` 的 `orient()` 早就把"方向无法判定"标成了
`direction_decided: False` 并打了 WARN，注释还写着"下游按 `direction_decided`
决定是否把这个方法算进共识"—— 而**它一个消费者都没有**：那些方法的拟时序
照样进 `cv_names` → 进一致性统计 → 进共识。`rho_vs_reference` 是 nan
（该方法退化成常数列）或恰为 0 时，`flipped=None` 意味着方向**没被校正**，
共识方向就是随机的，而下游"沿轨迹变化的基因""模块"全建在那个随机方向上
—— 产物里完全看不出来。**注释陈述的行为与代码实际行为对不上，是 E-58②
"看起来在算、其实没在算"的又一种形态。**

两道剔除理由不同、**分开记**（`cv_excluded.direction_reference` =
"定义上相关" / `direction_undecided` = "方向未知"），合并成一个数会让读者
以为剔的是同一类东西。并且**删掉"至少留一个"的静默回退**
（`cv_names if cv_names else names`）—— 全部方向都未判定时，那个回退恰好
把刚剔除的方法**全部放回共识**，而且 `mean_rho` 是 nan、`nan < 0.3` 为假，
连"一致性偏低"那条限制都不会加。现在这种情况判 `no_consensus`
（**已登记进 `STEP_ABORT_VALUES`**，否则会静默归成 `skip`）。

## 5. 每个"结论"都要带上它的适用范围

单细胞分析里，方法学限定不是免责声明，是**结论的适用范围**。
以下字段是必需的：

| 产物 | 必需字段 | 为什么 |
|---|---|---|
| `celltype_annotation.csv` | `score_margin`、`margin_state`、`assignment_confident` | margin 小时 assignment 不该被当结论；`margin_state` 把"不确定"与"算不出来"分开 |
| `trajectory_status.json` | `root_selection`、`limitations` | 拟时序方向完全依赖根的选择 |
| `communication_status.json` | `method`、`limitations` | 共表达不等于通讯 |
| `grn_status.json` | `method`、`limitations` | 没有 motif 剪枝就不是 SCENIC |

**报一个类型名 / 一个 p 值而不报它的不确定性，等于把不确定性藏起来。**

**并且"不确定"与"算不出来"要分开报（E-66 第二版同族，规则 4 同一条）。**
`assignment_confident` 是三态：`true`（`margin_state: ok`）/ `false`
（`low_margin`，两条路真的分不开）/ **`null`（`single_celltype` —— 没有
第二名可比较，或 `margin_undefined` —— 分数含 nan/inf）**。后两种不是
"不确定"，是**"这个指标在这里不适用"** —— 排查方向是补签名基因，
不是去比对两条注释路。把三种处境压成一个 `false`，日志就会说
「N/M 个簇 margin<0.05（assignment 不确定）」而其中可能一个簇的 margin
都没算出来。

## 6. 绘图 API 会随版本变，写防御性代码

实测踩过：
- `sc.pl.highly_variable_genes(return_fig=True)` 在 scanpy 1.12 被移除 → TypeError
- `ax.boxplot(labels=...)` 在 matplotlib 3.9 改名 `tick_labels`，3.11 直接 TypeError

所以：`requirements.txt` 钉上限；绘图尽量用当前版本支持的写法；
版本敏感的调用加回退。

## 7. 不要在循环里往 `adata.obs` 反复插列

`sc.tl.score_genes` 是往 `obs` 插列。循环 200 次会让 DataFrame
严重碎片化 —— pandas 刷一屏 `PerformanceWarning`，而且越来越慢。
用 `.pop()` 即时取走，obs 宽度不增长。

## 8. 置换检验要先把用到的列抽出来

`06_communication.py` 最初在 13714 基因的全矩阵上做 200 次置换 × 673 组合，
跑 4.5 分钟。抽成 66 列的配体/受体子矩阵后是 **16 秒**，结果完全一致。

**优化后必须验证结果一致**，不能只看"跑得更快了"。

## 9. `check_artifact_paths.mjs` 的路径归一化

artifact 路径形如 `data/<dataset_id>/raw.h5ad`，脚本里是
`data_dir / "raw.h5ad"`。两边归一到 **basename** 再比。
第一版只剥了 `data/` 前缀、没处理数据集 id 那一段，全是误报。

`data_dir / "cache"` 是下载/解压缓存目录，不是产物 —— 用"有没有
扩展名"判断跳过。

## 10. 数据源下载必须带 User-Agent

10x 的 CDN（`cf.10xgenomics.com`）对 `Python-urllib/3.x` 直接返回
**403 Forbidden**，而同样的 URL 用浏览器请求是 200。
报错信息只有 `HTTP Error 403: Forbidden`，看不出是 UA 的问题，
**很容易误判成"链接失效了"**。

## 11. CI 用 Python 3.12，不用最新的

`harmonypy` 等包在 Python 3.13/3.14 上还没有 wheel，要从源码编译，
且编译时常因缺 numpy 头文件失败（本地实测过）。
3.12 上这些包都有现成 wheel。

## 12. 每轮 CI 只修日志里明确显示的问题

不要凭猜测改代码。Part 1 有过一次把浮点末位分叉错误归因、白跑一轮的教训。
先读日志，再改，一次只改日志支持的那一处。

## 13. 图幅按毫米，宽度夹在标准栏宽内

参考规范：K-Dense `scientific-visualization` skill（样式文件已 vendored 到
`assets/publication.mplstyle`）。

**期刊栏宽是按毫米规定的**，英寸是排版软件内部单位。写英寸时"这图多宽"
要靠换算才知道，写毫米时一眼能对上投稿要求。

| 常量 | 值 | 用途 |
|---|---|---|
| `W_SINGLE` | 89 mm | 单栏 |
| `W_ONE_HALF` | 136 mm | 一栏半 |
| `W_DOUBLE` | 183 mm | 双栏（通栏）|

- **宽度随类别数增长的图必须夹住**：`min(W_DOUBLE, max(W_ONE_HALF, ...))`。
  无上限增长会画出装不进任何期刊一页的图 —— 实测修之前最宽的
  `domain_markers_dotplot` 是 370 mm。
- `save_fig()` 默认**不再用 tight bbox**。`bbox_inches="tight"` 会**改变物理
  输出尺寸**，让上面的毫米约定失效。溢出改由 `_content_overflow()` 检测，并
  **落盘成 `figure_overflow.json`**（不再是只打 WARN —— 见规则 25）。
- **`set_seed()` 末尾会 `apply_style()`** —— rcParams 在**图创建时**就被读取，
  在 `save_fig()` 里设样式已经太晚。
- 图上文字一律英文：matplotlib 自带的 DejaVu Sans 没有中文字形。

**两处与参考规范的偏差**（已写在 `.mplstyle` 文件头）：

1. `figure.constrained_layout.use: True` 全局开启。参考规范要求逐图 opt-in，
   但本仓库有 20+ 张图、8 个脚本，逐处改容易漏。
2. `font.sans-serif: DejaVu Sans, Arial, Helvetica`。参考规范首选 Arial，
   但 Ubuntu CI 上没有 Arial 而 Windows 上有 —— 会导致 CI 与本地渲染出
   **不同的字形**，破坏可复现性。DejaVu Sans 随 matplotlib 分发，处处一致。

**constrained layout 不会自动折行长标题。** 实测 `domains_on_he` 的单行
suptitle 超出 183 mm 宽 2.9%，被静默裁掉（`savefig.bbox: standard` 下文件照样
生成、`check_figures.mjs` 也照样报"有墨"）。长标题必须自己换行。

## 14. 可选步骤崩溃必须让验收变红

`main_analysis.py` 里可选步骤的 `failed` **不等于**"正确地跳过"。
`not_configured` / `not_applicable` / `disabled` / `not_run` 是设计如此；
`failed` 是崩了。

第一版把 `failed` 也归进 `ok=True`，后果是实测 `07_grn` 因漏 import 崩溃
而验收 40 项全绿 —— 而 `grn_status.json` 是**上一轮的**旧文件。

两条一起改才有效：

1. `failed` 记为不通过（可见但不阻断 job，因为步骤本身是可选的）
2. **每步开跑前先删掉自己的状态文件** —— 否则崩溃时旧文件冒充本轮结果

## 15. 静态检查要挡住"未定义名字"

**规则：** 静态检查必须做**逐作用域分析**，不能靠写死的名字名单 —— `py_compile` 只编译，看不出未定义名字。

**失败形态（E-61）：** 旧 `tools/check_py_syntax.mjs` 名单只有 **10 个名字**，而 `common.py` 导出 **74 个**，没列进去的照旧漏到 CI。实测 `spatial-pipeline-skill/scripts/03_spatial_domains.py` 漏 import `spot_radius_plot_units`，本地报"全部通过"，CI 跑 **25 分钟**到 H&E 叠图段才 `NameError`。**根因不是名单短了，而是判据的输入域与它要防的缺陷不匹配。** 现在 `tools/check_py_names.py` 用标准库 `symtable` 判三件事：未定义名字（引用却无绑定 → 必 `NameError`）、幽灵 import（`from common import X` 而 common 没有）、不安全构造（`import *` / `globals()` / `exec` / `eval` / `vars` / `locals` → **判红退出**）。

**三个坑：** ①不安全构造必须用 AST 判 —— 第一版扫子串时**检查器把自己的源码判红**（需豁免时同行写 `# py-names: unsafe-ok —— <理由>`，故意要写一句话）；②行号必须按作用域定位 —— 第一版把 `plt` 指到**另一个函数**的同名变量，**指向错的行号比不指行号更糟**；③提示要扣掉 common 自己 import 的名字（否则会写出 `from common import np`）。

**两仓该文件逐字节相同**（SHA256 `90D2ECEF11170306406C4FA356FD241189A2DFF5DC912ED7C060B8ADC783990E`，12530 字节），改一侧必须同步并比对哈希。标定：正向零命中；geo（纯 R）**不适用且不判红**；传错目录判红；反向逐类注入各自只响自己那条。

**15.2：** 触发动作是**手工清理死 import** —— 扫描器正确报出 `plot_marker_dotplot`，删除时把同一行相邻的 `spot_radius_plot_units`（**有**调用点）一起删了。**扫描器没错，是删除这一步错了**：删 import 要逐名核对，扫描器输出是**名单**不是**待删行号**。

> **"没有 Python 文件"是"不适用"，不是"失败"。**

台账：governance/15_ERROR_LEDGER.md E-61
完整原文（含全部实测证据与表格）：references/agents-detail.md#原规则-15

## 16. 每轮运行必须留下可追溯的运行清单（模块零）

参考规范：三大部分整合文档的「模块零」（§0.2–§0.4）。

**没有运行清单的分析结果不是结果。** 半年后拿到一份 `cell_communication.csv`，
如果不知道当时装的是哪个版本的 liana、输入的 h5ad 是哪个哈希、随机种子是多少，
那份 CSV 就**无法被复现，也无法被质疑** —— 而不可质疑的结论没有价值。

产物是 `results/<dataset_id>/run_manifest.json`，由 `common.py` 的清单层写出
（`init_manifest` / `capture_versions` / `probe_r_packages` / `record_input` /
`record_params` / `record_decision` / `record_human_review` /
`record_cross_language` / `manifest_summary`）。几条不能省的约定：
未装的工具记 `null`、不能省略键（"查过了没装"与"没查"是两件事）；
人工复核默认 `pending`，**不算失败但必须可见**；`init_manifest` 必须清掉上一轮；
清单在步骤跑完之后才登记输入；`pip freeze` 不用 subprocess 抓（沙箱下 EPERM）；
**R 包必须问 R**（`probe_r_packages()`）—— `importlib` 对 R 包**原理上**永远返回
"没有"，那个 `None` 会被读成"查过了，装不上"（E-41）。

台账：governance/15_ERROR_LEDGER.md 模块零（§0.2–§0.4）/ E-41
完整原文（含九个函数的字段表与三种失败原因）：references/agents-detail.md#原规则-16

## 17. 虚拟敲除 / 过表达是保留框架，重点是"没做什么"（§1.7 / §1.8）

**规则：** §1.7 / §1.8 是**保留框架**，落地在 `08_virtual_perturbation.py`。点名工具"装不了"必须给**实测原因**，不能只写"没装"。

| 工具 | 状态 | 原因 |
|---|---|---|
| `scTenifoldKnk` | ✅ **已接入** | R/CRAN 包（`scripts/lib/tenifold_knk.R` + CI 的 setup-r 层）|
| `PerturbNet` | ❌ 装不了 | PyPI 上 0.0.2/0.0.3 钉 `requires_python='<3.8,>=3.7'`、0.0.3b0/b1 钉 `'<3.11,>=3.10'` —— **没有任何一版支持 CI 的 3.12** |
| `RegVelo` | ❌ 装不了 | 要 spliced/unspliced 层（本数据没有），且拉入 `torch` + `scvi-tools` |

**E-41（探针查不到 ≠ 目标不存在）：** `probe_tools()` 曾对三个工具统一用 `importlib.util.find_spec(name)`，而它**只查 Python import 路径，对 R 包在原理上不可能返回 True** —— "我没装 R"被写成"这个工具不存在"，还被验收项固化成"**必须**声明不是 scTenifoldKnk"。**"不在 PyPI"那句理由是真的** —— 真理由 + 越界结论 = 一条读起来完全合理、实则自我封闭的记录。现在走 `_probe_r_package()`（`shutil.which("Rscript")` + `requireNamespace`），并**区分**"Rscript 不在 PATH"与"R 装了但包没装"。**每加一个 `find_spec`/`which`/`exists` 式判断，都要能回答"目标存在时它返回什么"**。`tools` 逐个记 `{available, reason}`：**"装不了"（有确切原因）和"没装"（疏忽）是两件事**。

**两个引擎并列，一张图只画一个：** `perturbation.engine` = `both`（默认）/ `first_order` / `tenifold`，**拼错直接报错、不静默退回**。两套结果**并列落盘**（`virtual_perturbation.csv` / `virtual_perturbation_tenifold.csv`），一致性**只比排名不比数值**（一阶是 z 分数的 L2 范数，tenifold 是流形欧氏距离，**量纲不同**）；两法共享 `07_grn` 上游，一致性高**不说明哪个对**。出图默认 tenifold（真方法），没跑出来退回一阶并写明；两法不画同一根 x 轴。`scTenifoldKnk` **没有细胞类型分辨率**，本仓没有按细胞类型重加权。

**空敲除（网络表达不了该扰动 ≠ 效应为 0）：** `transcriptomeWide` 敲除 = 把 WT 网络**那一行清零**（**包内源码** `R/scTenifoldKnk.R` L206–L207 `KO <- WT; KO[g, ] <- 0`）。该基因**出度本来就是 0** 时等于什么都没敲，`manifoldAlignment` 对两个相同网络返回的"距离"只剩浮点噪声（实测 **~1e-16**）——**不报错、不给 NA、不给零**，读表人只会得出"敲除无影响"，既有图门禁照过（它们查的是"有没有产出"不是"产出对不对"）。判据是**结构事实（出度）**不是"距离小"：R 侧 `zero_outdegree_targets()` 整行置 NA，meta 带 `empty_knockout_genes` + `target_outdegree`；Python 侧 `compress_tenifold_distances()` 用**那份结构证据**打标记（不从"整行 NA"反推），空敲除行报 `None` 而**不是 `0.0`**，`compare_engines()` **显式剔除并记账**而不是靠 `dropna()`。两侧各有自检：`Rscript scripts/lib/tenifold_knk.R --selftest` 与 `python tools/selftest_tenifold.py`（**不需要 R**，端到端跑 `run_tenifold_engine`，用假 Rscript）—— 抽函数时漏改的 `return` 旧变量名 `py_compile` 查不出，只会在 CI 上 tenifold 真跑通那一刻抛 `NameError`，而那时 R 已跑了 **6.5 分钟**。

**一阶近似：** `Δz_G = −z̄(G,C)`（敲除）/ `+overexpress_sd`；`Δz_Y = w(G,Y)·Δz_G`。`ko_magnitude = |Δz_G| × network_sensitivity` **混了两个因素必须分开报**：实测 `GATA2` magnitude **5.41** 而网络敏感度只有 **0.81**，`TBX21` 敏感度 **1.84**（全场最高）却只排第四；另落一列 `network_sensitivity = √(Σ_Y w²)`（与表达无关）。**不做多跳传播**（相关网络上多跳只放大噪声）。`PerturbNet` 始终没跑，`method` 必须仍否掉它。

**交接：** `perturbation.targets_csv` 指向 Part 1 的基因表，没有它时回退内部调控子并**必须在状态里写明回退**（`target_source` 必填）。Part 2 → Part 3 的 `clustered.h5ad` 契约三条：①`layers['counts']` 是**原始计数**（否则退回 `.X` log 值求参考谱 → **比例没有意义却仍像一组比例**）②`obs` 有细胞类型列（否则 KeyError）③基因集是**交集** —— 第 3 条是本仓实际情况（`clustered.h5ad` 来自 **2000 个高变基因的子集**）。契约**主动记进** `cluster_status.json` 的 `part3_reference`。`07_grn.py` 另写 `tf_regulon_edges.csv`（`tf, target, corr`）：**靶基因名与权重必须成对取**，分开写会在过滤 `corr > 0` 后错位而不报错。

> **名次会变，值不会** —— 当前 magnitude 最高是 `TYMS`（6.05 / 敏感度 0.94），而 5.41 / 0.81 / 1.84 三个值逐位未变；引用前看 `virtual_perturbation_status.json` 的 `top_by_effect`。

台账：governance/15_ERROR_LEDGER.md §1.7 / §1.8（含 E-41 探针教训）
完整原文（含全部实测证据与表格）：references/agents-detail.md#原规则-17

## 18. 笼统的 `except Exception` 会把代码 bug 记成环境问题

`try_celltypist()` 的第一版把 `celltypist.models_path` 当成
`pathlib.Path` 用了：

```python
path = ctm.models_path / info["model"]     # TypeError: str / str
```

而它其实是 **`str`**（celltypist **包内** `models.py:19` 是 `os.path.join(data_path, "models")`）。
这个 `TypeError` 被下面那个笼统的 `except Exception` 接住，于是状态里
写成：

> 模型 `Immune_All_Low.pkl` 拿不到（TypeError: ...）—— **CI 可能无外网**

**一个纯本地代码 / API 版本 bug 被记成了网络问题。** 下一个人会去查
runner 的出网策略、换镜像、加超时 —— 而真正要改的只有那一行。
（实测模型服务器 `celltypist.cog.sanger.ac.uk` 一直是可达的。）

所以**兜底的 `except Exception` 里必须把失败分类**，而不是把
`type(exc).__name__` 原样拼进一句预设的结论。现在拆成三类：

| 失败点 | status | 说明 |
|---|---|---|
| 路径构造 | `failed` | celltypist API 与调用方不匹配，**不是网络问题** |
| 下载抛异常 | `model_unavailable` | 服务器不可达或超时 |
| 下载没抛异常但文件仍不在 | `model_unavailable` | 名字不在模型清单里 |

第 3 类必须单独查：`download_models` 内部把每个模型的下载异常
**吞掉只打日志**（celltypist **包内** `models.py:512-517`），所以"下载失败"并不总是抛出来 ——
只看有没有异常会把"清单里没这个名字"漏成"下载成功"。

**规则：`except Exception` 里不要写结论，写事实。**
`{type(exc).__name__}: {exc}` 是事实；"CI 可能无外网"是猜测。
猜测写进产物就会被当成证据。

## 19. 点名工具跑了，就必须量化它与自建方法的一致性

和空间仓库同一条（那边写在 `spatial-pipeline-skill/AGENTS.md` 规则 20）。

**"跑通了"不是结论。** `main_analysis.py` 里 §2.4 CellTypist 与
§2.7 LIANA 的验收项分两层：

1. `status == "ok"` —— 工具跑成了没有
2. `compared == True` —— **它与自建方法的差异被量化了没有**

第 2 层是容易被省掉的。实测 LIANA 与自建打分的 Spearman rho ≈ 0.13，
而 top25 重叠 24/25 —— **头部一致、中段排序差异大**。
只报"LIANA 跑通了"会把这两件事都藏起来。

工具跑不成时（`model_unavailable` / `package_missing` / `failed`），
验收项**可见但非阻断**（`required=False`）：那是环境问题不是分析错了，
但**绝不能混在绿字里** —— 规则 4 的同一条理由。

## 20. `set_seed()` 管不到"库内部的迭代求解器"

**规则：** 看到"我明明设了种子"时，先问**"这个函数有没有 `seed` 参数"** —— 没有的话设多少遍都没用。`set_seed()` 只覆盖用 numpy/random 的随机调用；**第三方库内部的迭代求解器、并行归约、GPU 内核都不在里面**。

**失败形态：** 同一 Python 3.12、同一批包版本的六轮 CI（`8f3f56c0` / `9af13b42` / `5b241a3` / `db778bb` / `f99eab1` / `e549477`）里，`dpt`（+0.5856）与 `palantir`（+0.4674）**六轮逐位相同**，上游聚类数（10）、分辨率扫描、CytoTRACE 也完全一致，**只有 scFates 一个在变**：ρ = +0.5644 / +0.5296 / +0.5328 / +0.5328 / +0.5328 / +0.5646，平均 ρ +0.6363 / +0.6254 / +0.6265×3 / +0.6353，拓扑从 A 轮的 **4 片段 / 6 milestone** 变成之后六轮的 **6 片段 / 8 milestone**。

**两次归因都被否证：** ①"多线程 BLAS 归约顺序"—— **C 轮否证**（钉住后 +0.5296 → +0.5328，没有收敛）；②"`pynndescent` 的 Numba 并行"—— **D/E/F 三轮否证**（0.5328 / 0.5328 / **0.5646**，最后一个几乎回到未钉时的 +0.5644）。**留下的线索：钉并行度让离散的拓扑稳住了，但连续量 ρ 还在漂，看起来有两个吸引子 ≈+0.5328 与 ≈+0.5645，残留随机源尚未定位。**

**读源码得到的真路径：** `scf.pp.diffusion` → `palantir.utils.run_diffusion_maps` → `compute_kernel(..., backend="scanpy")` → `scanpy.neighbors.Neighbors.compute_neighbors(..., method=None)` → scanpy 默认 `method="umap"`。**`scFates.pp.diffusion`（1.2.5 `preprocessing/diffusion.py`）签名里根本没有 `seed`**，内部调 `run_diffusion_maps` 时也不传，靠 palantir 默认 `seed=0`；`05_trajectory.py` 只把 seed 传给了 `sct.tree()` / `sct.pseudotime()`（这两处是对的），**扩散图那一步没地方传**。特征求解器本身是干净的（`eigs(..., v0=rng.random(...))`）；simpleppt 的 PPT 主曲线树（`ppt.py:210-213`）也 seed 了，只是把上游末位差异放大成可见的 ρ 变化。排查顺序：`grep -n seed <包的源码>` —— 看 `seed` 是**出现在签名里**还是**只出现在 docstring 里**。

**并行度三个独立旋钮**（缺一个都不够，但**不保证**解决）：`OMP_NUM_THREADS`/`OPENBLAS_NUM_THREADS`/`MKL_NUM_THREADS`（BLAS 归约顺序）、`OPENBLAS_CORETYPE=Haswell`（按宿主 CPU 型号分发 SIMD 内核，线程数管不到）、`NUMBA_NUM_THREADS=1`（Numba `prange` 线程数）。**这五个变量现在都设着，ρ 仍然在漂（D/E/F 三轮）** —— 设它们是对的，但不要以为设了就可复现。另外**"版本不同"不是万能借口**：版本相同也能对不上，别一看到数字变了就归因到版本上。

**结论怎么报：** scFates 的 ρ **不要当单一确定值报** —— 实测范围 **+0.5296 ~ +0.5646**，平均 ρ **+0.6254 ~ +0.6363**，`05_trajectory.py` 会把这段写进 `trajectory_status.json` 的 `reproducibility` 字段（**结论的适用范围要跟着产物走，不能只写在 AGENTS 里**）。`dpt` / `palantir` / `cytotrace` 六轮逐位相同，可以按确定值报。

> **钉并行度能让拓扑稳定，但不足以让 ρ 稳定。**

台账：见 governance/15_ERROR_LEDGER.md 中与随机性/不可复现相关的条目
完整原文（含全部实测证据与表格）：references/agents-detail.md#原规则-20

---

## 21. 文档改动走独立的 `docs_check.yml`

`scrna_analysis.yml` 有 `paths:` 过滤（只跑 `scripts/` `tools/` `assets/`
`requirements.txt`），**`*.md` 的改动不触发它** —— 所以文档里的死链接
以前**没有任何门禁能挡住**，CI 每次都是绿的。

实测就踩到了：`references/methods.md` 与 `references/troubleshooting.md`
让读者去看 `pseudotime.py` 第 254 行，**但没说那是哪个包的** ——
它是 scFates 包内的文件，不在本仓库里。

现在 `tools/check_doc_refs.mjs` + `.github/workflows/docs_check.yml` 兜住这一类，
**约 20 秒**、不装依赖、不跑分析。两个逃生舱写在工具文件头：

| 逃生舱 | 标记 | 为什么必须写 |
|---|---|---|
| 指向**已删**的文件 | 已删 / 已移除 / 不再存在 / 曾经 / 当时的 / deleted | 读者要能区分"历史"和"笔误" |
| 指向**第三方包**源码 | 包内 / 上游 / 源码 / site-packages / 该包 | 读者要知道去哪个包里翻 |

**两个逃生舱都故意做成"要写一句话"的** —— 静默豁免会让检查退化成没有检查。

> **为什么文档不并进主 workflow：** 那样改一个错别字要跑完整的单细胞流水线
> （约 25 分钟），而且会被主流水线的偶发失败牵连 —— 文档改动因无关原因判红，
> 反而让"每次推送 CI 必须绿"这条失效。

## 22. 出图名要带代码坐标

参考规范：geo 仓库 AGENTS 规则 30（同一套约定，阶段号换成 `02`）。

**命名格式五个字段，用 `-` 连起来：**

```
<阶段>-<模块>-<图>-unit<单元>-<名称>
 02      03     01   unit1     umap-clusters
```

| 字段 | 取值 | 来源 |
|---|---|---|
| 阶段 | `02` | 单细胞 = 第二部分（geo `01` / spatial `03`）|
| 模块 | 两位数字 | **脚本文件名的前两位**（`03_cluster_annotate.py` → `03`）|
| 图 | 两位数字，从 `01` 起 | 该脚本内第几张图 |
| 单元 | `unit1` 起 | 同一张图里的功能单元 |
| 名称 | 小写连字符 slug | 图的内容 |

**为什么写成门禁而不是靠人记：** 图名和脚本序号是**两处**，而"图名里的
模块号写错了"没有任何东西能发现 —— 图照样生成、CI 照样绿、验收照样过，
只是读者按图名去 `scripts/` 里找代码时会**找错文件**。
`tools/check_fig_names.mjs` 查四条：格式合规、**模块号与脚本文件名一致**、
图号/单元号连续、全仓库无重名。它**只认字符串字面量**；经辅助函数传名的
调用会逐条列出（可见）但不判失败。

**图名改了，引用它的地方也要改**，否则验收项会查一个不存在的文件而静默
变成"永远 false"。实测 `main_analysis.py` 里 `virtual_perturbation_effect.png`
就是改名时漏掉的那个 —— **而 `check_doc_refs.mjs` 管不到这一类**：
它查的是源码/配置文件名，图名是**运行产物**、本地不存在。
改图名时手工 grep 一遍 `\.png` / `\.pdf`。

> 三个仓库的 `check_fig_names.mjs` 是**同一份文件**（geo 那份只多识别
> R 侧 `save_pdf` 的调用写法）。阶段号表写在文件头的 `PART_BY_REPO`，
> 仓库目录名认不出时它直接报错退出，不会静默放行。

## 23. 出图三条补充约定（评审 3.1/3.6/3.8 实测）

**23.1 退化分布必须显式处理。** 恒定值指标（如 pbmc3k 的 `pct_counts_hb`，
几乎全部细胞为 0）画小提琴会塌成一根裸竖线 + 顶端横线，
**看起来像渲染失败**。IQR=0 的面板改画 strip 散点，标题注明
`(no variance: N/M cells non-zero)` —— 少数非零点反而可见。

**23.2 跨图的同一分类变量必须同色同序。** cluster 在 `umap_clusters` 与
`pseudotime_umap` 第 3 面板里用不同调色板（tab20 vs 默认循环）时，
同一簇跨图变色、读者无法对照。簇色统一走**数值序逐簇 scatter +
PAL_CYCLE**（与 `axes.prop_cycle` 同源），需要就加显式图例。

**23.3 原始列名/变量名不能直接当图上标签。** `total_counts`、
`pct_counts_mt` 是数据结构泄露。统一映射成人类可读英文
（"Total counts per cell"、"Mitochondrial fraction (%)"），
映射表写在用它的脚本里（`01_qc.py` 的 `QC_LABELS`），三仓共用语义。

**23.4 行级 topN 做分类轴前先去重。** `reg.head(n)` 是行级排名，
同一 TF 的多行会让热图 x 轴出现重复列（实测 TBX21/TBX21、IRF1/IRF1）。
按 `drop_duplicates(subset=["tf"])` 去重再取 N。

> 散点标注防撞（`tf_specificity_scatter`）用"按 y 排序 + 上下交替偏移"
> 的纯绘图参数法；`adjustText` 不在依赖里，不要临时引入。

## 24. 图例一律图框外右侧、纵向排列（用户约定 v2，2026-09-23）

姊妹项目 `spatial-pipeline-skill/AGENTS.md` 规则 26 与
`geo-normal-pipeline-skill/AGENTS.md` 规则 31 是同一条约定；
本仓库此前**漏写了这条规则**（代码其实已经合规）——
补门禁时才发现文档缺口，所以补在这里。

**图例不能画在图框（panel）里面，也不能放在顶部** —— 框内会压住数据点，
顶部会把主图压扁变形（实测校准图 / UMAP / 去卷积条形图被压得很扁）。

- **matplotlib**：`fig.legend(loc="outside right center", ncol=1)`。
  `loc="outside ..."` **只对 `fig.legend()` 有效**，传给 `ax.legend()` 会报
  `ValueError: 'outside' option ... only works for figure legends`
  （实测 spatial run 35749568552 因此崩了整个 job）。
  所有 axes 级图例必须改成 `fig.legend(...)`。
- **`ncol=1` 强制纵向单列** —— 多图例时默认可能横排，必须显式指定。
- constrained layout 会自动为框外图例让出空间；**改完要亲读**确认没被裁掉
  （`savefig.bbox: standard` 下溢出是静默裁，见规则 13 同类问题）。

**门禁**：`node tools/check_legend_convention.mjs`
（三仓同一份，静态扫源码；CI 在"静态检查（不装依赖）"那一步跑）。

> **为什么需要门禁而不是靠记：** `check_py_syntax.mjs` 只查未定义名字，
> 看不见图例位置。一张图例压在数据点上的图，在
> "文件存在 / 有墨迹 / 图名合规 / 配色合规 / 图幅合规"眼里**全都是合格的** ——
> 这正是工作区治理层错误台账（`governance/15_ERROR_LEDGER.md`，**不在本仓库内**）E-06「门禁本身有盲区」的又一例。

## 25. 四条"图没了 / 图被裁了却全绿"的补强（Q-26，2026-09-25）

**规则：** E-48 与 E-49 两条**互相独立**的缺陷暴露了验收层与门禁层**四个各自独立的盲区**，四条补强各堵一个，缺一条那类缺陷就能再犯一次。

**25.1 嵌套 `status` 的 `failed` 必须判红（验收层）。** `status` 不只在顶层 —— `qc_status.json` / `integration_status.json` / `cluster_status.json` / `grn_status.json` 里都有**嵌套** `status`，而旧验收层只看顶层 `d.get("status")`。实测真 artifact `scrna-results-55`：顶层 `{ok: 7, not_configured: 1}`，**嵌套 `{ok: 5, failed: 1, not_applied: 1, not_done: 2}`** —— 唯一那条 `failed` 是 `grn_status.json → regulon_vs_pseudotime.status`，而它顶层是 `ok`：**顶层全绿、里面已经崩了**。判据**只把 `failed` / `error` / `fail` 判红**（`not_applied` / `not_done` 可见不阻断）。实现 `_iter_nested_status(obj, path=())`，扫**全部** `*status.json`（不只可选步骤那 5 个）。

**25.2 图验收必须条件化且覆盖全部出图脚本（验收层）。** `REQUIRED_FIGURES` 原来只有 **14 条**，**不含 `02-07-*`（三张 GRN 图）与 `02-08-01`** —— E-48 让 5 张图从未产出而验收层不知道它们该存在。补齐后还有第二个问题：`trajectory.enabled: false` 时 `02-05-*` 一张都不该有。改成由图名第 2 段反查所属步骤、从 `read_state(cfg)` 取 status 与 required：可选且未跑成 → `required: False`，否则 → `required: True`。**兜底是 25.1 那条独立扫描。**

**25.3 溢出检测必须落盘，不能只打 WARN（产物级）。** E-49 的形态：`_content_overflow()` **正确检测到** `width_overflow_frac=0.3523`、**正确打了 WARN**，然后**没有任何人读** —— 标题超宽 **35%**，`savefig.bbox: standard` 下被静默裁掉，图照样生成、门禁照样绿。**"检测到了"不等于"有人会知道"。** 现在 `_record_figure_overflow(cfg, name, bad)` 落盘到 `results/<dataset_id>/figure_overflow.json`，验收层读它并**判红**，detail 写明"**会被静默裁掉**"。三条约束：①跑前必须删旧文件 ②**只累积、不覆盖** ③记录可被修复。

**25.4 门禁要查"内容贴边"，不能只查"有没有墨"（门禁层）。** `check_figures.mjs` 原来只查"有没有墨 / 是不是糊死"——**被裁掉的图照样有墨**，这正是它漏掉 E-49 的原因。裁切的物理后果是**墨迹延伸到画布边缘**，量非背景像素外接框到左右边的距离，阈值 `EDGE_MIN_PX = 3`（标定：scrna 34 张 + spatial 76 张，三种候选阈值**都只命中 E-49 那一张**，零误伤）。两条刻意取舍：**只判左右不判上下**、**暗底图跳过**（整幅都是墨，判了全是假阳性）。

**25.5 上线首跑就抓到一条新缺陷。** 补强推送后第一次 CI（run `36140048953`，commit `e2d7a2c`）**验收判红**：`02-07-01-unit5-tf-smarca4-trend` `width_overflow_frac = 0.0371` —— `scripts/07_grn.py` 单 TF 面板标题是**单行**而面板宽只有 `W_SINGLE = 89 mm`（本地标定 `SMARCA4 +0.0371` **与 CI 逐位吻合**）。这批面板此前**从未渲染过**（Q-24 的三元素 `figsize` 让整段在 `plt.subplots` 就抛异常）—— **旧缺陷一直在掩盖新缺陷**。修法是**折行**（第二行宽度**与 TF 名无关**，固定 2.4% 余量）。

> **宽度受限的小面板，标题折行，不要删字；修好一个缺陷后要重跑一遍全部检查。**

台账：governance/15_ERROR_LEDGER.md Q-26（含 E-48 / E-49 / E-51，均属 `savefig.bbox: standard` 下静默裁切同族）
完整原文（含全部实测证据与表格）：references/agents-detail.md#原规则-25

## 26. native 崩溃绕过 `except`：`sc.pp.scrublet` 偶发 SIGSEGV（E-52，2026-09-25）

**规则：** native 崩溃不能靠 `except` 兜底 —— 凡有 C 扩展参与的关键步骤，先问"它要是崩了，我连日志都没有吗"。**给 native 崩溃留取证路径，比急着改代码值钱。**

**失败形态：** CI 偶发 **exit 139**，日志停在 `01_qc.py:226` 之后一行 `Segmentation fault (core dumped)`，**没有任何 Python traceback** —— `except Exception` 接不住 native SIGSEGV，`qc_status.json` 里连失败原因都没有。**全量 60 个 run 里 8 个 `run_attempt>1`，全部是"att1 失败 + rerun 成功"**，崩溃位置逐字相同，跨 5 天、跨 7 个 commit —— **"rerun 能过"把它掩盖了整整 5 天。**

**根因（`PYTHONFAULTHANDLER: 1` 的真栈定下来，run `36150495910`；此前只能靠猜）：** OpenBLAS 0.3.34 的 `dgemm_kernel_HASWELL` 栈越界。栈：`scipy/.../_interface.py:1118 in _matmat` → `_svds.py:511 in svds` → `sklearn/.../_pca.py:440 in fit` → `scanpy/.../_scrublet/pipeline.py:84 in pca` → `01_qc.py:79 in run_scrublet`，**一帧 numba 都没有**。上游 **OpenBLAS #6026**：Haswell 内核把 k 维打包进固定 **`0x7080`** 字节栈缓冲区、每步写 **96 字节**且**对 k 无上界**，**k ≳ 319** 时越界；**偶发是因为 level-3 blocking 按宿主 cache 拓扑选取**。**附带风险：越界若没碰到返回地址，进程正常返回而数值被污染 —— "没崩"不等于"算对了"。**

**处置：** 钉 `numpy>=2.1,<2.5.2`（**上游修复 0.3.35 尚未进入任何 numpy wheel**，只能钉、不能等）。复验 run `36153515064`：attempt 1 一次过，artifact `versions.numpy = 2.5.1` —— 两处对上。

> **"rerun 能过"不等于"没有问题"（重试率是可观测指标）；"本地复现不了"不等于"没问题"（k 按宿主 cache 拓扑选，设计如此）。不钉上界 = 把稳定性交给上游的发布节奏。**

台账：governance/15_ERROR_LEDGER.md E-52（任务行 governance/02_TASKLIST.md V-02）
完整原文（含全部实测证据与表格）：references/agents-detail.md#原规则-26

## 27. 步骤"没做成"必须传出来：返回值要接住、顶层 `status` 也要扫（E-56，2026-09-26）

**规则：** 步骤"没做成"必须传出来 —— **返回值要接住、顶层 `status` 也要扫**。凡"让失败可见"的修法，必须覆盖失败的所有承载位置：返回值、顶层字段、嵌套字段、日志、退出码。

**失败形态（E-56）：** `main_analysis.py` 的嵌套失败扫描写成 `if p and str(v).lower() in NESTED_FAILED_VALUES`，而 `_iter_nested_status(d)` 在**顶层**调用时 `path` 是**空元组**（falsy）—— **顶层 `status` 一个都扫不到**。本仓有 6 处以"写顶层非 ok 状态 + 正常 `return`"的方式失败：`scripts/05_trajectory.py:275-279` `bad_root`、`:313-317` `missing_clusters`、`:363-370` `insufficient_methods`、`scripts/07_grn.py:352-355` `no_tfs_in_data`、`scripts/08_virtual_perturbation.py:901-911` `no_candidates`、`scripts/04_pseudobulk_de.py:178-180` `no_celltype_testable` —— 结果这一步什么也没产出，而 `state.json` 记 `ok`、`acceptance.json` 记绿；且这四个词连 `NESTED_FAILED_VALUES = ("failed", "error", "fail")` 都不在。**第二处同病根：** `run_all` 把 `fn(cfg)` 的返回值**直接丢弃**、写死 `record_step(cfg, sid, "ok", ...)`，而步骤的失败路径是 `write_json(状态文件, status)` + `log_warn(...)` + `return status`（**不抛异常** → `except` 分支根本不进）；同一病根还在 **9 个步骤脚本各自的 `__main__` 块**里 —— "编排器跑"和"单独跑"两条路径**都**把失败记成成功。**这是 E-48 的另一半**：E-48 只覆盖"写进嵌套字段的失败"，**字段写在没人扫的位置**与**结果压根没被记下来**这两个变体照旧漏。

**处置：** `scripts/lib/common.py` 分两个**分开的**取值集合 —— `STEP_ABORT_VALUES`（判红，含通用 `failed`/`error`/`fail` + 上表六处 + 07 的 `no_regulons`、06 的 `no_pairs_in_data`/`no_signal`、04 的 `no_usable_pseudobulk`/`column_missing`）与 `STEP_SKIP_VALUES`（可见不阻断）；`result_status_of(res) -> str`（dict 带 `status` 键则取之，否则 `"ok"`）；`classify_step_result(v) -> str` 给 `"ok"`/`"abort"`/`"skip"`（**未登记取值归 `"skip"`**，只把确知失败的判红）；`record_step(..., result_status: str = None)` 让 `entry["result_status"]`（步骤自述"做成了没有"）与 `entry["status"]`（编排器记"有没有崩"）**分开存** —— `status="ok"` + `result_status="bad_root"` 正是要显形的组合。`main_analysis.py`：`res = fn(cfg)` → `record_step(..., result_status=result_status_of(res))`，`except` 分支加 `result_status="failed"`，嵌套扫描去掉 `if p` 改 `key = ".".join(p) + ".status" if p else "status"`（照抄 `spatial-pipeline-skill/scripts/main_analysis.py:380-398`，那边本来就是对的）；9 个步骤脚本的 `__main__` 块同步改造。

> **空元组/空字符串/空列表都是 falsy —— 凡"空值有语义"（顶层路径、根节点、默认分支）处，判据必须显式写 `if p else` 而不是 `if p`。顶层与嵌套是两个语义域，取值域重叠不等于可以共用判据（`tenifold.status` 可为 `timeout`/`no_candidates` 而内置引擎照样出结果，合并会让每个 job 都红，常年假红会训练人忽略告警）。标定必须用本次 CI 的 artifact，不能用工作区里可能陈旧的副本（`figure_review/scrna-pbmc3k` 内含 E-48 遗留的 `grn_status.json → regulon_vs_pseudotime.status = "failed"`，用它做正向标定会把真判据误判成误报）。**

台账：governance/15_ERROR_LEDGER.md E-56（任务行 governance/02_TASKLIST.md R-03；标定脚本 `D:\tmp\_s1\calib_e56.py`：正向 + 反向 19 类注入 + 回退版对照）
完整原文（含全部实测证据与表格）：references/agents-detail.md#原规则-27

## 28. 九条"看起来在算、其实没在算"的缺陷（E-58，2026-09-26）

**规则：** 九条"看起来在算、其实没在算"的缺陷，共同点是**代码在跑、产物齐全、状态全绿，但那个量取不到它声称要取的信息**。复核审计严重项 S2–S10 时 **7 条读码确认成立、1 条（S8）被实测否证、1 条（S3）根因比审计写的更具体**。

**28.1 "恒为 0 / 恒为真"的量比没有这个量更糟。** `scripts/03_cluster_annotate.py` 的 `compare_annotations()` 原来做 `str(own[c]).lower() == str(maj[c]).lower()`：marker 词表来自 `assets/celltype_markers.yml`（`T_cell`/`Platelet`），CellTypist 词表来自模型（`Tcm/Naive helper T cells`/`Megakaryocytes/platelets`）—— **两套词表没有任何一个字符串相等**，于是 `n_agree_exact: 0 / agreement_frac: 0.0` 恒成立；而逐簇看 `B_cell`↔`B cells`、`Platelet`↔`Megakaryocytes/platelets`、`Monocyte`↔`Classical monocytes` 明显一致。**一个恒为 0 的量看起来像一个结论（"两条路完全不一致"），实际只是词表不相交。** 凡两套词表/坐标系/口径要对齐处，先跑一遍看它**能不能取到非平凡值**。处置：新增 `assets/celltype_mapping.yml` + `load_celltype_mapping()`；没映射的簇**不进分母**、单列 `unmapped`；字段改名带 `mapped`（`n_mapped`/`n_unmapped`/`n_agree_mapped`/`agreement_frac_mapped`）并**删除旧字段名**；消费侧分三层判（没对比 → 红；有对比但 `n_mapped == 0` → 红；有映射 → 报真一致率）。

**28.2 注释写着正确做法、下一行做了相反的事（三处）。** `scripts/01_qc.py:78` 注释"用细胞数反推期望双细胞率"，代码 `clip(5000/n*0.01, 0.05, 0.10)` 在 n>1000 时**恒被下界截断到 0.05**（那个行为从不发生）；`scripts/00_fetch.py:147` 注释"判断整个 `X` 是不是计数"，代码 `X[:min(200, X.shape[0]), :]` 只抽前 200 行；`scripts/04_pseudobulk_de.py` 注释"必须用 counts layer"，代码 `layers["counts"] if "counts" in layers else adata.X` 静默退回 log 值。后果**产物里完全看不见**：scrublet 阈值偏保守（pbmc3k n=2652 实测新式 **0.02122** vs 旧式 **0.05000**，**高估一倍以上**）；计数校验在按样本拼接时可能给出 `is_counts: True` 而整体不是计数 —— 而这是下游全部方法学的前提；喂 log 值给 DESeq2 不报错，只给错的离散度估计。写注释时问：**"这行代码真的会走到我说的那条路吗？"** 处置：S9 改 10x 经验式 `clip(0.008*n/1000.0, 0.01, 0.10)` 并把**旧公式的值一并记进 status**（`expected_doublet_rate_old_rule`）；S10 改为**抽两批且抽法必须不同**（等距抽样作主判据 + 中段连续抽样交叉核对），两批占比差 >1 个百分点即判 `sampling_consistent: False`；S6 无 counts 层时返回 `counts_source: "missing"`，调用侧判成独立状态 **`missing_counts`**。

**28.3 检验家庭不能被任何"预过滤"缩小。** `scripts/06_communication.py:321` 原来 `if obs_score <= 0: continue` —— 零分组合根本不进 `rows`，`multipletests` 的分母成了"打分 > 0 的组合"。实测 pbmc3k 分母 **649 vs 1134**，**检验家庭缩小约 43%**，`p_adj_bh` 系统性偏小、`n_significant_bh` 上偏 —— 方向与紧邻注释担心的"校正组合数多导致假阳性"**相反**。处置：零分组合仍进 `rows`、`p_value` 记 **1.0**（零分观测本就无法被任何置换超越），另加 `tested` 布尔列。**分母必须是"被评估过的全部假设"。**

**28.4 方向/符号不能靠数据行序决定。** `scripts/04_pseudobulk_de.py:169` 原来 `str(sub["group"].unique()[0])` 当分子 —— 取的是**第一次出现的取值**，而 `sub` 顺序来自 `meta.groupby(...)` 即 obs 行序：**重排细胞顺序会让全部 `log2FoldChange` 变号，而 `padj` 一个都不变**（方向翻转符号对称），产物看起来完全正常。同类：`05_trajectory.py:439-444` 的共识拟时序把 `cytotrace` 算进去了，而 `direction_source == "cytotrace_fallback"` 时 `cytotrace` 正是方向参考（`:417`）—— 同文件 `:414-425` 算一致性时**正确排除了参考**，算共识却没有，**参考方法既定了方向又参与共识，等于自己给自己投票**；产物证据 `trajectory_direction.csv` 四个方法原始 rho **全为负、全部被翻转**。处置：S5 分组水平按 `sorted()` **字典序**定死，分子/分母拆成显式变量，实际方向写进 `status["contrast_used"]` + `contrast_rule`；S7 共识改用 `consensus_names = cv_names if cv_names else names`，并**同时算含参考方法的共识**、把两者 rho 写进 `consensus_vs_with_reference_rho`。

**28.5 审计报告也会错：落修法前先证伪。** 审计称 `scripts/07_grn.py:164` 的 `corr[top[:len(targets)]].mean()` 与 `targets` 错位，**实测证明不会**：`top = np.argsort(corr)[::-1][:N_TARGETS]` 按 corr **降序**，`corr[k] > 0` 过滤必然保留 `top` 的一个**前缀**，`top[:len(targets)]` 与 `[k for k in top if corr[k] > 0]` **恒等** —— 两万次随机对拍（含 `-inf` 自排除、含并列值、含三种尺度）**零次不等（`0/19944`）**。仍改成 `np.mean([w for _, w in pairs])`（等价更钝，不依赖"降序 → 前缀"），并在源码注释写明**审计前提不成立**。**"读码推断"不等于"跑一遍对拍"。**

**28.6 marker 基因宇宙错位（S2）。** `scripts/03_cluster_annotate.py:71` 用 `adata.var_names` 过滤签名基因，而**同一函数下一段用 `use_raw=True` 打分** —— `adata` 只剩 2000 HVG，CD3D/CD8A/CD14 这类**真实存在、只是没进 HVG** 的 marker 被判"缺失"，签名被削到无法区分（实测簇 1/簇 3 的 `T_cell` 与 `CD4_T` 分数逐位相同、margin 恰为 **0.0**）；**同文件 `:385` 的 dotplot 已经用的是 `adata.raw.var_names`**。处置：改用 `raw_names = set(adata.raw.var_names) if adata.raw is not None else hvg_names`，诊断拆成 `missing_markers`（真的没有）vs `not_in_hvg`（有、没进 HVG）—— 合并会让"签名被 HVG 削弱"和"这批数据没测到"看起来一样。

> **方向/符号来自显式排序，不来自数据行序（`unique()[0]`、`head(1)`、`groupby` 首元素都不行），且要把实际用到的方向写进产物。**

台账：governance/15_ERROR_LEDGER.md E-58（任务行 governance/02_TASKLIST.md R-03；标定脚本 `D:\tmp\_s2\calib_s2_s10.py`：每条配"回退版必须抓不到"对照，34 项全过）
完整原文（含全部实测证据与表格）：references/agents-detail.md#原规则-28

## 29. 门禁要带内建自检，自检必须被反向标定，且必须接进 CI 与 pre-push（E-62 / E-63，2026-09-26）

姊妹项目 `spatial-pipeline-skill/AGENTS.md` 规则 30 是同一批经验的另一半。

**没有自检的门禁只能证明"它没报错"，不能证明"它检查了"。** 两条实测：

| 台账 | 门禁 | 缺陷形态 |
|---|---|---|
| E-62 | `tools/check_legend_convention.mjs` | Python 侧没抹注释 → 注释里一个 `fig.legend(` 让括号配平**一路吞到文件尾**，其后所有真调用一个都没查，门禁照样打绿 |
| E-63 | `tools/check_figures.mjs` | WARN 落盘分支从落地起**一次都没执行过**，里面有两个必崩的错（`INK_FAIL_MIN` 未定义、报告路径用了循环变量 `dir`）|

两条的共同点：**假阴性** —— 门禁的失败方式不是"报错"，而是"什么都不报"，而"什么都没发现"与"检查通过了"在输出上完全一样。**假阴性比假阳性危险得多**：假阳性会被人骂着修掉，假阴性会被当成绿。

**29.1 自检要调真代码，不能自己重写一遍逻辑。** E-63 的自检第一版有 9 个用例，**用例 8 自己另写了一遍路径拼接**；反向标定把缺陷注回去（`outPath` 改回 `join(dir, ...)`）→ **自检仍然通过（exit=0）**。**抽函数**才解决：`checkDirs()` / `buildWarnReport()` / `writeWarnReport()` 三个纯函数，`main()` 与自检**都调它们**；抽完再注一次缺陷 → `ReferenceError: dir is not defined`、exit=1。**自检里重实现一遍被测逻辑，等于没测。**

**29.2 反向标定：逐个把原缺陷注回去，确认自检真的会红。** **正向通过证明不了任何事** —— 一个永远返回 True 的用例在干净产物上也是绿的。三处注入中 `outPath` 那条第一版**不红**（假自检），抽函数后才红；**尾斜杠那条是修完才发现的第三个缺陷**（`dirname("a/b/figures/")` 把最后一段当文件名，报告落进 `figures/` 里与图混在一起）。

**29.3 自检必须接进 CI 与 pre-push —— 没人跑的自检是同一类缺陷。** 现在三处都接：两个 `*_analysis.yml` 的「静态检查（不装依赖）」那一步、以及 `governance/hooks/pre-push.mjs` 的 `SELFTESTS` 表（日志标签 `（自检）` —— **必须区分"带镜像目录"与"带 `--selftest`"**）。

**29.4 判据的"通过数"必须能看见 0（E-64）。** `check_doc_refs.mjs` 补判据 C（跨仓引用）时，主体写在 `if (!isA && !isB) continue` **之后** —— 跨仓 token 正是"两条判据都不进"的那一类，于是**判据 C 的分支永远走不到**；报告那行又是 `if (nCheckedC)` 守卫的，计数恒 0 时**连打印都不打印**。**①新判据要放在早退分支之前**（死代码不报错、只是永远不执行）；**②计数行不能加 `if (n)` 守卫** —— **0 也是信息**。

**29.5 门禁的"检查范围"要和"它守护的动作"对齐（E-65）。** `governance/hooks/pre-push.mjs` 的 `changesOf(repo)` 原来只读 `git status --porcelain`（**只含未提交改动**），而 pre-push 唯一被调用的时刻就是"已经提交、还没推送" —— 那时 porcelain 为空，三仓全走 `无改动，跳过`，**[4/6] 静态门禁段一条都没跑**，而打印的是 `PRE-PUSH 通过（0 条提醒）。可以 push。` **这不是边角，是主路径。** 修法是取并集：再加 `git log --name-only --pretty=format: @{upstream}..HEAD`。

> **一个门禁段被整段跳过时，不能打印"通过"** —— `无改动，跳过` 用的是 `ok()`（绿勾），它和"查过了没问题"在输出上一样。**正确的提交顺序是：改完 → 跑 pre-push → 提交 → push。**
> **接了线但从不失败的检查，与没接线是一样的。** 每加一条自检，都要问"我怎样让它红一次"。

**两仓 `tools/check_figures.mjs` 与 `tools/check_legend_convention.mjs` 各自必须逐字节相同**（前者只在本仓与 spatial 仓之间，后者三仓同一份）。改一侧必须同步并比对 SHA256。

## 30. 「写出来了」不等于「有人读」：三种形态与三个守卫（E-69，2026-09-26）

**规则：** 一个字段 / 一条判据的价值不在于它被算出来，而在于**有人消费它**；而"算不出来"和"算出来很小"必须在产物里**长得不一样**。姊妹项目 `spatial-pipeline-skill/AGENTS.md` 规则 31 是同一次收口的另一半。

- **Form A：`nan` 参与比较会静默变成 `False`。** `nan > x` / `nan < x` 都是 `False` 且不报错 —— "**算不出来**"被当成"**算出来很小**"。本仓现场（已修）：`05_trajectory.py` 的 `off` 在 `cv_names` 只有 1 个方法时为空 → `mean_rho` / `min_rho` 记 `nan` → `mean_rho < 0.3` 那条限制**静默不触发**、日志打出 `+nan`、JSON 里写 `NaN`（非法字面量）。修法是**抽纯函数** `method_correlation_stats(cv_names, cmat) -> (mean_rho, min_rho, state, note)`，`state ∈ {ok, single_method, undefined}` —— **"只有一个方法"和"相关全是 nan"是两种处境、排查方向不同，不能压成一个 `nan`**；再配 `common.finite_round(x, n)`（nan / inf / None → `None`，**`0.0` 要保留、判空用 `is not None`**，不能用真值判断）。**抽出来才能被标定脚本直接调**（内联只能靠跑整条流水线，25 分钟一轮）。产出端落 `method_correlation_state` + 原因文案，验收层 `main_analysis.py` 拿它参与判红 —— **只把状态打进 detail 是装饰，判据本身必须由它参与**。
- **Form B：只写不读的状态字段。** 字段写进 JSON 而验收层没有消费者时，**坏值与"字段不存在"长得一模一样**。spatial 侧扫出 15 个、补了 12 条探针（`opt=False` = 无条件写、缺失判红「**产生端不再写了**」；`opt=True` = 条件写、缺失 PASS；`fn` 返回 `None` ⇒ 不适用）。本仓同族：`manifest_summary` 的 `n_human_review` / `human_review_confirmed`。
- **Form C：恒真判据。** `manifest:human_review` 的 `ok` 曾**写死 `True`** —— "**一个都没登记**"被写成"**全部已确认**"。修法：`_n_hr > 0` 才可能为真；默认 `pending` **不算失败但必须可见**。

**三个实现坑：** ①**消费端标定看不见产生端** —— 补的判据必须是**源码级**的；②`_code_only`（剥 COMMENT+STRING）**不能查字典键名**（键名本身就是 STRING token），要用只剥 COMMENT 的 `_no_comment`；③`tokenize` 的 token **不自带分隔符** —— `"".join(t.string)` 得到的既是 `"key":value`（所以查相邻 token 要按**无空格**形式写），也会把 `is not None` 拼成 `isnotNone`（所以查多 token 表达式**必然匹配不上**）。要查后者就按 token 的 `start`/`end` 坐标**补回原有空白**再匹配。

> **正向标定（干净产物）看不出这三形；反向标定必须逐个把原缺陷注回去，且判红必须伴随非空 FAIL 摘要。** 本仓这次是正向 37 项 + 反向 6 类注入（调用点退回裸 `nan` / 去掉 `is not None` 守卫 / 产出端退回裸 `round` / 去掉状态字段 / `finite_round` 不再挡 nan / 消费端不再用状态判红）—— **第一轮有 3 类注回去却全绿**，说明那 3 条判据当时压根不存在，补上才 6/6。

台账：governance/15_ERROR_LEDGER.md E-69（任务行 governance/02_TASKLIST.md R-03）
完整原文（含全部实测证据与表格）：references/agents-detail.md#原规则-30
