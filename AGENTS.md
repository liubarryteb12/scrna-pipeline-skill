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

## 5. 每个"结论"都要带上它的适用范围

单细胞分析里，方法学限定不是免责声明，是**结论的适用范围**。
以下字段是必需的：

| 产物 | 必需字段 | 为什么 |
|---|---|---|
| `celltype_annotation.csv` | `score_margin`、`assignment_confident` | margin 小时 assignment 不该被当结论 |
| `trajectory_status.json` | `root_selection`、`limitations` | 拟时序方向完全依赖根的选择 |
| `communication_status.json` | `method`、`limitations` | 共表达不等于通讯 |
| `grn_status.json` | `method`、`limitations` | 没有 motif 剪枝就不是 SCENIC |

**报一个类型名 / 一个 p 值而不报它的不确定性，等于把不确定性藏起来。**

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

`py_compile` **只做编译，看不出未定义名字**。漏 import 一个 `W_SINGLE`
时它照样报"语法通过"，要等运行时才炸 —— 实测因此白跑一整轮流水线。

`tools/check_py_syntax.mjs` 现在会检查：用到的 `W_SINGLE` / `W_ONE_HALF` /
`W_DOUBLE` / `mm` / `PAL` / `PAL_CYCLE` / `apply_style` 是否都 import 了。

名单是**写死的**，不是"common 导出的所有名字"。后者会把函数参数名当成用法
（`def verify_alignment(adata, log_info=None)` 的 `log_info`），要正确处理
得做作用域分析 —— 那是重写一个 linter。同时会剥掉注释和字符串再扫，
避免"名字只出现在注释里"的误报。

## 16. 每轮运行必须留下可追溯的运行清单（模块零）

参考规范：三大部分整合文档的「模块零」（§0.2–§0.4）。

**没有运行清单的分析结果不是结果。** 半年后拿到一份
`cell_communication.csv`，如果不知道当时装的是哪个版本的 liana、
输入的 h5ad 是哪个哈希、随机种子是多少，那份 CSV 就**无法被复现，
也无法被质疑** —— 而不可质疑的结论没有价值。

`common.py` 提供清单层，产物是 `results/<dataset_id>/run_manifest.json`：

| 函数 | 记录什么 | 规范条款 |
|---|---|---|
| `init_manifest` | 开新一轮（**清掉上一轮**） | — |
| `capture_versions` | 全量已装 Python 包 + `KEY_PACKAGES` + `R_KEY_PACKAGES` 逐个 | §0.3 |
| `probe_r_packages` | 起一次 `Rscript` 问 R 这些包装了没有 | §0.3 |
| `record_input` | 输入文件的 sha256 | §0.4 |
| `record_params` | 全部参数**含 seed** | §0.3 |
| `record_decision` | Agent 决策链（问题/结论/证据） | §0.4 |
| `record_human_review` | 人工复核节点及状态 | §0.4 |
| `record_cross_language` | 跨语言转换前后维度与**丢失字段** | §0.2 |
| `manifest_summary` | 供验收用的摘要 | — |

**几条不能省的约定：**

1. **未装的工具要记成 `null`，不能省略键。** `key_versions` 里
   `"liana": null` 和"没有 liana 这个键"是两件事：前者是"查过了，没装"，
   后者是"没查"。省略会让读者分不清。
   **但"查过了"的前提是问对了地方** —— R 包走 `R_KEY_PACKAGES` +
   `probe_r_packages()`，见第 6 条。
2. **人工复核未确认不算失败。** `human_review` 默认就是 `pending`，
   判成 FAIL 会让每个 job 都红，反而没人看。但必须**可见**。
3. **`init_manifest` 必须清掉上一轮。** 和规则 14 同一个道理：
   上一轮的清单留在那里冒充本轮，比没有清单更糟。
4. **清单在步骤跑完之后才登记输入。** 可选步骤这轮有没有产物，
   跑完才知道；在开头登记会把"上轮残留"记成本轮输入。
5. **`pip freeze` 不用 subprocess 抓。** 沙箱下管道捕获会 EPERM，
   用 `importlib.metadata` 枚举。**这条只对 Python 包成立** —— 见第 6 条。
6. **R 包必须问 R，`importlib` 永远问不出来。** `importlib.metadata` 与
   `importlib.util.find_spec` 查的都是 Python 的发行版数据库 / 模块查找器，
   对 R 包**原理上**永远返回"没有"。那个 `None` 会被读成"查过了，装不上"，
   而它真正的含义是"**问错了地方**"。

   实测的后果：K-01b 起 CI 装了 R、`scTenifoldKnk` 真跑了 6 分钟，而
   `run_manifest.json`（**复现依据**）仍记 `scTenifoldKnk: null` ——
   于是它与 `virtual_perturbation_status.json` 的
   `tools.scTenifoldKnk.available = true` 对**同一个工具**给出相反结论。

   所以拆成两条通道，且**失败原因分开写**：

   | 通道 | 实现 | 清单字段 |
   |---|---|---|
   | Python 包 | `importlib.metadata` 枚举（不起子进程） | `key_versions` |
   | R 包 | `probe_r_packages()`：一次 `Rscript -e requireNamespace` | `key_versions` + `r` |

   `r` 那段（`{rscript, r_version, packages, reason}`）是必要的：
   `key_versions` 是**平铺**的 name→version，读不出"哪些是 R 包、
   R 是什么版本、R 到底有没有装"。

   三种失败各有独立 `reason`，**不能混成一句**：包名非法（拒绝把任意
   字符串拼进 `Rscript -e`）/ `Rscript` 不在 PATH（本机没 CI 没装 R）/
   退出码非 0。"没装 R"与"R 装了但包没装"必须区分开，
   否则排查方向会被带偏。

   `08_virtual_perturbation.py` 的 `_probe_r_package()` **转调**
   `probe_r_packages()`，不再自己拼一遍 —— 同一件事只留一份代码
   （抄两份时验证的往往只是副本，见 E-41 与规则 15 的教训）。

`run_manifest.json` 落在 `results/` 下，随 artifact 一起上传 ——
**它必须和结果同时可及**，否则追溯链是断的。

## 17. 虚拟敲除 / 过表达是保留框架，重点是"没做什么"（§1.7 / §1.8）

参考规范：三大部分整合文档 §1.7 / §1.8。规范把这两节标为**保留框架**、
主语言 Python，候选靶基因由 Part 1 产出。落地在 `08_virtual_perturbation.py`。

**规范点名的三个工具，两个装不了、一个已经真跑了（理由都是实测的，不是"没装"）：**

| 工具 | 状态 | 原因 |
|---|---|---|
| `scTenifoldKnk` | ✅ **已接入** | R/CRAN 包（`scripts/lib/tenifold_knk.R` + CI 的 setup-r 层）。**"不在 PyPI"曾经被读成"用不了"** —— 见下 |
| `PerturbNet` | ❌ 装不了 | PyPI 上 0.0.2/0.0.3 钉 `requires_python='<3.8,>=3.7'`，0.0.3b0/b1 钉 `'<3.11,>=3.10'` —— **没有任何一版支持 CI 的 3.12** |
| `RegVelo` | ❌ 装不了 | 能装，但要 RNA velocity 的 spliced/unspliced 层（本数据没有），且拉入 `torch` + `scvi-tools` |

> **E-41 的教训（探针查不到 ≠ 目标不存在）。** 早先 `probe_tools()` 对三个工具
> 统一用 `importlib.util.find_spec(name)`，而 `find_spec` **只查 Python import
> 路径，对 R 包在原理上不可能返回 True**。于是"我没装 R"被写成"这个工具不存在"，
> 还被一条验收项固化成"**必须**声明不是 scTenifoldKnk"。**理由那句"不在 PyPI"
> 是真的** —— 真理由 + 越界结论 = 一条读起来完全合理、实则自我封闭的记录。
> 现在 `scTenifoldKnk` 走 `_probe_r_package()`（`shutil.which("Rscript")` +
> `requireNamespace`），并**区分**"Rscript 不在 PATH（CI 没装 R）"与
> "R 装了但包没装"。**每加一个 `find_spec`/`which`/`exists` 式判断，都要能
> 回答"目标存在时它返回什么"**；答不上来就说明这条判断测的是探针自己。

**"装不了"和"没装"是两件事。** 前者有确切原因、要写进状态文件；
后者是疏忽。所以 `tools` 里逐个记 `{available, reason}`，
验收项的判据是**每个不可用的都有原因**，而不是"全都不可用" ——
将来某个工具能装了，这条应该自动变 PASS 而不是 FAIL。

### 两个引擎并列，一张图只画一个

`perturbation.engine` 取 `both`（默认）/ `first_order` / `tenifold`。
**拼错不静默退回默认** —— 那会让"我配了 tenifold"变成一句假话，所以直接报错。

两套结果**并列落盘**（`virtual_perturbation.csv` 与
`virtual_perturbation_tenifold.csv`），一致性只比排名、不比数值
（一阶是 z 分数的 L2 范数，tenifold 是流形欧氏距离，**量纲不同**）。
一致性高**不说明哪个对**：两法共享同一个上游（`07_grn` 的共表达边），
可能只是共享了同一个偏差。

**出图默认画 tenifold（真方法），它没跑出来时退回一阶**，标题写明是哪个。
两法不画在同一根 x 轴上 —— 那会让人以为可以直接比大小。

**`scTenifoldKnk` 没有细胞类型分辨率**：整份数据只建一个网络，输出是
「扰动基因 × 网络基因」的全局距离。**本仓库没有**把距离按细胞类型重新加权 ——
那会造出一个既不是 scTenifoldKnk、也不是本仓库一阶近似的新方法，
然后借它的名字发出去。

### 空敲除：网络表达不了该扰动 ≠ 效应为 0

`scTenifoldKnk` 的 `transcriptomeWide` 模式敲除一个基因的方式是把 WT 网络的
**那一行清零**（**scTenifoldKnk 包内源码** `R/scTenifoldKnk.R` L206–L207
`KO <- WT; KO[g, ] <- 0`）。
若该基因**出度本来就是 0**，清一行全 0 的行等于什么都没敲，`KO` 与 `WT`
逐位相同，`manifoldAlignment` 对两个完全一样的网络对齐，返回的"距离"
只剩浮点噪声（实测 ~1e-16）。

**它不报错、不给 NA、不给零** —— 给出一排看起来完全正常的数，
读表的人只会得出"敲除这个基因没有影响"。所有既有图门禁都会通过，
因为它们查的是"有没有产出"不是"产出对不对"。

判据是**结构事实（出度）**，不是"距离是不是很小"—— 后者会把一个真实
但微弱的扰动也误判掉。R 侧 `zero_outdegree_targets()` 算出度并整行置 NA，
meta 里带 `empty_knockout_genes` + `target_outdegree`；Python 侧
`compress_tenifold_distances()` 用**那份结构证据**打标记
（而不是从"整行都是 NA"反推 —— 那会把 R 侧别的 NA 原因也误标成空敲除），
空敲除行报 `None` 而**不是 `0.0`**，并在 `compare_engines()` 里
**显式剔除并记账**，而不是靠 `dropna()` 悄悄少几个点。

**两侧各有自检**：`Rscript scripts/lib/tenifold_knk.R --selftest` 验出度计算，
`python tools/selftest_tenifold.py`（**不需要 R**）验 Python 侧有没有把
那份证据正确消费成 `None` / 标记列，并且**端到端跑一遍
`run_tenifold_engine`（用假的 Rscript）** —— 抽出函数时最容易漏改的
`return` 里的旧变量名，`py_compile` 查不出来，只会在 CI 上 tenifold
真跑通那一刻抛 `NameError`，而那时 R 已经跑了 6.5 分钟。

### 一阶近似：能报什么，不能报什么

用 `07_grn.py` 的调控子边表做**一阶、单跳**线性传播：

    Δz_G = −z̄(G,C)        （敲除）   /   +overexpress_sd（过表达）
    Δz_Y = w(G,Y) · Δz_G   对 G 的每个直接靶基因

**`ko_magnitude` 里混了两个因素，必须分开报。** 它是
`|Δz_G| × network_sensitivity`，所以排在前面的基因往往只是
**在该细胞类型里表达高**，而不是在网络里被连得紧 —— 实测同一份结果，
`GATA2` 的 magnitude 是 5.41（排第二）而网络敏感度只有 0.81，
`TBX21` 的敏感度 1.84（全场最高）却只排第四。所以额外落一列
`network_sensitivity = √(Σ_Y w²)`（与表达无关）。

> **名次会变，值不会。** 早先写的是"`GATA2` 排第一、`TBX21` 排第三"；
> 当前产物里 magnitude 最高的是 `TYMS`（6.05，敏感度 0.94），
> GATA2 退到第二、TBX21 退到第四 —— 而 **5.41 / 0.81 / 1.84 三个值逐位未变**。
> 引用前看 `virtual_perturbation_status.json` 的 `top_by_effect`，不要从文档抄。

**不做多跳传播。** 在相关网络上做多跳会放大噪声，看起来像"网络效应"，
其实只是把相关系数乘了几遍。

**`PerturbNet` 始终没跑，`method` 必须仍然否掉它。** 状态里的 `method`
显式写出"**不是 PerturbNet**"，验收项检查这个否定词在不在。
（原来那句"也不是 scTenifoldKnk"已随 K-01b 删除 —— 那条判据在真跑通
R 引擎之后会**把正确的结果判红**。现在的判据是"声明的方法与状态里的
引擎一致"：记了 `engines_used` 就必须逐个在方法串里被点名，
没跑 tenifold 时必须写明**为什么没跑** —— 否则"没跑"和"跑了没结果"
长得一样。）

### 跨部分交接只走 CSV（§0.2）

`perturbation.targets_csv` 指向 Part 1 交接的基因表（一列 `gene`，可选一列
`logFC`），相对路径按 `data_dir` 解析。**没有它时回退到内部调控子，
并且必须在状态里写明回退了** —— 此时 `signature_alignment` 为空是预期的，
不是 bug。回退不报错、只是把"Part 1 的疾病签名"换成了"本步自己的调控子"，
下游解读完全变了，所以 `target_source` 是必填字段。

### 本仓库还有第二条交接：Part 2 → Part 3 的参考 h5ad

上面那条是 **Part 1 → Part 2**（进）。还有一条**出**的：
Part 3（空间转录组）的 `deconvolution.reference: h5ad` 要的就是本仓库
`03_cluster_annotate.py` 写出的 **`clustered.h5ad`** ——
带细胞类型标签的单细胞参考。

**契约有三条，前两条不满足会让 Part 3 静默算错或 KeyError：**

| # | 要求 | 不满足时 Part 3 的表现 |
|---|---|---|
| 1 | `layers['counts']` 是**原始计数** | 退回用 `.X`（log 值）求参考谱 → **解出的比例没有意义，而它看起来仍像一组比例** |
| 2 | `obs` 里有细胞类型列（本仓库是 `celltype`）| 配 `celltype_key` 时 KeyError |
| 3 | 基因集是**交集** | 静默接受 —— 参考若是 HVG 子集，参考谱覆盖的基因就变少 |

**第 3 条最容易被忽略，而且它是本仓库的实际情况**：`clustered.h5ad` 来自
`integrated.h5ad`，那是**2000 个高变基因的子集**，不是全基因集。
Part 3 因此只在共同基因上建参考谱。

**交接的两端都要能被核对。** 所以产出侧把契约**主动记进**
`cluster_status.json` 的 `part3_reference`（含 `contract` 与 `limitations`），
验收项 `§0.2 Part 3 参考导出的契约` 据此检查 ——
只检查消费侧（Part 3 记没记）的话，产出侧悄悄丢掉 counts 层
**不会被任何人发现**。

### 候选基因的完整边表要落盘

`07_grn.py` 的 `tf_regulons.csv` 里 `top_targets` 只有前 12 个靶基因、
**且没有权重**，虚拟扰动需要全部靶基因及其相关系数。
所以 07 额外写 `tf_regulon_edges.csv`（`tf, target, corr`）。
**靶基因名与权重必须成对取** —— 分开写（一个列表取名字、另一个按位置取
相关系数）在过滤 `corr > 0` 之后会错位，而错位不报错，
只会给每个靶基因配上一个别人的相关系数。


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

姊妹项目 `geo-normal-pipeline-skill` 的规则 12 记了同一类问题（那边是极小的
P 值把浮点末位放大到可见）。单细胞这边实测的表现是**拓扑变了**。

**证据（同一 Python 3.12、同一批包版本的五轮 CI）：**

| 轮次 | commit | 并行度设置 | dpt | palantir | **scfates** | scFates 拓扑 | 平均 ρ |
|---|---|---|---|---|---|---|---|
| A | `8f3f56c0` | 无钉 | +0.5856 | +0.4674 | **+0.5644** | 4 片段 / 6 milestone | +0.6363 |
| B | `9af13b42`（**只改注释**） | 无钉 | +0.5856 | +0.4674 | **+0.5296** | 6 片段 / 8 milestone | +0.6254 |
| C | `5b241a3` | +OMP/OB/MKL/CORETYPE | +0.5856 | +0.4674 | **+0.5328** | 6 片段 / 8 milestone | +0.6265 |
| D | `db778bb` | +`NUMBA_NUM_THREADS=1` | +0.5856 | +0.4674 | **+0.5328** | 6 片段 / 8 milestone | +0.6265 |
| E | `f99eab1` | 同上 | +0.5856 | +0.4674 | **+0.5328** | 6 片段 / 8 milestone | +0.6265 |
| F | `e549477` | 同上 | +0.5856 | +0.4674 | **+0.5646** | 6 片段 / 8 milestone | +0.6353 |

`dpt` 与 `palantir` **六轮逐位相同**，上游的聚类数（10）、分辨率扫描、
CytoTRACE 也完全一致。**只有 scFates 一个在变。**

### 20.1 两次归因都被否证了

**第一次：多线程 BLAS 归约顺序。** 把 geo 规则 12 那两组变量搬了过来 ——
**C 轮否证**：钉住之后 scfates 从 +0.5296 变成 +0.5328，没有收敛。

**第二次：`pynndescent` 的 Numba 并行。** 读源码发现扩散图那一步经过
scanpy 默认的 `method="umap"`，于是补了 `NUMBA_NUM_THREADS=1` ——
**D/E/F 三轮否证**：0.5328 / 0.5328 / **0.5646**。
最后那个值几乎回到 A 轮未钉时的 +0.5644。

**两次否证留下的线索：**

- **拓扑稳住了。** A 轮是 4 片段/6 milestone，B 轮之后**六轮全部是
  6 片段/8 milestone** —— 钉并行度确实让**离散的**结构稳定了。
- **但连续量（ρ）还在漂**，而且看起来有两个吸引子：`≈+0.5328` 与
  `≈+0.5645`。**残留随机源尚未定位。**

**这一条比"我修好了"更有用：钉并行度能让拓扑稳定，但不足以让 ρ 稳定。**
第三次动手之前先读日志 —— 本仓库规则 12 就是为这个写的。

### 20.2 真正的路径（读源码得到）

`05_trajectory.py` 的 `compute_scfates()` 走的是：

```
scf.pp.diffusion(device="cpu")
  -> palantir.utils.run_diffusion_maps(data_df, n_components=10, knn=30, alpha=0)
       -> compute_kernel(..., backend="scanpy")            # palantir 默认后端
            -> scanpy.neighbors.Neighbors(temp)
                 .compute_neighbors(n_neighbors=30, n_pcs=0, method=None)
                      -> scanpy 默认 method="umap"
       -> diffusion_maps_from_kernel(kernel, n_components, seed=0)
            -> eigs(T, ..., v0=rng.random(...))            # v0 来自 seeded RNG，没问题
```

**两处关键事实：**

1. **扩散图那一步没有 seed 可传。** `scFates.pp.diffusion`（1.2.5
   `preprocessing/diffusion.py`）**签名里根本没有 `seed`**，内部调
   `run_diffusion_maps` 时也不传 —— 靠 palantir 的默认 `seed=0`。
   `05_trajectory.py` 把 seed 传给了 `sct.tree()` 和 `sct.pseudotime()`
   （这两处是对的），**但扩散图那一步没地方传**。
2. **特征求解器是干净的**（`eigs(..., v0=rng.random(...))`，`v0` 来自
   `np.random.default_rng(0)`）。
3. simpleppt 的 PPT 主曲线树（`ppt.py:210-213`，`np.random.seed(seed)`
   之后 `np.random.choice` 取初始节点）本身是 seed 了的 ——
   它只是把上游的末位差异放大成可见的 ρ 变化。

### 20.3 三条要记住的

1. **看到"我明明设了种子"时，先问"这个函数有没有 seed 参数"。**
   没有的话设多少遍都没用。`set_seed()` 只覆盖用 numpy/random 的随机调用；
   第三方库内部的迭代求解器、并行归约、GPU 内核都不在里面。
   排查顺序：`grep -n seed <包的源码>` —— 看 `seed` 是**出现在签名里**
   还是**只出现在 docstring 里**。
2. **并行度有三个独立的旋钮，缺一个都不够**（但它们**不保证**解决）：

   | 旋钮 | 管什么 |
   |---|---|
   | `OMP_NUM_THREADS` / `OPENBLAS_NUM_THREADS` / `MKL_NUM_THREADS` | BLAS 归约顺序 |
   | `OPENBLAS_CORETYPE` | OpenBLAS 按宿主 CPU 型号分发 SIMD 内核（线程数管不到） |
   | **`NUMBA_NUM_THREADS`** | **Numba `prange` 的线程数** |

   ```yaml
   env:
     OMP_NUM_THREADS: 1
     OPENBLAS_NUM_THREADS: 1
     MKL_NUM_THREADS: 1
     OPENBLAS_CORETYPE: Haswell
     NUMBA_NUM_THREADS: 1
   ```

   **这五个变量现在都设着，但 scFates 的 ρ 仍然在漂（D/E/F 三轮）。**
   设它们是对的（拓扑因此稳住了），**但不要以为设了就可复现。**
3. **"版本不同"不是万能借口。** 版本相同也能对不上。别一看到数字变了
   就归因到版本上 —— 那会掩盖真正的回归。先比对 `key_versions`，
   相同就去看**哪些量没变**。

### 20.4 结论怎么报

**scFates 的 ρ 不要当单一确定值报。** 实测范围 **+0.5296 ~ +0.5646**，
平均 ρ 因此是 **+0.6254 ~ +0.6363**。`05_trajectory.py` 会把这段写进
`trajectory_status.json` 的 `reproducibility` 字段 ——
**结论的适用范围要跟着产物走，不能只写在 AGENTS 里。**

`dpt` / `palantir` / `cytotrace` 六轮逐位相同，可以按确定值报。

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

起因是台账 E-48 与 E-49 两条**互相独立**的缺陷，它们暴露了验收层与门禁层
**四个各自独立的盲区**。四条补强分别堵一个，缺一条那类缺陷就还能再犯一次。

### 25.1 嵌套 `status` 的 `failed` 必须判红（验收层）

`status` 不只在顶层。`qc_status.json` / `integration_status.json` /
`cluster_status.json` / `grn_status.json` 里都有**嵌套**的 `status` 字段，
而旧验收层只看顶层 `d.get("status")`。

实测（真 artifact `scrna-results-55`）：顶层分布 `{ok: 7, not_configured: 1}`，
**嵌套分布 `{ok: 5, failed: 1, not_applied: 1, not_done: 2}`** ——
唯一那条 `failed` 是 `grn_status.json → regulon_vs_pseudotime.status = failed`
（`ValueError: Invalid unit 6.800000000000001 in 'figsize'`，E-48 的后果），
而它的顶层 `status` 是 `ok`。**顶层全绿、里面已经崩了。**

判据必须**只把 `failed` / `error` / `fail` 判红**：

| 嵌套值 | 含义 | 处置 |
|---|---|---|
| `failed` / `error` / `fail` | 崩了 | **判红（required）** |
| `not_applied` / `not_done` | 设计如此地没做 | 可见不阻断 |
| `ok` | 正常 | 放行 |

误判的代价是双向的：把 `not_applied`（单样本不做整合）判红会让每个 job 都红，
把 `failed` 放行则正是 E-48。实现是 `_iter_nested_status(obj, path=())`
递归产出 `(路径, 值, 同级 reason)`，扫**全部** `*status.json`
（不只可选步骤那 5 个 —— 上面四个文件里三个是必需步骤的）。

### 25.2 图验收必须条件化，且要覆盖全部出图脚本（验收层）

`REQUIRED_FIGURES` 原来只有 14 条，**不含 `02-07-*`（三张 GRN 图）与
`02-08-01`** —— E-48 让 5 张图从未产出，而验收层根本不知道它们该存在。

补齐之后还有第二个问题：**图不能无条件要求产出**。`trajectory.enabled: false`
时 `02-05-*` 一张都不该有，旧写法一关轨迹就必红。所以改成由图名第 2 段
（`02-07-01` 的 `07`）反查所属步骤，从 `read_state(cfg)` 取该步骤的
status 与 required：

- 所属步骤**可选且未跑成** → `required: False`，detail 写明
  `步骤 grn = not_configured，本轮不要求产出`
- 否则 → `required: True`，缺失时 detail 带 `**缺失**（步骤 grn = ok）`

**兜底是 25.1 那条独立扫描** —— 步骤崩了的时候，即使图项因条件化放行，
嵌套 `failed` 仍会把它判红。

### 25.3 溢出检测必须落盘，不能只打 WARN（产物级）

E-49 的形态：`_content_overflow()` **正确检测到了**
`width_overflow_frac=0.3523`、**正确打了 WARN**，然后**没有任何人读** ——
标题超宽 35%，在 `savefig.bbox: standard` 下被静默裁掉，图照样生成、
门禁照样绿。

**"检测到了"不等于"有人会知道"。** 现在 `save_fig()` 里的
`_record_figure_overflow(cfg, name, bad)` 把溢出**落盘**到
`results/<dataset_id>/figure_overflow.json`（`{figures: {名: bad}, n_overflow}`），
验收层读它并**判红（required）**，detail 写明"**会被静默裁掉**"。

三条配套约束：

1. **与状态文件同理，跑前必须删旧文件**（`_ovf_stale.unlink()`）——
   否则图修好了旧记录还在，验收把已修好的图判红。
2. **只累积、不覆盖** —— 同一次运行多张图溢出要全部记下。
3. **记录可被修复** —— 标题改短后不应再新增记录（已做双向验证）。

### 25.4 门禁要查"内容贴边"，不能只查"有没有墨"（门禁层）

`check_figures.mjs` 原来只查"有没有墨 / 是不是糊死"。
**被裁掉的图照样有墨** —— 这正是它当年漏掉 E-49 的原因。

裁切的物理后果是**墨迹延伸到画布边缘**，所以量非背景像素外接框到左右边的
距离。阈值 `EDGE_MIN_PX = 3`，实测标定（scrna 34 张 + spatial 76 张）：

- scrna 边距分布 `{0:1, 10:9, 11:2, 12:6, 13:13, 14:2, 41:1}` ——
  `L<3 或 R<3`、`L==0 或 R==0`、`L<6 或 R<6` **都只命中 `02-08-01` 一张**
  （就是 E-49 那张），零误伤。
- spatial 76 张全部通过，零误伤。

两条刻意的取舍：

1. **只判左右，不判上下。** 实测 `02-05-05-unit1-pseudotime-by-cluster`
   上边距 = 0 —— 那是布局取舍，不是裁切。单看某一边会误伤。
2. **暗底图跳过。** 整幅都是"墨"，外接框必满幅，判了全是假阳性。

> 本文件与 `spatial-pipeline-skill/tools/check_figures.mjs` **必须逐字相同**
> （两仓同一份，见该文件头注释）；改一侧必须同步另一侧并比对哈希。

### 25.5 上线首跑就抓到一条新缺陷 —— 修好一个缺陷会暴露被它掩盖的下一个

四条补强推送后第一次 CI（run `36140048953`，commit `e2d7a2c`）**验收判红**：

```
[FAIL] 图 02-07-01-unit5-tf-smarca4-trend 内容超出画布  **会被静默裁掉**：{'width_overflow_frac': 0.0371}
验收失败：1 项必需检查未通过
```

`scripts/07_grn.py` 的单 TF 面板标题是**单行**，而面板宽度只有
`W_SINGLE = 89 mm`：

```python
ax_t.set_title(f"{tf} activity along pseudotime (mean ± 1 SD per bin)")   # 旧
```

本地用仓库自己的 `apply_style` + `_content_overflow` 逐名标定，得到
`FOSL2 +0.0093 / ATF4 −0.0019 / E2F2 −0.0021 / SMARCA4 +0.0371`
（**与 CI 逐位吻合**）；假想的长 TF 名 `SMARCAD1 +0.0485`。即**在阈值附近
徘徊、随 TF 名长度而变** —— 换一份数据必然再踩。

**这批面板此前从未渲染过**：Q-24 的三元素 `figsize` 让整段在 `plt.subplots`
就抛异常，4 张面板一张都没画出来。所以这是它们的**首次渲染**，
缺陷与图同时诞生 —— **旧缺陷一直在掩盖新缺陷**。

修法是**折行**而不是删字：

```python
ax_t.set_title(f"{tf} activity along pseudotime\n(mean ± 1 SD per bin)")
```

第二行宽度**与 TF 名无关**，本地对全部 8 个名字（含 3 个假想长名）标定
**一律 −0.0238**（固定 2.4% 余量）。对照候选：去掉尾巴只剩 −0.0143
（余量太薄）、去掉 `along` 在 `SMARCAD1` 上仍 `+0.0228` 判红。

**两条规则：**

1. **宽度受限的小面板，标题折行，不要删字。** 删字只是把阈值往下挪一点，
   输入一长又踩到；折行让第二行宽度与输入无关，是**结构性**修法。
2. **修好一个缺陷后要重跑一遍全部检查。** "图重新出现"这个动作本身必须
   触发一次完整复检 —— 不能只看"图现在有了"，因为它第一次出现时可能就带着
   一个此前没人有机会看见的毛病（本条就是）。

台账：`governance/15_ERROR_LEDGER.md` E-51（E-49 的同族，都是
`savefig.bbox: standard` 下静默裁切）。


## 26. native 崩溃绕过 `except`：`sc.pp.scrublet` 偶发 SIGSEGV（E-52，2026-09-25）

CI 的「跑流水线」步骤偶发 **exit 139**。日志停在 `01_qc.py:226` 的
`过滤: …` 之后一行 `Segmentation fault (core dumped)`，**没有任何
Python traceback**。

**它是长期存在的，不是哪一轮引入的。** 全量 60 个 run 里有 8 个
`run_attempt>1`，**全部**是"att1 失败 + rerun 成功"形态，最早一个是
`06977d0`（2026-09-20T11:34:50Z）。**这 8 个的 att1 崩溃位置逐字相同**，
跨 5 天、跨 7 个 commit。**"rerun 能过"把它掩盖了整整 5 天。**

### 26.1 崩溃点怎么锁死的

`01_qc.py:226` 的 `过滤:` 是崩溃前最后一行，下一句是 L229
`db = run_scrublet(adata, cfg)`，里面 L79 `sc.pp.scrublet(...)`。
`scripts/lib/common.py` 第 62 行的 `print(..., flush=True)` **逐行 flush**，
所以"日志最后一行"就是"代码走到哪一行"，不存在缓冲区丢日志。

**但 `01_qc.py:85` 的 `except Exception` 接不住它** —— native SIGSEGV
直接杀进程，`qc_status.json` 里连一行失败原因都没有。

### 26.2 根因：OpenBLAS 0.3.34 的 `dgemm_kernel_HASWELL` 栈越界

**这一步是靠 `PYTHONFAULTHANDLER: 1` 拿到的真栈定下来的**（run
`36150495910`）—— 在此之前只能靠猜。栈的落点（逐字，自下往上读）：

```
scipy/sparse/linalg/_interface.py:1118  in _matmat      ← return self.A @ X
scipy/sparse/linalg/_interface.py:451   in _shared_matmat
scipy/sparse/linalg/_interface.py:491   in matmat
scipy/sparse/linalg/_eigen/_svds.py:511 in svds         ← Av = X_matmat(eigvec)
scipy/_lib/_util.py:306                 in wrapper
sklearn/decomposition/_pca.py:740       in _fit_truncated
sklearn/decomposition/_pca.py:542       in _fit
sklearn/decomposition/_pca.py:440       in fit
scanpy/preprocessing/_scrublet/pipeline.py:84  in pca
scanpy/preprocessing/_scrublet/__init__.py:448 in _scrublet_call_doublets
scanpy/preprocessing/_scrublet/__init__.py:243 in _run_scrublet
scanpy/preprocessing/_scrublet/__init__.py:298 in scrublet
scripts/01_qc.py:79                     in run_scrublet
```

**栈里一帧 numba 都没有。** 崩溃是一次**稠密 GEMM**：`_scrublet/pipeline.py`
第 79 行先把稀疏矩阵 `.toarray()`，第 84 行交给
`PCA(svd_solver="arpack")`，sklearn 走 `svds` → ARPACK 求特征向量后做
`X_matmat(eigvec)`，而 `self.A` 已是**稠密 ndarray**，于是落到
`dgemm_kernel_HASWELL`。

上游 **OpenBLAS #6026**（2026-09-11 报、09-13 关）的形态与本仓**逐项吻合**：

- 环境原文即 `OpenBLAS 0.3.34, as bundled in the scipy-openblas64 wheel
  shipped with numpy 2.5.2 and 2.5.3` + `OPENBLAS_CORETYPE=Haswell` +
  GitHub 托管 `ubuntu-latest` —— 正是本仓 job env 凑齐的组合。
- **机制**：Haswell 内核把 k 维分块打包进一个固定 `0x7080` 字节的**栈
  缓冲区**，每步写 96 字节且**对 k 无上界**；k ≳ 319 时越界覆盖
  callee-saved 寄存器与返回地址槽，最后那条 `ret` 直接 SIGSEGV。
- **为什么偶发、为什么 rerun 能过**：level-3 的 blocking（决定 k 的大小）
  **在运行期按宿主 cache 拓扑选取**。GitHub 的 Azure fleet 混着不同型号
  CPU —— 落到 k 会变大的机器上就崩，落到 k ≤ 256 的机器上就过。
  **同一 commit、同一批包版本、相隔 8 分钟的两个 run 一成一败**，只能
  这样解释。
- 上游实测**受影响的机器上约三次崩一次** —— 与本仓观测到的 rerun 频率一致。

**附带风险：越界若只覆盖了保存寄存器、没碰到返回地址，进程会正常返回而
数值被污染。** 即"没崩"不等于"算对了"。

**各 numpy wheel 自带的 OpenBLAS**（读 `numpy.libs/libscipy_openblas64_*.so`
里的版本串得到）：

| numpy | 自带 OpenBLAS | 状态 |
|---|---|---|
| 2.4.6 | 0.3.31.188.0 | 正常 |
| 2.5.0 / 2.5.1 | 0.3.33.112.0 | **最后一个正常版** |
| 2.5.2 | 0.3.34.0.0 | 回归 |
| **2.5.3**（本仓此前浮动到的） | **0.3.34.106.0** | 回归 |

### 26.3 处置

**修复是钉 numpy 上界**（`requirements.txt`）：
`numpy>=2.1,<2.5.2` —— 下界 2.1 是 `anndata 0.13` 的硬要求，上界躲开
OpenBLAS 0.3.34。**上游修复（0.3.35）尚未进入任何 numpy wheel**
（PyPI 上最新仍是 2.5.3），所以只能钉，不能等。

**上线复验（run `36153515064` / commit `f19eff9`，2026-09-25T15:20Z）**：
**attempt 1 一次过** —— 这是本条 60 个 run 里**第一次** att1 就成功
（此前 8 个 rerun 全是 att1 失败）。日志 segfault 命中 **0**；CI 实际装的
是 `numpy-2.5.1-cp312-cp312-manylinux_2_27_x86_64…whl`，artifact
`run_manifest.json` 的 `versions.numpy = 2.5.1` —— **两处对上，说明钉上界
真的生效，而不是"这次恰好没崩"**。artifact `scrna-results-61` 亲验：
`acceptance.json` = **75 项全过 / 0 失败**、`figure_overflow.json`
不存在、`figures/` 39 PNG + 39 PDF（`02-07-01-unit1..5` 全在）。

> **"一次过"本身不是证明。** #6026 是**宿主相关**的偶发缺陷 —— 单轮成功
> 只说明"这一轮的宿主 k ≤ 256"。真正的证据是**版本已经落到正常侧**
> （2.5.1 → OpenBLAS 0.3.33.112.0），所以这一轮的成败**不再是随机的**。
> 要再确认，观测口径是**后续若干轮的 rerun 率应降到 0**，而不是
> "再看一轮是不是绿的"。

**`PYTHONFAULTHANDLER: 1` 继续留着**，理由变了：它这次兑现了价值（没有它
根因定不下来），而 numpy 上界只覆盖 OpenBLAS 这一个来源。

**`NUMBA_THREADING_LAYER: workqueue` 保留，但它的立项理由已被推翻。**
它是按"numba 并行运行期"的假设加的，而栈证明 numba 不在路径上 ——
带上它照崩（run `36150495910`）。保留只剩一条防御性理由：numba 的 `omp`
线程层会加载一份 libgomp，而 scikit-learn 自己 bundle 了一份
（`sklearn.utils._openmp_helpers`，实测在崩溃轮的 218 个 extension
modules 里），**同进程两份 OpenMP 运行时**是上游有记录的崩溃形态。

### 26.4 五条规则

1. **native 崩溃不能靠 `except` 兜底。** `except Exception` 只覆盖 Python
   层异常；凡是有 C 扩展参与的关键步骤，都要问一句"它要是崩了，我会不会
   连日志都没有"。
2. **"rerun 能过"不等于"没有问题"。** 重试成功是掩盖，不是修复。本条被
   掩盖了 5 天、8 次 rerun —— **应当把重试率本身当成一项可观测指标**。
3. **给 native 崩溃留一条取证路径，比急着改代码值钱。** 从"日志停在
   `过滤:`"到"OpenBLAS 栈越界"之间隔了一整轮错误的归因；`faulthandler`
   一条 env 就把它终结了。**改代码之前先确认自己有没有能力观测。**
4. **依赖的二进制缺陷也是本仓的缺陷。** 这次崩的不是本仓任何一行代码，
   而是 numpy wheel 里打包的 OpenBLAS。**依赖"没钉上界"就是把自己的
   稳定性交给上游的发布节奏** —— `requirements.txt` 开头的告诫因此
   不只是"API 会变"，也包括"二进制会坏"。
5. **同一份二进制在不同宿主机上行为可以不同。** "本地复现不了"不等于
   "没问题"：k 的大小取决于宿主 cache 拓扑，这是**设计如此**的行为。
   遇到宿主相关的间歇失败，先找"哪个量在运行期按机器选"。

台账：`governance/15_ERROR_LEDGER.md` E-52；任务行 `governance/02_TASKLIST.md` V-02。

## 27. 步骤"没做成"必须传出来：返回值要接住、顶层 `status` 也要扫（E-56，2026-09-26）

### 27.1 现象：一条已存在但从未生效的判据

`main_analysis.py` 的嵌套失败扫描本意是"任何 `*status.json` 里写了失败就判红"，
但写成 `if p and str(v).lower() in NESTED_FAILED_VALUES`。`_iter_nested_status(d)`
在**顶层**调用时 `path` 是**空元组** —— falsy —— 于是**顶层 `status` 一个都扫不到**。

而本仓有 6 处代码正是以"写顶层非 ok 状态 + 正常 `return`"的方式失败：

| 位置 | 取值 |
|---|---|
| `scripts/05_trajectory.py:275-279` | `bad_root` |
| `scripts/05_trajectory.py:313-317` | `missing_clusters` |
| `scripts/05_trajectory.py:363-370` | `insufficient_methods` |
| `scripts/07_grn.py:352-355` | `no_tfs_in_data` |
| `scripts/08_virtual_perturbation.py:901-911` | `no_candidates` |
| `scripts/04_pseudobulk_de.py:178-180` | `no_celltype_testable` |

**结果：这一步什么也没产出，`state.json` 记 `ok`、`acceptance.json` 记绿。**
而且 `insufficient_methods` / `bad_root` / `missing_clusters` / `no_candidates`
这四个词**连 `NESTED_FAILED_VALUES = ("failed", "error", "fail")` 都不在** ——
即使 `if p` 修好，顶层这四个也照样漏。

### 27.2 第二处同病根：返回值被丢弃

`run_all` 的步骤循环把 `fn(cfg)` 的返回值**直接丢弃**，写死
`record_step(cfg, sid, "ok", ...)`。步骤函数的"跑完了、但结果是『没做成』"
路径是 `write_json(状态文件, status)` + `log_warn(...)` + `return status`，
**不抛异常** → `except` 分支根本不进。

同一病根还在 **9 个步骤脚本各自的 `__main__` 块**里（独立运行时同样丢弃
返回值、无条件记 `ok`）—— 即"编排器跑"和"单独跑"两条路径**都**把失败记成成功。

**这是 E-48 的另一半。** E-48 当时的修法是"让嵌套失败可见"，但只覆盖了
"写进嵌套字段的失败"：写进**顶层**的失败与"返回值里的失败"都没被覆盖。
E-48 的教训是"异常被降级成一个没人读的字段"，本条是它的两个变体：
**字段写在没人扫的位置**，和**结果压根没被记下来**。

### 27.3 处置

`scripts/lib/common.py`：

- 两个**分开的**取值集合：`STEP_ABORT_VALUES`（判红）与 `STEP_SKIP_VALUES`
  （可见不阻断，只作文档用途）。前者含通用 `failed`/`error`/`fail` 加上表
  六处的具体取值，以及 07 的 `no_regulons`、06 的 `no_pairs_in_data`/`no_signal`、
  04 的 `no_usable_pseudobulk`/`column_missing`。
- `result_status_of(res) -> str`：dict 且带 `status` 键则取之，否则 `"ok"`。
- `classify_step_result(v) -> str`：`"ok"` / `"abort"` / `"skip"`。
  **未登记的取值归 `"skip"`** —— 只把确知是失败的判红（假阳性比假阴性更危险）。
- `record_step(..., result_status: str = None)`：`entry["result_status"]` 与
  `entry["status"]` **分开存**。后者是编排器记的"有没有崩"，前者是步骤自述的
  "做成了没有" —— `status="ok"` + `result_status="bad_root"` 正是要显形的组合。
  默认 `None` → 不写这个键，向后兼容。

`scripts/main_analysis.py`：步骤循环 `res = fn(cfg)` → `record_step(...,
result_status=result_status_of(res))`；`except` 分支加 `result_status="failed"`；
步骤验收循环新增"步骤 X 自述状态"检查；嵌套扫描去掉 `if p`，改
`key = ".".join(p) + ".status" if p else "status"`（照抄姊妹仓 spatial 的写法，
`spatial-pipeline-skill/scripts/main_analysis.py:380-398` 本来就是对的）。
9 个步骤脚本的 `__main__` 块同步改造。

### 27.4 四条规则

1. **"让失败可见"的修法要覆盖失败的所有承载位置。** 修完要问："这个失败还可能
   被记在哪里？"—— 返回值、顶层字段、嵌套字段、日志、退出码。E-48 只修了
   嵌套字段一处，另外两处照旧漏了。
2. **`if p` 在"空路径 = 顶层"的语义下必然滤掉顶层。** 空元组/空字符串/空列表
   都是 falsy，而写的时候看起来完全合理。凡"空值有语义"的地方（顶层路径、
   根节点、默认分支），判据必须显式写 `if p else` 而不是 `if p`。
3. **顶层与嵌套是两个语义域，取值域重叠不等于可以共用判据。** `tenifold.status`
   可以是 `timeout`/`no_candidates` **而内置引擎照样出了结果（顶层 `ok`）**；
   `annotation.celltypist` 的 `model_unavailable`、`liana` 的 `api_not_found`
   同理。合并两个集合会让**每个 job 都红**，而常年假红会训练人忽略告警。
4. **标定必须用本次 CI 的 artifact，不能用工作区里可能陈旧的副本。**
   `figure_review/scrna-pbmc3k` 是旧 run 的产物，内含 E-48 遗留的
   `grn_status.json → regulon_vs_pseudotime.status = "failed"` ——
   用它做正向标定会把真判据误判成误报。

台账：`governance/15_ERROR_LEDGER.md` E-56；任务行 `governance/02_TASKLIST.md` R-03。
标定脚本：`D:\tmp\_s1\calib_e56.py`（正向 + 反向 19 类注入 + 回退版对照）。




