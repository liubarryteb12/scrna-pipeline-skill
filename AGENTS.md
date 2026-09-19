# 单细胞流水线的工程规则

姊妹项目 `geo-brca-microarray-skill/AGENTS.md` 的规则在这里同样适用。
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
  输出尺寸**，让上面的毫米约定失效。溢出改由 `_content_overflow()` 检测并告警。
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
| `capture_versions` | 全量已装包 + `KEY_PACKAGES` 逐个 | §0.3 |
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
2. **人工复核未确认不算失败。** `human_review` 默认就是 `pending`，
   判成 FAIL 会让每个 job 都红，反而没人看。但必须**可见**。
3. **`init_manifest` 必须清掉上一轮。** 和规则 14 同一个道理：
   上一轮的清单留在那里冒充本轮，比没有清单更糟。
4. **清单在步骤跑完之后才登记输入。** 可选步骤这轮有没有产物，
   跑完才知道；在开头登记会把"上轮残留"记成本轮输入。
5. **`pip freeze` 不用 subprocess 抓。** 沙箱下管道捕获会 EPERM，
   用 `importlib.metadata` 枚举。

`run_manifest.json` 落在 `results/` 下，随 artifact 一起上传 ——
**它必须和结果同时可及**，否则追溯链是断的。
