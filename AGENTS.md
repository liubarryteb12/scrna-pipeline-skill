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
