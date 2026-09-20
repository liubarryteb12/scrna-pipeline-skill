# 故障排查

按"报错信息"索引。每条都记了**为什么**，因为知道原因才能判断
换个数据会不会再犯。

> **"每轮运行在什么条件下跑出来的"看 `run_manifest.json`**，
> 结构见 [`module0.md`](module0.md)。排查"数值和上次不一样"时先看它的
> `key_versions` 与 `seed`。
>
> **但"版本不同"不是唯一解释，别把它当成万能借口。** 实测同一份代码、
> **同一个 Python 3.12 + 同一批包版本**的两轮 CI，轨迹一致性给过
> `+0.6363` 和 `+0.6254` —— 那两轮之间只隔了一个**只改注释**的 commit。
> 所以顺序是：**先比对版本 → 版本相同再比对 `key_versions` 里有没有
> BLAS 相关项 → 都相同就怀疑多线程归约**（见下面"数值不稳定"一节）。
> **先比对版本，再怀疑代码**这条仍然对，但"版本一样就一定是回归"是错的。

---

## 下载相关

### `HTTP Error 403: Forbidden` 下载 10x 数据

**原因**：10x 的 CDN（`cf.10xgenomics.com`）对 `Python-urllib/3.x`
的 User-Agent 直接返回 403。同一个 URL 用浏览器或 PowerShell 请求是 200。

**容易误判**：报错信息只有 `403 Forbidden`，看起来像"链接失效了"，
于是去换数据集 —— 方向完全错了。

**已修**：`00_fetch.py` 的 `download()` 带了浏览器 UA。

### `下载不完整: N/M 字节`

网络中断。`download()` 会重试 3 次，并把内容写到 `.part` 再改名 ——
所以中断不会留下半个文件被下次当成"已缓存"。

### 解压后 `找不到 matrix.mtx`

10x 的 tar 里目录名带版本后缀（`filtered_gene_bc_matrices/hg19/`），
不同版本结构不同。`read_10x_tar()` 会搜含 `matrix.mtx` 的目录，
不写死路径。若仍找不到，说明这个 tar 不是 10x 计数矩阵格式。

---

## 输入校验

### `输入矩阵不是原始整数计数，拒绝继续`

`00_fetch.py` 的硬门禁。**这是故意的，不要绕过。**

把已经 normalize 或 log 过的矩阵当计数喂进去：
- QC 指标（`total_counts`）失去意义
- `seurat_v3` 的负二项拟合前提被破坏
- 拟bulk 的 DESeq2 负二项模型不成立

**而所有图和数字看起来都正常。**

修法：用原始计数矩阵（10x 的 `matrix.mtx` / `filtered_feature_bc_matrix`）。

### `细胞数 N < 50`

过滤后剩不到 50 个细胞。检查 `qc` 段的阈值是不是对本数据集太严 ——
高线粒体组织（心肌、肝）用 `max_pct_mt: 10` 会滤掉大部分细胞。

---

## 依赖相关

### `ModuleNotFoundError: No module named 'skmisc'`

`seurat_v3` 口味的 HVG 需要 `scikit-misc`。流水线会**降级到 `seurat`
口味并记录**（`hvg_fallback` 字段），不会失败。

但要注意：两种口味选出的 HVG 集合不同，**下游 PCA/聚类/marker
全部跟着变**。所以降级不是"没影响"。

修法：`pip install scikit-misc`。

### 双细胞检测 `status: package_missing`

scrublet 需要 `scikit-image`。注意它的**导入名是 `skimage`**，
不是 `scikit-image` —— 用错名字做检查会得到"没装"的假结论。

修法：`pip install scikit-image`。

### `harmonypy` 装不上（编译失败）

Python 3.13/3.14 上 `harmonypy` 还没有 wheel，要从源码编译，
常因缺 numpy 头文件失败。

**已修**：CI 用 Python 3.12，这些包都有现成 wheel。

### `pip install harmonypy` 失败但 `integration.method: none`

不影响 —— 单样本数据不需要 harmony，代码不会 import 它。
多供体数据才需要。

---

## scanpy / matplotlib API

### `TypeError: highly_variable_genes() got an unexpected keyword argument 'return_fig'`

scanpy 1.12 移除了 `sc.pl.highly_variable_genes` 的 `return_fig`。
它画到当前 figure 上，用 `plt.gcf()` 取。

**已修**，并且 `requirements.txt` 钉了 `scanpy<1.13`。

### `TypeError: Axes.boxplot() got an unexpected keyword argument 'labels'`

matplotlib 3.9 把 `labels` 改名为 `tick_labels`，3.11 起旧名字直接报错。

**已修**：`05_trajectory.py` 先试 `tick_labels`，TypeError 时回退 `labels`。

### `PerformanceWarning: DataFrame is highly fragmented`

`sc.tl.score_genes` 往 `adata.obs` 插列。循环 200 多次后 DataFrame
严重碎片化 —— 警告刷屏且越来越慢。

**已修**：`07_grn.py` 用 `adata.obs.pop(name)` 即时取走。

---

## 分析结果相关

### 拟bulk DE 记 `not_configured`

`design.sample_key` / `design.group_key` 为空。**这是正确行为，不是 bug。**

单细胞数据的样本量单位是供体/病人，不是细胞。没有这两列就无法构造
拟bulk，**也不应退化成细胞层面的检验**（那会给出 `p = 1e-300` 的假阳性）。

修法：在 config 的 `design` 段填上 obs 里实际存在的列名。

### `design.sample_key='X' 不在 obs 里`

列名拼错，或者数据本身没有样本信息（单样本 10x 数据通常没有）。
报错信息里会列出实际的 obs 列名。

### 细胞通讯只有 3/38 对可用

**这是 HVG 子集造成的。** 配体/受体大多是低表达的细胞因子/趋化因子，
几乎不进 HVG。

实测：HVG 子集（2000 基因）3 对可用；全基因集（13714 基因）27 对。

**已修**：`06` 和 `07` 强制走 `adata.raw`，`raw` 为空时直接报错。

### 细胞类型注释的 `score_margin` 是 0.0

第一名和第二名打分完全相同。**最常见的原因是签名基因完全重叠。**

实测：pbmc3k 里 `T_cell`（CD3D/CD3E/CD3G/TRAC/IL7R/LTB）与
`CD4_T`（CD4/IL7R/CCR7/LTB/TCF7）共享 IL7R 和 LTB，且 CD4 在
数据里表达低，导致两者打分一模一样。

**这不是 bug，是打分法的固有局限。** `assignment_confident: false`
如实标出。要区分需要看 CD4 的具体表达或换有监督方法。

### 拟时序的根看起来不对

现在的自动选根是"CytoTRACE GCS 最高的细胞"，仍**是统计判据不是生物学
判据**。（旧版用"PAGA 连通度最高"，在 PBMC3k 上挑中了浆细胞 —— 终末
分化 —— 当根；连通度高只说明它在图上居中，而居中既可能是祖细胞，
也可能是终末态。改用 GCS 后 DPT 与 GCS 的相关从 +0.25 升到 +0.59。）

修法：在 config 里设 `trajectory.root_cluster` 指定根。
`trajectory_status.json` 的 `root_selection` 会记录根是怎么来的。

### 各方法的拟时序符号对不上

**这是正常现象，不是 bug。** 实测四种方法的原始拟时序两两相关从
−0.50 到 +0.86 —— 符号都不一样。拟时序的符号是任意的：DPT 从根出发，
根选早或选晚，整条轴就反过来。

本流水线按「值越大越晚」统一方向后才比较。`trajectory_direction.csv`
记录每个方法是否被翻转。**未配置 marker 时退回 CytoTRACE 作方向参考，
产物里 `direction_source` 会写明**。

要生物学方向，在 config 里填 `trajectory.early_markers` / `late_markers`。

### `scFates` 的 `pseudotime()` 崩在 `IndexError: index 0 is out of bounds`

`pd.Series(adata.uns["graph"]["milestones"]) == t][0]` 报 IndexError。

**原因：漏了 `tl.cleanup()`。** 官方顺序是
`pp.diffusion → tl.tree → tl.cleanup → tl.root → tl.pseudotime`。
不 cleanup 时 `map_cells` 要读 `graph['milestones']`，而它直到
`pseudotime.py` 第 254 行才被写入 —— 调用在第 87 行，顺序矛盾。

### `scFates` 报 `rpy2 installation is necessary`

`tl.test_association` / `tl.test_fork` 需要 rpy2 + R + mgcv。
这与"云端不装 R"的设计冲突，所以**本流水线不用这两个函数** ——
沿轨迹的基因分析是自己实现的（`trajectory_genes.csv` /
`trajectory_modules.csv`）。

### `scFates` 的 `method='epg'` 报 `TypeError: object of type 'int' has no len()`

`elpigraph.utils.getPseudotime` → `networkx.add_edges_from`。
这是 elpigraph 与当前 networkx 的版本不兼容。用默认的
`method='ppt'`（本流水线的选择）。

### 轨迹报"没有方向性"

**这是正确行为。** 没有 RNA 速率（spliced/unspliced）时拟时序只能说
相似度排序，不能说分化方向。`trajectory_status.json` 的 `scvelo`
字段记 `not_done` 并说明原因。要速率需要从 Cell Ranger 的 `velocyto`
或 `--include-introns` 输出重新开始。

### `regulon×拟时序` 里几乎全部"显著"

n 大时正常。本数据 n=2652，|rho| 只要约 0.06 就能过 BH<0.05，
201 个里 190 个"显著"。**看效应量不看显著个数** ——
`tf_activity_vs_pseudotime.csv` 有 rho，产物里也写了效应量分布
（本数据 |rho| 中位 0.316，48 个 >0.5）。

---

## CI 相关

### job 绿了但 artifact 里缺文件

**这是 Part 1 踩过的真坑。** 验收清单检查的是 runner 工作目录里的文件，
artifact 的 `path:` 是**完全独立的列表**。两者不一致时验收通过、
job 绿、但用户下载后缺文件。

**已修**：`tools/check_artifact_paths.mjs` 交叉核对两边，
并在 CI 的静态检查阶段跑（装依赖之前）。

### 改了 `*.md` 但 CI 没触发

workflow 有 `paths:` 过滤，只在这些路径变更时跑：
`scripts/**`、`tools/**`、`assets/**`、`requirements.txt`、workflow 本身。

**这是故意的** —— 改文档不该烧 CI 时间。想强制跑用 `workflow_dispatch`。

### `check_figures.mjs` 报某张图是空白

**空白图是最难发现的问题**：job 绿、文件存在、大小正常，只有打开看
才知道没画出东西。

常见原因：
- 数据被上游全过滤掉了（散点为空）
- 图例占了整个画布
- 颜色映射把数据映射成了背景色

看 `results/<id>/figures/<name>.png` 确认。

---

## 本地环境

### 本地 `import scanpy` 失败

本地不装也能工作 —— 推云端看日志是设计好的路径。若要本地跑：

```bash
pip install -r requirements.txt
```

Python 版本用 3.10-3.12。3.13+ 上部分包没有 wheel。

### 中文字体

图里的标题和标签用英文，所以不依赖中文字体。
JSON 和 CSV 里的中文是 UTF-8，读的时候注意编码
（Windows PowerShell 需要 `$env:PYTHONIOENCODING="utf-8"`）。
