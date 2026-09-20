# 方法学与文献依据

本文记录每个步骤的方法选择及其依据，以及**该方法的适用范围**。
最后一条比前一条重要 —— 一个方法报出的数字，离开它的前提就没有意义。

---

## 0. 为什么单细胞分析的方法学限定必须写进产物

单细胞流程里绝大多数错误**不会报错**。HVG 顺序反了不报错、在细胞层面
跑检验不报错、用 HVG 子集做通讯分析不报错、拟时序根选错不报错。
它们都只是给出一个**看起来正常但错误**的结果。

所以本流水线的约定是：每个结论都带上 `limitations` 字段。
那不是免责声明，是结论的适用范围。

---

## 1. 质量控制

### 硬阈值过滤

`min_genes` / `max_genes` / `max_pct_mt` / `min_cells`。

线粒体比例上限的依据：破损细胞/死细胞的胞质 RNA 流失，线粒体转录本
相对富集。阈值随组织差异很大 —— PBMC 常用 5-10%，心肌、肝等
高线粒体组织可以到 20-30%。**用同一套阈值跨组织是错的。**

### 双细胞检测：scrublet

Wolock SL, Lopez R, Klein AM. "Scrublet: Computational Identification of
Cell Doublets in Single-Cell Transcriptomic Data." *Cell Systems* 2019.
doi:10.1016/j.cels.2018.11.005

原理：模拟双细胞（随机取两个细胞的表达谱相加），在 PCA 空间里看每个
真实细胞的最近邻里模拟双细胞的比例，得到 doublet score。

**局限**：
- 依赖 `scikit-image` 做邻居图；缺它时流水线记 `package_missing` 而非静默跳过
- 期望双细胞率是超参。本流水线按细胞数反推并夹在 5%-10%
- 同型双细胞（同种细胞的两个）检测不到 —— 它们的表达谱和单细胞一样

### ambient RNA：**不做，并如实记录**

胞外游离 RNA（"汤"）污染的校正需要：
- **SoupX**（R）—— 需要**空液滴**（raw 矩阵里的空液滴 barcode）
- **CellBender**（Python）—— 需要 GPU 或极长的 CPU 训练

GitHub 托管 runner **没有 GPU**（官方文档：Free/Pro/Team/Enterprise 的
"Maximum concurrent GPU jobs" 都是 "Not applicable"；GPU 只在 larger
runners 上有，另行计费）。所以本流水线默认 `qc.ambient_rna: none`，
在 `qc_status.json` 里记 `not_done` 加原因。

**不用代理指标冒充校正结果。** `simple` 模式只报血红蛋白基因的
背景水平，并明确标注 `heuristic_only` + "这不是 ambient RNA 校正"。

---

## 2. 标准化与高变基因

### 标准化

`normalize_total(target_sum=1e4)` + `log1p`。

`target_sum` 的选择会影响下游的差异表达幅度，但对聚类结构影响很小。
1e4 是 scanpy 的惯例值。

### 高变基因：seurat_v3

Stuart T, et al. "Comprehensive Integration of Single-Cell Data."
*Cell* 2019. doi:10.1016/j.cell.2019.05.031

`seurat_v3` 用负二项分布拟合基因的均值-方差关系，取标准化方差最高的
基因。**它吃原始计数，必须在 normalize 之前算。**

另两种口味（`seurat`、`cell_ranger`）相反，吃 log 后的数据。
**顺序反了不报错，只是给出错的 HVG。**

依赖 `scikit-misc` 的 loess 拟合。缺它时降级到 `seurat` 口味，
并记录 `hvg_flavor_requested` vs `hvg_flavor_used`。

---

## 3. 降维、整合、聚类

### PCA

`arpack` 求解器，`random_state` 固定。

### 批次整合

- `harmony`：Korsunsky I, et al. "Fast, sensitive and accurate integration
  of single-cell data with Harmony." *Nature Methods* 2019.
  doi:10.1038/s41592-019-0619-0
- `combat`：Johnson WE, Li C, Rabinovic A. "Adjusting batch effects in
  microarray expression data using empirical Bayes methods."
  *Biostatistics* 2007. doi:10.1093/biostatistics/kxj037

**单样本数据不整合是正确的选择，不是偷懒。** 多供体数据必须整合，
否则聚类被供体差异主导 —— 那是技术差异，不是细胞类型差异。

`integration_status.json` 会区分 `not_applied`（不需要）与
`not_configured`（需要但没配），两者的含义完全不同。

### Leiden 聚类

Traag VA, Waltman L, van Eck NJ. "From Louvain to Leiden: guaranteeing
well-connected communities." *Scientific Reports* 2019.
doi:10.1038/s41598-019-41695-z

用 `flavor="igraph"`（scanpy 1.10+ 的推荐路径）。

**分辨率决定簇的粒度，而下游所有分析都建立在这个粒度上。**
本流水线额外跑一个分辨率扫描（0.2→2.0）并把"多少簇随分辨率怎么变"
落盘，让"选 1.0"变成有依据的决定。

---

## 4. 细胞类型注释

### 打分法及其固有局限

`sc.tl.score_genes`：对每个细胞的签名基因做相对表达打分（减去表达
水平匹配的随机对照基因集），再按簇取均值。

**这不是注释，是打分提示。** 局限：
1. 签名不全 —— 稀有类型拿低分后被归到"最接近的"，而不是"未知"
2. 组织特异 —— 内置签名是外周血/免疫为主的
3. **状态与类型混淆** —— 激活/耗竭的 T 细胞仍拿高 T 细胞分
4. **签名重叠** —— 实测 pbmc3k 里 `T_cell` 与 `CD4_T` 的 marker
   完全重叠（IL7R、LTB 两边都有），打分差 **0.0**

所以每个 assignment 都带 `score_margin` 与 `assignment_confident`。
**margin 小的 assignment 不该被当成结论** —— 只报类型名等于把
不确定性藏起来。

### 第二条独立证据：CellTypist（文档 §2.4 点名，**本流水线已实现**）

marker 打分是**无监督**的：它只知道"这个簇里这几个基因高"。
CellTypist 用预训练的 logistic 回归模型（免疫图谱）给每个细胞打分，
是**独立的第二条证据**。两者不一致的簇才值得人工看。

- 模型：`Immune_All_Low.pkl`（可在 config 里改 `analysis.celltypist_model`）
- 输入：**`adata.raw` 全基因集**，不是 HVG 子集 —— 模型依赖的大部分基因
  不在 HVG 里，喂 HVG 会让标签退化成噪声
- 一致性的量化在**簇层面**（marker 给每簇一个标签，CellTypist 给每个细胞
  一个标签，先按簇取众数再比），并报每簇的 `celltypist_purity`
- 跑不成时如实记 `package_missing` / `model_unavailable` / `failed`，
  **marker 打分仍然是主结果**，不会因为第二条证据缺失而消失

**一个踩过的坑（值得记下来）：** 第一版把 `celltypist.models_path` 当成
`pathlib.Path` 用了（`models_path / model_name`），而它其实是 **`str`**
（`models.py:19` 是 `os.path.join(...)`）—— 抛 `TypeError`。
那个 `TypeError` 被笼统的 `except Exception` 接住，于是状态里写成
**「模型拿不到 —— CI 可能无外网」**。**一个纯本地代码 bug 被记成了网络问题**，
下一个人会去查 runner 的出网策略。现在失败被拆成三类分别记原因：
路径构造失败 → `failed`；下载抛异常 / 下载后文件仍不存在 →
`model_unavailable`。

### 更严格的做法（本流水线未实现）

- 有参考数据集时用 **label transfer**（Seurat 的 `FindTransferAnchors`
  或 scanpy 的 `sc.tl.ingest`）
- 用 **scANVI** 这类需要自己训练的有监督模型（CellTypist 用的是现成模型）
- 人工核对每个簇的 marker 后再定名

---

## 5. 拟bulk 差异表达

### 为什么不能在细胞层面检验

这是单细胞差异分析里最严重的系统性错误。

把每个细胞当一个独立样本，是把"细胞数"当成了"样本量"。
10000 个细胞来自 3 个供体，自由度是 **2** 不是 9997。
检验会给出 `p = 1e-300` 级别的"极显著"结果 —— **假阳性看起来
比真信号还显著**。

### 拟bulk + DESeq2

Love MI, Huber W, Anders S. "Moderated estimation of fold change and
dispersion for RNA-seq data with DESeq2." *Genome Biology* 2014.
doi:10.1186/s13059-014-0550-8

实现用 `pydeseq2`：Muzellec B, Teleńczuk M, Cabeli V, Andreux M.
"PyDESeq2: a python package for bulk RNA-seq differential expression
analysis." *Bioinformatics* 2023. doi:10.1093/bioinformatics/btad547

流程：按 (样本 × 分组 × 细胞类型) 加总**原始计数** → 负二项 GLM。

**必须用原始计数，不能用 normalize/log 后的 X。** 把 log 值喂进
负二项模型不会报错，只会给出错的离散度估计。

### 门槛

- 每个 (样本, 细胞类型) 组合至少 10 个细胞 —— 更少的话拟bulk 是噪声
- 每组至少 2 个生物学重复 —— 1 个重复算不出组内方差
- 不满足的细胞类型**记入 `skipped_celltypes` 并说明原因**

---

## 6. 轨迹推断

### 为什么必须至少两种方法

单一方法的拟时序是一个**一维坐标**，把高维状态压成一条线；不同算法
压出来的线可以完全不同。实测（PBMC3k，2652 细胞）四种方法的**原始**
拟时序两两 Spearman 相关：

| | dpt | palantir | scFates | CytoTRACE |
|---|---|---|---|---|
| dpt | 1.000 | **−0.389** | **+0.855** | +0.250 |
| palantir | −0.389 | 1.000 | **−0.500** | **−0.467** |
| scFates | +0.855 | −0.500 | 1.000 | +0.534 |
| CytoTRACE | +0.250 | −0.467 | +0.534 | 1.000 |

**符号都不一样。** 不做交叉验证就报"轨迹"，报的是某个算法的一次输出，
不是数据里的结构。

### 四种方法

| 方法 | 实现 | 原始约定 |
|---|---|---|
| DPT | scanpy `tl.diffmap` + `tl.dpt` | 0 = 根 = 早 |
| Palantir | `palantir.core.run_palantir`（扩散图 + Markov 链） | 0 = 起始细胞 = 早 |
| scFates | `scFates.tl.tree/pseudotime`（主曲线树，Slingshot 的 Python 移植） | 0 = 根 = 早 |
| CytoTRACE | **本仓库自行实现**的 GCS 统计量 | 1 = 分化潜能高 = 早 |

Wolf FA, Hamey FK, Plass M, et al. "PAGA: graph abstraction reconciles
clustering with trajectory inference through a topology preserving map of
single cells." *Genome Biology* 2019. doi:10.1186/s13059-019-1663-x

Haghverdi L, Büttner M, Wolf FA, Buettner F, Theis FJ. "Diffusion
pseudotime robustly reconstructs lineage branching."
*Nature Methods* 2016. doi:10.1038/nmeth.3971

Setty M, Kiseliovas V, Levine J, et al. "Characterization of cell fate
probabilities in single-cell data with Palantir."
*Nature Biotechnology* 2019. doi:10.1038/s41587-019-0068-4

Street K, Risso D, Fletcher RB, et al. "Slingshot: cell lineage and
pseudotime inference for single-cell transcriptomics."
*BMC Genomics* 2018. doi:10.1186/s12864-018-4772-0
（Python 端用 scFates 实现，Faure L, Soldatov R, Kharchenko PV,
Adameyko I. *Bioinformatics* 2023. doi:10.1093/bioinformatics/btac746）

Gulati GS, Sikandar SS, Wesche DJ, et al. "Single-cell transcriptional
diversity is a hallmark of developmental potential."
*Science* 2020. doi:10.1126/science.aax0249
（CytoTRACE 原论文。**本仓库实现的是它的核心 GCS 统计量，不是该包** ——
CytoTRACE2 不在 PyPI 上）

### 方向校正：拟时序的符号是任意的

DPT 从根出发，根选在早期还是晚期，整条轴就反过来。所以本流水线：

1. 每个方法各自算出拟时序（原始符号，不假设谁对）
2. 用一个**方向参考**判断每个方法是否需要翻转
3. 统一约定成「**值越大越晚**」后再互相比较

方向参考优先级：配置的 `trajectory.early_markers` / `late_markers`
（生物学判据）→ 退回 CytoTRACE 分化潜能分（**记进 `direction_source`**）。

**方向参考不能算进一致性统计。** 退回模式下参考就是 CytoTRACE 本身，
它与参考的相关恒为 ±1 —— 那是定义不是证据。第一版把它算进去了，
一致性被抬到 0.587；排除后是 **0.636**（DPT/Palantir/scFates 三方），
这才是真实的一致性。

### 选根：从 PAGA 连通度改成 CytoTRACE

原来的自动选根是"PAGA 连通度最高的簇"。在 PBMC3k 上它挑中了**浆细胞**
（终末分化）当根 —— 连通度高只说明"在图上居中"，而居中既可能是
祖细胞也可能是终末态。

改用 CytoTRACE GCS 最高的细胞作根后，DPT 与 GCS 的相关从 **+0.25
升到 +0.59**。

**但这仍然是统计判据，不是生物学判据。** 要下方向性结论仍需人工用
已知的早期/晚期 marker 复核。

### 沿轨迹的基因与模块

- `trajectory_genes.csv`：与共识拟时序的 Spearman 相关（表达基因全集）
- `trajectory_modules.csv`：按拟时序分箱 → 每个基因的平滑表达曲线 →
  行 z-score → k-means 聚成 6 个模块。**行 z-score 是为了让模块反映
  "形状"而不是"表达量"**，否则高表达基因会主导聚类
- 分支点：scFates 的 `seg` 片段 + 每片段的拟时序中位数

### scFates 的两个坑（实测）

1. **`tl.pseudotime()` 之前必须先 `tl.cleanup()`。** 官方顺序是
   `diffusion → tree → cleanup → root → pseudotime`。漏掉 cleanup 时
   `pseudotime()` 必定崩在
   `pd.Series(uns['graph']['milestones']) == t][0]` 的 IndexError ——
   因为 `map_cells` 读 milestones 时它还没被写入（scFates 包内的
   `pseudotime.py` 第 254 行才写入，而调用在第 87 行）。
2. **`tl.test_association` / `tl.test_fork` 需要 rpy2 + R + mgcv。**
   这与"云端不装 R"的设计冲突，所以**没有用**这两个函数 ——
   沿轨迹的基因分析是本仓库自己实现的。
   另外 `method='epg'` 与当前 networkx 不兼容
   （`add_edges_from` 收到 int 就 TypeError），所以用 `method='ppt'`。

### 适用范围（关键）

- PAGA 给出的是**簇间连通性**，不是分化方向
- 拟时序是**一维坐标**，分支过程（一个祖先进两种细胞）会被压成
  "先后"，而实际是"并列"
- **没有 RNA 速率（spliced/unspliced）时不能下方向性结论。**
  本流水线只有计数矩阵，`scvelo` 记 `not_done`。要速率需要从
  Cell Ranger 的 `velocyto` 或 `--include-intrins` 输出重新开始
- **方法间一致性只是内部一致性。** 几种方法都错向同一个伪轨迹时，
  它们依然彼此高度相关。**一致不等于正确**
- `regulon × 拟时序` 的显著性在 n 大时很廉价：本数据 n=2652，
  |rho| 只要约 0.06 就能过 BH<0.05（201 个里 190 个"显著"）。
  要看效应量分布（本数据 |rho| 中位 0.316，48 个 >0.5）

---

## 7. 细胞通讯

### 本流水线的方法

数据库驱动的配体-受体共表达打分：
`配体在发送簇的平均表达 × 受体在接收簇的平均表达`，
配合簇标签置换检验（200 次）+ BH 多重检验校正。

### 主工具：LIANA（文档 §2.7 点名，**本流水线已跑**）

- **LIANA**：Dimitrov D, et al. "Comparison of methods and resources for
  cell–cell communication inference from single-cell RNA-Seq data."
  *Nature Communications* 2022. doi:10.1038/s41467-022-30755-0

LIANA 把多个通讯打分方法（CellPhoneDB / NATMI / Connectome / logFC 等）
聚合成一个 **consensus rank aggregate**。本流水线跑它的
`rank_aggregate`，用的是**全基因集**（`adata.raw`），与自建打分同一套
簇标签。**两者用同一份输入，所以差异只来自方法本身** —— 这才是有信息量的
对照。产物里报共同组合数、Spearman rho、top25 重叠。

**与自建打分的差异必须报，不能只报"跑通了"。** 实测 rho 并不高
（≈0.13），top25 却重叠 24/25 —— 这说明**头部一致、中段排序差异大**。
只报一个数会把这两件事混成一件。

### 其他标准工具（本流水线未用）

- **CellChat**：Jin S, et al. "Inference and analysis of cell-cell
  communication using CellChat." *Nature Communications* 2021.
  doi:10.1038/s41467-021-21246-9 —— R 包，CI 不装 R
- **CellPhoneDB v2**：Efremova M, Vento-Tormo M, Teichmann SA, Vento-Tormo R.
  "CellPhoneDB: inferring cell–cell communication from combined expression
  of multi-subunit ligand–receptor complexes." *Nature Protocols* 2020.
  doi:10.1038/s41596-020-0292-x —— LIANA 内部已调用它的打分方法

### 适用范围（关键）

1. **共表达不等于通讯。** 没有空间信息时，只能说两类细胞**分别**
   表达了配体和受体，不代表它们在组织里相邻
2. 表达量是稳态丰度，**不等于蛋白水平，也不等于分泌量**
3. **自建打分是启发式，不是 consensus rank aggregate** —— 所以主结论
   以 LIANA 为准，自建打分作为可复算的对照
4. **内置库只有 38 对**（免疫为主），远少于 CellChatDB 的数千对。
   LIANA 用的是它自带的 `consensus` 资源，对数多得多；
   两个来源的"没找到显著通讯"都是**假阴性**，不是真的没有通讯

---

## 8. 转录因子调控

### 本流水线的方法

Pearson 相关取每个 TF 的 top 30 共表达靶基因 → AUCell 式活性打分
（scanpy `score_genes`，校正表达水平与测序深度）。

### SCENIC（本流水线未实现）

Van de Sande B, et al. "A scalable SCENIC workflow for single-cell gene
regulatory network analysis." *Nature Protocols* 2020.
doi:10.1038/s41596-020-0336-2

Aibar S, et al. "SCENIC: single-cell regulatory network inference and
clustering." *Nature Methods* 2017. doi:10.1038/nmeth.4463

SCENIC 的流程是 GRNBoost2/GENIE3（梯度提升/随机森林）推断 TF-靶基因
关系，**然后用 cisTarget 做 motif 富集剪枝** —— 只有靶基因启动子区
真的有该 TF 的结合 motif 才保留。

**剪枝是 SCENIC 的关键一步**，它把"共表达"收紧成"可能有直接调控"。
本流水线**没有 motif 剪枝**（需要 cisTarget 的排名数据库，几百 MB），
所以报的是**共表达模块，用 TF 命名而已**。

### 适用范围（关键）

一个 TF 和一组基因共表达，可能是它调控它们，也可能是：
- 它们被同一个上游因子调控
- **它们只是同一种细胞类型的标志物（最常见）**
- 拷贝数变异导致的共表达

所以每个调控子报 `cluster_specificity` 与 `frac_cells_expressing_tf` ——
让"这个调控子是不是只是细胞类型的代理"可以被看出来。

---

## 9. 参考的综述

- Luecken MD, Theis FJ. "Current best practices in single-cell RNA-seq
  analysis: a tutorial." *Molecular Systems Biology* 2019.
  doi:10.15252/msb.20188746
- Heumos L, et al. "Best practices for single-cell analysis across
  modalities." *Nature Reviews Genetics* 2023.
  doi:10.1038/s41576-023-00586-w
- Squair JW, et al. "Confronting false discoveries in single-cell
  differential expression." *Nature Communications* 2021.
  doi:10.1038/s41467-021-25960-2
  （**拟bulk 那一节的直接依据** —— 该文系统展示了细胞层面检验的
  假阳性问题）
