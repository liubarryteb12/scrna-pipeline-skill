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

### 更严格的做法（本流水线未实现）

- 有参考数据集时用 **label transfer**（Seurat 的 `FindTransferAnchors`
  或 scanpy 的 `sc.tl.ingest`）
- 用 **scANVI / CellTypist** 这类有监督模型
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

### PAGA

Wolf FA, Hamey FK, Plass M, et al. "PAGA: graph abstraction reconciles
clustering with trajectory inference through a topology preserving map of
single cells." *Genome Biology* 2019. doi:10.1186/s13059-019-1663-x

### 扩散拟时序（DPT）

Haghverdi L, Büttner M, Wolf FA, Buettner F, Theis FJ. "Diffusion
pseudotime robustly reconstructs lineage branching."
*Nature Methods* 2016. doi:10.1038/nmeth.3971

### 适用范围（关键）

- PAGA 给出的是**簇间连通性**，不是分化方向
- DPT 的方向**完全依赖根的选择**。根选错则整条轨迹反向。
  本流水线的自动选根是"连通度最高"，那是**启发式**：
  连通度高只说明它在图上居中，而居中既可能是祖细胞，也可能是
  被各种中间态包围的终末态
- 拟时序是**一维坐标**，分支过程（一个祖先进两种细胞）会被压成
  "先后"，而实际是"并列"
- **没有 RNA 速率（spliced/unspliced）时不能下方向性结论。**
  本流水线只有计数矩阵，所以报的是相似度排序，不是分化方向。
  要方向性需要 `scvelo`（需 spliced/unspliced 定量）或 `CellRank`

---

## 7. 细胞通讯

### 本流水线的方法

数据库驱动的配体-受体共表达打分：
`配体在发送簇的平均表达 × 受体在接收簇的平均表达`，
配合簇标签置换检验（200 次）+ BH 多重检验校正。

### 标准工具（本流水线未用）

- **LIANA**：Dimitrov D, et al. "Comparison of methods and resources for
  cell–cell communication inference from single-cell RNA-Seq data."
  *Nature Communications* 2022. doi:10.1038/s41467-022-30755-0
- **CellChat**：Jin S, et al. "Inference and analysis of cell-cell
  communication using CellChat." *Nature Communications* 2021.
  doi:10.1038/s41467-021-21246-9
- **CellPhoneDB v2**：Efremova M, Vento-Tormo M, Teichmann SA, Vento-Tormo R.
  "CellPhoneDB: inferring cell–cell communication from combined expression
  of multi-subunit ligand–receptor complexes." *Nature Protocols* 2020.
  doi:10.1038/s41596-020-0292-x

### 适用范围（关键）

1. **共表达不等于通讯。** 没有空间信息时，只能说两类细胞**分别**
   表达了配体和受体，不代表它们在组织里相邻
2. 表达量是稳态丰度，**不等于蛋白水平，也不等于分泌量**
3. 打分是启发式，**不是 LIANA 的 consensus rank aggregate**
4. **内置库只有 38 对**（免疫为主），远少于 CellChatDB 的数千对。
   覆盖不全时"没找到显著通讯"是**假阴性**，不是真的没有通讯

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
