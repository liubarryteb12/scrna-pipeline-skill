# 线2 / 单细胞转录组 + 虚拟敲除线

> 术语见 `GLOSSARY.md`，隔离约束见 `docs/isolation.md`。
> 本目录只放线2 代码。线1、线3 代码不得出现在这里，本目录脚本也不得引用另外两条线。

## 输入

已上数数据 GSE254707（BJD 2024 银屑病表皮 pDC 研究）：**Smart-seq2 孔板单细胞**，
不是 10x 目录结构。算力节点 `raw_data/` 三件：

- `GSE254707_RAW.tar`：4988 个按细胞的 `GSM*_*.txt.gz`（外层纯 tar 不压，内层 gz），
  每个文件 = 基因符号 + 计数两列（CRLF 行尾，末尾 5 行 featureCounts `__` 汇总行须排除）；
- `GPL21697_family.soft.gz`：平台注释，`^SAMPLE` 块带每细胞元数据
  （donor / diagnosis / region / class / nCount / nFeature / percent.ercc / 作者 celltype）；
- `grouping.txt`：组名清单（LS / NLS / HC），非样本到组映射。

区域码惯例：PP = 皮损（LS，2161 细胞）、PN = 非皮损（NLS，934）、H = 健康对照（HC，1893）；
8 个供体（H061/H062/H066 健康，P185-P190 银屑病）；207 个 doublet（作者用于 CIM-seq，
`class` 字段标注，celltype_paper=NA）。步骤 1 `ingest_counts` 负责把上述形态统一为
cell x gene 稀疏矩阵（`scrna_counts.h5ad`）。

## 线内步骤（顺序）

1. 读入与格式统一
2. QC（按细胞与基因双阈值过滤）
3. 归一化 / 高变基因
4. 去批次 / 整合
5. 降维与聚类
6. 细胞类型注释
7. 拟时序
8. 细胞间通讯
9. 虚拟敲除（末端模块）
10. DR genes 汇总输出

## 工具分工（线内混用）

| 语言 | 负责 |
|---|---|
| Python | Scanpy 流程、拟时序（scanpy DPT，Monocle3 运行环境缺位）、细胞间通讯（liana consensus，运行环境 R 无 CellChat/CellPhoneDB）、虚拟敲除（scTenifoldpy 0.5.1 tenifoldknk，DR genes） |
| R | Seurat 侧注释（如后续启用） |

主语言与版本、随机种子记到 `envs/manifest.txt`，分析中途不更换。

## 虚拟敲除约束

- 虚拟敲除在本线末端执行，输入必须是本线 `results/` 内的中间产物。
- 产出叫 **DR genes**，不是 DEG，不得与差异分析概念混用（见 `GLOSSARY.md`）。
- 敲除对象来自本线注释或本线 hub 基因，**不得引用线1 的 GEO 交集基因作输入**。

## 产出去向

算力节点 `<BIO_DATA_ROOT>/projects/<项目>/line2_scrna/{raw_data,processed,results,figures}/`
（脚本 `bin/line2_scrna/` 与环境声明 `envs/line2_scrna/` 全局一份，所有项目共用）。
任务台账 JSON 落 `<BIO_DATA_ROOT>/system_logs/`，`line_id` = `line2_scrna`，`project_id` = 所属项目。
**产物回搬（用户拍板 2026-10-13）**：processed/results 产物与 system_logs 台账
必须通过 `tools/ops/07_pull_outputs.ps1` 回搬至本地 `outputs/` 目录存档（ADB 通道 +
逐文件字节校验，同名默认跳过、-Force 覆盖），算力节点保留运行态正本。

## 本目录结构（本地代码本体，与算力节点侧数据目录一一对应）

本地只放代码；数据、产物、环境全部在算力节点 `<BIO_DATA_ROOT>/projects/<项目>/line2_scrna/`（代码与环境在顶层全局一份）。

```
src/line2_scrna/
├── README.md            # 本文件
├── scripts/             # 线内分析脚本（.R / .py 可混用，但只在本线内）
└── checkpoint/          # 桥上检查点脚本，命名 checkpoint_*.R / checkpoint_*.py
```

桥接检查点（Python 出 `.h5ad` 给 R 读之类）必须校验四项：基因名、矩阵方向、细胞名、embedding。
约束与校验规则见 `docs/isolation.md`。

## 状态

步骤 1 `ingest_counts` 已落地（2026-10-13）：4988 个按细胞计数文件 + GPL soft 元数据 →
cell x gene 稀疏矩阵 `scrna_counts.h5ad`（33213 基因，含 102 个 ERCC；nnz 1824 万，密度 0.11）。
自校验全过：GSM 集合 tar 与 soft 一致、基因表跨文件一致、`__` 汇总行不泄漏、
区域码与诊断自洽、与 soft 元数据对账 1% 容差内（782/4988 细胞三项精确一致，
残差为作者再定量差异，如实报告不静默，见踩坑台账 #103）。界面验收 exit 0、台账 completed/passed。

步骤 2 `qc_cells` 已落地（2026-10-13）：论文固定阈值 nCount_RNA>50000 / nFeature_RNA>1500 /
ERCC<15%（不随结果调整）+ 基因过滤（保留细胞中表达细胞数 >=3）。实测：前两项 0 失败
（deposited 数据 min 53020/1525 恰在阈值上，作者已预过滤），ERCC 超限 216 个为实际剔除项；
**207 个 doublet 处置拍板：排除出主分析矩阵**（CIM-seq 混合图谱、celltype_paper=NA，
混入下游会制造虚假中间态），全部 4988 个细胞指标与去留原因落 `qc_cell_metrics.tsv` 可审计。
产物 `scrna_qc.h5ad`（4587 细胞 × 25063 基因，LS 1913 / NLS 834 / HC 1840，8 供体全保留）+
`qc_report.txt`。界面验收 exit 0、台账 completed/passed。

步骤 3 `normalize_hvg` 已落地（2026-10-13）：ERCC spike-in 剔除 94 个（在归一化前，
防扭曲库大小与 HVG 选择）→ scanpy normalize_total(1e4) + log1p（与 Seurat LogNormalize
数值等价，原始整数计数存 layers['counts']）→ HVG seurat flavor n_top_genes=2000
（dispersion 口径，全数据集选取未分 batch）。自校验：归一化后库大小全部精确 1e4、
X 全有限非负、HVG 精确 2000 且 0 ERCC 混入、obs 19 列不丢。top HVG 为
CCL22/IL17A/MMP12/TNF/CXCL13 等银屑病免疫/角化 marker，生物学合理。
产物 `scrna_norm.h5ad`（4587 细胞 × 24969 基因）+ `hvg_genes.tsv` + `norm_report.txt`。
界面验收 exit 0、台账 completed/passed。

环境声明已修正：`envs/manifest.txt` 原误抄 line1 的 R 4.3.3，本线主线是 Scanpy 流程
（README 工具分工），修为 Python 3.12.3（运行环境实测值）；line2 无历史任务，不存在中途换语言问题。

步骤 4 `integrate_harmony` 已落地（2026-10-13）：HVG 子集（2000 基因，沿用步骤 3 选择）→
scale(max_value=10) → PCA(50) → Harmony 校正（batch = donor 8 供体，harmonypy 原生调用，
random_state 用引擎注入 seed；scanpy 1.12.4 的 `sce.pp.harmony_integrate` 包装器对 Z_corr
做 .T 与运行环境 harmonypy 形状不兼容，已绕开，见踩坑台账 #104）。实测：PC1 方差比 9.77%、
cum30 36.36%；donor silhouette 前 -0.1798 / 后 -0.1917（前后基本持平且为负，说明 donor
不是主导方差来源，如实报告不调参凑数）；8 供体、组计数（LS 1913/NLS 834/HC 1840）不变。
**配对设计混淆警告已写入报告**：healthy 供体只贡献 HC、银屑病供体贡献 LS+NLS，
batch=donor 与 group 部分共线，harmony 校正会顺带削去部分组间信号，解释时注意。
产物 `scrna_integrated.h5ad`（4587 细胞 × 2000 HVG，X=scaled、layers['counts'] 原始整数、
obsm 含 X_pca 与 X_pca_harmony 双嵌入）+ `integrate_report.txt`。界面验收 exit 0、
台账 completed/passed。

步骤 5 `cluster_umap` 已落地（2026-10-13）：在 X_pca_harmony（harmony 校正后嵌入）
上 neighbors(n_neighbors=15, n_pcs=30) → UMAP（random_state=引擎 seed）→
Leiden（igraph flavor、resolution=1.0，口径固定不调参）。实测 21 个簇
（size 17-684），簇 x group 交叉显示 LS 主导簇（0/1/2/6/9 等）与 HC 主导簇
（10/11/12/14 等）分离清晰、混合簇居中；簇 x donor 交叉显示无大簇被单一供体
独占（>=100 细胞簇单一供体占比 <=50% 自校验，防 harmony 失效回退；小簇不参与，
样本波动会误伤，坑 #105）。产物 `scrna_clustered.h5ad`（4587 × 2000，
obsm 三嵌入 X_pca/X_pca_harmony/X_umap、obsp 邻居图、obs['leiden']）+ `cluster_report.txt`
（簇大小 + 簇 x group + 簇 x donor crosstab）。界面验收 exit 0、台账 completed/passed。

步骤 6 `annotate_celltypes` 已落地（2026-10-13）：先出证据后下结论，对 21 个簇做
wilcoxon vs 其余簇 top20 marker 基因，再与 16 类文献常见皮肤/免疫 marker 集
（大类粒度，14 类+unknown，不比作者注释细）比对打分（hits/min(5,len)，阈值 0.4，
不过记 unknown 不硬猜），并与作者注释 27 个细分类映射大类做细胞级对照审计
（总体一致率 0.8129，不一致簇逐簇列出供人工复核，不静默抹平）。
实测 21 簇注释：簇 19=pDC 20/20（marker TCF4/IL3RA/IRF8/GZMB，正中论文主题）、
簇 7=melanocyte 7/7 marker（387/389）、簇 15=progenitor_keratinocyte 56/56、
簇 17=channel_keratinocyte 58/62、簇 4=T_cell（作者 Treg 93/104，粒度差异已注明）；
簇 12（167 细胞，作者 Activated Transition keratinocytes 161/167）与簇 18
（122 细胞）marker 打分不过阈值记 unknown，不硬猜。
log2fc inf 420 处为组内全零表达基因的 wilcoxon 已知行为，置 NaN 不进结论。
产物 `scrna_annotated.h5ad`（4587 × 2000，obs 加 celltype / celltype_paper_major 两列、
uns['annotation'] 含方法/marker 集/映射表/簇注释/一致率、uns['cluster_markers'] DataFrame
420 行=21 簇 x top20）+ `cluster_markers.tsv` + `annotate_report.txt`。
界面验收 exit 0（台账 a69f2a4b completed/passed），files 面板 15 产物可见，巡检 30 条全合法。

步骤 7 `pseudotime` 已落地（2026-10-13）：拟时序路线因运行环境 R 侧 monocle3/monocle/
slingshot/tradeSeq 全缺位、cellrank 依赖链过重，**选 scanpy 内置 DPT（零新依赖）**：
在 harmony 邻居图上 diffmap(n_comps=15) → 根选 progenitor_keratinocyte 簇内
DC 空间离质心最远细胞为 iroot → sc.tl.dpt(n_dcs=10, n_branchings=0)。
**不可达细胞 dpt=inf 是 scanpy 正常输出不是 bug**：melanocyte/macrophage_mDC/
Langerhans_cell/pDC 四类非表皮谱系细胞（773 个）与根在扩散图上不连通，
自校验用豁免集管理（谱系细胞出现 inf 才 fail），报告注明免疫细胞 dpt
无生物学拟时序意义、只解读表皮谱系。实测主链单调合理：progenitor（median 0.131）
→ suprabasal/activated/mitotic/basal（0.586-0.614）→ channel（0.945）。
产物 `scrna_pseudotime.h5ad`（4587 × 2000，obs['dpt_pseudotime']、obsm['X_diffmap']、
uns['iroot']、uns['pseudotime'] 含 per_celltype/per_cluster 统计）+ `pseudotime_report.txt`。
界面验收 exit 0（台账 423da81f completed/passed），files 面板 17 产物可见，巡检 31 条全合法。
**工具分工表更新**：拟时序从 R 侧 Monocle3 改为 Python 侧 scanpy DPT。

步骤 8 `cell_communication` 已落地（2026-10-13）：装 liana 1.10.0 走
consensus 资源聚合（4620 对 LR，含 CellPhoneDB/CellTalkDB/ConnectomeDB2020
等汇总），`li.mt.rank_aggregate` groupby=celltype（13 大类）、expr_prop=0.1、
min_cells=5、n_perms=1000、seed 注入，RRA 聚合。
**输入矩阵口径（原型实测定稿）**：表达矩阵取步骤 3 norm 全基因产物（24969
基因）而非步骤 7 的 2000 HVG 矩阵，consensus 的 2016 个 LR 基因在全基因
矩阵覆盖 81.3%，在 HVG 只覆盖 21.9%，HVG 按全局方差选取与 LR 相关性弱，
用它做通讯分析会系统性漏掉 78% 的配体受体；标签取步骤 6 annotated 的
obs['celltype']，按索引对齐。两个 rank 分工写进报告：specificity_rank 衡量
LR 对相对其他组的特异性（生物学结论以它为准），magnitude_rank 衡量表达强度
（B2M-CD1A/VIM-CD44/DSG1-DSC3 等管家 pair 会霸占 magnitude 榜但特异性低）。
unknown 类（289 细胞，混合类）保留在结果中但不作生物学结论。
实测 top：LC→T_cell HLA-DQB2-LAG3（抗原呈递）、NK_ILC/T_cell→LC
CSF2-CSF1R、channel→progenitor RSPO1-LGR5 与 FGF7-FGFR1（表皮干细胞龛）、
pDC→pDC BST2-LILRA4（自分泌）、macrophage_mDC→pDC CD14-TLR9，均符合
皮肤免疫/银屑病文献常识。n_perms=100 与 1000 的 top20 完全一致（原型验证）。
产物 `scrna_cellcomm.h5ad`（4587 × 2000，obs 全列保留、uns['cell_communication']
含口径参数/覆盖率/top20/per-pair top1 两个 DataFrame）+ `cellcomm_interactions.tsv`
（57027 行=13x13 对 x 过滤后 LR，含自分泌对 4953 个）+ `cellcomm_report.txt`。
界面验收 exit 0（台账 8a96a231 completed/passed），files 面板 20 产物可见，
巡检 32 条全合法。**工具分工表更新**：细胞间通讯从 R 侧 CellChat 改为
Python 侧 liana（运行环境 R 无 CellChat/CellPhoneDB）。
**连带修复（坑 #107）**：装 liana 时 pip 把 pandas 3.0.6 降级到 2.3.3，
anndata 0.13 写 obs 字符串列因 allow_write_nullable_strings 开关未定而崩，
7 个写 h5ad 脚本统一显式关掉该开关，annotate_celltypes 回归重跑产物与
降级前逐字节一致（一致率 0.8129 不变）。

步骤 9 `virtual_knockout` 已落地（2026-10-13）：装 scTenifoldpy 0.5.1（安装前
dry-run 证实只新增 typer/rich/requests 等包，numpy/pandas/scanpy/anndata 核心栈
未动）走 tenifoldknk 路线。**敲除口径（原型实测定稿）**：群体取 Psoriasis
2747 细胞（避开步骤 4 donor-group 共线警告；pDC 仅 20 个细胞不够建网），
敲除基因取 TNF（本线 top HVG、Psoriasis 表达占比 24.5% 稳过 QC 闸；pDC marker
TCF4 表达 6.55% 贴闸，原型对比后选 TNF），表达矩阵取 HVG 2000 上的 counts。
**API 结构（源码级确认）**：`scTenifoldKnk` 类不接受顶层 random_state，必须走
`nc_kws={"random_state": seed}` + `td_kws={"random_state": seed}`；`d_regulation`
步骤自动把敲除基因排除出期望计算（与 R 一致）。管线 QC（sc_QC+cpm_norm）剔
216 细胞/775 基因后建 1225 基因网络，adj p<0.05 的 4 个非敲除 DR 基因
TMEM45A/S100A9/S100A2/S100A8（角化 + 警报素，TNF 敲除的银屑病下游讲得通）；
换 seed top20 重合 19/20。产物 `ko_dr_genes.tsv`（DR 全表 + is_knockout_gene
列）、`scrna_knockout.h5ad`（2747 × 2000，obs['tenifold_qc_kept'] 2531/2747、
var 六列 tenifold_*、uns['virtual_knockout'] 含完整口径与 warnings）、
`knockout_report.txt`（口径说明 5 条：DR≠DEG、敲除基因排名解释、HVG 基因空间
边界、计算模拟非实验验证、种子稳定性）。界面验收 exit 0（台账 5db1ad45
completed/passed），files 面板三产物可见，巡检 33 条全合法。
**工具分工表更新**：虚拟敲除走 Python 侧 scTenifoldpy 0.5.1。

步骤 10 `summarize_dr` 已落地（2026-10-13）：线内末端汇总模块，读步骤 9 两产物
（*knockout*.h5ad 的 uns['virtual_knockout'] 口径 + uns['annotation']['marker_sets']
16 类 + obs['celltype'] + layers['counts']，*ko_dr*.tsv DR 全表），不重算敲除、
不改上游口径。汇总口径：显著 DR 基因 = adj p<0.05 且非敲除基因，按 p-value 升序
rank；每基因输出 direction（FC>1 up / FC<1 down）、marker_celltypes（命中步骤 6
marker 集的细胞类型）、top_expr_celltype / top_expr_prop（counts>0 细胞占比最高
的类，unknown 混合类排除）。自校验：tsv 列名/敲除基因唯一/adj p 单调非降、
uns 口径键完整、显著数与 uns n_significant_excluding_ko 一致、rank 连续、
汇总基因集合 == 显著集、**uns seed 与引擎注入 seed 必须一致**（汇总复述的是
上游那一次敲除，种子对不上就是两份口径）。实测：1225 DR 基因中显著非敲除
4 个（TMEM45A/S100A9/S100A2/S100A8，全 up，3/4 命中 activated_keratinocyte
marker，top_expr 全为角化类 0.972-1.000，与步骤 9 结论一致），敲除基因 TNF
按 p-value 排名第 1。产物 `dr_genes_summary.tsv`（12 列）+ `dr_summary_report.txt`
（口径复述 + 统计 + 明细 + 5 条局限说明）。界面验收 exit 0（台账 48bb85cb
completed/passed），files 面板两产物可见，巡检 34 条全合法。**线2 十步清单
至此全部完成**（步骤 1-10，界面验收全 exit 0、台账全 completed/passed）。

---

## aims 线2 下游清单补齐（2026-10-13 晚，步骤 14-18）

用户拍板"节点能做的做完、图要能目视审核"后，补齐 aims 线2 下游十项里缺的
5 项（原步骤 1-10 只覆盖了 QC/注释/拟时序/通讯/虚拟敲除）。

**步骤 14 `scrna_composition`（细胞类型比例/组成分析）**：比例单位 = donor
（每 donor 内 13 类占比谱，行和=1）。组间检验 diagnosis 两级（Healthy 3 donor
vs Psoriasis 5 donor）+ region 次级；统计口径与 line3 步骤 9 同源（教科书 KW
含 tie 校正 + 2 组精确置换枚举 C(8,3)=56 无蒙特卡洛误差 + eta2）。
**关键设计混淆如实记录**：healthy 供体只贡献 region=H，psoriasis 供体同贡献
LS/NLS，region 级差异不得解读为区域效应。实测 diagnosis 显著 3/13：
basal_keratinocyte（Healthy 0.371 → Psoriasis 0.089，eta2 0.83）、
suprabasal_keratinocyte（0.289 → 0.083，eta2 0.77）、unknown（0.105 → 0.032）；
T_cell 0.0035 → 0.245、macrophage_mDC 0.0018 → 0.104（免疫浸润方向正确）。
产物 `scrna_composition.h5ad`（obs 加 celltype_donor_prop/donor_n_cells，
uns['composition'] 含比例矩阵/计数矩阵/group_tests_json）+ 
`celltype_composition.tsv`（104 行 = 8 donor × 13 类）+ `composition_report.txt`。

**步骤 15 `deg_groups`（组间差异表达）**：逐 celltype 内做 diagnosis 两组
wilcoxon（scanpy rank_genes_groups + BH）。**表达矩阵口径踩坑**：不能用
adata.X（步骤 4 已把它换成 scale z-score，负值致 log2 出 nan，实测 268/2000
基因 nan），也不能喂原始 counts（expm1 overflow），改为从 layers['counts']
重算 log 归一化（normalize_total(1e4)+log1p 复现步骤 3）。DEG 判据固定
adj p<0.05 且 |log2fc|>1。min 20 cells/组门槛，pDC/Langerhans/NK_ILC/T_cell
（Healthy 组 <20）跳过。实测 324 个显著 DEG：mitotic 166、basal 63、
unknown 40、suprabasal 34、melanocyte 21、channel/progenitor 0。
产物 `scrna_deg.h5ad` + `deg_all_by_celltype.tsv`（14000 行）+ 
`deg_significant_by_celltype.tsv`（324 行）+ `deg_report.txt`。

**步骤 16 `subcluster`（目标细胞类型亚群再聚类）**：目标 T_cell（585 细胞），
子集内重算 counts→lognorm→HVG→scale→PCA(30)→neighbors(15)→Leiden(0.5)→UMAP。
实测 6 个亚群（size 275/103/90/69/33/15）。产物 `scrna_subcluster.h5ad`
（585 × 1664）+ `subcluster_markers.tsv`（300 行）+ 
`subcluster_diagnosis_crosstab.tsv` + `scrna_deg_subclustered.h5ad`（原对象
贴回 subcluster 列）+ `subcluster_report.txt`。

**步骤 17 `enrichment`（功能富集 GO/KEGG/Reactome）**：gseapy 1.3.1，
基因集来自 Enrichr（缓存到 `../references/gmt_cache/`，离线复用）：
GO_Biological_Process_2023（5406 集）/ KEGG_2021_Human（320）/ 
Reactome_Pathways_2024（2100）。ORA（显著 DEG 作前景，超几何）+ prerank GSEA
（全基因 log2FC 排序，nperm=1000）。实测 ORA 2480 行、174 个显著 term
（脂肪酸代谢/转录调控等）；GSEA 16030 行。产物 `enrichment_ora.tsv` + 
`enrichment_gsea.tsv` + `enrichment_report.txt`。
**出图（`enrichment_figures`，2026-10-14 补）**：2 张 PNG 落 `line2_scrna/figures/`：
`l2_11_ora_top.png`（ORA top20，条形=-log10 校正 p，按 celltype 分色，
标注 odds ratio + 校正 p）、`l2_12_gsea_top.png`（GSEA top20，条形=NES 红激活/蓝抑制，
标注 NES + FDR q 值）。

**步骤 18 `tf_activity`（转录因子调控分析）**：TF 靶基因库 ORA
（ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X 92 集 + ChEA_2022 684 集），
显著 TF = 其 ChIP 靶基因在 DEG 中富集；附加 TF 自身 log2FC。
**SCENIC 不可行如实记录**（pySCENIC + cisTarget 人类库 >1GB，运行环境内存/算力
不支持，待更大设备）。产物 `tf_activity.tsv`（2653 行）+ `tf_report.txt`。

**出图（`scrna_figures` 步骤）**：10 张 PNG 落 `line2_scrna/figures/`
（不混进 processed）：UMAP×4（celltype/diagnosis/donor/leiden）、
组成堆叠条/箱线、marker top5 热图、DEG 火山（红上调/蓝下调/灰非显著 + 图例）、
DEG 计数条、PCA by diagnosis（95% 置信椭圆）。前端"结果图件浏览"面板可
内联渲染 + 点击放大。

**本线仍未做（节点生态不适配，待设备）**：无（线2 十项已补齐）。

**图件规范**：见 `docs/figure_spec.md`；三线图件统一落 `<line>/figures/`。

---

## aims 线2 TF 调控精度升级（2026-10-14，步骤 19 `scenic_tf`）

**背景**：步骤 18 `tf_activity` 用"TF 靶基因库 ORA"近似 TF 调控活性（显著 TF =
其 ChIP 靶基因在 DEG 中富集）。本步骤升级为**标准 SCENIC**（GRNBoost2 GRN 推断
-> cisTarget hg38 motif 富集 -> AUCell 活性打分），**取代步骤 18 作为 TF 调控分析
的主路线**（tf_activity 产物保留为对照）。

**口径与算力**：基因空间 = 表达广度 top 5000（env `SCENIC_N_GENES` 可调）；
GRNBoost2 n_estimators=5000（SCENIC 默认，不降树数）；dask 线程模式
（Linux 运行环境内 Nanny 无法 fork worker）。全基因（20990 个 >=10 细胞基因）在节点端
GRNBoost2 约需 90h，不可行，报告如实记录。实测全量 2480.7s（41 min）：
309 TF、450907 GRN 边、**129 regulons x 4587 cells**。

**结果**：top regulon 生物学合理，MXI1/SOX4 -> melanocyte、NFYB/TFDP1/ILF2 -> mitotic、
ING4/USF2/CHD2/HSF1/RELA -> pDC（论文主题 pDC）、NFKB1 -> macrophage_mDC、
HLF/ARID5B -> T_cell、MAFB -> suprabasal_keratinocyte。

**产物**（落 processed/）：`scenic_adjacencies.tsv`（GRN 全表）/
`scenic_regulons.csv`（cisTarget 富集后 regulon）/
`scenic_auc_mtx.tsv`（细胞 x regulon AUCell 活性）/
`scenic_regulon_by_celltype.tsv`（celltype x regulon 平均活性）/
`scenic_tf_summary.tsv`（每 regulon 靶基因数/首选 celltype/特异性 z）/
`scenic_report.txt`。

**出图**（`scenic_figures` 步骤）：2 张 PNG 落 `line2_scrna/figures/`
（top30 regulon x celltype 热图、top20 regulon 特异性条形图）。

**环境变更**：运行环境 dask/distributed 2026.8.0 -> 2024.8.0（pyscenic 0.12.1 兼容，
只动 dask 栈，核心科学栈零变更）；cisTarget hg38 库 722MB 常驻
`public_references/cistarget/`。详见 `ledger/踩坑台账.md` #123。

**验证**：回读 11 项全 PASS；台账 completed/passed；04 巡检 PASS；产物 07 回搬
`outputs/`。

