# 出图规范（aims.txt 第 3 条落地）

> 版本：v2.3（2026-11-08，组合图退役：生产默认只出单图）
> 适用范围：三条线所有分析步骤的图件产物（PCA/UMAP/SpatialFeaturePlot/PPI/DR 网络等）。
> 原则：算力节点资源受限，图件必须"可读、可复现、可归档"。

## 零、组图退役裁定（v2.3）

- 生产母版默认只输出「一张图一个论点」的单图，不再生成多面板组合图。
- 原 G1 到 G5 组图契约自 v2.3 起作废，历史文档中的组合图条款不再作为出图要求。
- 已退役的生产图件：line1 `l1_32_deg_overview.png`，line3 `l3_<species>_16_domain_panel.png`。
- 组合模板 `fig_l3_multipanel`、`fig_l3_combinationpanel`、`fig_l3_sankeybubbleplot`
  一律保持 draft/pending，不得进入 Ready，也不得带 preview 或 render_evidence。
- 同一论点的分面（按样本或分组展开）仍按单图处理；不得把不同统计主题拼成一张组合图。

## 一、图件落点与命名

| 项 | 规则 |
|---|---|
| 目录 | **一律 `<line>/figures/`**（工作台图廊与交付包只认这个目录）；`processed/` 只放数据表与检查点，**不得放图件** |
| 命名 | `<line前缀>_<序号>_<fig语义>.png`，例：`l1_02_pca.png`、`l3_human_10_dr_gsea_top.png`；同语义的步骤级图与汇总图用后缀区分（如 `l1_02_pca_step3.png` 为步骤级） |
| 扩展名 | **静态图统一 PNG**（dpi=300，`bbox_inches="tight"`）；**禁止 PDF/SVG 作为交付图件**；交互式不放服务器（零浏览器面） |
| DPI 声明 | **必须显式声明，不许吃默认值**：R `png()/jpeg()/tiff()` 必须带 `res=300`（R 默认 72dpi，漏写即废图）；Python 必须在 rcParams 声明 `savefig.dpi=300` 或 `savefig(dpi=)`。源级判据 C8/C9，见 §七配套工具链 |
| 单位 | 尺寸按 12x8 英寸基准缩放；图内字号 >= 8pt（缩放与传输后仍可读） |
| 元数据 | 每图必须能在同步骤 report 里找到生成它的口径（种子、输入产物、关键参数、文件名） |
| 回搬 | 图件落 `figures/` 后由任务完成自动回搬（aims §7.3），本地 `outputs/<项目>/<线>/figures/` |

## 二、图件分类与必备元素

### 1. 嵌入类（PCA/UMAP/harmony）
- 必备：点=观测（cell/spot），颜色=分类标签（celltype/domain/leiden），图例可读
- 同一论点需要按物种/样本/分组展开时，用分面单图承载，不跨样本混色
  （domain 身份样本私有，坑 #114 同源）；不同统计主题不得拼成组合图
- 标题注明：嵌入来自哪一步产物（例 "X_pca_harmony from step4 integrate_harmony"）

### 2. 空间类（SpatialFeaturePlot/域图/DR 空间映射）
- 必备：逐样本画（跨样本像素坐标不可比，坑 #114），一样本一图
- 叠层信息：spot/域边界可选，捕获区位置用 array_row/col 或样本内像素坐标，标题注明坐标系
- 颜色：连续量用单色系 ramp（viridis 默认），禁彩虹色系（jet）
- 目标基因定位/DR 映射图：同时给出"表达值"与"该域/邻域统计"两个视图

### 3. 统计类（箱线/火山/热图/富集条形）
- 箱线：须显示个体点（n 小时必须），标注 n
- 火山图/散点：阈值线注明（adj p<0.05、|log2FC|>1），轴标签带基因 ID 类型（symbol/Ensembl）
- 热图：行=基因，列=观测，基因名用 symbol；聚类树可选但须与数值口径一致
- 富集条形：NES/富集比 + 校正 p 双标注，条形按显著性排序

### 4. 网络类（PPI/DR 调控网络/LR 网络）
- 节点=基因/细胞类型，边权标注（权重/相关性/表达比）
- 大网络只画 top N（默认 50 节点、100 边）并在标题注明截断规则
- 布局：spring 布局固定随机种子（可复现）

## 三、样式基线（三线统一）

```python
import matplotlib
matplotlib.use("Agg")            # 算力节点无显示面，禁止交互后端
import matplotlib.pyplot as plt

plt.rcParams.update({
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.grid": False,          # 默认无网格，按需开
    "figure.autolayout": True,
    "savefig.bbox": "tight",
    "image.cmap": "viridis",
})
```

- 分辨率统一 300 DPI（满足 Web 预览 + 发表级双场景；Web 端 CSS 缩放显示，放大清晰）

- 配色：分类标签用 seaborn 定性调色板（colorblind 友好）；连续量单色 ramp
- **二分类（Control vs Treat / WT vs Disease）统一标准色**：Control/WT = `#3B4992`（冷蓝），Treat/Disease = `#EE0000`（暖红）；多分类仍用 seaborn 色盲友好色板与 `RdBu_r`
- **虚拟敲除产出图表标题强制带角标**：`[Virtual KO - DR Genes]` 标识在 caption 首段或图注条首，防与 DEG 混淆
- 每脚本出图后 `plt.close("all")`，防算力节点内存泄漏
- 图内中文禁用（字体缺失风险），标签一律英文

## 三之二、图注（caption）与图例（legend）位置

图注是图的说明文字（含标题语义与关键口径），**必须置于坐标轴框（axes）之外**，
不得压住数据点或图例。

**图注默认位置一律为图幅左下角**（2026-11-06 用户裁定：「每个图的介绍默认写在图的左下角，
并不要影响主图展示」）。实现方式是把文字锚在 `y<0`（轴框下方、画布之外），
由 `savefig(bbox_inches="tight")` 向下扩边容纳，因此**不占用、也不挤压主图数据区**。

| 对象 | 位置 | 排布 | 硬性约束 |
|---|---|---|---|
| 图注 caption | **图幅左下角（轴框外，y<0）** | 纵向多行，一条一行 | 不压数据不压图例 |
| 图例 legend | **轴框外右侧（全部图幅方向）** | 纵向单列 | 条目不重叠；必须留够 subplots_adjust(right) |

图注与图例位置均不按图幅横竖分支（2026-11-07 用户定稿）：
- 旧规范把横版图图例放轴框下方，导致部分图图例仍在主图区内或与图注重叠；
- 统一改为轴框外右侧纵向排列后，图注独占左下角，互不干扰；
- 多组图例（size legend + colorbar 等）按 slot 依次右移，各自纵向，互不叠加。

统一实现（三线脚本各内嵌同名函数，落盘时由 `bbox_inches="tight"` 自动扩边包含）：

```python
def caption(fig, *parts):
    """图注默认置于图幅左下角、轴框之外（y<0），靠 savefig(bbox="tight")
    向下扩边容纳，因此不占用、也不挤压主图数据区（2026-11-06 用户裁定）。
    不再按图幅横竖分支选位置：旧规范把竖版图注放右上，会与图例列抢空间。"""
    parts = [str(p) for p in parts if p not in (None, "")]
    if not parts:
        return
    text = "\n".join(parts)
    fig.text(0.01, -0.02, text, ha="left", va="top", fontsize=8, linespacing=1.5)
```

- 调用处不再用 `ax.set_title()` / `fig.suptitle()`，一律改 `caption(fig, 条目1, 条目2, ...)`
- 条目按语义拆分（如 `"line1: DEG volcano"`、`"adj p<0.05 (FDR)"`、`"n=2000 genes"`），
  无论横竖版一律**一条一行**（横版曾用 `|` 拼成一行，v2.2 起废弃：长条目横拼会超出图幅右缘）
- 字号固定 8pt（图内正文 9pt 或 8pt，图注不喧宾夺主）
- **图内文字不得互相重叠**：点密集的散点/嵌入图，样本或类别标签必须自动避让
  （用 `adjustText` 的 `adjust_text(texts, ax=ax, expand=(1.2, 1.5), arrowprops=...)`），
  不得靠缩小字号（字号下限 8pt 不可破）或删标签来"解决"重叠。
- 图例（legend）是数据映射，不属图注。**任何用颜色编码分类的图都必须有图例**（散点/柱状/箱线/堆叠），
  否则颜色不可读。连续量的图例是 colorbar，且必须带 `label` 说明它是什么（如 `expr`、`KO DR expr`、
  `STRING combined score`），光有刻度数字不算可读。
- 位置规则：**图例一律移到轴框外**（v2.2 起由 `legend_out()` 统一接管，见 §三之三），
  不再允许轴内 `loc="best"`：轴内任何位置在大类别数/主图占满绘图区时都会压数据。
  **源级也不得写轴内 `loc`**（`loc="upper right"` / `"lower right"` / `"best"` 等）：
  `layout_fix` 只是兜底，不能作为"先写违规再靠 guard 纠正"的理由。
  - 所有图幅方向：轴框外**右侧**，纵向单列（`bbox_to_anchor=(1.02, 1.0)`，
    `ncol=1`）；
  - **图例右侧必须配合 `subplots_adjust(right=...)` 预留空间**（典型 `right=0.80`，
    条目多或标签长时更小如 `right=0.75`），严禁 `right>=0.95` 后再放右侧图例
    （2026-11-07 实测：`right=0.98` + `bbox_to_anchor=(1.02, 1.0)` 导致图例
    被裁掉或叠在轴框上，是全部「图例不在右侧」缺陷的共同根因）；
  - 条目 > 12 条时仍按列分组，但**列方向不变**（仍在轴框右侧水平展开，
    每列 ≤ 12 行纵向排列），不搬回下方。
  - 分面/多子图图：**不做图例搬迁**（各分面有自己的图例，统一外移会互相压扁），
    分面内图例按原样保留。
  - `frameon=False`（去边框），落盘靠 `savefig.bbox="tight"` 自动扩边包含图例。
  - 颜色数必须够：分类数超过 `tab20` 的 20 色时要拼接更大色板（`tab20`+`tab20b`+`tab20c` = 60 色），
    否则循环撞色后同色不同类，图例失效。
- 唯一豁免：颜色不是分类的唯一标识时（如柱状图 x 轴标签已写明分组、单色柱状图），可不加图例。

## 三之三、落盘前统一收口（layout guard）

**问题背景**（2026-10-18 用户审图反馈）：单个脚本逐图手调版式会漏，
问题集中在四类：① 长标签（基因名/细胞类型名）超出轴框；② 阈值线标注放中间被点云压住；
③ 图例图注留在轴框内压数据；④ 网络类节点标签在 spring 布局下互相压字。

**收口机制**：每个出图脚本内嵌同一份 helper 块并调用 `_install_layout_guard()`，
它接管 `Figure.savefig`，脚本里**任何** `fig.savefig(...)` 落盘前自动先走 `layout_fix()`。
新增出图代码不需要记得手工调用，版式规则自动生效（脚本独立部署，不能引共享模块，
故按"同一份代码逐脚本内嵌"的方式复制，改一处必须同步 17 份）。

| helper | 作用 | 关键约束 |
|---|---|---|
| `wrap_label(text, width=18)` | 长标签折行：优先按 `_` 边界折，仍超宽再按字符折 | **只用于轴内标注**（条形末端、柱顶）；每次调用点显式包裹 |
| `fit_text_edges(fig, ax)` | 轴内文本右缘/顶缘/左缘/下缘越出轴框时自动扩 `xlim`/`ylim` 收进框内 | 按像素迭代，至少留 6px 边距（比审计容差 3px 更严）；扩一轮无改善即停手 |
| `legend_out(ax, fig, leg)` | 图例外置：所有图幅方向统一轴框外右侧单列（>16 条按列分组） | 只在单坐标轴图生效；`frameon=False` |
| `layout_fix(fig)` | 落盘前总入口：单轴图做图例搬迁 + 文本收边 | 多子图（分面）图只做文本收边，不搬图例 |

**硬规则与例外**：

- **tick 标签一律不换行**。tick 标签本就在轴框外，`savefig.bbox="tight"` 会把画布扩到装下它；
  强行换行反而让多行标签在密集类目轴上行距不够、互相压字（2026-10-18 实测：
  `l1_13_gsea_top` 的 20 行 GO term 换行后行间压字，比不换行更差，已回退）。
  超长 tick 标签按原样交给 matplotlib 的截断/旋转处理。
- 阈值线标注位置固定：**p 阈值（`adj p=0.05` / `p=0.05`）贴左侧轴框内缘**
  （`ha="left"`，x 取 `xlo + 0.02*(xhi-xlo)`）；**log2FC 竖直标注置轴框顶部**
  （`va="top"`，`ha="right"`，`rotation=90`）。放中间会被点云/条形遮挡。
- 网络类节点标签必须自动避让：`nx.draw_networkx_labels()` 之后追加
  `adjust_text(ax.texts, ax=ax, expand=(1.4, 1.6), arrowprops=dict(arrowstyle="-", color="0.75", lw=0.4))`。
  spring 布局每次运行都可能不同，不加避让就是"偶发压字"（不可复现的间歇性 FAIL）。
- 轴内数值标注在条形/柱体之外时（如瀑布图条末挂 `基因 值`），x 轴需按文本宽度留位，
  由 `fit_text_edges` 自动扩 `xlim` 兜底。

**机器验收**（几何审计，不靠肉眼看图）：

```bash
# 单脚本（用当前 cwd 的数据重渲染，默认不落盘，只判定；FAIL 退出码 1）
python3 tools/reporting/audit_figure_layout.py --script <脚本路径> --cwd <数据目录> --runner bash
# 加 --write 才会真落盘；--runner auto 有 bash 时自动走真脚本（多块脚本/环境变量都对得上）
```

判定项（容差 3px）：A 轴内标注 ⊆ 轴框；B 图例在轴框外；C 图注不与轴框相交；
D 标注/图例条目互不重叠（重叠 > 6px² 判违规）；E tick 标签越出画布 → INFO 不计 FAIL
（`tight` 裁切会容纳）。**负对照纪律**：任何"全 PASS"结论必须先在回退版本上证明会红
（删掉 `wrap_label` / `fit_text_edges` 后必须复现 A 项 FAIL）。

## 四、算力约束（与 aims.txt 第 4 节对齐）

- 出图进程同计算任务：`taskset -c 4,5,6 nice -n 10 ionice -c 2 -n 4`
- 大矩阵图（全基因热图）禁止；只画 top N，全量数值落 tsv 由图外工具按需渲染
- 单图渲染超 300s（5 分钟）视为异常，回查是否误画全矩阵

  预算口径（2026-11-06 用户裁定）：「部份生图慢可以容忍，单图生成不超 5min 不算问题。」
  因此从 60s 放宽到 300s。放宽的是**时间**，不是**信息量底线**：
  上面那条「大矩阵图禁止」按可读性执行，与耗时无关，不因预算放宽而失效。
  实践含义：遇到慢图先把 render_ms 记下来判读，确属结构错误（画了上万行无意义矩阵）
  才改画法；只为提速而删减图件内容不算合规优化。

## 五、归档与回搬

- 图件随步骤产物一起过 07 回搬到本地 `outputs/`
- report 里登记：图文件名、对应口径、生成耗时
- PNG 全部字节校验（与 h5ad/tsv 同规则）

## 六、当前落地状态

三线共 21 个出图脚本（含 `pca_matrix` / `deg_limma` / `wgcna_modules` / `celltype_annotation` 等步骤级 R 出图），PNG 全部独立落 `<line>/figures/`，统一走本规范：

| 线 | 脚本 | 图件 |
|---|---|---|
| line1_bulk | bulk_figures | 9 张（相关性/PCA/火山×2/MA×2/topvar/QC×2） |
| line1_bulk | pca_matrix / deg_limma | 2 张步骤级 PNG（`l1_02_pca_step3` / `l1_03_volcano_limma`） |
| line1_bulk | wgcna_modules | 4 张（`l1_28_wgcna_softpower` / `l1_29_wgcna_module_trait` / `l1_30_wgcna_sample_tree` / `l1_31_wgcna_mm_gs`） |
| line1_bulk | ppi_string / key_genes / ml_model / external_validation | 6 张（PPI 网络+hub/交集/ROC/重要性/外部 ROC） |
| line1_bulk | enrichment_figures | 4 张（GO BP/CC/MF GSEA / 模块 ORA / KEGG ORA / KEGG GSEA 各 top20） |
| line1_bulk | ml_hub_analysis | 10 张（LASSO 路径/CV、RFE 曲线、单因素、筛选交集、RF+GBDT 重要性、SHAP、SHAP 瀑布、聚类热图、相关网络） |
| line2_scrna | scrna_figures | 28 张（UMAP×4/组成×14/marker/DEG 火山×7/计数/PCA） |
| line2_scrna | scenic_figures | 2 张（regulon 活性热图/特异性条形） |
| line2_scrna | enrichment_figures | 2 张（ORA / GSEA 各 top20） |
| line2_scrna | target_gene_localization | 1 张（目标基因表达定位 dotplot） |
| line2_scrna | gsva_pathways | 1 张（GSVA 通路活性热图 top30） |
| line2_scrna | dr_post_analysis | 3 张（DR GSEA / DR 调控网络 / DR∩DEG∩marker 交集） |
| line3_spatial | spatial_figures | 每物种 11 张（域图/目标基因/KO DR 逐样本 + 域统计 + Moran/DR/UMAP/组成） |
| line3_spatial | dr_post_analysis | 每物种 2 张（DR GSEA / DR∩域 marker 交集） |
| line3_spatial | enrichment_figures | 每物种 2 张（域 marker ORA / GSEA 各 top20） |

图件数：line1 **35** / line2 37 / line3 69，三线合计 **141**。
（line1 由 29 → 33：pca_matrix / deg_limma / wgcna_modules 原先把 PDF 写在 `processed/`，
2026-10-19 按本规范改为 300 DPI PNG 落 `figures/`；再由 33 → 35：WGCNA 补
`l1_30_wgcna_sample_tree`（样本离群检测）与 `l1_31_wgcna_mm_gs`（kME/GS；`l1_27` 已被
`cox_analysis` 森林图占用，故顺延到 31）。）

**已知限制**：base R graphics（`png()` 设备）无法套用 ggplot 三模板，图幅按 4:3 基准自定；
ggplot 图件（`ggsave`）由 `Rprofile` 垫片自动套模板。

### 6.1 SCI 期刊发表格式核对表（本规范的站位）

SCI 期刊对图件的通用要求与本规范的对应关系。**本规范覆盖期刊硬性要求中的格式与可读性部分**，
投稿前的期刊专属要求（格式说明页的字数限制、文件类型、单双栏宽度）由交付方按目标期刊核对。

| 期刊通用要求 | 本规范条款 | 状态 |
|---|---|---|
| 位图 >=300 dpi（TIFF/EPS/PDF 优先，PNG 需 >=300） | §一 300 DPI PNG | 覆盖 |
| 矢量格式供出版（AI 可编辑） | §一 PDF/SVG 仅作可选副本 | 覆盖 |
| 图内字号可读（期刊通常要求 >=6pt，实际建议 >=8pt） | §一字号 >=8pt | 覆盖（严于期刊） |
| 图例可读且不遮挡数据 | §三之二 图例外置 | 覆盖 |
| 轴标签含单位与基因 ID 类型 | §二.3 / §三之二 | 覆盖 |
| 统计方法、阈值、样本量标注 | §二.3 / §三之三 | 覆盖 |
| 图注（figure legend）为独立成段的文字，逐条并列 | §三之二 caption 一条一行 | 覆盖 |
| 色盲可辨配色，不用 jet/rainbow | §三 配色基线 | 覆盖 |
| 颜色不是唯一区分手段（纹理/形状补充） | §三之二 唯一豁免 | 部分覆盖 |
| 图件命名与正文引用一致 | §一 命名规范 | 覆盖 |
| 图片宽度符合期刊栏宽（单栏 ~85mm / 双栏 ~170mm） | §一 12x8 英寸基准 | 间接覆盖 |
| 一张图一个论点，不在单图堆叠多面板 | §三之二 单图原则 | 覆盖 |
| 图内不含中文 | §三 图内中文禁用 | 覆盖 |
| 每张图可在方法学部分复现（软件、参数、种子） | §一元数据登记台账 | 覆盖 |
| **投稿前需人工核对**：目标期刊的 `author guidelines` 格式页、版权协议、原始数据提交要求 | 无 | **不覆盖，交付方负责** |

**机器验收与评分**（见 §七）：`tools/reporting/audit_figure_layout.py`（几何 A–E 判定）、
`cn_debug/figure_score.py`（五维百分制，>=75 达标）、
`cn_debug/figure_code_audit.py`（源码级静态扫描）。


**2026-10-14 全量核查修复**（对照本规范逐条，见 `ledger/操作台账.md` 续三十）：

- 图内字号统一 `>= 8pt`：原多处 5/6/7pt 的 tick 标签、图例、数值标注全部提至 8pt，
  图幅相应放大以免重叠（如 topvar 热图 50 基因名 5pt→8pt，figsize 加高至 14in）。
- rcParams 基线补全：8 个脚本全部显式声明 `axes.grid=False`、`image.cmap=viridis`
  （原部分脚本缺项）。
- 火山 / MA 图阈值线在图内加文字标注（`p=0.05`、`log2FC=±1`），不只靠 caption。
- 嵌入类（UMAP/PCA）caption 注明来源步骤（如 `X_umap from step5 cluster_umap`）。
- 空间类 caption 注明坐标系（`sample-local obsm['spatial']`）。
- 目标基因 / KO DR 映射图补「域统计」第二视图（`*_by_domain` 条形图），
  与逐样本空间图配对，满足规范的双视图要求。
- 网络类（PPI）边权映射为宽度 + 颜色并附 colorbar，边数限 top100（节点 top30）。
- marker 热图改为行=基因、列=观测（原为行=细胞类型，违反热图方向规范）。
- line2 图名统一加 `l2_` 前缀，与 line1（`l1_`）/ line3（`l3_`）风格一致。

**2026-10-14 二次核查修复（图例缺失专项，见 `ledger/踩坑台账.md` #126/#127）**：

- line3 `01_domains` 逐样本空间域图：补域图例（原只有颜色无说明），图例置轴框外右侧。
- line3 `02_target_gene` / `03_ko_dr` 逐样本图：colorbar 补 `label`（`expr` / `KO DR expr`）。
- line3 `07_domain_composition`：补域图例；域数达 57 > 20 色，色板改为
  `tab20+tab20b+tab20c`（60 色）防撞色，图例置轴框外右侧两列。
- line1 `05a/05b` QC 图：x 轴是样本名、组别只能靠颜色区分，补组图例（LS/NLS）。
- line2 `06_composition_box`：补图例（两组箱体色 + 个体点三项）。
- 规范图例条升级为硬约束：任何颜色编码分类的图都必须有图例；colorbar 必须带 label；
  分类数超 20 色必须扩色板。

**2026-10-14 三次核查修复（全量 102 张逐张目视审查，见 `ledger/踩坑台账.md` #128）**：

- **审查方式**：102 张全部逐张读图（非抽查），4 个审查单元并行，每条线/每张图给出 PASS 或具体问题。
- line1 `02_pca`：样本标签在点密集处互相重叠 → 引入 `adjustText` 自动避让（新增依赖 adjustText 1.4.0，仅新增自身）。
- line1 `05a/05b`：条形顶格时轴内无空白，group 图例压住最右条形 → 图例移到轴框外右侧。
- line2 `04_umap_leiden`：21 类 > tab20 的 20 色循环撞色 → 色板扩为 60 色（同 line3 方案）。
- line2 `06_composition_box_*`（13 张）：Psoriasis 组箱体高时图例压数据 → 图例统一移到轴框外右侧。
- 误报复核：`l2_07_marker_top5` / `l2_scenic_01` 的 colorbar 实际带 label `z`，`l1_05b` 图例在空白区，三者原本合规，未改动。
- 复核：修复后对受影响的 line1（15）+ line2（30）再次全量逐张复核，全 PASS。

**2026-10-14 富集出图补齐（见 `ledger/操作台账.md` 续三十五）**：

- 新增 `src/line1_bulk/scripts/enrichment_figures`：`l1_13_gsea_top`（GSEA GO BP top20）、
  `l1_14_module_ora_top`（WGCNA 模块 ORA top20）、`l1_15_kegg_top`（KEGG top20）。
- 新增 `src/line2_scrna/scripts/enrichment_figures`：`l2_11_ora_top`（逐 celltype ORA top20）、
  `l2_12_gsea_top`（逐 celltype prerank GSEA top20）。
- 口径落实：条形按显著性排序（主键校正 p；并列时 line1 GSEA 按 |NES|、ORA/KEGG 按原始 p，
  line2 GSEA 按 |NES|，避免并列值退化成原始顺序）；GSEA 条形 = NES（红激活/蓝抑制）+ 图例，
  ORA/KEGG 条形 = -log10 校正 p；**NES/富集比 + 校正 p 双标注**（GSEA 标 NES + adj p/FDR q，
  ORA/KEGG 标 GeneRatio 或 odds ratio + 校正 p）；图注置于轴框外；图例不压数据。
- 图件数：line1 15→18，line2 30→32，三线合计 102→107。

**2026-10-14 高级分析补图（见 `ledger/操作台账.md` 续三十六）**：

- 新增 line2 `target_gene_localization`：`l2_13_target_gene_dotplot`（点大小=%表达，颜色=平均表达，
  colorbar 带 label `mean expr (log-norm)`，尺寸图例 `% expressing`）。
- 新增 line2 `gsva_pathways`：`l2_14_gsva_heatmap`（top30 变异通路，行=通路列=细胞类型，
  colorbar 带 label `GSVA z-score`）。
- 新增 line2 `dr_post_analysis`：`l2_15_dr_gsea_top`（DR prerank GSEA，条形=NES + 图例，
  标 NES+FDR q）、`l2_16_dr_regulatory_network`（SCENIC regulon TF→DR 边，蓝=TF 红=DR 基因，
  图例在轴框外，spring 布局固定 seed）、`l2_17_dr_overlap`（DR/DEG/marker 交集规模）。
- 新增 line3 `dr_post_analysis`：`l3_{mouse,human}_10_dr_gsea_top`（小鼠/人源分库 GSEA）、
  `l3_{mouse,human}_11_dr_overlap`（DR∩域 marker 交集）。
- 图件数：line2 32→37，line3 每物种 11→13，三线合计 107→116。
- 缺陷修复：line2 DR 调控网络原 0 边（`scenic_regulons.csv` 的 `TargetGenes` 为 pyscenic
  ctx 原生 repr，含 `np.float64(...)` 调用节点，`ast.literal_eval` 全解析失败）→ 改正则提取
  基因名，修复后 59 边 / 29 TF / 9 DR 靶。

**2026-10-14 line1 hub 基因 ML 模块（见 `ledger/操作台账.md` 续三十七）**：

- 新增 `src/line1_bulk/scripts/ml_hub_analysis`：第一类筛选（LASSO/RFE/单因素）+ 第二类
  重要性评估（RF/GBDT/SHAP/聚类/相关网络），9 张图 `l1_16`–`l1_24`。
- `l1_16_lasso_path`：glmnet 系数路径，x 轴反转对齐 glmnet 约定，选中基因红色 + 图例。
- `l1_17_lasso_cv`：二项偏差 vs log(lambda)，lambda.min/1se 虚线，error bar ±1 SE。
  （原顶部 nzero 双轴已移除：glmnet 在 n=16 共线特征下路径 df 非单调，双轴会误导；
  改为图注写明 lambda.min 选中基因数。）
- `l1_18_rfe_curve`：选中基因数 vs 4 折交叉验证 AUC（LOO 单样本折无法算 AUC，故用 4 折）。
- `l1_19_univariate`：逐基因单因素 logistic，条形 -log10(p)，红/蓝区分 OR 方向，虚线 p=0.05。
- `l1_20_hub_overlap`：LASSO/RFE/单因素三类及其两两交集、真 3-way、实际 core 集。
- `l1_21_importance`：RF + GBDT 重要性分组条形（不减少基因数）。
- `l1_22_shap`：mean|SHAP|（TreeExplainer on RF）。
- `l1_23_cluster_heatmap`：行=基因列=样本热图 + 树状图 + 分组色条，图注标 ARI。
- `l1_24_corr_network`：共表达网络，边色=Spearman r（RdBu_r + colorbar），边宽∝|r|。
- 图件数：line1 18→27，三线合计 116→125。

**2026-10-14 功能富集三线补齐（见 `ledger/操作台账.md` 续三十八）**：

- 核查发现三类缺口：GO CC/MF 三线全缺（line1 `gseGO`/`enrichGO` 只传 ont="BP"，
  line2 库列表只有 GO_BP）；line1 缺 KEGG/Reactome GSEA；line3 完全没有功能富集步骤
  （aims 线3 未要求，用户拍板"三线全补齐"）。
- line1 `enrichment_bulk`：GSEA 改 GO BP/CC/MF 三本体（结果合并落 `enrichment_gsea_go.tsv`，
  含 ont 列；BP 3453 / CC 529 / MF 794 terms）+ 新增 KEGG GSEA（`enrichment_gsea_kegg.tsv`，
  285 terms，143 显著）；WGCNA 模块 ORA 改三本体（1262 行，BP 1010 / CC 111 / MF 141）。
- line2 `enrichment`：库列表加 `GO_Cellular_Component_2023` / `GO_Molecular_Function_2023`
  （5 库并行 ORA + GSEA）。
- line3 新增 `enrichment` + `enrichment_figures`：逐域 ORA（前景 = 每域 adj p<0.05 marker
  最多 top100，不足 10 退 top50）+ 逐域 prerank GSEA（logfoldchanges 排序）；
  mouse 用 KEGG_2019_Mouse/WikiPathways_2019_Mouse + GO 人库大写化匹配，
  human 用 KEGG_2021_Human/Reactome + GO 人库；产物 `enrichment_{ora,gsea}_{mouse,human}.tsv`
  （mouse ORA 62396 行 / human 29736 行）+ 4 张图 `l3_{mouse,human}_12/13`。
- 出图口径：l1_13 按本体分色（蓝 BP / 绿 CC / 紫 MF）+ 图例，y 标签带本体前缀；
  l1_14 y 标签带本体前缀；新增 l1_25 KEGG GSEA。
- 图件数：line1 27→28，line3 61→65，三线合计 125→130。

**2026-10-14 line1 ML 配色优化 + SHAP 瀑布图（见 `ledger/操作台账.md` 续三十九）**：

- 用户新分析方案对照核查：核心流程（DEG→富集→PPI→hub→LASSO→RF）与进阶分析
  （共表达网络/SHAP/聚类）已全落地，缺的实现项只有 SHAP 瀑布图；其余为配色差异性问题。
- 新增 `l1_26_shap_waterfall`：单样本预测分解（选 SHAP 绝对值总和最大的样本），
  条形 = base value + 逐基因 SHAP 累积（f(X)=E[f(X)]+sum(SHAP)），红=推向 LS / 蓝=推向 NLS，
  条上标注基因名 + 贡献值。坑：shap 0.53 的 `expected_value` 是 per-class ndarray，
  直接 `float()` 抛 "only 0-dimensional arrays..."，需 `arr[-1]` 取 LS 类。
- 配色优化（差异性突出）：
  - `l1_19_univariate`：原红/蓝按 OR 方向分色（5 个显著基因混在蓝堆里不可辨）
    → 改为红=显著（p<0.05 加星加粗）、蓝深浅 ~ p 值渐变（越接近 0.05 越深）。
  - `l1_20_hub_overlap`：原 8 柱 8 色语义混乱 → 三层语义配色（深蓝=单方法、
    浅蓝=交集、橙=下游使用的 core）+ 数值标签 core 加粗金色 + 图例。
  - `l1_21_importance`：红蓝对比强化（#2166AC/#B2182B）+ 每条形附数值标注 +
    x 轴留 22% 标签空间。
- 全量复核 line1 ML 10 张 + 富集 4 张（逐张读图）：全部合规（图注轴框外、图例不压数据、
  字号 >= 8pt、colorbar 带 label）。
- 图件数：line1 28→29，三线合计 130→131。

**2026-10-18 版式统一收口（见 `ledger/踩坑台账.md` #140–#143）**：

用户审图反馈四类版式问题，本轮按 §三之二 / §三之三 收口，17 个出图脚本全部内嵌同一份
canonical helper（`caption`/`wrap_label`/`fit_text_edges`/`legend_out`/`layout_fix`/`_install_layout_guard`）：

- **长标签溢出图框**：`l1_26_shap_waterfall`、`l2_scenic_02_top_regulons` 等条形图的
  条形末端/柱顶挂标（基因名、细胞类型名）超出轴框 → `wrap_label(width=18)` 折行 +
  `fit_text_edges` 四边扩轴兜底。
- **p 阈值标注被点云遮挡**：`l2_08_deg_volcano_*`（7 张）与 line1 火山（fdr/nominal）的
  `adj p=0.05` 标注由中间移到**左侧轴框内缘**；`log2FC=±1` 竖直标注移到轴框**顶部**。
  顺带修正 line1 火山口径文案（`adj p=0.05` / `p=0.05` 随 fdr/nominal 分支取）。
- **图例图注位置**：图注固定在图幅左下角；图例由 `legend_out()` 统一外置到轴框右侧
  （所有图幅方向单列纵向，>16 条按列分组，`frameon=False`）。
- **网络节点标签压字**：`l2_16_dr_regulatory_network` 与 `l1_07_ppi_network` 的
  spring 布局节点标签加 `adjust_text` 自动避让（后者的 STRING API 每次重取导致布局变化，
  属**间歇性**压字风险，靠避让而非固定种子兜住）。
- **新增机器验收**：`tools/reporting/audit_figure_layout.py`（接管 `Figure.savefig` 做 A–E 几何判定，
  容差 3px，默认不落盘）。真数据全量审计 **135 张图 fail=0**（line1 29 / line2 37 / line3 69）。
  负对照已验证：回退 `wrap_label`+`fit_text_edges` 后 `l2_scenic_02` 必复现 A 项 FAIL。
- **踩坑回退**：tick 标签换行方案（走 `set_*ticklabels`）实测在 `l1_13_gsea_top` 的密集类目轴上
  反而压字，整体回退，并在 §三之三 写明"tick 标签一律不换行"的理由。
- 图件数不变（line1 29 / line2 37 / line3 69，三线合计 135），本轮为版式重渲染。

**2026-10-19 全量核查与打分系统上线（见 `ledger/踩坑台账.md` #180、#181–#183）**：

用户提出五项要求：规范预核查、代码+图片全维度合规、质量打分系统（>=75 达标）、
WGCNA 分组逻辑专项、统一错误台账。本轮全部落地：

- **规范补齐**：新增 §6.1 SCI 期刊格式核对表、§七 出图质量打分系统、§八 分组与比对契约。
- **代码修正**：源码级静态扫描（`cn_debug/figure_code_audit.py`，C1–C7 项）查出
  **8 处轴内 `loc="best"/"upper right"/"lower right"` 图例**与 1 处 `cex=0.65`（<8pt），
  分布 bulk_figures / enrichment_figures / ml_hub_analysis / dr_post_analysis / scrna_figures /
  cox_analysis，已全部改为轴框外并统一 8pt 下限。C4 项最初把 ggplot `geom_point(size=)`、
  gseapy `min_size=` 误报为字号，已收紧为只查 matplotlib `fontsize=` 与 R `cex=`。
- **图件校验**：`cn_debug/figure_audit_scan.py` 对 30 张本地图件逐张判定
  （格式 / DPI / 可读面积 / 边缘留白 / 非空白），**30/30 全部 PASS**。
- **打分系统**：`cn_debug/figure_score.py`，五类 20 项满分 100，实测 **30 张均 100 分**。
- **WGCNA 分组专项**：见 §八 与踩坑 #181。


## 七、出图质量打分系统（五维百分制，>=75 达标）

**为什么要有分数**：肉眼审图会漂移，且"合规/不合规"的二值判断无法表达"能交付但不理想"。
打分把规范条款量化为可比较、可跟踪的数值，作为**交付前门禁**与**回归对照基线**。

执行入口：`python3 cn_debug/figure_score.py --dir <本地图件目录> [--csv out.csv] [--json out.json]`

**权重与细则（满分 100）**：

| 维度 | 分值 | 细则 |
|---|---|---|
| **A 格式合规性** | 30 | A1 扩展名 PNG 5 · A2 DPI==300 8 · A3 最短边>=1200px 5 · A4 文件>20KB 4 · A5 非空白(std>5) 4 · A6 边缘无裁切(edge_ink<5%) 4 |
| **B 参数准确性** | 20 | B1 图内标注⊆轴框 6 · B2 图例在轴框外 5 · B3 图注不与轴框相交 5 · B4 标注/图例互不重叠 4（由 `tools/reporting/audit_figure_layout.py` 的 A–E 几何判定提供） |
| **C 标注规范性** | 20 | C1 `figures.js` caption 登记命中 6 · C2 命名符合 `l<n>_[前缀_]<序号>[_分号]_语义.png` 4 · C3 图件与 `STEP_CATALOG` 声明一致 6 · C4 台账/report 登记生成口径 4 |
| **D 视觉可读性** | 20 | D1 无低 DPI 4 · D2 无 <8pt 字号 6 · D3 图例未用轴内 `loc` 6 · D4 未用 PDF/SVG 作交付件 4 |
| **E 逻辑正确性** | 10 | E1 分组/方向来自 `project.json` contrasts 或报告注明的派生规则 5 · E2 相关值域 [-1,1] 且方向与 case/ctrl 定义一致 5 |

**判定标准**：

- `>= 75` → **PASS**，可进交付。
- `60–74` → **FIX**，需按缺失项整改后重跑。
- `< 60` → **FAIL**，不得进交付。

**一票否决**（无论总分，直接 FAIL）：A1（非 PNG）、A2（DPI≠300）、B1（图内标注出框）、
C2（命名不合规）。这四项任意一项失败说明图件不可用或不可归档，得分再高也没有意义。

**三态判决与外部判决导入（2026-10-19 起，FE-012 关闭）**：

B1–B4（几何项）与 D2/D3（源级字号/图例项）**不再允许默认满分**。它们的唯一合法来源是
机器判决，必须显式导入；导不到就记 **UNVERIFIED**，"没测"不等于"通过"：

```
python3 cn_debug/figure_score.py --dir <图件目录> \
  --audit-csv cn_debug/figure_audit_line1.csv \      # tools/reporting/audit_figure_layout.py --csv 的逐图判决（B1–B4）
  --code-csv  cn_debug/code_audit.csv \              # figure_code_audit.py --csv 的源码扫描（D2/D3）
  --ledger-dir outputs/system_logs                   # 台账反查"图→出图脚本"，供 D2/D3 归因
```

- 判决优先级：一票否决 FAIL > 存在未验证项 **UNVERIFIED**（不得进交付）> 按分数线 PASS/FIX/FAIL。
- 退出码：FAIL 或 UNVERIFIED 均返回 1，即**未验证的图件不能算交付**。
- 未验证项在输出里带 `?` 后缀（如 `B1?`），与"已验证不通过"区分：前者要补测，后者要整改。

**已知边界**：几何审计靠 patch `matplotlib.Figure.savefig`，**只覆盖 Python 出图**；
R `png()` 设备的图件（`pca_matrix` / `deg_limma` / `wgcna_modules` 等）跑出来的 `figs=0`
不代表合规，必须走人工审查登记，或补 R 侧审计（方案待定，见 `workspace/handoff/archive/2026-10-19_figure_polish_plan.md` §待决策）。

**实测基线（2026-10-19，PRJ_E2E_GSE277173 / line1_bulk，外部判决导入后）**：

```
python3 cn_debug/figure_score.py --dir outputs/PRJ_E2E_GSE277173/line1_bulk/figures \
  --audit-csv cn_debug/figure_audit_line1.csv --code-csv cn_debug/code_audit.csv \
  --ledger-dir outputs/system_logs
# score: figs=33 PASS=24 FIX=0 FAIL=0 UNVERIFIED=9 avg=94.4 threshold=75
```

9 个 UNVERIFIED 的来源与处置（**不接受"看着没问题就改 PASS"**）：
6 张 R 图（`l1_02_pca_step3` / `l1_03_volcano_limma` / `l1_28`~`l1_31`，几何审计不覆盖 R，
走 `ledger/figure_manual_review.md` 人工清单）、2 张 PPI 图（`ppi_string` 每次重取 STRING API
不可复现，禁跑审计）、1 张陈旧副本（`l1_27_wgcna_mm_gs.png`，远端已无此图，见 FE-017）。

**配套工具链**（三者互补，缺一不可）：

| 工具 | 层级 | 检查什么 |
|---|---|---|
| `cn_debug/figure_code_audit.py` | 源码 | C1 轴内图例 / C2 禁用后端 / C3 低 DPI / C4 低字号 / C5 图件落错目录 / C6 PDF 交付 / C7 图内中文 / **C7 R `png()` 未声明 `res=`** / **C8 Python `savefig` 未声明 dpi** |
| `tools/reporting/audit_figure_layout.py` | 渲染过程 | A 标注出框 / B 图例位置 / C 图注相交 / D 重叠 / E tick 越界（INFO） |
| `cn_debug/figure_audit_scan.py` + `figure_score.py` | 交付物 | DPI / 面积 / 边缘 / 空白 / 命名 / caption 登记，汇总为百分制得分 |

---

## 八、分组与比对契约（"vs" 的语义）

**硬规则**：**`vs` 表示两个组别进行组间比对，不是把两个组合并为一个组。**

1. **方向必须显式**：任何按分组拆分的图与统计，必须能说清 case（实验组）与 ctrl（对照组）分别是谁。
   case/ctrl 只能来自 `project.json` 的 `grouping_design.contrasts`；没有该契约时按
   `grouping_design.groups[].is_reference` 或字典序派生，且**必须在报告 footnote 里写明派生了**。
2. **禁止硬编码组名**：出图与统计代码不得写死 `"LS"` / `"NLS"` 或任何具体组名；
   组名从 `group_map.tsv` 与 project.json 读。写死组名会让换项目时静默算错方向。
3. **多组分列不合并**：有三组及以上时，性状/颜色/分面按组分列，`vs` 逐对计算，
   不得把两组压成一个合成标签。
4. **图注必须带方向**：涉及组间比对的图，图注写明 `contrast <case> vs <ctrl>`（例
   `l1_29` 写 `contrast LS_vs_NLS (case=LS vs ctrl=NLS)`）。
5. **二分类配色固定**：ctrl = `#3B4992`（冷蓝），case = `#EE0000`（暖红）；分组数 > 2 时才用定性色板。
   颜色不得与 case/ctrl 角色脱钩（不得出现"红色代表对照组"）。

**WGCNA 模块-性状相关图的落地要求**（2026-10-19 专项，踩坑 #181）：

- 性状矩阵按分组展开：主对比对一列 `case=1 / ctrl=0`，另可为每组各加一列指示变量；
  **不得把多个组压成单一 0/1 列**。
- 模块与每一列分别算 Pearson 相关；报告与图注只引用**主对比对**那一列，并写明依据。
- `cor > 0` 的含义必须写成"该模块特征在 case 组相对升高"，`cor < 0` 为"相对降低"；
  方向随 case/ctrl 定义，不固定在某个组名上。
- 图注固定三要素：对比对名、case/ctrl 名、trait 来源（project.json 或派生的）。

---
## 比例尺的物理长度怎么验（2026-11-07 定）

比例尺是图上唯一的物理尺度标尺，错 5% 就是结论级错误。判据与量测方法如下。

**标定口径（2026-11-07 实测改判）**：以 **Visium 晶格 pitch 100 µm** 为准，
`um_per_px = 100 / 相邻 spot 中心距(px)`。不要把 `spot_diameter_fullres` 当基准：
它是 Space Ranger 对光斑直径的逐图估计，本数据集 5 个样本在 87.06~90.19 px 间波动，
而最近邻恒为 132.0 数据单位。用 `65/d_px` 会得到 0.7207~0.7466 的分散值，
使标注 500 µm 的条实际画出 507.9~525.5 µm（系统性偏长 1.5~5.1%）。
按晶格标定后 500 µm = 660.0 数据单位 = 恰好 5 个晶格间距，与平台规格严格一致。

**两份实现必须同算法**：出图母版的 `_spot_um_per_px` 与 `ops/spatial_calibration._first_shell_pitch`
各有一份 pitch 计算。两处都只取「前 6 个邻居全部落在第一壳层（<1.2 倍中位）」的内部点再取中位，
否则边缘 spot（最近邻跳到第二壳层，约 1.73 倍）会把 pitch 抬高，导致图上的 500 µm
与代码里的 200 µm 在不同切片上分叉。合成晶格自测见 `tools/reporting/spatial_pitch_selftest.py`。

**绘制时机**：比例尺必须登记后置，在 `fit_text_edges` 等所有会改 `xlim/ylim` 的布局动作
全部收敛之后再画，画完还要再收边一轮并确认轴限不再变化。数据坐标系下的图元，
轴限一变长度就变，先画后调布局必然失真。

**两级验证，各有分工**：

1. **图内自检（防布局漂移）**：画完读**真实 artist 的端点**（`line.get_xdata()`），用
   `transData` 换算显示像素，再除以每数据单位像素数得到实际微米数，误差 >3% 直接失败。
   禁止用 `bar` 重算端点：`dxb - dxa` 恒等于 `bar` 的像素长度，那是恒等式，永远通过。
   这一级的边界：只防「画完被后续布局挤压」，**不防「一开始常量就算错」**。
2. **成品图像素核验（防渲染链路失真）**：`tools/reporting/figure_scalebar_png_measure.py` 有两种用法。
   - `measure(png)`：以晶格 pitch（Visium 相邻 spot 恒为 100 µm）为绝对尺反推微米数。
     **只适用于 spot 阵列清晰的图型**（如 01_domains）；点云过密或配色复杂的图型上，
     x 投影自相关会挑到错误的峰，量出来会系统性偏大 5~15%，不能当判据用。
   - `cross_check_dir(dir)`：把 manifest 的 `bar_px` 与每张 PNG 的端刻度间距逐张比对。
     这条不依赖晶格，因此对**所有**图型成立，证明「画出来的条长 == 代码记录的条长」，
     即 dpi、bbox、布局重排对比例尺零失真。实测 33 张图全部一致（19/19 + 14/14）。

**量测口径的两个坑（第 73 轮踩过）**：

- 比例尺长度取**两端刻度中心之差**，不取横线暗段长度：线宽 1.8pt 加抗锯齿会在两端各留
  约 2-3 像素的晕，133 px 的暗段对应的真实条长是 128 px，直接用暗段长会把比例尺算大 4%。
- 六角晶格沿 x 方向投影的周期是 `pitch × cos30°`，不是 `pitch`。把投影周期 22 px 当成
  真实 pitch，会让反推的长度偏小 15%，进而把 500 µm 误判成 470 µm。
