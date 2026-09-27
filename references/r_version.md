# scrna R 版立项（R-05，2026-09-26）

> **本文档是 scrna 仓 R 版的实现契约。** R 版与 Python 版**并存**，用户按配置二选一；
> 两版共用同一份 `assets/*.yml` 配置、同一套产物路径与文件名、同一套图名 ——
> 所以下游的 `governance/manuscript_map.yml` / `governance/figure_captions_draft.csv` /
> `governance/SCI_HANDOFF.md` **对两版都成立，一个字不改**。
> spatial 仓的对应文档 R-16 另写，那边的方法学差异（SpaGCN → BayesSpace/Giotto 等）
> 逐条标注。
>
> **关于「引用了还不存在的文件」**：这是**立项**文档，引用的 R 侧产物
> （R 入口脚本、R 公共库、R 侧 CI workflow、spatial 仓对应文档）都是
> R-06..R-16 要新建的东西，写此文时磁盘上还没有。按本仓门禁的规矩
>（不指向不存在的文件），这些**计划中的**名字一律不加反引号、不加路径前缀，
> 让引用检查器分得清「已存在的引用」与「计划中的占位」；
> 等 R-06..R-16 落地后再补成真实引用。

---

## 1. 结论摘要

1. **并存、二选一**：`scripts/main_analysis.py`（Python）与 R 版入口脚本（R）都存在，
   由配置顶层新键 pipeline.language: python | r 选择；**一次运行只出一种语言的产物**，
   不混装。产物目录与文件名相同。
2. **配置零拷贝**：R 版用 R 的 yaml 包读同一份 `assets/config.pbmc3k.yml`
   （158 行，顶层键 `dataset_id / source / design / qc / norm / reduce / integration /
   trajectory / communication / grn / perturbation / analysis`）与 6 个资源表
   （`datasets / celltype_mapping / celltype_markers / ligand_receptor / qc_gene_sets / tf_list`）。
3. **图名集合相等**：R 版产出的 39 个图名必须与 Python 版**完全一致**
   （新增门禁判据 `figures:parity`，见 §6）。因此 caption CSV、图组文档、
   governance 仓根的 SCI_HANDOFF.md 的对接端口对两版是同一个。
4. **9 步里 1 步零差异、5 步同方法不同实现、3 步真方法差异**（§3 的等级表）。
   差异最大的两步（01 双细胞检测、05 拟时序）都**不改变结论形态**，
   只改变"哪条路径得到的数"—— 状态 JSON 里照旧写明用的哪条路径。
5. **08 是纯 R 既有资产**：Python 版本来就是 shell 出去调 `scripts/lib/tenifold_knk.R`
   （567 行）。R 版直接调用同一份文件 → **零差异**，且 `tools/check_r_syntax.mjs`
   现有的 R↔Python 参数/产物对齐检查（L277-361）原样适用。
6. **工作量量级**：Python 版 9 步脚本共 **5654 行** + `main_analysis.py` **1270 行** +
   `lib/common.py` **1589 行** = 8513 行。R 版目标量级相当（geo 仓的 R 公共库（≈1380 行）
   可作结构模板）。

---

## 2. 已核实的事实（立项依据，全部 grep 到源码）

### 2.1 Python 版现状

| 项 | 事实 | 出处 |
|---|---|---|
| 步骤表 | 9 步：`fetch / qc / integrate / cluster_annotate / pseudobulk_de / trajectory / communication / grn / virtual_perturbation` | `scripts/main_analysis.py:38-49` |
| 状态文件 | 8 键 `STEP_STATUS_FILES`；`fetch` 是唯一豁免（`_STEPS_WITHOUT_STATUS_FILE = {"fetch"}`，写的是 `data_dir/dataset_info.json`） | `main_analysis.py:55-64`、`:80` |
| 启动期对齐 | `{s[0] for s in STEPS} - set(STEP_STATUS_FILES) - _STEPS_WITHOUT_STATUS_FILE` 非空即 `raise RuntimeError` | `main_analysis.py:82-89` |
| 人工复核 | 6 项 `HUMAN_REVIEW_NODES`，默认 `pending` 不自动转 `confirmed` | `main_analysis.py:95-105` |
| 输入登记 | 5 项 `INPUT_FILES`：`raw.h5ad / dataset_info.json / qc_filtered.h5ad / integrated.h5ad / clustered.h5ad` | `main_analysis.py:108-114` |
| 必需产物 | `REQUIRED_FILES`：`state.json / qc_status.json / integration_status.json / cluster_status.json / markers_all.csv / cluster_resolution_scan.csv / celltype_annotation.csv / celltype_scores.csv / …` | `main_analysis.py:117+` |
| 语言字段 | `init_manifest(cfg, language: str = "python")` → `run_manifest.json` 已带 `language` | `scripts/lib/common.py:637` |
| 跨语言桥 | `record_cross_language()`：只走 CSV 桥接；对象级工具限 `zellkonverter` / `anndata2ri`，**禁用 sceasy** | `common.py:929` |
| R 探针 | `probe_r_packages()`（`Rscript -e 'cat(...packageVersion(p)...)')`，包名白名单 `^[A-Za-z][A-Za-z0-9.]*$` | `common.py:686` |
| 既有 R 资产 | `scripts/lib/tenifold_knk.R` 567 行，CI 里 `Rscript scripts/lib/tenifold_knk.R --selftest` | `.github/workflows/scrna_analysis.yml:218` |
| 图名门禁 | `tools/check_fig_names.mjs:151-153` 的过滤器**已经接受 `.R`**：`/^\d+_.*\.(R\|py)$/`；跨脚本重名判红在 L188-191 | `tools/check_fig_names.mjs` |
| R 静态门禁 | `tools/check_r_syntax.mjs`（15181 字节）：walk `scripts/` 收 `.R`（L38-45），已含 R↔Python 参数/产物名对齐（`walkPy`、`rParamNames`、`rOutputNames`，L277-361） | `tools/check_r_syntax.mjs` |
| 图库 | 39 张（PNG+PDF 各 39）：`02-01`×6、`02-02`×2、`02-03`×12、`02-05`×10、`02-06`×1、`02-07`×7、`02-08`×1 | `figure_review/INDEX.txt` |
| 脚本行数 | `00_fetch` 356 / `01_qc` 427 / `02_integrate` 288 / `03_cluster_annotate` 815 / `04_pseudobulk_de` 321 / `05_trajectory` 1259 / `06_communication` 558 / `07_grn` 498 / `08_virtual_perturbation` 1132（合计 **5654**） | `wc -l scripts/*.py` |

### 2.2 各步实际方法（grep 到源码，防止立项文档写错方法）

| 步 | 实测方法 | 出处 |
|---|---|---|
| 00 fetch | `urllib.request` 下载（**UA 伪装** —— Python-urllib 直连 403 的实测教训）+ `sc.read_10x_mtx` / `sc.read_10x_h5` / `ad.read_h5ad` + `validate_counts` | `scripts/00_fetch.py:22,49-91,144,148-155,158` |
| 01 qc | `sc.pp.calculate_qc_metrics` + **`sc.pp.scrublet`** + `sc.pp.filter_cells/genes` | `scripts/01_qc.py` |
| 02 integrate | `normalize_total / log1p / highly_variable_genes / combat` + `sc.tl.pca` + **`sc.external.pp.harmony_integrate`（harmonypy）** | `scripts/02_integrate.py` |
| 03 cluster_annotate | `sc.pp.neighbors` + `sc.tl.umap / leiden / rank_genes_groups / score_genes` + `sc.get.rank_genes_groups_df` | `scripts/03_cluster_annotate.py` |
| 04 pseudobulk_de | **pydeseq2**（DESeq2 的 Python 移植） | `scripts/04_pseudobulk_de.py` |
| 05 trajectory | 自实现 CytoTRACE GCS（`compute_cytotrace_gcs`）+ `compute_dpt` + **`compute_palantir`（扩散图 + Markov）** + **`compute_scfates`（`scFates` 的 **PyPI Python 移植**，"Slingshot 的 Python 移植"）**；`stable_methods = ["dpt","palantir","cytotrace"]` | `scripts/05_trajectory.py:60,109,120,149,443-456,1201` |
| 06 communication | **LIANA**（`import liana as li`，`LIANA_N_PERMS = 100`，`mt.rank_aggregate`）为主；自建共表达打分为退回。**CellChat 只是 LIANA 的 consensus 资源之一，不是独立工具** | `scripts/06_communication.py:7-15,105,108-140` |
| 07 grn | **共表达推断，不是 SCENIC**（docstring 明说无 cisTarget motif 剪枝）；TF 列表来自 `assets/tf_list.yml`；调控子活性 = **AUCell 式**（`sc.tl.score_genes`）；每个调控子报 `n_cells_expressing_tf` 与 `cluster_specificity` | `scripts/07_grn.py:5-24,53,146` |
| 08 virtual_perturbation | scTenifoldKnk —— Python 侧 shell 出去跑 `scripts/lib/tenifold_knk.R` | `scripts/08_virtual_perturbation.py` |

### 2.3 可参考代码（用户自己的 R 代码，`可参考代码/`）

| 目录 | 对应步骤 | 里面实际用的包 |
|---|---|---|
| `95.单细胞：数据读取与质控` | 00 / 01 | `10X数据读取.R` / `H5读取.R` / `1.单矩阵读入.R`：Seurat、hdf5r、Matrix、tidyverse |
| `96.单细胞：鉴定高变基因` | 02 | Seurat、scater、tidyverse |
| `97.单细胞：PCA降维及去批次` | 02 | Seurat、harmony、PCA（tidyverse） |
| `99.单细胞：自动注释` | 03 | SingleR、celldex、clustree、monocle3、人工注释（Seurat）、BiocParallel |
| `101.单细胞：富集分析` | 03/04 后续 | clusterProfiler、org.Hs.eg.db、DOSE、GSVA、UCell、irGSEA |
| `104.单细胞：拟时序分析` | 05 | monocle3、celldex、tricycle、assertthat |
| （贯穿） | 出图 | SCpubr、ggplot2、patchwork、cowplot、ggsci、randomcoloR、future |

> 这些文件的头部都带同一大块 `library()` 清单（模板头），**不能**据此断言"每个文件都用了
> 所有这些包"；以各目录对应的步骤语义为准。

---

## 3. 技术选型表

**差异等级**：
- **D 零差异** —— 同一份代码/同一个工具，连数值都应逐位一致。
- **A 同方法不同实现** —— 方法名一样（DESeq2 / ComBat / AUCell），数值可能有末位差
  （规则 12 的口径：可声称结构与量级可复现，不声称逐字节）。
- **B 换实现同思路** —— 等价操作（`normalize_total` ↔ `NormalizeData`），
  数值量级可比但不逐位。
- **C 真方法差异** —— 换了算法/统计量，**必须在状态 JSON 与文档里标注**，
  并写对结论形态的影响。

| 步 | Python 实测 | R 选型 | 等级 | 必须标注的差异 |
|---|---|---|---|---|
| 00 fetch | `urllib`（UA 伪装）+ `sc.read_10x_mtx/h5`/`read_h5ad` | `download.file()`（**UA 问题在 R 侧要重新踩一遍**，用 `options(HTTPUserAgent=...)`）；`Seurat::Read10X()` / `Read10X_h5()`（hdf5r）；GEO 源走 `GEOquery::getGEO()`；`validate_counts` 等价实现 | B | 读进来的矩阵**坐标约定**（`rownames=基因符号`）必须与 Python 版一致，否则后面全部错位 |
| 01 qc | `calculate_qc_metrics` + **scrublet** + filter | `Seurat::PercentageFeatureSet(pattern="^MT-")` + **`scDblFinder::scDblFinder`** + `scater::isOutlier`（MAD）+ 手工阈值 | **C** | 双细胞检测器不同：scrublet 是"模拟双体 + 分类器"，scDblFinder 是"人工双体 + 建模"；**检出数量与阈值口径会不同** → `qc_status.json` 里 `n_doublets` 等字段两版不可直接比 |
| 02 integrate | `normalize_total/log1p/HVG/ComBat/PCA` + harmonypy | `NormalizeData / FindVariableFeatures / ScaleData / RunPCA` + `sva::ComBat` + `harmony::RunHarmony` | A | 归一化的 scale factor 缺省值（`target_sum=1e4`）两侧一致；HVG 选择准则（seurat_v3 vs dispersion）**必须对齐配置**，否则 HVG 集不同 → PCA/聚类全不同 |
| 03 cluster_annotate | neighbors + umap + leiden + rank_genes_groups + score_genes | `FindNeighbors` + `RunUMAP` + `FindClusters(algorithm=4)`（Leiden）+ `FindAllMarkers`（Wilcoxon）+ `AddModuleScore` + `SingleR`+`celldex` + `clustree` | A/B | `AddModuleScore` 与 `sc.tl.score_genes` 同源（scanpy 是照 Seurat 移植的），但 bin 数缺省值不同（25 vs 100）→ 必须在配置里显式对齐 |
| 04 pseudobulk_de | **pydeseq2** | **`DESeq2`**（母实现） | **A** | 拟合算法同（负二项 GLM + Wald/LRT），数值可能有末位差；`sizeFactors`/`design` 公式必须同构 |
| 05 trajectory | 自实现 GCS + DPT + **Palantir** + **scFates（Python 移植）** | 同一套**自实现 GCS（逐行转写公式，保证两侧同值）** + `destiny::dpt` + **Palantir 无 R 实现 → 丢弃** + **`scFates`（CRAN 原版）** | **C** | ① Palantir 丢弃后方法数从 4 → 3，`trajectory_status.json` 的 `method_correlation` 仍成立（≥2 方法即可）；`stable_methods` 名单变为 `["dpt","scfates","cytotrace"]`。② scFates 从"Python 移植"换成"CRAN 原版"，**理论上更接近原始实现**，但数值会有差异。③ **`scFates` 是否在 CRAN 上可装，R-07 实测确认**；装不上则退回 `monocle3` 或 `slingshot` 并在文档标注 |
| 06 communication | LIANA（`rank_aggregate`，100 perms）+ 自建共表达退回 | `liana`（**GitHub saezlab/liana**，非 CRAN/Bioc）+ 同一套自建退回 | A/B | `liana` 的 R 版是 GitHub 包，CI 里要 `remotes::install_github`（有版本漂移风险）；装不上就退回自建 —— 与 Python 的退回逻辑**同构** |
| 07 grn | **共表达推断**（无 motif 剪枝）+ AUCell 式打分 | **同样不做 SCENIC**：`cor()` 共表达矩阵 + `AddModuleScore`（AUCell 式）+ 同样的 `n_cells_expressing_tf` / `cluster_specificity` | **A** | **不许"升级"成真 SCENIC**：加了 motif 剪枝后"共表达模块"就变成"可能有直接调控"，Python 版 docstring L15-21 写的所有诚实性说明会全部作废，两版结论不可比。**方法保持一致是硬约束** |
| 08 virtual_perturbation | shell → `scripts/lib/tenifold_knk.R` | **直接调同一份 `lib/tenifold_knk.R`** | **D** | 零差异；`check_r_syntax.mjs` 的 R↔Python 对齐检查原样适用 |

> **为什么 07 坚持"不升级"**：R 版的价值是"用户可选实现"，不是"换个方法重算"。
> 一旦 R 版用了 SCENIC，`tf_regulons.csv` 里同一行在两版里的**含义**都不一样，
> 任何"两版对照"都会变成假对照。方法升级要单独立项（S-01 那类），混进语言移植里
> 是最糟糕的做法。

---

## 4. 命名契约与共存

### 4.1 脚本与入口

```
scripts/
  00_fetch.R            ← 与 Python 同号同语义
  01_qc.R
  ...
  08_virtual_perturbation.R
  main_analysis.R       ← 入口；结构复刻 main_analysis.py
  lib/
    common.R            ← 基础设施（geo 仓 common.R ≈1380 行作结构模板）
    tenifold_knk.R      ← 既有资产，两版共用，不许复制出第二份
  r_deps.R              ← CI 安装脚本（版本清单写死在这里，不靠 DESCRIPTION）
```

- **`check_fig_names.mjs:151-153` 已接受 `.R`** —— 所以 R 脚本放进 `scripts/` 后会被
  同一个门禁扫到。这是**好事**：R 版图名自动进图名门禁。
- **代价**：跨脚本重名判红（`check_fig_names.mjs:188-191`）会把 R/Python 同名图判成
  重名 → 必须把判据改成**按语言分域**（§6.1）。

### 4.2 语言选择接口

```yaml
# assets/config.pbmc3k.yml 顶层新增（缺省 python，两版都认）
pipeline:
  language: python      # python | r
```

- Python 侧：`scripts/main_analysis.py` 读到 r 就**立即报错退出**（不许静默改跑 Python）。
- R 侧：R 版入口脚本读到 python 同样报错退出。
- `use.sh` 不动（它只管把 skill 注册进 agent 的发现路径）。

### 4.3 图名（§1 的第 3 条）

- R 版必须产出**与 Python 版完全相同的 39 个图名**。
- 唯一允许的差异要写成 `# R_FIG_DIFF: <图名> — <一句话原因>`（脚本内注释标记），
  门禁看到这个标记才放行 —— 与 `DYNAMIC_FIG_BASES` / `CONDITIONAL_FIGURES` 的豁免
  哲学一致：**豁免必须显式、带原因、可被 grep**。

---

## 5. 产物兼容边界

| 层 | 兼容性 | 说明 |
|---|---|---|
| 目录结构 | **逐字节同** | `data/`、`results/<dataset_id>/figures/`、`run_manifest.json`、`acceptance.json` |
| 配置 | **同一份文件** | `assets/*.yml`，R 用 `yaml::read_yaml()` |
| 图名 | **集合相等** | 39 个名字一致（§6.1 判据） |
| 图幅 | **同判据** | `W_SINGLE=89 / W_ONE_HALF=136 / W_DOUBLE=183 mm`（R 版公共库复刻 `mm()` 语义）；`tools/check_fig_sizes.mjs` 不区分语言 |
| 图例 | **同判据** | 框外右侧、单列（`check_legend_convention.mjs` 已支持 R 侧 —— geo 仓就是这么用的） |
| 调色 | **同判据** | `PAL[...]`，`check_palette.mjs` 需双语（R-09） |
| **中间对象** | **不跨语言** | Python 走 `.h5ad`，R 走 `.rds`。`INPUT_FILES`/`REQUIRED_FILES` 的**语义标签**（中文说明）必须一致，**文件名扩展允许不同** |
| 状态 JSON | **schema 同** | `qc_status.json` 等 8 个文件的**键集合**必须一致（R-08 的判据），值可因方法差异不同 |
| 验收 | **check id 集合同** | `acceptance.json` 的 `checks[].id` 集合必须一致（R-17 判据）—— 这样下游对两版是同一个接口 |

> **为什么不走 `.h5ad`**（`zellkonverter`/`anndata`）：那会拉进 basilisk/conda 运行时，
> CI 装环境又慢又脆，而且**用户根本不需要两版互换中间对象** —— 需要对照时走 CSV
> （`record_cross_language` 的既定约定）。所以 R 中间对象用 `.rds`，文档里写明
> "两版产物目录同名但中间对象格式不同，**对照请用 CSV 产物**"。

---

## 6. 门禁改动清单

### 6.1 `tools/check_fig_names.mjs`（R-09）
1. 重名判定按**语言分域**：`seen` 的键从 `<name>` 改成 `<lang>:<name>`（语言由文件扩展名推出）。
2. 新增判据 **`figures:parity`**：Python 版声明的图名集合与 R 版的必须**集合相等**；
   不等时，每一条差异都必须能在对应 `.R`/`.py` 文件里找到 `# R_FIG_DIFF:` 标记，
   否则判红（列出缺失/多余的名字）。
3. `PART`（阶段号）两版相同（都是 `02`）。

### 6.2 `tools/check_r_syntax.mjs`（R-09，扩展既有能力）
- 现有 R↔Python 参数/产物对齐（L277-361）已覆盖 tenifold 桥；扩展到
  R 版入口脚本 ↔ `scripts/main_analysis.py`：
  - `STEPS` 的 9 个 id 一致；
  - `STEP_STATUS_FILES` 的 8 键一致；
  - `HUMAN_REVIEW_NODES` 的 6 项一致；
  - `INPUT_FILES`/`REQUIRED_FILES` 的**条数**与**中文说明**一致（扩展名允许不同）。

### 6.3 `tools/check_palette.mjs`（R-09）
- 现在只认 Python 侧的 `PAL[...]`；要支持 R 侧 `PAL <- c(...)`。
- **硬约束：geo 仓已有可参考实现**（geo 是 R 仓，它的 palette 门禁就是 R 侧的）。

### 6.4 其余
- `check_figures.mjs` / `check_fig_sizes.mjs` / `check_fig_dpi.mjs` /
  `check_legend_convention.mjs` / `check_artifact_paths.mjs`：**不用改**（按目录扫）。
- `check_py_syntax.mjs`：不扫 R，不动。
- 新增一个 R 依赖清单门禁？**不新增** —— 版本清单放 scripts 目录下的 r_deps.R，
  由 CI 的 Rscript 调它打印版本，门禁只比对它输出的文本与清单一致
  （复用 geo 仓的版本探针做法）。

---

## 7. 验收层复刻要求（R-17 的规格）

R 版入口脚本的验收层必须复刻 `scripts/main_analysis.py` 的 `run_acceptance()`：

1. `chk(cid, kind, ok, detail, severity="required")` 的**调用位置参数含义保持一致**
   （第 2 个是 `kind`，不是 `severity` —— 这是本仓踩过的真实坑）。
2. `acceptance.json` 的 `checks[].id` 集合与 Python 版**完全一致**（静态判据比对两边的
   `chk("...")` 字面量）。
3. `n_checks` 数值可能不同（方法差异导致某些检查分支不同）—— **允许**，
   但**差异必须能被 `R_FIG_DIFF`/`R_CHK_DIFF` 类标记解释**（同 §6.1 的哲学）。
4. `HUMAN_REVIEW_NODES` 默认 `pending`，不许自动转 `confirmed`。

---

## 8. CI 设计（R-10 的规格）

**直接照抄 geo 仓已验证的模式**（`geo-normal-pipeline-skill/.github/workflows/geo_analysis.yml:89-164`）：

```yaml
- uses: r-lib/actions/setup-r@v2
  with: { r-version: 'release', use-public-rspm: true }
- uses: actions/cache@v4
  with:
    path: ${{ env.R_LIBS_USER }}
    key: rlib-${{ runner.os }}-scrna-r-v1     # 递增 v1/v2/... 才会重装
- uses: r-lib/actions/setup-r-dependencies@v2
  with: { cache-version: 1, extra-packages: <清单> }
- run: Rscript -e 'for (p in c("Seurat","harmony","scDblFinder",...)) cat(...)'
- run: Rscript scripts/main_analysis.R --config assets/config.pbmc3k.yml
- run: node tools/check_figures.mjs results/pbmc3k/figures
- run: node tools/check_fig_names.mjs           # 含 figures:parity
- run: node tools/check_r_syntax.mjs
- run: node tools/check_palette.mjs
```

- **为什么不用 renv**：geo 仓用 `setup-r-dependencies` + 手动 `actions/cache` 已经稳定
  跑了多轮（注释里写明"无 DESCRIPTION/lockfile 时它自带的缓存不生效"）。
  沿用已验证路径，不引入新变量。
- **新 workflow 文件**（R-10 建，写此文时还没有）：.github/workflows/scrna_r_analysis.yml，
  `paths` 过滤 `scripts/**/*.R` + `tools/**` + `assets/**` + 自身，
  另加 `workflow_dispatch`。**Python 文件的改动不触发 R job**（各跑各的，
  免得每次 Python 提交都白跑 10 分钟 R 环境）。
- artifact 名：scrna-r-results-<n>。

---

## 9. 任务分解（挂进 `governance/02_TASKLIST.md`）

| 编号 | 内容 | 依赖 |
|---|---|---|
| R-05 | 立项文档（本文档） | — |
| R-06 | scripts/lib/ 下的 R 公共库：`mm()/W_*/apply_style/save_fig/write_json/read_yaml/load_config/ensure_dirs/log_info/record_step/finite_round/probe_r_packages/init_manifest/record_cross_language` + 状态文件 schema 与 Python 对齐 | R-05 |
| R-07 | 9 个步骤脚本 `00-08_*.R`（含方法学差异标注，按 §3 的表逐行实现） | R-06 |
| R-08 | scripts 目录下的 R 版入口脚本：STEPS/状态文件/HUMAN_REVIEW/输入输出登记 + 验收层 | R-06, R-07 |
| R-09 | 门禁：`check_fig_names.mjs` 按语言分域 + `figures:parity`；`check_r_syntax.mjs` 扩到 main 层；`check_palette.mjs` 双语 | R-07 |
| R-10 | CI：R 侧 workflow 文件 + scripts 目录下的 r_deps.R | R-08 |
| R-17 | R 版验收（跑通 + 产物比对 + 亲读像素） | R-10 |

> spatial 仓对应 R-11..R-16（spatial 仓的立项文档 R-16 另写），其中
> **SpaGCN → BayesSpace/Giotto 的方法学差异表**是 spatial 侧的硬要求。

---

## 10. R1 纪律与本地可做的验证

- **本地不执行 R/Python 分析代码**（R1）—— 所有分析只在 GitHub Actions 上跑。
- 本地**可以**做（既定口径）：
  - 静态门禁（Node）：`check_r_syntax.mjs` / `check_fig_names.mjs` / `check_palette.mjs`
  - R 语法静态检查 + `py_compile` 等价的 R 语法 lint
  - 合成数据的 matplotlib 出图（验证 `save_fig`/`fit_fig_to_*` 等 helper）
  - 模块 import + 直接函数调用
- **R 侧本地没有 R 解释器** —— 所以 R 代码的第一次真实执行发生在 CI。
  这意味着 R-07 的每一步都必须在 CI 里逐步验证（先 `--selftest`，再真跑），
  不能一次写完九步再一起丢给 CI。

---

## 11. 未决问题（R-07/R-09 要回来拍板的）

1. **`scFates` 是否在 CRAN 上可装**（§3 05 行）—— **已核实（2026-09-26，R-07）：CRAN 上没有 `scFates`（HTTP 404，带 Chrome UA 复核非拦截），也没有任何同名 R 包 —— "CRAN 原版 scFates" 这个说法本身不成立**（Python `scFates` 从未发布 R 版）。按预案退回 Bioconductor：`slingshot`（HTTP 200）作主曲线树实现 + `TSCAN`（HTTP 200）备选；`destiny`（HTTP 200）承接 DPT。**裁决：05 行的 R 版方法集 = 自实现 GCS（与 Python 逐行同式）+ `destiny::DiffusionMap`+DPT + `slingshot`（主曲线拟时序，对应 scfates 槽位）；`stable_methods = ["dpt","slingshot","cytotrace"]`。** 方法学等级仍记 C（实现全换）；slingshot 的可复现性**没有六轮 CI 证据**，`reproducibility` 段的 stable/unstable 名单与证据文字必须改写、不能照抄 Python 版的 scfates 段。
2. **`liana` R 版从 GitHub 装**的版本漂移（§3 06 行）—— 是否要钉 commit。
3. **`AddModuleScore` 与 `score_genes` 的 bin 数**（§3 03 行）—— 配置里显式对齐，
   还是文档里标注差异。
4. **HVG 准则**（§3 02 行）—— `seurat_v3` 与 Seurat 的 `vst` 是否等价，
   需要一个实测对照。
5. `qc_status.json` 的双细胞字段两版数值不可直接比（§3 01 行）——
   验收层是否要给这条加 `kind="content"` 的降级标注。
6. **双版本图库的 `figure_review/` 镜像目录**：R 版产物与 Python 版产物同名，
   镜像按数据集存放会互相覆盖 —— 是按 `figure_review/<dataset>/` 分语言子目录，
   还是只镜像"当前语言"的产物（待 R-17 定）。
