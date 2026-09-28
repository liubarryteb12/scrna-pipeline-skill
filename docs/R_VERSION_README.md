# scRNA-seq 流水线 · R 语言版本介绍

> 生成日期：2026-09-28 · 对应 commit `65c236a` · **状态：代码全部落地，CI 排障收尾中**

单细胞 RNA-seq 分析流水线的 **R 语言版本**：与 Python 版**并存、二选一**，共用同一份
配置、同一套图名、同一个验收契约。用户可按语言偏好任选其一，产物互不依赖。

---

## 1. 这是什么

同一套 9 步单细胞分析（下载 → QC → 整合 → 聚类注释 → 拟bulk差异 → 轨迹 → 细胞通讯 →
转录因子调控 → 虚拟扰动），原来只有 Python/scanpy 实现；现在提供语义对齐的 R/Seurat
实现，特点：

- **配置零拷贝**：R 版用 R 的 `yaml` 包读**同一份** `assets/config.pbmc3k.yml`
  与 6 个资源表（datasets / celltype_mapping / celltype_markers / ligand_receptor /
  qc_gene_sets / tf_list），不需要维护两份配置。
- **一次运行只出一种语言的产物**：配置顶层 `pipeline.language: python | r` 选择；
  产物目录与文件名相同，不混装。
- **图名集合相等**：R 版产出的图名必须与 Python 版完全一致（静态门禁
  `figures:parity` + `# R_FIG_DIFF:` 豁免机制逐张对齐）→ caption CSV、叙事映射、
  SCI 稿件对接端口（SCI_HANDOFF）对两版是**同一个**。
- **方法学差异全部显式登记**：每一步按 D/A/B/C 差异等级标注（见 §4），
  状态 JSON 里照旧写明实际用的哪条路径，不冒充。

## 2. 代码结构（当前规模）

```
scrna-pipeline-skill/
  scripts/
    main_analysis.R        # 编排入口（~1200 行）：步骤表/验收层 23 组检查/退出码
    00_fetch.R             # 数据下载 + 读入（download.file + Read10X + rhdf5）
    01_qc.R                # QC 指标 + scDblFinder 双细胞检测 + 过滤（5 张小提琴图）
    02_integrate.R         # Normalize/HVG(vst)/PCA + ComBat/Harmony（Seurat 3/5 双兼容）
    03_cluster_annotate.R  # 邻居+UMAP+bluster Leiden 扫描 + FindAllMarkers + AddModuleScore
    04_pseudobulk_de.R     # 拟bulk 聚合 + DESeq2（母实现，与 pydeseq2 同法）
    05_trajectory.R        # GCS(自实现) + destiny DPT + slingshot + cytoTRACE(自实现)
    06_communication.R     # liana(GitHub 版) + 配体受体共表达打分
    07_grn.R               # 共表达推断(无 motif 剪枝) + AUCell 式打分
    08_virtual_perturbation.R  # TenifoldNet 虚拟扰动（与 Python 逐字节同款）
    lib/common.R           # 基础设施 ~1290 行：配置/日志/落盘/图/状态/审计
    lib/tenifold_knk.R     # 既有资产 567 行，两版共用（08 的引擎）
    r_deps.R               # CI 包安装脚本（18 包，唯一事实源）
  references/r_version.md  # R 版规格书（316 行，11 节）：立项依据/选型表/契约/门禁/CI
  .github/workflows/
    scrna_r_analysis.yml   # R 版 CI（geo 模式：setup-r + rspm + 手动缓存）
```

代码量：9 步脚本 + main_analysis.R ≈ **6050 行**，`lib/common.R` ≈ **1290 行**
（Python 版对照：9 步 5654 行 + main_analysis.py 1270 行 + common.py 1589 行）。

## 3. 运行方式

```bash
# 本地 / 服务器
Rscript scripts/main_analysis.R --config assets/config.pbmc3k.yml

# 或只跑部分步骤（逗号分隔）
Rscript scripts/main_analysis.R --config assets/config.pbmc3k.yml --steps 00_fetch,01_qc

# 自检（验收层自测，不跑分析）
Rscript scripts/main_analysis.R --selftest

# 云端（推荐）：推到 GitHub 即自动触发 R job，结果在 Actions artifact 下载
```

CI（`.github/workflows/scrna_r_analysis.yml`）：
- **触发隔离**：只认 `scripts/**/*.R`、tools、assets 与 workflow 自身的改动；
  Python 文件提交不会触发 R job，反之亦然。
- 环境：ubuntu-latest + R release + `{RSPM}` 二进制源 + 手动缓存 key
  `rlib-<os>-scrna-r-v1`；18 个包经 `r_deps.R` 固定版本安装；liana 走
  `remotes::install_github("saezlab/liana")`（`GITHUB_TOKEN` 显式注入）。
- 产物：`scrna-r-results-<n>` artifact（图 PNG/PDF + 全部状态 JSON + CSV + run manifest）。

## 4. 方法学差异（D/A/B/C 等级表，全表见 `references/r_version.md` §3）

| 步 | R 选型 | 等级 | 一句话差异 |
|---|---|---|---|
| 00 fetch | download.file + Read10X/hdf5r + GEOquery | B | 读入坐标约定必须与 Python 一致 |
| 01 qc | **scDblFinder**（Python: scrublet） | **C** | 双细胞检测器不同，检出数不可直接比 |
| 02 integrate | Seurat Normalize/FVF/PCA + sva::ComBat + harmony | A | HVG 准则必须配置对齐 |
| 03 cluster | Seurat + bluster Leiden + AddModuleScore + SingleR | A/B | bin 数已显式对齐；AddModuleScore 空分位行为差异已记录 |
| 04 pseudobulk | **DESeq2**（pydeseq2 的母实现） | A | 同负二项 GLM，末位差可接受 |
| 05 trajectory | destiny DPT + **slingshot**（Palantir 无 R 实现→丢弃）+ 自实现 GCS/cytoTRACE | **C** | 方法 4→3，≥2 方法一致性检查仍成立 |
| 06 communication | liana（GitHub saezlab/liana） | A/B | 版本漂移风险；退回逻辑与 Python 同构 |
| 07 grn | cor() 共表达 + AUCell 式（**不升级成 SCENIC**） | **A** | 升级会作废两版可比性，硬约束 |
| 08 perturbation | 同一份 `lib/tenifold_knk.R` | **D** | 零差异，连数值都应逐位一致 |

> **为什么 07 坚持"不升级"**：R 版的价值是"用户可选实现"，不是"换个方法重算"。
> 换了 SCENIC，同一行 TF 调控在两版里的**含义**都不一样，两版对照就是假对照。

## 5. 质量保障（与 Python 版同款纪律）

**静态门禁（全部本地可跑，Node 实现）**：
- `tools/check_r_syntax.mjs` —— R 语法 + **R↔Python 参数/产物名对齐**检查
- `tools/check_fig_names.mjs` —— 出图命名（按语言分域，R/Python 同名不误报）
- `tools/check_palette.mjs` —— 调色板纪律（PAL 常量，无字面 hex）
- `tools/check_doc_refs.mjs` —— 文档引用逐条真实存在（当前 272 条全绿）

**验收层（`main_analysis.R` 内置，与 Python 版 23 组检查 id 完全一致）**：
- 步骤状态三态（PASS/FAIL/INFO）+ 嵌套状态扫描（`*status.json` 逐个判红）
- 图产物三桶（missing/waived/skipped_optional）+ 动态图名计数 + 超宽图检测
- 轨迹内容级检查（方法数≥2、方向来源、limitations≥4）+ CellTypist/LIANA 双层校验
- 诚实性语义：跳过必须给 reason、空敲除≠效应 0、HUMAN_REVIEW pending 不自动确认

**双向标定**：每个新检查都做过"注入缺陷→必须红 / 干净基线→必须绿"两向验证；
已抓出并在落地前修掉 9+ 处真缺陷（如 R `[[` 越界会崩掉整层验收）。

## 6. 当前状态（截至 2026-09-28，commit 65c236a）

**已落地**：9 步 R 脚本 + common.R + main_analysis.R（验收层 23 组检查）+ 4 静态门禁
全绿 + CI workflow + r_version.md 规格书 + README R 节 + AGENTS.md 修订。

**CI 排障进行中**（GitHub Actions 实跑逐层剥皮，每轮修复一个真实 API/语义缺陷，
已 25 轮）：fetch/QC/整合/聚类注释/**拟bulk DE**/**轨迹**/**细胞通讯**/**GRN** 全部跑通，
关键修复包括：Seurat 5 Assay5 槽位兼容、矩阵方向语义（细胞×基因 vs 基因×细胞）、
bluster NNGraphParam 真实签名、slingshot 矩阵直传、AddModuleScore 零表达基因过滤等。

**剩余收尾项**（详见 `governance/02_TASKLIST.md` R-10b 行）：
- 虚拟扰动（08）步骤 CI 首验证（其引擎 tenifold_knk.R selftest 已单独通过）
- slingshot 伪时序修复的 CI 确认（SCE 包装已弃用改矩阵直传）
- 验收层偶发红项清零（跨语言交接计数等展示性检查）

**已知的两版行为差异**（已记录于 `references/r_version.md` §11，不影响结论形态）：
- AddModuleScore 空分位严格报错 vs scanpy score_genes 静默容忍
- liana 版本漂移（GitHub 包，无版本 pin，已评估）
- 05 轨迹方法集差异（Palantir/scFates 无 R 实现 → dpt+slingshot 替换）

## 7. 文档与治理

- 规格书：`references/r_version.md`（R 版唯一事实源：立项核实、选型表、命名契约、
  产物兼容边界、门禁改动、验收层规格、CI 设计、R1 纪律、未决问题）
- 任务台账：`governance/02_TASKLIST.md` R-04…R-10 行（全程可追溯：
  每一轮 CI 排障的根因/修复 commit/证据）
- 缺陷台账：`governance/15_ERROR_LEDGER.md`（R 版排障沉淀的通用教训，
  如"长度恰好匹配时语义反转不报错比不匹配更危险"）
