# 单细胞转录组流水线

<p align="center">
  <a href="https://github.com/liubarryteb12/scrna-pipeline-skill/actions/workflows/scrna_analysis.yml"><img alt="scRNA pipeline (Python)" src="https://github.com/liubarryteb12/scrna-pipeline-skill/actions/workflows/scrna_analysis.yml/badge.svg?branch=main"></a>
  <a href="https://github.com/liubarryteb12/scrna-pipeline-skill/actions/workflows/scrna_r_analysis.yml"><img alt="scRNA pipeline (R)" src="https://github.com/liubarryteb12/scrna-pipeline-skill/actions/workflows/scrna_r_analysis.yml/badge.svg?branch=main"></a>
  <a href="https://github.com/liubarryteb12/scrna-pipeline-skill/actions/workflows/docs_check.yml"><img alt="Docs check" src="https://github.com/liubarryteb12/scrna-pipeline-skill/actions/workflows/docs_check.yml/badge.svg?branch=main"></a>
  <img alt="Python" src="https://img.shields.io/badge/Python-3.10--3.12-3776AB">
  <img alt="R" src="https://img.shields.io/badge/R-4.3%2B-276DC3">
</p>

> **L2 分析线：单细胞转录组。** 覆盖 QC、双细胞检测、整合、聚类注释、拟bulk 差异表达、轨迹、细胞通讯、转录因子调控与虚拟扰动，提供 Python 与 R 双版本，产物以 GitHub Actions artifact 形式交付。

## 三线体系

| 线 | 仓库 | 数据 | 核心问题 |
|---|---|---|---|
| L1 | [`geo-normal-pipeline-skill`](https://github.com/liubarryteb12/geo-normal-pipeline-skill) | bulk 芯片 | 差异表达、通路富集、共表达、预后签名 |
| **L2** | **本仓库：`scrna-pipeline-skill`** | 单细胞 | 细胞类型、轨迹、通讯、调控 |
| L3 | [`spatial-pipeline-skill`](https://github.com/liubarryteb12/spatial-pipeline-skill) | Visium 空间转录组 | 空间域、空间高变基因、邻域、空间通讯 |

三条线共用同一套工程约定：配置驱动、状态落盘、验收判据、静态门禁。本仓库负责把单细胞分析的样本单位、基因集与结论强度处理正确。

## 这个仓库是什么

- **框架，不是焊死的流水线。** 步骤、参数、产物与验收项都随分析需求变化。加一步、换一种方法、关掉某个可选步骤，都是预期用法。
- **样本量的单位是生物学重复。** 差异表达只做拟bulk，不在细胞层面跑检验。
- **双语言实现。** Python 版用于日常与 CI 主路径，R 版保留同一套分析语义，产物与图名集合对齐。
- **云端仓库是唯一真源。** 本 skill 不需要安装，也不依赖某台机器上的目录，按需从 GitHub 拉取即可。

## 分析流程

```text
10x / h5ad 数据
   ↓
00 取数与硬校验    → 拒绝非整数计数矩阵
   ↓
01 质控            → 硬阈值 + scrublet 双细胞检测 + ambient RNA 记录
   ↓
02 标准化与降维    → HVG(seurat_v3) → PCA → 批次整合
   ↓
03 聚类与注释      → Leiden + 分辨率扫描 + marker + 细胞类型打分(带 margin)
   ↓
04 拟bulk DE       → pydeseq2，样本量单位是生物学重复
   ↓
05 轨迹            → PAGA + 四种拟时序交叉验证，方向校正与根的选择记录
   ↓
06 细胞通讯        → 配体-受体 + 置换检验 + BH 校正
   ↓
07 转录因子调控    → 共表达推断 + AUCell 式活性打分
   ↓
08 虚拟扰动        → scTenifold 网络敲除，空敲除按结构事实施加标记
   ↓
验收清单 + artifact
```

可选步骤的失败不会让 CI 变红，但每一步都必须在状态文件中写下 `status` 与原因。`not_done` 表示明确记录“做不了以及为什么”，与“做了没问题”严格区分。

## Python / R 双版本

R 版与 Python 版并存，任选其一，产物互不依赖。

```bash
Rscript scripts/main_analysis.R --config assets/config.pbmc3k.yml
```

- CI 入口：`.github/workflows/scrna_r_analysis.yml`。R job 只认 `.R`、`tools`、`assets` 与本 workflow 的改动，Python 文件的提交不会触发 R job，反之亦然。
- 依赖清单：`scripts/r_deps.R`，共 18 个包；liana 走 GitHub。装包的唯一事实源是 workflow 的包名清单。
- 差异记录在 [`references/r_version.md`](references/r_version.md)，按 D 零差异、A 同方法不同实现、B 换实现同思路、C 真方法差异分级。示例：01 双细胞检测用 scDblFinder（Python 是 scrublet）；05 轨迹用 destiny DPT 与 slingshot（Palantir、scFates 无 R 实现）；08 虚拟扰动与 Python 共用 `scripts/lib/tenifold_knk.R`，逐字节同款。
- 图名集合与 Python 版相等，`R_FIG_DIFF` 豁免机制保证逐张对齐。R 版中间对象走 `.rds`，不跨语言读 `.h5ad`，对照表走 CSV。

## 快速开始

### Python 版

```bash
pip install -r requirements.txt
python scripts/main_analysis.py --config assets/config.pbmc3k.yml
```

### R 版

```bash
Rscript scripts/main_analysis.R --config assets/config.pbmc3k.yml
```

推送到 GitHub 后，workflow 会自动运行，结果在 Actions 页面的 artifact 里下载。

### 作为 Agent Skill 使用

```bash
./use.sh                          # 拉取或更新到 ~/.cache/dsh-skills/，打印路径
./use.sh --register               # 需要本机 agent 直接发现时加
./use.sh --ref v1.0               # 钉住某一版
SKILL=$(./use.sh --print-path)    # 只取路径，便于脚本使用
./use.sh --clean                  # 清理缓存副本并撤销注册
```

拉取后会校验 `SKILL.md` 存在。远端改名或换结构时会明确报错，不会安静地给一个空目录。

## 换成自己的数据

```bash
cp assets/config.pbmc3k.yml assets/config.my_data.yml
```

改这几处：

```yaml
dataset_id: my_data
source:
  kind: h5ad                       # 10x_tar | 10x_h5 | h5ad
  url: https://.../my_data.h5ad
design:
  sample_key: donor                # 哪一列是供体或病人
  group_key: condition             # 哪一列是分组
  reference_group: control         # 对照组取值
  batch_key: batch                 # 多供体时填
integration:
  method: harmony                  # 多供体时必须
```

然后运行 `python scripts/main_analysis.py --config assets/config.my_data.yml`。

## 三条最容易踩的坑

### 1. 样本量是生物学重复，不是细胞

在细胞层面跑检验会把细胞数当成样本量，10000 个细胞来自 3 个供体时自由度是 2，不是 9997，而 Wilcoxon 会给出 `p = 1e-300`。这类假阳性看起来极其显著。本流水线只做拟bulk。

### 2. 配体、受体、转录因子不在 HVG 里

实测 HVG 子集上 38 对配体-受体只有 3 对可用，全集上 27 对。用 HVG 做通讯分析会得到大量假阴性，而假阴性看起来和“真的没有”一样。需要特定基因的分析统一走 `adata.raw` 全基因集。

### 3. 注释和轨迹的结论比聚类弱

- 注释是打分提示。`celltype_annotation.csv` 的 `score_margin` 告诉你第一名与第二名差多少；margin 接近 0 时，那个 assignment 不该被当结论。`margin_state` 把“不确定”（`low_margin`）与“算不出来”（`single_celltype`、`margin_undefined`）分开，后者的 `assignment_confident` 是 `null` 而不是 `false`。
- 拟时序的符号是任意的。实测四种方法（DPT、Palantir、scFates、CytoTRACE）在 PBMC3k 上的原始拟时序两两相关从 -0.50 到 +0.86，符号都不一样。本流水线按“值越大越晚”统一方向后再互相比较，统一后的三方平均 rho 为 +0.636。
- 方向参考不能算进一致性统计。未配置 marker 时退回 CytoTRACE 作方向参考，它与参考的相关恒为 ±1，那是定义不是证据，必须排除。
- 没有 RNA 速率时，拟时序只能说相似度排序，不能说分化方向。`scvelo` 记 `not_done`，因为输入只有一套计数，没有 spliced 与 unspliced。
- `regulon × 拟时序` 的“显著个数”在 n 大时很廉价。n=2652 时 |rho| 只要约 0.06 就能过 BH<0.05，201 个里 190 个“显著”。要看效应量分布。

## 产物

`results/<dataset_id>/`：

| 文件 | 内容 |
|---|---|
| `acceptance.json` | 验收清单结论 |
| `state.json` | 每步状态与耗时 |
| `qc_status.json` | QC 结论，含双细胞与 ambient RNA 的实际状态 |
| `integration_status.json` | HVG 口味、PCA 方差、整合方法 |
| `cluster_status.json` | 簇数、分辨率扫描、注释与 margin |
| `markers_all.csv` | 每簇的 marker 基因 |
| `celltype_annotation.csv` | 注释、margin、置信标记 |
| `pseudobulk_de.csv` | 拟bulk 差异表达，需要分组信息 |
| `paga_connectivities.csv` | PAGA 连通性 |
| `pseudotime_by_cluster.csv` | 每簇拟时序分布 |
| `cell_communication.csv` | 配体-受体打分、置换 p 值、BH |
| `tf_regulons.csv` | 调控子、靶基因、簇特异性 |
| `figures/` | 19 张图，PNG + PDF，图名格式见 AGENTS 规则 22 |

## 质量门禁

```bash
node tools/check_py_syntax.mjs
python tools/check_py_names.py .
node tools/check_r_syntax.mjs
node tools/check_fig_names.mjs
node tools/check_palette.mjs
node tools/check_legend_convention.mjs
node tools/check_figures.mjs results/pbmc3k/figures
node tools/check_fig_sizes.mjs results/pbmc3k/figures
node tools/check_fig_dpi.mjs results/pbmc3k/figures
node tools/check_artifact_paths.mjs
node tools/check_doc_refs.mjs
```

R workflow 额外跑 `node tools/check_artifact_paths.mjs --lang r` 与 `Rscript scripts/lib/tenifold_knk.R --selftest`。CI 还会跑图例、空白图、文档引用的内建自检，避免“检查了但什么都没查”的假阴性。

## 目录结构

```text
assets/       配置、数据集登记表、marker 签名、配体受体库、TF 列表
scripts/      00-08 步骤脚本 + main_analysis 编排器，Python 与 R 双版
scripts/lib/  公共库，含 tenifold_knk.R
tools/        静态检查、拟bulk 自检、tenifold 自检
references/   methods.md、r_version.md、module0.md、troubleshooting.md
docs/         R_VERSION_README.md、figure_groups.md
AGENTS.md     工程规则与踩坑记录
SKILL.md      Agent Skill 元数据
```

## 文档导航

| 文档 | 内容 |
|---|---|
| [`AGENTS.md`](AGENTS.md) | 工程规则、验收判据、踩坑记录 |
| [`references/methods.md`](references/methods.md) | 方法学依据与参数 |
| [`references/r_version.md`](references/r_version.md) | R 版与 Python 版的逐条差异 |
| [`references/module0.md`](references/module0.md) | 模块零运行时规范 |
| [`references/troubleshooting.md`](references/troubleshooting.md) | 故障排查 |
| [`docs/R_VERSION_README.md`](docs/R_VERSION_README.md) | R 版运行说明 |

## 赞赏支持

如果这个仓库帮你节省了时间，欢迎请作者喝杯咖啡。赞赏完全自愿，不影响代码使用与后续更新。

<p align="center">
  <img src="assets/sponsor-qr.png" alt="伊霖的赞赏码" width="260">
</p>

> 感谢支持，我会尽快更新代码，维护代码。

## 授权

本仓库暂未附 LICENSE 文件。如需商用或二次分发，请先联系作者确认。
