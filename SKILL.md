---
name: scrna-pipeline-skill
description: Run an end-to-end single-cell RNA-seq analysis pipeline (QC, doublet detection, normalization, HVG, PCA, batch integration, Leiden clustering, cell-type annotation, marker genes, pseudobulk DESeq2 differential expression, cross-validated trajectory inference with DPT/Palantir/scFates/CytoTRACE plus direction correction and gene modules, ligand-receptor communication, TF regulon activity projected onto pseudotime) and ship it as a GitHub Actions workflow that uploads results as an artifact. Use when the user asks for single-cell analysis, scRNA-seq, 10x Genomics processing, scanpy clustering, cell type annotation, pseudobulk DE, trajectory inference, pseudotime, Monocle/Slingshot/Palantir/CytoTRACE alternatives in Python, RNA velocity questions, cell-cell communication, regulon or SCENIC-style analysis, or a reproducible cloud-run single-cell workflow.
---

# 单细胞转录组流水线

在 GitHub Actions 上跑完整单细胞分析，结果作为 artifact 下载。

## 这是一个**框架**，不是一条焊死的流水线

本仓库提供的是**生信分析的骨架与判据**：数据门禁、方法学约定、验收项、
产物清单。**具体跑什么由输入数据和配置文件决定**，步骤本身可增删 ——
加一步、换一种方法、关掉某个可选步骤，都是预期用法，不是"改坏了"。

所以「这个仓库能做什么」的答案在 `assets/config.*.yml` 和验收项里，
**不在目录结构里**。`scripts/` 中没有任何一处硬编码某个疾病或某个平台。

后续会有一个独立的「流水线编排模块」，让使用者挑选分析模块并串起来，
再与本仓库对接。**那部分不在本仓库职责范围内** —— 本仓库只负责把每一步做对。

## 怎么拿到它：云端仓库是唯一真源

本 skill **不需要"安装"**，也不依赖任何一台机器上的目录。真源是 GitHub 仓库：

    https://github.com/liubarryteb12/scrna-pipeline-skill

要用的时候从云端拉下来：

```bash
./use.sh                          # 拉取/更新到 ~/.cache/dsh-skills/，打印路径
./use.sh --register               # 需要本机 agent 直接发现它时才加
./use.sh --ref v1.0               # 钉住某一版
SKILL=$(./use.sh --print-path)    # 只取路径，便于脚本里用
./use.sh --clean                  # 清掉缓存副本并撤销注册
```

**拉取后会校验 `SKILL.md` 存在。** 远端改名或换结构时会明确报错，
而不是安静地给一个空目录 —— 实测过：`git clone` 失败时后面的步骤照样会跑，
最后就是靠这道校验拦住的。

> `use.sh` 里每个可能失败的步骤都显式 `|| die`，**不依赖 `set -e`**。
> 实测（bash 5.3）在 `resolved="$(pull)"` 这种「函数在命令替换里」的结构下，
> 函数内部的失败不一定会中止外层脚本。出错的路径必须自己说出来。

**与 `geo-normal-pipeline-skill` 是姊妹项目**，同一套工程约定：
云端跑、产物必上传、每个可选步骤的"没做"都要留记录。

## 什么时候用它

- "帮我分析这个单细胞数据集"
- "10x 数据怎么聚类/注释"
- "拟bulk 差异表达"
- "细胞通讯 / 轨迹 / 转录因子调控"
- 要一个**可复现、云端跑、带产物**的单细胞流程

## 快速开始

```bash
# 1. 用内置的健康 PBMC 数据集验证流程能跑
python scripts/main_analysis.py --config assets/config.pbmc3k.yml

# 2. 换成自己的数据：复制配置改 source 段
cp assets/config.pbmc3k.yml assets/config.my_data.yml
#    改 dataset_id 与 source.dataset（见 assets/datasets.yml）
python scripts/main_analysis.py --config assets/config.my_data.yml
```

结果在 `results/<dataset_id>/`，图在 `results/<dataset_id>/figures/`。

## 流水线步骤

| 步骤 | 脚本 | 产出 |
|---|---|---|
| 取数与硬校验 | `00_fetch.py` | `raw.h5ad`、`dataset_info.json` |
| 质控 | `01_qc.py` | `qc_filtered.h5ad`、`qc_status.json`、`qc_cells.csv` |
| 标准化与降维 | `02_integrate.py` | `integrated.h5ad`、`integration_status.json` |
| 聚类与注释 | `03_cluster_annotate.py` | `clustered.h5ad`、`markers_all.csv`、`celltype_annotation.csv` |
| 拟bulk 差异表达 | `04_pseudobulk_de.py` | `pseudobulk_de.csv`、`pseudobulk_status.json` |
| 轨迹推断 | `05_trajectory.py` | `paga_connectivities.csv`、`pseudotime_by_cluster.csv` |
| 细胞通讯 | `06_communication.py` | `cell_communication.csv` |
| 转录因子调控 | `07_grn.py` | `tf_regulons.csv`、`tf_activity_by_cluster.csv` |

编排器 `main_analysis.py` 跑完所有步骤后执行**验收清单**：
必需步骤是否成功、必需产物是否存在、图是否非空白。
任一项必需失败 → 退出码 1。

## 三条必须知道的事

### 1. 样本量的单位是生物学重复，不是细胞

在细胞层面直接跑检验（Wilcoxon 等）是把"细胞数"当成了"样本量"。
10000 个细胞来自 3 个供体，自由度是 2 不是 9997 —— 而检验会给出
`p = 1e-300`。**这是单细胞差异分析里最严重的系统性错误，而且它的
假阳性看起来极其显著。**

所以本流水线只做**拟bulk**：按 (样本 × 分组 × 细胞类型) 加总计数，
再用 pydeseq2 做负二项检验。没有 `design.sample_key` / `design.group_key`
时记 `not_configured` 并说明原因，**不退化成细胞层面的检验**。

### 2. 需要特定基因的分析必须用全基因集，不是 HVG

配体、受体、转录因子大多是**低表达**的，几乎不会进高变基因。
实测：在 HVG 子集（2000 基因）上，38 对配体-受体里只有 **3 对**可用；
在全集（13714 基因）上是 **27 对** —— 差 9 倍。

后果是"没找到显著通讯"变成**假阴性**，而假阴性看起来和"真的没有通讯"
一模一样。`06` 和 `07` 强制走 `adata.raw`，缺失时直接报错而不是
悄悄退回 HVG。

### 3. 注释和轨迹的结论强度不一样，产物里要能看出来

- **细胞类型注释**是打分提示，不是结论。每个簇的 assignment 带
  `score_margin`（第一名与第二名的差）。实测 pbmc3k 里簇 2 和 6 的
  margin 是 **0.0** —— 因为 `T_cell` 和 `CD4_T` 的 marker 完全重叠，
  打分法分不开。`assignment_confident: false` 如实标出。
- **轨迹**：没有 RNA 速率时拟时序**不能**说明方向性，只能说明相似度
  排序。`trajectory_status.json` 的 `limitations` 字段写明了这一点，
  以及根是怎么选的。

## 换成自己的数据

1. **数据集**：`assets/datasets.yml` 登记，或直接在 config 的 `source`
   段写 `kind` + `url`。支持 `10x_tar` / `10x_h5` / `h5ad`。
2. **分组**（要做拟bulk DE 才需要）：填 `design.sample_key`（哪列是
   供体/病人）、`design.group_key`（哪列是分组）、`design.reference_group`。
3. **批次**：多供体数据必须填 `design.batch_key` 并把
   `integration.method` 设为 `harmony`。单样本数据保持 `none` ——
   那是正确的选择，不是偷懒。
4. **细胞类型 marker**：`assets/celltype_markers.yml`。**内置那套是
   外周血/免疫为主的**，用到肿瘤、脑、肝等组织上会给出无意义的结果。

## 目录

```
assets/          配置、数据集登记表、marker 签名、配体受体库、TF 列表
scripts/         流水线步骤（00-07 + main_analysis）
scripts/lib/     公共库（日志、配置、出图、状态）
tools/           静态检查（语法、artifact 路径、图的空白检测）+ 自检
references/      方法学说明与故障排查
```

## 云端运行

推送到 GitHub 后 `.github/workflows/scrna_analysis.yml` 自动跑：
静态检查 → 装依赖 → 跑流水线 → 自检 → 检查图 → 上传 artifact。

**本地不需要装任何东西**，改代码推到云端看日志即可。
