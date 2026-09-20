# 单细胞转录组流水线

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

在 GitHub Actions 上跑的完整单细胞分析流程。云端运行，结果作为 artifact 下载。

**姊妹项目**：[`geo-normal-pipeline-skill`](https://github.com/liubarryteb12/geo-normal-pipeline-skill)
（GEO 芯片数据挖掘）。两者共用同一套工程约定。

---

## 它做什么

```
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
04 拟bulk DE       → pydeseq2（样本量单位是生物学重复）
05 轨迹            → PAGA + 四种拟时序交叉验证（方向校正 + 根的选择记录）
06 细胞通讯        → 配体-受体 + 置换检验 + BH 校正
07 转录因子调控    → 共表达推断 + AUCell 式活性打分
   ↓
验收清单 + artifact
```

## 快速开始

**本地验证流程**（需要 Python 3.10-3.12）：

```bash
pip install -r requirements.txt
python scripts/main_analysis.py --config assets/config.pbmc3k.yml
```

**云端运行**（推荐）：推到 GitHub，workflow 自动跑，结果在 Actions 页面的
artifact 里下载。

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
  sample_key: donor                # 哪一列是供体/病人
  group_key: condition             # 哪一列是分组
  reference_group: control         # 对照组取值
  batch_key: batch                 # 多供体时填
integration:
  method: harmony                  # 多供体时必须
```

然后 `python scripts/main_analysis.py --config assets/config.my_data.yml`。

## 三条最容易踩的坑

### 1. 样本量是生物学重复，不是细胞

在细胞层面跑检验会把"细胞数"当成"样本量"，给出 `p = 1e-300` 的
假阳性 —— 而且看起来极其显著。本流水线只做拟bulk。

### 2. 配体/受体/转录因子不在 HVG 里

实测：HVG 子集上 38 对配体-受体只有 3 对可用，全集上 27 对。
用 HVG 做通讯分析会得到大量假阴性，而假阴性看起来和"真的没有"一样。

### 3. 注释和轨迹的结论比聚类弱

- 注释是打分提示。`celltype_annotation.csv` 的 `score_margin` 告诉
  你第一名和第二名差多少；margin ≈ 0 时那个 assignment 不该被当结论。
- **拟时序的符号是任意的。** 实测四种方法（DPT / Palantir / scFates /
  CytoTRACE）在 PBMC3k 上的**原始**拟时序两两相关从 **−0.50 到 +0.86**
  —— 符号都不一样。本流水线按「值越大越晚」统一方向后才互相比较，
  统一后的三方平均 rho = **+0.636**。只跑一种方法时，报出来的
  "轨迹"是某个算法的一次输出，不是数据里的结构。
- **方向参考不能算进一致性统计。** 未配置 marker 时退回 CytoTRACE 作
  方向参考，它与参考的相关恒为 ±1 —— 那是定义不是证据，必须排除。
- 没有 RNA 速率时，拟时序只能说相似度排序，**不能**说分化方向。
  `scvelo` 记 `not_done`（输入只有一套计数，没有 spliced/unspliced）。
- `regulon × 拟时序` 的"显著个数"在 n 大时很廉价：n=2652 时 |rho| 只要
  约 0.06 就能过 BH<0.05（201 个里 190 个"显著"）。要看效应量分布。

## 产物

`results/<dataset_id>/`：

| 文件 | 内容 |
|---|---|
| `acceptance.json` | 验收清单结论 |
| `state.json` | 每步的状态与耗时 |
| `qc_status.json` | QC 结论（含双细胞、ambient RNA 的实际状态） |
| `integration_status.json` | HVG 口味、PCA 方差、整合方法 |
| `cluster_status.json` | 簇数、分辨率扫描、注释（含 margin） |
| `markers_all.csv` | 每簇的 marker 基因 |
| `celltype_annotation.csv` | 注释 + margin + 置信标记 |
| `pseudobulk_de.csv` | 拟bulk 差异表达（需要分组信息） |
| `paga_connectivities.csv` | PAGA 连通性 |
| `pseudotime_by_cluster.csv` | 每簇拟时序分布 |
| `cell_communication.csv` | 配体-受体打分 + 置换 p 值 + BH |
| `tf_regulons.csv` | 调控子 + 靶基因 + 簇特异性 |
| `figures/` | 14 张图（PNG + PDF） |

## 工程约定

见 [`AGENTS.md`](AGENTS.md)。要点：

- 可选步骤的"没做"必须和"做了没问题"长得不一样
- 每个结论都带上它的适用范围（`limitations` 字段）
- 静态检查在前、装依赖在后
- artifact 清单与脚本实际写出的文件由 `tools/check_artifact_paths.mjs` 交叉核对

## 目录

```
assets/      配置、数据集登记表、marker 签名、配体受体库、TF 列表
scripts/     00-07 步骤 + main_analysis 编排器
scripts/lib/ 公共库
tools/       静态检查 + 拟bulk 自检
references/  方法学与故障排查
```

## 引用

方法学参考见 [`references/methods.md`](references/methods.md)。
