# 单细胞转录组流水线

在 GitHub Actions 上跑的完整单细胞分析流程。云端运行，结果作为 artifact 下载。

**姊妹项目**：[`geo-brca-microarray-skill`](https://github.com/liubarryteb12/geo-brca-microarray-skill)
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
05 轨迹            → PAGA + DPT（带根的选择记录与方法学限定）
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
- 没有 RNA 速率时，拟时序只能说相似度排序，**不能**说分化方向。

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
