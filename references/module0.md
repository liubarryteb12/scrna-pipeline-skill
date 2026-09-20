# 模块零：运行清单

参考规范：三大部分整合文档「模块零：语言与运行时规范」。
姊妹项目 `geo-normal-pipeline-skill` / `spatial-pipeline-skill` 的
`references/module0.md` 是同一份文档，接口**同名同义**，
三部分的清单可以并排读。

---

## 0. 为什么要有这一层

**没有运行清单的分析结果不是结果。**

半年后拿到一份 `cell_communication.csv`，如果不知道当时装的是哪个版本的
liana、输入的 h5ad 是哪个哈希、随机种子是多少、置换跑了几次，
那份 CSV 就**无法被复现，也无法被质疑** —— 而不可质疑的结论没有价值。

这一层不产出任何生物学结论，它只回答一个问题：
**"这份结果是在什么条件下跑出来的？"**

产物：`results/<dataset_id>/run_manifest.json`

**为什么不塞进 `state.json`：** `state.json` 记的是"这一步跑没跑成"，
每步重写；manifest 记的是"本轮是在什么条件下跑出来的"，是证据，
写入后不该再变。混在一起会让后者被前者覆盖。

---

## 1. 规范条款与接口的对应

| 条款 | 要求 | 接口 | 落到 manifest 的哪个字段 |
|---|---|---|---|
| §0.3 | `pip freeze` **全量**输出 | `capture_versions()` | `versions` |
| §0.3 | 关键工具**逐个**记版本 | `capture_versions()` 的 `KEY_PACKAGES` | `key_versions` |
| §0.3 | 所有随机过程固定种子并记录 | `record_params()`（`params` 里含 `seed`）+ `m["seed"]` | `params`、`seed` |
| §0.4 | 输入数据哈希 | `record_input()` | `inputs` |
| §0.4 | 全部关键参数 | `record_params()` | `params` |
| §0.4 | Agent 决策链 | `record_decision()` | `decisions` |
| §0.4 | 人工干预记录 | `record_human_review()` | `human_review` |
| §0.2 | 跨语言转换前后维度、丢失字段 | `record_cross_language()` | `cross_language` |

---

## 2. 接口清单

| 函数 | 作用 | 备注 |
|---|---|---|
| `MANIFEST_NAME` | `"run_manifest.json"` | 三个仓库一致 |
| `KEY_PACKAGES` | 文档点名的关键工具 | 见 §4 |
| `NAMED_TOOLS` | 点名但用不了的工具 + 理由 | 见 §5 |
| `manifest_path(cfg)` | 清单路径 | |
| `read_manifest(cfg)` | 读清单；文件不存在返回 `{}` | |
| `init_manifest(cfg, language="python")` | **开新一轮，清掉上一轮** | 见 §3.1 |
| `capture_versions(cfg, ...)` | 全量已装包 + `KEY_PACKAGES` 逐个 | |
| `record_input(cfg, path, ...)` | 输入文件的 sha256 与字节数 | |
| `record_params(cfg, params)` | 参数（`update` 合并，可多次调用） | |
| `record_decision(cfg, node, q, a, evidence)` | 决策链 | |
| `record_human_review(cfg, node, required, status, note)` | 人工复核节点 | |
| `record_cross_language(cfg, src, dst, fmt, before, after, lost, ...)` | 跨语言转换 | |
| `manifest_summary(cfg)` | 供验收用的摘要 | 见 §3.3 |
| `probe_named_tools(log, only)` | 把 `NAMED_TOOLS` 整理成可写进状态 JSON 的表 | |
| `named_tools_note()` | 一句话说明为什么多数点名工具没用上 | |

---

## 3. 四条不能省的约定

### 3.1 `init_manifest` 必须清掉上一轮

上轮的清单留在那里冒充本轮，**比没有清单更糟** —— 它看起来是证据，
实际是过期证据。这和 AGENTS 规则 14（每步开跑前先删自己的状态文件）
是同一个道理。

### 3.2 输入哈希在**步骤跑完之后**才登记

可选步骤这轮有没有产物，跑完才知道；在开头登记会把"上轮残留"记成本轮输入。

### 3.3 `manifest_summary` 同时报两个缺失数

```
inputs_missing           所有缺失的输入（可见，供人看）
inputs_missing_required  只有 required=True 的缺失（供验收判 FAIL）
```

**第一版只有一个数，把可选输入缺失也算成失败**，结果每轮都红。
可选输入缺失是**正常状态**（Part 1 的交接表没配、外部队列不存在），
判成 FAIL 会让红灯失去意义。但它必须**可见** —— 所以两个数都要有。

### 3.4 未装的工具要记成 `None`，不能省略键

```json
"key_versions": { "liana": "1.4.0", "celltypist": null }
```

`"celltypist": null` 和"没有 celltypist 这个键"是两件事：前者是
"**查过了，没装**"，后者是"**没查**"。省略会让读者分不清。

### 3.5 `pip freeze` 不用 `subprocess` 抓

**沙箱下管道捕获会 EPERM**（`child_process` 默认 `stdio: 'pipe'`）。
用 `importlib.metadata` 枚举已装包。
包名归一化走 **PEP 503**（`_norm_pkg`）：`scikit-misc` / `scikit_misc` /
`Scikit.Misc` 是同一个包，不归一会让 `key_versions` 里出现查不到的键。

---

## 4. `KEY_PACKAGES`：三部分的清单

文档点名的工具**逐个列出**，装了的记版本、没装的记 `null`。

| 部分 | 数量 | 内容 |
|---|---|---|
| Part 1 | 15 | `GEOquery` `limma` `WGCNA` `clusterProfiler` `GSVA` `glmnet` `survival` `survminer` `timeROC` `rms` `STRINGdb` + `TRRUST` `ChEA3` + `scTenifoldKnk` `PerturbNet` `RegVelo` |
| **Part 2（本仓库）** | 25 | `scanpy` `anndata` `scvi-tools` `cellbender` `harmonypy` `scvelo` `celltypist` `pyscenic` `liana` `doubletdetection` `scrublet` + 拟时序 `palantir` `scfates` `cytotrace` + R 包 `monocle3` `slingshot` `cellchat` `soupx` `scdblfinder` + `scTenifoldKnk` `PerturbNet` `RegVelo` |
| Part 3 | 16 | `SpatialDE` `SpatialDE2` `spacexr` `BayesSpace` `SPARK-X` `cell2location` `STAGATE` `SpaGCN` `SpaceFlow` `stLearn` `ISORT` `Bering` `BOMS` + LIANA 等 |

**Part 2 的清单里为什么有 R 包（`monocle3` / `slingshot` / `cellchat` /
`soupx` / `scdblfinder`）：** 规范说 Part 2 是"Python (+R via rpy2)"，
但**本仓库 CI 里没有 R**，所以这条路径实际不存在。这些键**本来就该是
`null`** —— 但必须留着，否则读者不知道"是没查还是不该有"。
确切原因写在 `NAMED_TOOLS` 里（见 §5）。

---

## 5. `NAMED_TOOLS`：点名工具"为什么没用上"的登记

**这是本仓库最容易被误读的地方。** 每一步都有产物、都有图 ——
`02_integrate.py` 有归一化、`04_pseudobulk_de.py` 有差异表、
`05_trajectory.py` 有四条轨迹、`07_grn.py` 有调控子 ——
看起来"§2 做完了"，而文档点名的方法**多数没用上**。

所以把"为什么没用上"逐条落盘。**判据是"理由写了没有"，
不是"工具跑了没有"。** 将来某个工具能装了，验收应该依然 PASS
（理由变成"已装"），而不是因为 `available=False` 就变红 ——
那会把"如实记录"惩罚成失败。

四类 `kind`：

| kind | 含义 | 本仓库的例子 |
|---|---|---|
| `r_package` | R/Bioconductor 包，CI 无 rpy2 | `SoupX` `SCTransform` `scran` `DESeq2` `Monocle3` `CellChat` |
| `not_on_pypi` | 真包不在 PyPI | `CytoTRACE2` |
| `deps` / `needs_*` | PyPI 有真包，依赖链或资源跑不动 | `CellBender`（GPU）`scVelo`（缺 spliced/unspliced 层）`pySCENIC`（GB 级 motif 库）`SingleR`（缺参考数据集） |
| `name_taken` | **PyPI 上那个名字是另一个不相干的包** | `edgeR` `Slingshot` |

### `name_taken` 是最危险的一类

因为 `pip install` 会**成功**，装进来的是完全无关的东西。实测元数据：

```
pip install edgeR     -> "Redirect Microsoft Edge to your preferred browser"
                         （作者 Daniel Agans，github.com/phwelo/edger）
pip install slingshot -> "Index Migration for ElasticSearch"
                         （作者 Thierry Jossermoz）
```

`edgeR` 是 §2.5 点名的差异分析工具，`slingshot` 是 §2.6 点名的轨迹工具。
**这两个名字绝对不能进 `requirements.txt`。**
装不上会立刻报错，**装错了要到 import 或跑出结果才发现** ——
而那时"结果"可能已经在图上看着挺像回事了。

**反例：`SingleR` 不是顶名的。** PyPI 上的 `SingleR` 0.5.0 是
**BiocPy/singler**，R 那个算法的官方 Python 绑定（作者 Aaron Lun）。
**判断依据是 summary / author / project_urls，不是"名字存不存在"。**

---

## 6. 人工复核节点

| 节点 | required | 说明 |
|---|---|---|
| `celltype_labels` | TRUE | 细胞类型注释的最终标签 |
| `doublet_threshold` | TRUE | 双细胞判定的阈值 |
| `trajectory_root` | TRUE | 拟时序根的选择（方向完全依赖它） |
| `communication_pairs` | FALSE | 通讯配对的生物学合理性 |
| `virtual_perturbation_targets` | FALSE | 虚拟扰动靶基因的合理性（§1.7/§1.8） |

**默认 `pending`，不算失败。** 自动化流水线不能替人签字 ——
把未确认的节点记成已确认，等于把复核节点变成摆设。
但必须**可见**：验收里作为"可见但不阻断"的项列出。

---

## 7. 清单必须和结果同时可及

`run_manifest.json` 落在 `results/<dataset_id>/` 下，随 artifact 一起上传
（artifact 包含整个数据集目录）—— **它必须和结果同时可及，
否则追溯链是断的。** 只把清单写在 CI 日志里不行：日志会滚掉，文件不会。

---

## 8. 怎么读一份清单

```python
import json, pathlib
m = json.loads(pathlib.Path("results/pbmc3k/run_manifest.json").read_text())
m["language"]        # "python"
m["seed"]            # 随机种子
m["inputs"]          # 每个输入文件的 sha256 与字节数
m["key_versions"]    # 点名工具逐个的版本（null = 查过了没装）
m["params"]["named_tools"]   # 点名工具的缺口登记
m["decisions"]       # Agent 决策链：问题 / 结论 / 证据
m["human_review"]    # 人工复核节点及状态
m["cross_language"]  # 跨语言转换记录
```

**验收读的是 `manifest_summary(cfg)["inputs_missing_required"]`，不是
`inputs_missing`** —— 后者包含可选输入的缺失，那是正常状态。

---

## 9. 数值为什么会对不上（两个不同的原因）

清单会记下 `versions`，所以**"同一个 commit 在两台机器上给出不同数字"
是可以解释的**，不是玄学。但**"版本不同"只解释了一部分**，实测有两类：

### 9.1 版本差异（本地 vs CI）

| 量 | 本地（Python 3.14 / scipy 1.18） | CI（Python 3.12） |
|---|---|---|
| 四条轨迹的平均 Spearman ρ | `+0.6259` | `+0.6363` |
| Leiden 标签的编号 | 与 CI 不同（簇划分相同、编号顺序不同） | — |

**这是版本差异，不是回归。** 所以比对数值前先比对 `key_versions`。

### 9.2 库内部的并行 kNN（CI 内部，同一版本）

**同一个 Python 3.12、同一批包版本的多轮 CI 也会对不上。** 实测六轮：

| 轮次 | commit | 并行度 | dpt | palantir | **scfates** | 拓扑 | 平均 ρ |
|---|---|---|---|---|---|---|---|
| A | `8f3f56c0` | 无钉 | +0.5856 | +0.4674 | **+0.5644** | 4 片段 / 6 mil | +0.6363 |
| B | `9af13b42`（**只改注释**） | 无钉 | +0.5856 | +0.4674 | **+0.5296** | 6 / 8 | +0.6254 |
| C | `5b241a3` | +BLAS 线程/CORETYPE | +0.5856 | +0.4674 | **+0.5328** | 6 / 8 | +0.6265 |
| D | `db778bb` | +`NUMBA_NUM_THREADS` | +0.5856 | +0.4674 | **+0.5328** | 6 / 8 | +0.6265 |
| E | `f99eab1` | 同上 | +0.5856 | +0.4674 | **+0.5328** | 6 / 8 | +0.6265 |
| F | `e549477` | 同上 | +0.5856 | +0.4674 | **+0.5646** | 6 / 8 | +0.6353 |

**`dpt` 与 `palantir` 六轮逐位相同**，上游的聚类数（10）、分辨率扫描、
CytoTRACE 也完全一致。**只有 scFates 一个在变。**

> **两次归因都被否证了，记在这里比"我修好了"有用。**
>
> 1. **"多线程 BLAS 归约顺序"** —— 钉住 OMP/OPENBLAS/MKL + CORETYPE
>    之后（C 轮）scfates 仍从 +0.5296 变 +0.5328。
> 2. **"`pynndescent` 的 Numba 并行"** —— 补上 `NUMBA_NUM_THREADS=1`
>    之后（D/E/F 轮）三轮给 +0.5328 / +0.5328 / **+0.5646**，
>    最后那个几乎回到未钉时的 +0.5644。
>
> **但钉并行度不是白做**：离散的**拓扑稳住了** —— A 轮是 4 片段/6
> milestone，B 轮之后六轮全部是 6 片段/8 milestone。
> **它让结构稳定，不足以让 ρ 稳定。残留随机源尚未定位。**

**路径**（读源码得到）：

```
scf.pp.diffusion(device="cpu")            # 签名里没有 seed
  -> palantir.run_diffusion_maps(...)     # 靠 palantir 默认 seed=0
       -> compute_kernel(backend="scanpy")
            -> scanpy 默认 method="umap"
       -> diffusion_maps_from_kernel(..., seed=0)
            -> eigs(T, ..., v0=rng.random(...))     # 这步是干净的
```

- **特征求解器没问题**：`v0` 来自 `np.random.default_rng(0)`。
- **扩散图那一步没有 seed 可传**（`scFates.pp.diffusion` 1.2.5 签名里没有）。
- simpleppt 的 PPT 主曲线树（`ppt.py:210-213`）本身是 seed 了的，
  它只是把上游的末位差异放大成可见的 ρ 变化。

**处理**：两个 Python 仓库的 workflow 在 **job 级**钉住三组变量 ——
**它们是三个独立的旋钮，但实测不足以解决这个问题**：

| 旋钮 | 管什么 |
|---|---|
| `OMP_NUM_THREADS` / `OPENBLAS_NUM_THREADS` / `MKL_NUM_THREADS` | BLAS 归约顺序 |
| `OPENBLAS_CORETYPE=Haswell` | OpenBLAS 按宿主 CPU 型号分发 SIMD 内核（线程数管不到） |
| `NUMBA_NUM_THREADS=1` | Numba `prange` 的线程数 |

> **教训**：`set_seed()` 只覆盖"用 numpy/random 的随机调用"。
> 第三方库内部的迭代求解器、并行归约、并行 kNN 都不在里面。
> **看到"我明明设了种子"时，先问一句"这个函数有没有 seed 参数"** ——
> 没有的话，设多少遍都没用。
>
> **更重要的教训：设了环境变量不等于问题解决了。**
> 两次归因都看着很有道理（源码路径也读对了），但 CI 日志两次都否证了。
> **`dpt`/`palantir` 逐位相同这条线索，一开始就该让"BLAS"出局** ——
> 真是归约顺序的话，同一条代码路径上的其他方法也该漂。

### 9.3 报数的时候

`trajectory_status.json` 的 `reproducibility` 字段记的就是这一段：
哪些方法可以按定值报、哪些必须带范围。**别只写在文档里** ——
拿到产物的人不一定读得到这里。

---
