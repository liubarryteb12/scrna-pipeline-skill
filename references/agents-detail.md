# 长文规则原文归档（AGENTS.md 精简版对应的完整原文）

> 本文件是原本写在 `AGENTS.md` 里的长文规则叙事的**未删节归档**：`AGENTS.md` 对每条规则只保留一份精简版，
> 因为 harness 在 **65536 字节**处截断注入的指令文件，长文会让排在后面的规则**整条进不去**。
> **这里没有任何内容被取代** —— 需要完整的事故经过与实测证据时读本文件。
> 持久的事故记录仍然是 `governance/15_ERROR_LEDGER.md`（该仓库在工作区，位于 `..\governance\`）。

## 原规则 15. 静态检查要挡住"未定义名字"

`py_compile` **只做编译，看不出未定义名字**。漏 import 一个 `W_SINGLE`
时它照样报"语法通过"，要等运行时才炸 —— 实测因此白跑一整轮流水线。

`tools/check_py_syntax.mjs` 现在会检查：用到的 `W_SINGLE` / `W_ONE_HALF` /
`W_DOUBLE` / `mm` / `PAL` / `PAL_CYCLE` / `apply_style` 是否都 import 了。

名单是**写死的**，不是"common 导出的所有名字"。后者会把函数参数名当成用法
（`def verify_alignment(adata, log_info=None)` 的 `log_info`），要正确处理
得做作用域分析 —— 那是重写一个 linter。同时会剥掉注释和字符串再扫，
避免"名字只出现在注释里"的误报。

### 15.1 写死的名单只覆盖了它自己 —— 现在是真正的逐作用域分析（E-61）

上面那份名单**只有 10 个名字**，而 `common.py` 实际导出 **74 个**。
于是它挡住的全是"恰好被列进去的"，没列进去的照旧漏到 CI。
姊妹项目空间仓库实测踩到（2026-09-26）：`spatial-pipeline-skill/scripts/03_spatial_domains.py` 漏 import
`spot_radius_plot_units`，本地 `check_py_syntax.mjs` 报"全部通过"，
CI 跑 25 分钟到 H&E 叠图段才 `NameError`，后续步骤全没跑。

**根因不是"名单短了一点"，而是判据的输入域与它要防的缺陷不匹配** ——
漏 import 的可以是任意一个导出。把名单补全只是把同一种漂移往后推一次
（`common.py` 下次加函数又会漂）。

规则 15 当年放弃做作用域分析的理由是"那是重写一个 linter"。
**那个顾虑是对的，但作用域分析不必自己写**：标准库 `symtable` 就是
CPython 编译器的符号表，按作用域给出每个名字是 local / parameter /
imported / free（闭包）/ global，正是这里需要的。

现在 `tools/check_py_names.py` 做三件事：

| 检查 | 判据 |
|---|---|
| 未定义名字 | 某名字在某作用域**被引用**，却既非该作用域局部绑定、也非闭包自由变量、模块层也没有、也不是内置 → 运行时必然 `NameError` |
| 幽灵 import | `from common import X` 而 common 模块级没有 `X` |
| 不安全构造 | 出现 `import *` / `globals()` / `exec` / `eval` / `vars` / `locals` → **判红退出**，而不是假装通过 |

**三个实现上的坑（都实测踩过）：**

1. **不安全构造必须用 AST 判，不能扫子串。** 第一版扫裸子串，于是本文件
   自己的模式元组命中了它自己 —— **检查器把自己的源码判红**。
   "检查器要检查的东西"与"检查器描述自己要检查什么"在文本上无法区分，
   只有语法结构能区分（同规则 25.1：`stripComments` 把字符串换成 `""` 后
   把要检查的东西本身擦掉了）。确实需要时在同一行写
   `# py-names: unsafe-ok —— <理由>` 豁免（**故意做成要写一句话的**）。
2. **行号必须按作用域定位，不能全文件找首次出现。** 第一版把 `plt` 报成
   `common.py:1037` 之外的别处 —— 那是**另一个函数**里的同名变量，
   那个函数自己 import 了它，行号指向一处**没问题的代码**。
   **指向错的行号比不指行号更糟**：下一个人会去读一段正确的代码然后困惑。
3. **"common 提供这个名字"的提示要扣掉 common 自己 import 进来的。**
   `import numpy as np` 被删时报出「common 导出过 np」，而 `np` 出现在
   common 命名空间里只是因为 **common 自己也 import 了 numpy**。
   照着改会写出 `from common import np` —— 把 common 的内部依赖当接口用。
   两个集合分工：含 import 的用于**幽灵 import 判据**，只有 `def`/赋值的
   才用于**提示**。

**两仓的这个文件逐字节相同**（SHA256
`90D2ECEF11170306406C4FA356FD241189A2DFF5DC912ED7C060B8ADC783990E`，12530 字节）——
与 `check_figures.mjs` / `check_legend_convention.mjs` 同样的约定：
改一侧必须同步另一侧并比对哈希。

**标定：** 正向两仓零命中；geo（纯 R）判为**不适用且不判红**；传错目录判红；
反向逐类注入 —— 删 common import、删第三方 import（证明不限于 common 导出）、
幽灵 import、`import *`、`globals()`、`eval()` 各自只让它自己那条响；
逃生舱写了理由后放行。

> **"没有 Python 文件"是"不适用"，不是"失败"。** geo 是纯 R 仓库，跑到
> 这里必须放行 —— 判红会让 pre-push 在 geo 上永远红，而那条告警与 geo 的
> 改动毫无关系。但"仓库根目录传错了"必须与它长得不一样，所以分开检查。

**门禁接线：** CI 在两个 workflow 的「静态检查（不装依赖）」那一步跑；
`governance/hooks/pre-push.mjs` 的 `PY_GATES` 表（Python 门禁用解释器跑，
不是 `process.execPath`；本机没有 `python` 时记为**提醒**而不是通过 ——
"没跑成"和"跑过了没问题"必须长得不一样）。

### 15.2 手工删 import 会连带删掉同一行的活名字

E-61 的直接触发动作是**手工清理死 import**：死名字扫描器正确报出
`plot_marker_dotplot`，但执行删除时把同一 import 行里相邻的
`spot_radius_plot_units`（它**有**调用点，不在死名单里）一起删了。

**扫描器没错，是执行删除这一步错了。** 所以：删 import 时逐名核对
"这个名字在文件里还有引用吗"，不要按行删。死名字扫描器的输出是**名单**，
不是**待删行号**。

---

## 原规则 16. 每轮运行必须留下可追溯的运行清单（模块零）

参考规范：三大部分整合文档的「模块零」（§0.2–§0.4）。

**没有运行清单的分析结果不是结果。** 半年后拿到一份
`cell_communication.csv`，如果不知道当时装的是哪个版本的 liana、
输入的 h5ad 是哪个哈希、随机种子是多少，那份 CSV 就**无法被复现，
也无法被质疑** —— 而不可质疑的结论没有价值。

`common.py` 提供清单层，产物是 `results/<dataset_id>/run_manifest.json`：

| 函数 | 记录什么 | 规范条款 |
|---|---|---|
| `init_manifest` | 开新一轮（**清掉上一轮**） | — |
| `capture_versions` | 全量已装 Python 包 + `KEY_PACKAGES` + `R_KEY_PACKAGES` 逐个 | §0.3 |
| `probe_r_packages` | 起一次 `Rscript` 问 R 这些包装了没有 | §0.3 |
| `record_input` | 输入文件的 sha256 | §0.4 |
| `record_params` | 全部参数**含 seed** | §0.3 |
| `record_decision` | Agent 决策链（问题/结论/证据） | §0.4 |
| `record_human_review` | 人工复核节点及状态 | §0.4 |
| `record_cross_language` | 跨语言转换前后维度与**丢失字段** | §0.2 |
| `manifest_summary` | 供验收用的摘要 | — |

**几条不能省的约定：**

1. **未装的工具要记成 `null`，不能省略键。** `key_versions` 里
   `"liana": null` 和"没有 liana 这个键"是两件事：前者是"查过了，没装"，
   后者是"没查"。省略会让读者分不清。
   **但"查过了"的前提是问对了地方** —— R 包走 `R_KEY_PACKAGES` +
   `probe_r_packages()`，见第 6 条。
2. **人工复核未确认不算失败。** `human_review` 默认就是 `pending`，
   判成 FAIL 会让每个 job 都红，反而没人看。但必须**可见**。
3. **`init_manifest` 必须清掉上一轮。** 和规则 14 同一个道理：
   上一轮的清单留在那里冒充本轮，比没有清单更糟。
4. **清单在步骤跑完之后才登记输入。** 可选步骤这轮有没有产物，
   跑完才知道；在开头登记会把"上轮残留"记成本轮输入。
5. **`pip freeze` 不用 subprocess 抓。** 沙箱下管道捕获会 EPERM，
   用 `importlib.metadata` 枚举。**这条只对 Python 包成立** —— 见第 6 条。
6. **R 包必须问 R，`importlib` 永远问不出来。** `importlib.metadata` 与
   `importlib.util.find_spec` 查的都是 Python 的发行版数据库 / 模块查找器，
   对 R 包**原理上**永远返回"没有"。那个 `None` 会被读成"查过了，装不上"，
   而它真正的含义是"**问错了地方**"。

   实测的后果：K-01b 起 CI 装了 R、`scTenifoldKnk` 真跑了 6 分钟，而
   `run_manifest.json`（**复现依据**）仍记 `scTenifoldKnk: null` ——
   于是它与 `virtual_perturbation_status.json` 的
   `tools.scTenifoldKnk.available = true` 对**同一个工具**给出相反结论。

   所以拆成两条通道，且**失败原因分开写**：

   | 通道 | 实现 | 清单字段 |
   |---|---|---|
   | Python 包 | `importlib.metadata` 枚举（不起子进程） | `key_versions` |
   | R 包 | `probe_r_packages()`：一次 `Rscript -e requireNamespace` | `key_versions` + `r` |

   `r` 那段（`{rscript, r_version, packages, reason}`）是必要的：
   `key_versions` 是**平铺**的 name→version，读不出"哪些是 R 包、
   R 是什么版本、R 到底有没有装"。

   三种失败各有独立 `reason`，**不能混成一句**：包名非法（拒绝把任意
   字符串拼进 `Rscript -e`）/ `Rscript` 不在 PATH（本机没 CI 没装 R）/
   退出码非 0。"没装 R"与"R 装了但包没装"必须区分开，
   否则排查方向会被带偏。

   `08_virtual_perturbation.py` 的 `_probe_r_package()` **转调**
   `probe_r_packages()`，不再自己拼一遍 —— 同一件事只留一份代码
   （抄两份时验证的往往只是副本，见 E-41 与规则 15 的教训）。

`run_manifest.json` 落在 `results/` 下，随 artifact 一起上传 ——
**它必须和结果同时可及**，否则追溯链是断的。

> 本段原为 `AGENTS.md` 规则 16 的正文（2026-09-26 E-69 收口时移入本归档）。


---

## 原规则 17. 虚拟敲除 / 过表达是保留框架，重点是"没做什么"（§1.7 / §1.8）

参考规范：三大部分整合文档 §1.7 / §1.8。规范把这两节标为**保留框架**、
主语言 Python，候选靶基因由 Part 1 产出。落地在 `08_virtual_perturbation.py`。

**规范点名的三个工具，两个装不了、一个已经真跑了（理由都是实测的，不是"没装"）：**

| 工具 | 状态 | 原因 |
|---|---|---|
| `scTenifoldKnk` | ✅ **已接入** | R/CRAN 包（`scripts/lib/tenifold_knk.R` + CI 的 setup-r 层）。**"不在 PyPI"曾经被读成"用不了"** —— 见下 |
| `PerturbNet` | ❌ 装不了 | PyPI 上 0.0.2/0.0.3 钉 `requires_python='<3.8,>=3.7'`，0.0.3b0/b1 钉 `'<3.11,>=3.10'` —— **没有任何一版支持 CI 的 3.12** |
| `RegVelo` | ❌ 装不了 | 能装，但要 RNA velocity 的 spliced/unspliced 层（本数据没有），且拉入 `torch` + `scvi-tools` |

> **E-41 的教训（探针查不到 ≠ 目标不存在）。** 早先 `probe_tools()` 对三个工具
> 统一用 `importlib.util.find_spec(name)`，而 `find_spec` **只查 Python import
> 路径，对 R 包在原理上不可能返回 True**。于是"我没装 R"被写成"这个工具不存在"，
> 还被一条验收项固化成"**必须**声明不是 scTenifoldKnk"。**理由那句"不在 PyPI"
> 是真的** —— 真理由 + 越界结论 = 一条读起来完全合理、实则自我封闭的记录。
> 现在 `scTenifoldKnk` 走 `_probe_r_package()`（`shutil.which("Rscript")` +
> `requireNamespace`），并**区分**"Rscript 不在 PATH（CI 没装 R）"与
> "R 装了但包没装"。**每加一个 `find_spec`/`which`/`exists` 式判断，都要能
> 回答"目标存在时它返回什么"**；答不上来就说明这条判断测的是探针自己。

**"装不了"和"没装"是两件事。** 前者有确切原因、要写进状态文件；
后者是疏忽。所以 `tools` 里逐个记 `{available, reason}`，
验收项的判据是**每个不可用的都有原因**，而不是"全都不可用" ——
将来某个工具能装了，这条应该自动变 PASS 而不是 FAIL。

### 两个引擎并列，一张图只画一个

`perturbation.engine` 取 `both`（默认）/ `first_order` / `tenifold`。
**拼错不静默退回默认** —— 那会让"我配了 tenifold"变成一句假话，所以直接报错。

两套结果**并列落盘**（`virtual_perturbation.csv` 与
`virtual_perturbation_tenifold.csv`），一致性只比排名、不比数值
（一阶是 z 分数的 L2 范数，tenifold 是流形欧氏距离，**量纲不同**）。
一致性高**不说明哪个对**：两法共享同一个上游（`07_grn` 的共表达边），
可能只是共享了同一个偏差。

**出图默认画 tenifold（真方法），它没跑出来时退回一阶**，标题写明是哪个。
两法不画在同一根 x 轴上 —— 那会让人以为可以直接比大小。

**`scTenifoldKnk` 没有细胞类型分辨率**：整份数据只建一个网络，输出是
「扰动基因 × 网络基因」的全局距离。**本仓库没有**把距离按细胞类型重新加权 ——
那会造出一个既不是 scTenifoldKnk、也不是本仓库一阶近似的新方法，
然后借它的名字发出去。

### 空敲除：网络表达不了该扰动 ≠ 效应为 0

`scTenifoldKnk` 的 `transcriptomeWide` 模式敲除一个基因的方式是把 WT 网络的
**那一行清零**（**scTenifoldKnk 包内源码** `R/scTenifoldKnk.R` L206–L207
`KO <- WT; KO[g, ] <- 0`）。
若该基因**出度本来就是 0**，清一行全 0 的行等于什么都没敲，`KO` 与 `WT`
逐位相同，`manifoldAlignment` 对两个完全一样的网络对齐，返回的"距离"
只剩浮点噪声（实测 ~1e-16）。

**它不报错、不给 NA、不给零** —— 给出一排看起来完全正常的数，
读表的人只会得出"敲除这个基因没有影响"。所有既有图门禁都会通过，
因为它们查的是"有没有产出"不是"产出对不对"。

判据是**结构事实（出度）**，不是"距离是不是很小"—— 后者会把一个真实
但微弱的扰动也误判掉。R 侧 `zero_outdegree_targets()` 算出度并整行置 NA，
meta 里带 `empty_knockout_genes` + `target_outdegree`；Python 侧
`compress_tenifold_distances()` 用**那份结构证据**打标记
（而不是从"整行都是 NA"反推 —— 那会把 R 侧别的 NA 原因也误标成空敲除），
空敲除行报 `None` 而**不是 `0.0`**，并在 `compare_engines()` 里
**显式剔除并记账**，而不是靠 `dropna()` 悄悄少几个点。

**两侧各有自检**：`Rscript scripts/lib/tenifold_knk.R --selftest` 验出度计算，
`python tools/selftest_tenifold.py`（**不需要 R**）验 Python 侧有没有把
那份证据正确消费成 `None` / 标记列，并且**端到端跑一遍
`run_tenifold_engine`（用假的 Rscript）** —— 抽出函数时最容易漏改的
`return` 里的旧变量名，`py_compile` 查不出来，只会在 CI 上 tenifold
真跑通那一刻抛 `NameError`，而那时 R 已经跑了 6.5 分钟。

### 一阶近似：能报什么，不能报什么

用 `07_grn.py` 的调控子边表做**一阶、单跳**线性传播：

    Δz_G = −z̄(G,C)        （敲除）   /   +overexpress_sd（过表达）
    Δz_Y = w(G,Y) · Δz_G   对 G 的每个直接靶基因

**`ko_magnitude` 里混了两个因素，必须分开报。** 它是
`|Δz_G| × network_sensitivity`，所以排在前面的基因往往只是
**在该细胞类型里表达高**，而不是在网络里被连得紧 —— 实测同一份结果，
`GATA2` 的 magnitude 是 5.41（排第二）而网络敏感度只有 0.81，
`TBX21` 的敏感度 1.84（全场最高）却只排第四。所以额外落一列
`network_sensitivity = √(Σ_Y w²)`（与表达无关）。

> **名次会变，值不会。** 早先写的是"`GATA2` 排第一、`TBX21` 排第三"；
> 当前产物里 magnitude 最高的是 `TYMS`（6.05，敏感度 0.94），
> GATA2 退到第二、TBX21 退到第四 —— 而 **5.41 / 0.81 / 1.84 三个值逐位未变**。
> 引用前看 `virtual_perturbation_status.json` 的 `top_by_effect`，不要从文档抄。

**不做多跳传播。** 在相关网络上做多跳会放大噪声，看起来像"网络效应"，
其实只是把相关系数乘了几遍。

**`PerturbNet` 始终没跑，`method` 必须仍然否掉它。** 状态里的 `method`
显式写出"**不是 PerturbNet**"，验收项检查这个否定词在不在。
（原来那句"也不是 scTenifoldKnk"已随 K-01b 删除 —— 那条判据在真跑通
R 引擎之后会**把正确的结果判红**。现在的判据是"声明的方法与状态里的
引擎一致"：记了 `engines_used` 就必须逐个在方法串里被点名，
没跑 tenifold 时必须写明**为什么没跑** —— 否则"没跑"和"跑了没结果"
长得一样。）

### 跨部分交接只走 CSV（§0.2）

`perturbation.targets_csv` 指向 Part 1 交接的基因表（一列 `gene`，可选一列
`logFC`），相对路径按 `data_dir` 解析。**没有它时回退到内部调控子，
并且必须在状态里写明回退了** —— 此时 `signature_alignment` 为空是预期的，
不是 bug。回退不报错、只是把"Part 1 的疾病签名"换成了"本步自己的调控子"，
下游解读完全变了，所以 `target_source` 是必填字段。

### 本仓库还有第二条交接：Part 2 → Part 3 的参考 h5ad

上面那条是 **Part 1 → Part 2**（进）。还有一条**出**的：
Part 3（空间转录组）的 `deconvolution.reference: h5ad` 要的就是本仓库
`03_cluster_annotate.py` 写出的 **`clustered.h5ad`** ——
带细胞类型标签的单细胞参考。

**契约有三条，前两条不满足会让 Part 3 静默算错或 KeyError：**

| # | 要求 | 不满足时 Part 3 的表现 |
|---|---|---|
| 1 | `layers['counts']` 是**原始计数** | 退回用 `.X`（log 值）求参考谱 → **解出的比例没有意义，而它看起来仍像一组比例** |
| 2 | `obs` 里有细胞类型列（本仓库是 `celltype`）| 配 `celltype_key` 时 KeyError |
| 3 | 基因集是**交集** | 静默接受 —— 参考若是 HVG 子集，参考谱覆盖的基因就变少 |

**第 3 条最容易被忽略，而且它是本仓库的实际情况**：`clustered.h5ad` 来自
`integrated.h5ad`，那是**2000 个高变基因的子集**，不是全基因集。
Part 3 因此只在共同基因上建参考谱。

**交接的两端都要能被核对。** 所以产出侧把契约**主动记进**
`cluster_status.json` 的 `part3_reference`（含 `contract` 与 `limitations`），
验收项 `§0.2 Part 3 参考导出的契约` 据此检查 ——
只检查消费侧（Part 3 记没记）的话，产出侧悄悄丢掉 counts 层
**不会被任何人发现**。

### 候选基因的完整边表要落盘

`07_grn.py` 的 `tf_regulons.csv` 里 `top_targets` 只有前 12 个靶基因、
**且没有权重**，虚拟扰动需要全部靶基因及其相关系数。
所以 07 额外写 `tf_regulon_edges.csv`（`tf, target, corr`）。
**靶基因名与权重必须成对取** —— 分开写（一个列表取名字、另一个按位置取
相关系数）在过滤 `corr > 0` 之后会错位，而错位不报错，
只会给每个靶基因配上一个别人的相关系数。

---

## 原规则 20. `set_seed()` 管不到"库内部的迭代求解器"

姊妹项目 `geo-normal-pipeline-skill` 的规则 12 记了同一类问题（那边是极小的
P 值把浮点末位放大到可见）。单细胞这边实测的表现是**拓扑变了**。

**证据（同一 Python 3.12、同一批包版本的五轮 CI）：**

| 轮次 | commit | 并行度设置 | dpt | palantir | **scfates** | scFates 拓扑 | 平均 ρ |
|---|---|---|---|---|---|---|---|
| A | `8f3f56c0` | 无钉 | +0.5856 | +0.4674 | **+0.5644** | 4 片段 / 6 milestone | +0.6363 |
| B | `9af13b42`（**只改注释**） | 无钉 | +0.5856 | +0.4674 | **+0.5296** | 6 片段 / 8 milestone | +0.6254 |
| C | `5b241a3` | +OMP/OB/MKL/CORETYPE | +0.5856 | +0.4674 | **+0.5328** | 6 片段 / 8 milestone | +0.6265 |
| D | `db778bb` | +`NUMBA_NUM_THREADS=1` | +0.5856 | +0.4674 | **+0.5328** | 6 片段 / 8 milestone | +0.6265 |
| E | `f99eab1` | 同上 | +0.5856 | +0.4674 | **+0.5328** | 6 片段 / 8 milestone | +0.6265 |
| F | `e549477` | 同上 | +0.5856 | +0.4674 | **+0.5646** | 6 片段 / 8 milestone | +0.6353 |

`dpt` 与 `palantir` **六轮逐位相同**，上游的聚类数（10）、分辨率扫描、
CytoTRACE 也完全一致。**只有 scFates 一个在变。**

### 20.1 两次归因都被否证了

**第一次：多线程 BLAS 归约顺序。** 把 geo 规则 12 那两组变量搬了过来 ——
**C 轮否证**：钉住之后 scfates 从 +0.5296 变成 +0.5328，没有收敛。

**第二次：`pynndescent` 的 Numba 并行。** 读源码发现扩散图那一步经过
scanpy 默认的 `method="umap"`，于是补了 `NUMBA_NUM_THREADS=1` ——
**D/E/F 三轮否证**：0.5328 / 0.5328 / **0.5646**。
最后那个值几乎回到 A 轮未钉时的 +0.5644。

**两次否证留下的线索：**

- **拓扑稳住了。** A 轮是 4 片段/6 milestone，B 轮之后**六轮全部是
  6 片段/8 milestone** —— 钉并行度确实让**离散的**结构稳定了。
- **但连续量（ρ）还在漂**，而且看起来有两个吸引子：`≈+0.5328` 与
  `≈+0.5645`。**残留随机源尚未定位。**

**这一条比"我修好了"更有用：钉并行度能让拓扑稳定，但不足以让 ρ 稳定。**
第三次动手之前先读日志 —— 本仓库规则 12 就是为这个写的。

### 20.2 真正的路径（读源码得到）

`05_trajectory.py` 的 `compute_scfates()` 走的是：

```
scf.pp.diffusion(device="cpu")
  -> palantir.utils.run_diffusion_maps(data_df, n_components=10, knn=30, alpha=0)
       -> compute_kernel(..., backend="scanpy")            # palantir 默认后端
            -> scanpy.neighbors.Neighbors(temp)
                 .compute_neighbors(n_neighbors=30, n_pcs=0, method=None)
                      -> scanpy 默认 method="umap"
       -> diffusion_maps_from_kernel(kernel, n_components, seed=0)
            -> eigs(T, ..., v0=rng.random(...))            # v0 来自 seeded RNG，没问题
```

**两处关键事实：**

1. **扩散图那一步没有 seed 可传。** `scFates.pp.diffusion`（1.2.5
   `preprocessing/diffusion.py`）**签名里根本没有 `seed`**，内部调
   `run_diffusion_maps` 时也不传 —— 靠 palantir 的默认 `seed=0`。
   `05_trajectory.py` 把 seed 传给了 `sct.tree()` 和 `sct.pseudotime()`
   （这两处是对的），**但扩散图那一步没地方传**。
2. **特征求解器是干净的**（`eigs(..., v0=rng.random(...))`，`v0` 来自
   `np.random.default_rng(0)`）。
3. simpleppt 的 PPT 主曲线树（`ppt.py:210-213`，`np.random.seed(seed)`
   之后 `np.random.choice` 取初始节点）本身是 seed 了的 ——
   它只是把上游的末位差异放大成可见的 ρ 变化。

### 20.3 三条要记住的

1. **看到"我明明设了种子"时，先问"这个函数有没有 seed 参数"。**
   没有的话设多少遍都没用。`set_seed()` 只覆盖用 numpy/random 的随机调用；
   第三方库内部的迭代求解器、并行归约、GPU 内核都不在里面。
   排查顺序：`grep -n seed <包的源码>` —— 看 `seed` 是**出现在签名里**
   还是**只出现在 docstring 里**。
2. **并行度有三个独立的旋钮，缺一个都不够**（但它们**不保证**解决）：

   | 旋钮 | 管什么 |
   |---|---|
   | `OMP_NUM_THREADS` / `OPENBLAS_NUM_THREADS` / `MKL_NUM_THREADS` | BLAS 归约顺序 |
   | `OPENBLAS_CORETYPE` | OpenBLAS 按宿主 CPU 型号分发 SIMD 内核（线程数管不到） |
   | **`NUMBA_NUM_THREADS`** | **Numba `prange` 的线程数** |

   ```yaml
   env:
     OMP_NUM_THREADS: 1
     OPENBLAS_NUM_THREADS: 1
     MKL_NUM_THREADS: 1
     OPENBLAS_CORETYPE: Haswell
     NUMBA_NUM_THREADS: 1
   ```

   **这五个变量现在都设着，但 scFates 的 ρ 仍然在漂（D/E/F 三轮）。**
   设它们是对的（拓扑因此稳住了），**但不要以为设了就可复现。**
3. **"版本不同"不是万能借口。** 版本相同也能对不上。别一看到数字变了
   就归因到版本上 —— 那会掩盖真正的回归。先比对 `key_versions`，
   相同就去看**哪些量没变**。

### 20.4 结论怎么报

**scFates 的 ρ 不要当单一确定值报。** 实测范围 **+0.5296 ~ +0.5646**，
平均 ρ 因此是 **+0.6254 ~ +0.6363**。`05_trajectory.py` 会把这段写进
`trajectory_status.json` 的 `reproducibility` 字段 ——
**结论的适用范围要跟着产物走，不能只写在 AGENTS 里。**

`dpt` / `palantir` / `cytotrace` 六轮逐位相同，可以按确定值报。

---

## 原规则 25. 四条"图没了 / 图被裁了却全绿"的补强（Q-26，2026-09-25）

起因是台账 E-48 与 E-49 两条**互相独立**的缺陷，它们暴露了验收层与门禁层
**四个各自独立的盲区**。四条补强分别堵一个，缺一条那类缺陷就还能再犯一次。

### 25.1 嵌套 `status` 的 `failed` 必须判红（验收层）

`status` 不只在顶层。`qc_status.json` / `integration_status.json` /
`cluster_status.json` / `grn_status.json` 里都有**嵌套**的 `status` 字段，
而旧验收层只看顶层 `d.get("status")`。

实测（真 artifact `scrna-results-55`）：顶层分布 `{ok: 7, not_configured: 1}`，
**嵌套分布 `{ok: 5, failed: 1, not_applied: 1, not_done: 2}`** ——
唯一那条 `failed` 是 `grn_status.json → regulon_vs_pseudotime.status = failed`
（`ValueError: Invalid unit 6.800000000000001 in 'figsize'`，E-48 的后果），
而它的顶层 `status` 是 `ok`。**顶层全绿、里面已经崩了。**

判据必须**只把 `failed` / `error` / `fail` 判红**：

| 嵌套值 | 含义 | 处置 |
|---|---|---|
| `failed` / `error` / `fail` | 崩了 | **判红（required）** |
| `not_applied` / `not_done` | 设计如此地没做 | 可见不阻断 |
| `ok` | 正常 | 放行 |

误判的代价是双向的：把 `not_applied`（单样本不做整合）判红会让每个 job 都红，
把 `failed` 放行则正是 E-48。实现是 `_iter_nested_status(obj, path=())`
递归产出 `(路径, 值, 同级 reason)`，扫**全部** `*status.json`
（不只可选步骤那 5 个 —— 上面四个文件里三个是必需步骤的）。

### 25.2 图验收必须条件化，且要覆盖全部出图脚本（验收层）

`REQUIRED_FIGURES` 原来只有 14 条，**不含 `02-07-*`（三张 GRN 图）与
`02-08-01`** —— E-48 让 5 张图从未产出，而验收层根本不知道它们该存在。

补齐之后还有第二个问题：**图不能无条件要求产出**。`trajectory.enabled: false`
时 `02-05-*` 一张都不该有，旧写法一关轨迹就必红。所以改成由图名第 2 段
（`02-07-01` 的 `07`）反查所属步骤，从 `read_state(cfg)` 取该步骤的
status 与 required：

- 所属步骤**可选且未跑成** → `required: False`，detail 写明
  `步骤 grn = not_configured，本轮不要求产出`
- 否则 → `required: True`，缺失时 detail 带 `**缺失**（步骤 grn = ok）`

**兜底是 25.1 那条独立扫描** —— 步骤崩了的时候，即使图项因条件化放行，
嵌套 `failed` 仍会把它判红。

### 25.3 溢出检测必须落盘，不能只打 WARN（产物级）

E-49 的形态：`_content_overflow()` **正确检测到了**
`width_overflow_frac=0.3523`、**正确打了 WARN**，然后**没有任何人读** ——
标题超宽 35%，在 `savefig.bbox: standard` 下被静默裁掉，图照样生成、
门禁照样绿。

**"检测到了"不等于"有人会知道"。** 现在 `save_fig()` 里的
`_record_figure_overflow(cfg, name, bad)` 把溢出**落盘**到
`results/<dataset_id>/figure_overflow.json`（`{figures: {名: bad}, n_overflow}`），
验收层读它并**判红（required）**，detail 写明"**会被静默裁掉**"。

三条配套约束：

1. **与状态文件同理，跑前必须删旧文件**（`_ovf_stale.unlink()`）——
   否则图修好了旧记录还在，验收把已修好的图判红。
2. **只累积、不覆盖** —— 同一次运行多张图溢出要全部记下。
3. **记录可被修复** —— 标题改短后不应再新增记录（已做双向验证）。

### 25.4 门禁要查"内容贴边"，不能只查"有没有墨"（门禁层）

`check_figures.mjs` 原来只查"有没有墨 / 是不是糊死"。
**被裁掉的图照样有墨** —— 这正是它当年漏掉 E-49 的原因。

裁切的物理后果是**墨迹延伸到画布边缘**，所以量非背景像素外接框到左右边的
距离。阈值 `EDGE_MIN_PX = 3`，实测标定（scrna 34 张 + spatial 76 张）：

- scrna 边距分布 `{0:1, 10:9, 11:2, 12:6, 13:13, 14:2, 41:1}` ——
  `L<3 或 R<3`、`L==0 或 R==0`、`L<6 或 R<6` **都只命中 `02-08-01` 一张**
  （就是 E-49 那张），零误伤。
- spatial 76 张全部通过，零误伤。

两条刻意的取舍：

1. **只判左右，不判上下。** 实测 `02-05-05-unit1-pseudotime-by-cluster`
   上边距 = 0 —— 那是布局取舍，不是裁切。单看某一边会误伤。
2. **暗底图跳过。** 整幅都是"墨"，外接框必满幅，判了全是假阳性。

> 本文件与 `spatial-pipeline-skill/tools/check_figures.mjs` **必须逐字相同**
> （两仓同一份，见该文件头注释）；改一侧必须同步另一侧并比对哈希。

### 25.5 上线首跑就抓到一条新缺陷 —— 修好一个缺陷会暴露被它掩盖的下一个

四条补强推送后第一次 CI（run `36140048953`，commit `e2d7a2c`）**验收判红**：

```
[FAIL] 图 02-07-01-unit5-tf-smarca4-trend 内容超出画布  **会被静默裁掉**：{'width_overflow_frac': 0.0371}
验收失败：1 项必需检查未通过
```

`scripts/07_grn.py` 的单 TF 面板标题是**单行**，而面板宽度只有
`W_SINGLE = 89 mm`：

```python
ax_t.set_title(f"{tf} activity along pseudotime (mean ± 1 SD per bin)")   # 旧
```

本地用仓库自己的 `apply_style` + `_content_overflow` 逐名标定，得到
`FOSL2 +0.0093 / ATF4 −0.0019 / E2F2 −0.0021 / SMARCA4 +0.0371`
（**与 CI 逐位吻合**）；假想的长 TF 名 `SMARCAD1 +0.0485`。即**在阈值附近
徘徊、随 TF 名长度而变** —— 换一份数据必然再踩。

**这批面板此前从未渲染过**：Q-24 的三元素 `figsize` 让整段在 `plt.subplots`
就抛异常，4 张面板一张都没画出来。所以这是它们的**首次渲染**，
缺陷与图同时诞生 —— **旧缺陷一直在掩盖新缺陷**。

修法是**折行**而不是删字：

```python
ax_t.set_title(f"{tf} activity along pseudotime\n(mean ± 1 SD per bin)")
```

第二行宽度**与 TF 名无关**，本地对全部 8 个名字（含 3 个假想长名）标定
**一律 −0.0238**（固定 2.4% 余量）。对照候选：去掉尾巴只剩 −0.0143
（余量太薄）、去掉 `along` 在 `SMARCAD1` 上仍 `+0.0228` 判红。

**两条规则：**

1. **宽度受限的小面板，标题折行，不要删字。** 删字只是把阈值往下挪一点，
   输入一长又踩到；折行让第二行宽度与输入无关，是**结构性**修法。
2. **修好一个缺陷后要重跑一遍全部检查。** "图重新出现"这个动作本身必须
   触发一次完整复检 —— 不能只看"图现在有了"，因为它第一次出现时可能就带着
   一个此前没人有机会看见的毛病（本条就是）。

台账：`governance/15_ERROR_LEDGER.md` E-51（E-49 的同族，都是
`savefig.bbox: standard` 下静默裁切）。

---

## 原规则 26. native 崩溃绕过 `except`：`sc.pp.scrublet` 偶发 SIGSEGV（E-52，2026-09-25）

CI 的「跑流水线」步骤偶发 **exit 139**。日志停在 `01_qc.py:226` 的
`过滤: …` 之后一行 `Segmentation fault (core dumped)`，**没有任何
Python traceback**。

**它是长期存在的，不是哪一轮引入的。** 全量 60 个 run 里有 8 个
`run_attempt>1`，**全部**是"att1 失败 + rerun 成功"形态，最早一个是
`06977d0`（2026-09-20T11:34:50Z）。**这 8 个的 att1 崩溃位置逐字相同**，
跨 5 天、跨 7 个 commit。**"rerun 能过"把它掩盖了整整 5 天。**

### 26.1 崩溃点怎么锁死的

`01_qc.py:226` 的 `过滤:` 是崩溃前最后一行，下一句是 L229
`db = run_scrublet(adata, cfg)`，里面 L79 `sc.pp.scrublet(...)`。
`scripts/lib/common.py` 第 62 行的 `print(..., flush=True)` **逐行 flush**，
所以"日志最后一行"就是"代码走到哪一行"，不存在缓冲区丢日志。

**但 `01_qc.py:85` 的 `except Exception` 接不住它** —— native SIGSEGV
直接杀进程，`qc_status.json` 里连一行失败原因都没有。

### 26.2 根因：OpenBLAS 0.3.34 的 `dgemm_kernel_HASWELL` 栈越界

**这一步是靠 `PYTHONFAULTHANDLER: 1` 拿到的真栈定下来的**（run
`36150495910`）—— 在此之前只能靠猜。栈的落点（逐字，自下往上读）：

```
scipy/sparse/linalg/_interface.py:1118  in _matmat      ← return self.A @ X
scipy/sparse/linalg/_interface.py:451   in _shared_matmat
scipy/sparse/linalg/_interface.py:491   in matmat
scipy/sparse/linalg/_eigen/_svds.py:511 in svds         ← Av = X_matmat(eigvec)
scipy/_lib/_util.py:306                 in wrapper
sklearn/decomposition/_pca.py:740       in _fit_truncated
sklearn/decomposition/_pca.py:542       in _fit
sklearn/decomposition/_pca.py:440       in fit
scanpy/preprocessing/_scrublet/pipeline.py:84  in pca
scanpy/preprocessing/_scrublet/__init__.py:448 in _scrublet_call_doublets
scanpy/preprocessing/_scrublet/__init__.py:243 in _run_scrublet
scanpy/preprocessing/_scrublet/__init__.py:298 in scrublet
scripts/01_qc.py:79                     in run_scrublet
```

**栈里一帧 numba 都没有。** 崩溃是一次**稠密 GEMM**：`_scrublet/pipeline.py`
第 79 行先把稀疏矩阵 `.toarray()`，第 84 行交给
`PCA(svd_solver="arpack")`，sklearn 走 `svds` → ARPACK 求特征向量后做
`X_matmat(eigvec)`，而 `self.A` 已是**稠密 ndarray**，于是落到
`dgemm_kernel_HASWELL`。

上游 **OpenBLAS #6026**（2026-09-11 报、09-13 关）的形态与本仓**逐项吻合**：

- 环境原文即 `OpenBLAS 0.3.34, as bundled in the scipy-openblas64 wheel
  shipped with numpy 2.5.2 and 2.5.3` + `OPENBLAS_CORETYPE=Haswell` +
  GitHub 托管 `ubuntu-latest` —— 正是本仓 job env 凑齐的组合。
- **机制**：Haswell 内核把 k 维分块打包进一个固定 `0x7080` 字节的**栈
  缓冲区**，每步写 96 字节且**对 k 无上界**；k ≳ 319 时越界覆盖
  callee-saved 寄存器与返回地址槽，最后那条 `ret` 直接 SIGSEGV。
- **为什么偶发、为什么 rerun 能过**：level-3 的 blocking（决定 k 的大小）
  **在运行期按宿主 cache 拓扑选取**。GitHub 的 Azure fleet 混着不同型号
  CPU —— 落到 k 会变大的机器上就崩，落到 k ≤ 256 的机器上就过。
  **同一 commit、同一批包版本、相隔 8 分钟的两个 run 一成一败**，只能
  这样解释。
- 上游实测**受影响的机器上约三次崩一次** —— 与本仓观测到的 rerun 频率一致。

**附带风险：越界若只覆盖了保存寄存器、没碰到返回地址，进程会正常返回而
数值被污染。** 即"没崩"不等于"算对了"。

**各 numpy wheel 自带的 OpenBLAS**（读 `numpy.libs/libscipy_openblas64_*.so`
里的版本串得到）：

| numpy | 自带 OpenBLAS | 状态 |
|---|---|---|
| 2.4.6 | 0.3.31.188.0 | 正常 |
| 2.5.0 / 2.5.1 | 0.3.33.112.0 | **最后一个正常版** |
| 2.5.2 | 0.3.34.0.0 | 回归 |
| **2.5.3**（本仓此前浮动到的） | **0.3.34.106.0** | 回归 |

### 26.3 处置

**修复是钉 numpy 上界**（`requirements.txt`）：
`numpy>=2.1,<2.5.2` —— 下界 2.1 是 `anndata 0.13` 的硬要求，上界躲开
OpenBLAS 0.3.34。**上游修复（0.3.35）尚未进入任何 numpy wheel**
（PyPI 上最新仍是 2.5.3），所以只能钉，不能等。

**上线复验（run `36153515064` / commit `f19eff9`，2026-09-25T15:20Z）**：
**attempt 1 一次过** —— 这是本条 60 个 run 里**第一次** att1 就成功
（此前 8 个 rerun 全是 att1 失败）。日志 segfault 命中 **0**；CI 实际装的
是 `numpy-2.5.1-cp312-cp312-manylinux_2_27_x86_64…whl`，artifact
`run_manifest.json` 的 `versions.numpy = 2.5.1` —— **两处对上，说明钉上界
真的生效，而不是"这次恰好没崩"**。artifact `scrna-results-61` 亲验：
`acceptance.json` = **75 项全过 / 0 失败**、`figure_overflow.json`
不存在、`figures/` 39 PNG + 39 PDF（`02-07-01-unit1..5` 全在）。

> **"一次过"本身不是证明。** #6026 是**宿主相关**的偶发缺陷 —— 单轮成功
> 只说明"这一轮的宿主 k ≤ 256"。真正的证据是**版本已经落到正常侧**
> （2.5.1 → OpenBLAS 0.3.33.112.0），所以这一轮的成败**不再是随机的**。
> 要再确认，观测口径是**后续若干轮的 rerun 率应降到 0**，而不是
> "再看一轮是不是绿的"。

**`PYTHONFAULTHANDLER: 1` 继续留着**，理由变了：它这次兑现了价值（没有它
根因定不下来），而 numpy 上界只覆盖 OpenBLAS 这一个来源。

**`NUMBA_THREADING_LAYER: workqueue` 保留，但它的立项理由已被推翻。**
它是按"numba 并行运行期"的假设加的，而栈证明 numba 不在路径上 ——
带上它照崩（run `36150495910`）。保留只剩一条防御性理由：numba 的 `omp`
线程层会加载一份 libgomp，而 scikit-learn 自己 bundle 了一份
（`sklearn.utils._openmp_helpers`，实测在崩溃轮的 218 个 extension
modules 里），**同进程两份 OpenMP 运行时**是上游有记录的崩溃形态。

### 26.4 五条规则

1. **native 崩溃不能靠 `except` 兜底。** `except Exception` 只覆盖 Python
   层异常；凡是有 C 扩展参与的关键步骤，都要问一句"它要是崩了，我会不会
   连日志都没有"。
2. **"rerun 能过"不等于"没有问题"。** 重试成功是掩盖，不是修复。本条被
   掩盖了 5 天、8 次 rerun —— **应当把重试率本身当成一项可观测指标**。
3. **给 native 崩溃留一条取证路径，比急着改代码值钱。** 从"日志停在
   `过滤:`"到"OpenBLAS 栈越界"之间隔了一整轮错误的归因；`faulthandler`
   一条 env 就把它终结了。**改代码之前先确认自己有没有能力观测。**
4. **依赖的二进制缺陷也是本仓的缺陷。** 这次崩的不是本仓任何一行代码，
   而是 numpy wheel 里打包的 OpenBLAS。**依赖"没钉上界"就是把自己的
   稳定性交给上游的发布节奏** —— `requirements.txt` 开头的告诫因此
   不只是"API 会变"，也包括"二进制会坏"。
5. **同一份二进制在不同宿主机上行为可以不同。** "本地复现不了"不等于
   "没问题"：k 的大小取决于宿主 cache 拓扑，这是**设计如此**的行为。
   遇到宿主相关的间歇失败，先找"哪个量在运行期按机器选"。

台账：`governance/15_ERROR_LEDGER.md` E-52；任务行 `governance/02_TASKLIST.md` V-02。

---

## 原规则 27. 步骤"没做成"必须传出来：返回值要接住、顶层 `status` 也要扫（E-56，2026-09-26）

### 27.1 现象：一条已存在但从未生效的判据

`main_analysis.py` 的嵌套失败扫描本意是"任何 `*status.json` 里写了失败就判红"，
但写成 `if p and str(v).lower() in NESTED_FAILED_VALUES`。`_iter_nested_status(d)`
在**顶层**调用时 `path` 是**空元组** —— falsy —— 于是**顶层 `status` 一个都扫不到**。

而本仓有 6 处代码正是以"写顶层非 ok 状态 + 正常 `return`"的方式失败：

| 位置 | 取值 |
|---|---|
| `scripts/05_trajectory.py:275-279` | `bad_root` |
| `scripts/05_trajectory.py:313-317` | `missing_clusters` |
| `scripts/05_trajectory.py:363-370` | `insufficient_methods` |
| `scripts/07_grn.py:352-355` | `no_tfs_in_data` |
| `scripts/08_virtual_perturbation.py:901-911` | `no_candidates` |
| `scripts/04_pseudobulk_de.py:178-180` | `no_celltype_testable` |

**结果：这一步什么也没产出，`state.json` 记 `ok`、`acceptance.json` 记绿。**
而且 `insufficient_methods` / `bad_root` / `missing_clusters` / `no_candidates`
这四个词**连 `NESTED_FAILED_VALUES = ("failed", "error", "fail")` 都不在** ——
即使 `if p` 修好，顶层这四个也照样漏。

### 27.2 第二处同病根：返回值被丢弃

`run_all` 的步骤循环把 `fn(cfg)` 的返回值**直接丢弃**，写死
`record_step(cfg, sid, "ok", ...)`。步骤函数的"跑完了、但结果是『没做成』"
路径是 `write_json(状态文件, status)` + `log_warn(...)` + `return status`，
**不抛异常** → `except` 分支根本不进。

同一病根还在 **9 个步骤脚本各自的 `__main__` 块**里（独立运行时同样丢弃
返回值、无条件记 `ok`）—— 即"编排器跑"和"单独跑"两条路径**都**把失败记成成功。

**这是 E-48 的另一半。** E-48 当时的修法是"让嵌套失败可见"，但只覆盖了
"写进嵌套字段的失败"：写进**顶层**的失败与"返回值里的失败"都没被覆盖。
E-48 的教训是"异常被降级成一个没人读的字段"，本条是它的两个变体：
**字段写在没人扫的位置**，和**结果压根没被记下来**。

### 27.3 处置

`scripts/lib/common.py`：

- 两个**分开的**取值集合：`STEP_ABORT_VALUES`（判红）与 `STEP_SKIP_VALUES`
  （可见不阻断，只作文档用途）。前者含通用 `failed`/`error`/`fail` 加上表
  六处的具体取值，以及 07 的 `no_regulons`、06 的 `no_pairs_in_data`/`no_signal`、
  04 的 `no_usable_pseudobulk`/`column_missing`。
- `result_status_of(res) -> str`：dict 且带 `status` 键则取之，否则 `"ok"`。
- `classify_step_result(v) -> str`：`"ok"` / `"abort"` / `"skip"`。
  **未登记的取值归 `"skip"`** —— 只把确知是失败的判红（假阳性比假阴性更危险）。
- `record_step(..., result_status: str = None)`：`entry["result_status"]` 与
  `entry["status"]` **分开存**。后者是编排器记的"有没有崩"，前者是步骤自述的
  "做成了没有" —— `status="ok"` + `result_status="bad_root"` 正是要显形的组合。
  默认 `None` → 不写这个键，向后兼容。

`scripts/main_analysis.py`：步骤循环 `res = fn(cfg)` → `record_step(...,
result_status=result_status_of(res))`；`except` 分支加 `result_status="failed"`；
步骤验收循环新增"步骤 X 自述状态"检查；嵌套扫描去掉 `if p`，改
`key = ".".join(p) + ".status" if p else "status"`（照抄姊妹仓 spatial 的写法，
`spatial-pipeline-skill/scripts/main_analysis.py:380-398` 本来就是对的）。
9 个步骤脚本的 `__main__` 块同步改造。

### 27.4 四条规则

1. **"让失败可见"的修法要覆盖失败的所有承载位置。** 修完要问："这个失败还可能
   被记在哪里？"—— 返回值、顶层字段、嵌套字段、日志、退出码。E-48 只修了
   嵌套字段一处，另外两处照旧漏了。
2. **`if p` 在"空路径 = 顶层"的语义下必然滤掉顶层。** 空元组/空字符串/空列表
   都是 falsy，而写的时候看起来完全合理。凡"空值有语义"的地方（顶层路径、
   根节点、默认分支），判据必须显式写 `if p else` 而不是 `if p`。
3. **顶层与嵌套是两个语义域，取值域重叠不等于可以共用判据。** `tenifold.status`
   可以是 `timeout`/`no_candidates` **而内置引擎照样出了结果（顶层 `ok`）**；
   `annotation.celltypist` 的 `model_unavailable`、`liana` 的 `api_not_found`
   同理。合并两个集合会让**每个 job 都红**，而常年假红会训练人忽略告警。
4. **标定必须用本次 CI 的 artifact，不能用工作区里可能陈旧的副本。**
   `figure_review/scrna-pbmc3k` 是旧 run 的产物，内含 E-48 遗留的
   `grn_status.json → regulon_vs_pseudotime.status = "failed"` ——
   用它做正向标定会把真判据误判成误报。

台账：`governance/15_ERROR_LEDGER.md` E-56；任务行 `governance/02_TASKLIST.md` R-03。
标定脚本：`D:\tmp\_s1\calib_e56.py`（正向 + 反向 19 类注入 + 回退版对照）。

---

## 原规则 28. 九条"看起来在算、其实没在算"的缺陷（E-58，2026-09-26）

R-03 第二批复核的是审计严重项 S2–S10。九条里 **7 条读码确认成立、1 条（S8）
被实测否证、1 条（S3）根因比审计写的更具体**。共同点是：**代码在跑、产物齐全、
状态全绿，但那个量取不到它声称要取的信息。**

### 28.1 "恒为 0 / 恒为真"的量比没有这个量更糟

`scripts/03_cluster_annotate.py` 的 `compare_annotations()` 原来做
`str(own[c]).lower() == str(maj[c]).lower()`。marker 词表来自
`assets/celltype_markers.yml`（`T_cell`/`Platelet`），CellTypist 词表来自模型
（`Tcm/Naive helper T cells`/`Megakaryocytes/platelets`）—— **两套词表没有任何
一个字符串相等**，于是 `n_agree_exact: 0 / agreement_frac: 0.0` 恒成立。而逐簇看
`B_cell`↔`B cells`、`Platelet`↔`Megakaryocytes/platelets`、`Monocyte`↔
`Classical monocytes` 明显一致。

**一个恒为 0 的量看起来像一个结论（"两条路完全不一致"），实际只是词表不相交。**
凡两套词表、两套坐标系、两个口径要对齐的地方，先跑一遍**看它到底能不能取到
非平凡值**（同 Q-27 的 `figures:dynamic` 按前缀计数恒真）。

处置：新增 `assets/celltype_mapping.yml` + `load_celltype_mapping()`；没有映射的
簇**不进分母**、单独列 `unmapped`（猜一个映射等于把"我没定义"伪装成"不一致"）；
字段改名带 `mapped`（`n_mapped`/`n_unmapped`/`n_agree_mapped`/
`agreement_frac_mapped`），**旧字段名一并删除**，避免读者混用两个分母。消费侧
`scripts/main_analysis.py` 分三层判：没对比 → 红；有对比但 `n_mapped == 0` → 红
（映射表缺条目或词表已变，正是旧缺陷的形态）；有映射 → 报真的一致率。

### 28.2 注释写着正确做法、下一行做了相反的事（三处）

| 位置 | 注释说 | 代码做 |
|---|---|---|
| `scripts/01_qc.py:78` | "用细胞数反推期望双细胞率" | `clip(5000/n*0.01, 0.05, 0.10)` 在 n>1000 时**恒被下界截断到 0.05**，那个行为从不发生 |
| `scripts/00_fetch.py:147` | "判断整个 `X` 是不是计数" | `X[:min(200, X.shape[0]), :]` 只抽前 200 行 |
| `scripts/04_pseudobulk_de.py` | "必须用 counts layer，不能用 `X`" | `layers["counts"] if "counts" in layers else adata.X` 静默退回 log 值 |

三处的共同后果都是**产物里完全看不见**：scrublet 阈值偏保守、`predicted_doublet`
偏少（pbmc3k n=2652 实测新式 0.02122 vs 旧式 0.05000，**高估一倍以上**）；计数校验
在数据按样本拼接时可能给出 `is_counts: True` 而整体不是计数 —— 而这是下游全部
方法学的前提；喂 log 值给 DESeq2 不报错，只给错的离散度估计。

写注释时问一句：**"这行代码真的会走到我说的那条路吗？"** 尤其 `clip` 的上下界、
`if/else` 的两个分支、`[:200]` 这类抽样切片。

处置：S9 改 10x 经验式 `clip(0.008*n/1000.0, 0.01, 0.10)` 并把**旧公式的值一并
记进 status**（`expected_doublet_rate_old_rule`）—— 让读者看见差距，而不是相信
一个注释；S10 改为**抽两批且抽法必须不同**（等距抽样作主判据 + 中段连续抽样作
交叉核对），两批占比差 >1 个百分点即判 `sampling_consistent: False`；S6 没有
counts 层时返回 `counts_source: "missing"` 而非退回，调用侧判成独立状态
**`missing_counts`**（上游契约被破坏，与"这批数据不适合做拟bulk"必须长得不一样）。

### 28.3 检验家庭不能被任何"预过滤"缩小

`scripts/06_communication.py:321` 原来 `if obs_score <= 0: continue` —— 零分组合
根本不进 `rows`，于是 `multipletests` 的分母是"打分 > 0 的组合"而非"全部被评估的
组合"。实测 pbmc3k 上分母 649 vs 1134，**检验家庭缩小约 43%**，`p_adj_bh` 系统性
偏小、`n_significant_bh` 上偏 —— 方向恰好与紧邻注释所担心的"校正组合数多导致的
假阳性"**相反**。

处置：零分组合仍进 `rows`，`p_value` 记 **1.0**（零分的观测值本就无法被任何置换
超越，1.0 是它在检验家庭里的正确取值），另加 `tested` 布尔列区分"做过置换"与
"恒为 0"。凡多重检验，**分母必须是"被评估过的全部假设"**。

### 28.4 方向/符号不能靠数据行序决定

`scripts/04_pseudobulk_de.py:169` 原来写 `str(sub["group"].unique()[0])` 当分子 ——
取的是"**第一次出现的取值**"，而 `sub` 的顺序来自 `meta.groupby(...)`，即 obs
行序。**重排细胞顺序会让全部 `log2FoldChange` 变号，而 `padj` 一个都不变**（对比
方向翻转是符号对称的），产物看起来完全正常。

同类：`05_trajectory.py:439-444` 的共识拟时序把 `cytotrace` 也算进去了，而
`direction_source == "cytotrace_fallback"` 时 **`cytotrace` 正是方向参考**（`:417`）
—— 同一文件在 `:414-425` 算方法间一致性时**正确地把参考排除了**，算共识时却没有，
**参考方法既定了方向、又参与共识，等于自己给自己投票**。产物证据：
`trajectory_direction.csv` 四个方法的原始 rho **全为负、全部被翻转**，即"方向修正"
在做全部的工作。

处置：S5 分组水平按 `sorted()` **字典序**定死（可复现、与行序无关），分子/分母拆成
显式变量，实际对比方向写进 `status["contrast_used"]` + `contrast_rule`；S7 共识改用
`consensus_names = cv_names if cv_names else names`（与一致性统计同口径），并**同时
算一个含参考方法的共识**供对比、把两者 rho 写进 `consensus_vs_with_reference_rho`
—— "差多少"本身就是信息。

### 28.5 审计报告也会错：落修法前先证伪

审计称 `scripts/07_grn.py:164` 的 `corr[top[:len(targets)]].mean()` 与 `targets`
错位（`top` 是位置数组、`targets` 是名字列表）。**实测证明不会**：
`top = np.argsort(corr)[::-1][:N_TARGETS]` 是**按 corr 降序**排的，所以
`corr[k] > 0` 这个过滤必然保留 `top` 的一个**前缀** —— `top[:len(targets)]` 与
`[k for k in top if corr[k] > 0]` **恒等**。两万次随机对拍（含 `-inf` 自排除、
含并列值、含三种尺度）**零次不等**（`0/19944`）。

处置：仍改成 `np.mean([w for _, w in pairs])`（等价、更钝 —— 不依赖"降序 → 前缀"
这个推理，以后若有人把排序改成升序或改成按 p 值选，这里不会跟着错），并在源码
注释里**写明审计前提不成立**，防止后人照审计报告改回去。

**"读码推断"不等于"跑一遍对拍"**，尤其当论断依赖"某个数组的排序性质"时。

### 28.6 marker 基因宇宙错位（S2）

`scripts/03_cluster_annotate.py:71` 原来用 `adata.var_names` 过滤签名基因，而
**同一函数下一段用 `use_raw=True` 打分** —— `adata` 此时只剩 2000 HVG，
CD3D/CD8A/CD14 这类**在数据里真实存在、只是没进 HVG** 的 marker 被判成"缺失"，
签名被削到无法区分（实测簇 1/簇 3 的 `T_cell` 与 `CD4_T` 分数逐位相同、margin
恰为 0.0）。**同一个文件 `:385` 的 dotplot 已经用的是 `adata.raw.var_names`。**

处置：改用 `raw_names = set(adata.raw.var_names) if adata.raw is not None else
hvg_names`，并把诊断拆成两个字段：`missing_markers`（数据里真的没有）vs
`not_in_hvg`（有、只是没进 HVG）—— **合并成一个字段会让"签名被 HVG 削弱"和
"这批数据没测到"看起来一样**。

### 28.7 四条规则

1. **"恒为 0 / 恒为真"的量比没有这个量更糟** —— 它看起来像结论。新写一个
   比例、一致率、覆盖率时，先跑一遍确认它**能取到非平凡值**。
2. **注释陈述的行为必须与代码实际行为对得上** —— 逐条检查 `clip` 的边界、
   `if/else` 的两个分支、抽样切片的范围。
3. **多重检验的分母是"被评估过的全部假设"** —— 任何 `continue` 掉一批样本
   再算 FDR 的写法，等价于偷偷改分母。
4. **方向/符号来自显式排序，不来自数据行序** —— `unique()[0]`、`head(1)`、
   `groupby` 首元素都不行；且要把实际用到的方向写进产物。

台账：`governance/15_ERROR_LEDGER.md` E-58；任务行 `governance/02_TASKLIST.md` R-03。
标定脚本：`D:\tmp\_s2\calib_s2_s10.py`（每条都配"回退版必须抓不到"的对照，34 项全过）。

---

## 原规则 29. 门禁要带内建自检，自检必须被反向标定，且必须接进 CI 与 pre-push（E-62 / E-63，2026-09-26）

姊妹项目 `spatial-pipeline-skill/AGENTS.md` 规则 30 是同一批经验的另一半。

**没有自检的门禁只能证明"它没报错"，不能证明"它检查了"。** 两条实测：

| 台账 | 门禁 | 缺陷形态 |
|---|---|---|
| E-62 | `tools/check_legend_convention.mjs` | Python 侧没抹注释 → 注释里一个 `fig.legend(` 让括号配平**一路吞到文件尾**，其后所有真调用一个都没查，门禁照样打绿 |
| E-63 | `tools/check_figures.mjs` | WARN 落盘分支从落地起**一次都没执行过**，里面有两个必崩的错（`INK_FAIL_MIN` 未定义、报告路径用了循环变量 `dir`）|

两条的共同点：**假阴性**。门禁的失败方式不是"报错"，而是"什么都不报" ——
而"什么都没发现"与"检查通过了"在输出上完全一样。**假阴性比假阳性危险得多**：
假阳性会被人骂着修掉，假阴性会被当成绿。

### 29.1 自检要调真代码，不能自己重写一遍逻辑

E-63 的自检第一版有 9 个用例，**用例 8 自己另写了一遍路径拼接**
（`join(dirname(join(dir, "figures")), "warn_report.json")`），没调真代码。
反向标定把缺陷注回去（`outPath` 改回 `join(dir, ...)`）→ **自检仍然通过（exit=0）**。

**抽函数**才解决：`checkDirs()` / `buildWarnReport()` / `writeWarnReport()`
三个纯函数（不 print 不 exit），`main()` 与自检**都调它们**。抽完再注一次缺陷 →
`ReferenceError: dir is not defined`、exit=1。

> 同 E-62 的「检查器要检查的东西，与检查器描述自己要检查什么，在纯文本上
> 无法区分」是同一个坑的两种形态：**自检里重实现一遍被测逻辑，等于没测。**

`check_figures.mjs` 的 `--selftest` 自带零依赖 PNG 编码器（`crc32` +
`encodePng` + `deflateSync`）合成用例，9 个用例里 3 条是 E-63 回归，
输出 `自检通过（9 个用例，含 3 条 E-63 回归）`。

### 29.2 反向标定：逐个把原缺陷注回去，确认自检真的会红

**正向通过证明不了任何事** —— 一个永远返回 True 的用例在干净产物上也是绿的。

| 注入 | 期望 | 实测 |
|---|---|---|
| `inkFailBlank: MIN_INK` → `INK_FAIL_MIN` | 回归 1/2 红 | ✅ `ReferenceError` |
| `outPath` 改回 `join(dir, ...)` | 回归 2 红 | ✅ 第一版**不红**（假自检）→ 抽函数后红 |
| 尾斜杠处理删掉 | 回归 3 红 | ✅ |

**尾斜杠那条是修完才发现的第三个缺陷**：`dirname("a/b/figures/")` 给出
`a/b/figures`（尾斜杠把最后一段当成文件名），报告落进 `figures/` 里与图混在一起。
修法 `String(dirs[0]).replace(/[\\/]+$/, "")` 后再 `dirname`。

### 29.3 自检必须接进 CI 与 pre-push —— 没人跑的自检是同一类缺陷

写了 `--selftest` 却只在本地手敲，等于又造了一个"从未执行过的分支"。
现在三处都接：

| 位置 | 内容 |
|---|---|
| `.github/workflows/scrna_analysis.yml`（本仓）/ `spatial-pipeline-skill/.github/workflows/spatial_analysis.yml` | 「静态检查（不装依赖）」那一步跑 `check_legend_convention.mjs --selftest` + `check_figures.mjs --selftest` |
| `geo-normal-pipeline-skill/.github/workflows/geo_analysis.yml` | 同一步跑 `check_legend_convention.mjs --selftest`（geo 无 `check_figures.mjs`）|
| `governance/hooks/pre-push.mjs` | 新增 `SELFTESTS` 表，在**所有**静态门禁之后跑，日志标签是 `（自检）` |

**日志标签必须区分"带镜像目录"与"带 `--selftest`"** —— 两者都走 `extraArgs`，
但一个是拿真实产物判、一个是拿合成用例判，长得一样就没法排查。

### 29.4 判据的"通过数"必须能看见 0 —— 否则死代码与"没有这类输入"无法区分（E-64）

`check_doc_refs.mjs` 补判据 C（跨仓引用）时，主体写在
`if (!isA && !isB) continue` **之后** —— 而跨仓 token 正是"两条判据都不进"
的那一类，于是**判据 C 的分支永远走不到**。更糟的是报告那行是
`if (nCheckedC) console.log(...)` 守卫的：计数恒 0 时**连打印都不打印**，
输出看起来与本仓没有跨仓引用**完全一样**。三仓复跑全绿、毫无异常迹象。

**两条规则：**

1. **新判据要放在早退分支之前** —— 分流顺序错了分支就是死代码，
   而**死代码不报错、只是永远不执行**。
2. **计数行不能加 `if (n)` 守卫** —— **0 也是信息**。"检查了 0 条"与
   "根本没检查"必须在输出上长得不一样。

### 29.5 门禁的"检查范围"要和"它守护的动作"对齐（E-65）

`governance/hooks/pre-push.mjs` 的 `changesOf(repo)` 原来只读
`git status --porcelain`（**只含未提交改动**）。而 pre-push 是 `git push`
的钩子，**它唯一被调用的时刻就是"已经提交、还没推送"** —— 那时 porcelain
为空，三仓全走 `无改动，跳过`，**[4/6] 静态门禁段一条都没跑**，
而打印的是 `PRE-PUSH 通过（0 条提醒）。可以 push。`

**这不是边角，是主路径**：正常情况下它每次都在空转，给出虚假的安心。
修法是取并集 —— 除未提交改动外，再加
`git log --name-only --pretty=format: @{upstream}..HEAD`（**已提交未推送**）。

> **一个门禁段被整段跳过时，不能打印"通过"。** `无改动，跳过` 用的是
> `ok()`（绿勾），它和"查过了没问题"在输出上一样。写门禁前先跑一次
> "什么都没改"的路径 —— 如果它空转时也说通过，那它有改动时说的通过
> 也不可信。**正确的提交顺序是：改完 → 跑 pre-push → 提交 → push。**

**`SELFTESTS` 的接线也做了反向标定**：把 `check_figures.mjs` 自检里
"全白判红"用例的条件改成 `false` → pre-push 输出
`✗ spatial-pipeline-skill tools/check_figures.mjs 未通过`、
`PRE-PUSH 未通过：1 项判红 —— 禁止 push。` —— **说明接线真的会拦，
而不是只在日志里多打一行 `✓`。**

> **接了线但从不失败的检查，与没接线是一样的。** 每加一条自检，都要问
> "我怎样让它红一次"——答不上来就说明它现在是个装饰。

**两仓 `tools/check_figures.mjs` 与 `tools/check_legend_convention.mjs`
各自必须逐字节相同**（`check_figures.mjs` 只在本仓与 spatial 仓之间，
`check_legend_convention.mjs` 三仓同一份）。改一侧必须同步并比对 SHA256。

---

## 原规则 30. 「写出来了」不等于「有人读」：三种形态与三个守卫（E-69，2026-09-26）

> 这是 `AGENTS.md` 规则 30 的**未删节原文**。`AGENTS.md` 只保留一份精简版（因为 harness 在 65536 字节处截断注入的指令文件）。

**规则句：** 一个字段 / 一条判据的价值不在于它被算出来，而在于**有人消费它**；而"算不出来"和"算出来很小"必须在产物里**长得不一样**。

E-68 是同族第一次自查（规则 29 那批门禁的连带产物），E-69 是拿它的判据**回头扫全仓**抓到的第二次 —— **第三次跨仓抓到东西**（spatial 侧 15 个只写不读的字段、一个恒真的人工复核判据、以及验收层看不见的产生端；scrna 侧同形一处）。

---

### 30.1 Form A：`nan` 参与比较会静默变成 `False`

`nan > x` / `nan < x` / `nan == x` 全是 `False`，而且**不报错、不打日志、不抛异常**。后果是"**算不出来**"和"**算出来很小**"在产物里长得一模一样 —— 状态字段会一本正经地报告一个**由 nan 推出的结论**。

三处现场（spatial 侧，`spatial-pipeline-skill/AGENTS.md` 规则 31 的同一次收口）：

| 现场 | 缺陷 | 假结论 |
|---|---|---|
| `spatial-pipeline-skill/scripts/05_deconvolution.py` 重建误差 | `np.nanmean` / `np.nanpercentile` / `np.nanmax` 在全 nan 时返回 nan 并继续算；`frac_unreliable` 的分母用 `len(errors)` 而不是有限值个数 | "重建误差的中位数是 nan" 被当成一个可比较的量 |
| `spatial-pipeline-skill/scripts/07_spatial_communication.py` z 分数 | `null_sd = float(np.nanstd(perm_means)) or 1e-9` —— 全 nan 时 `nan` 是 truthy，`or` **不触发** → `z = nan` 写进 top5；`denom` 全零时 `100.0/1e-9 = 1e11` 的**假放大** | "这个配体-受体对显著富集" |
| `spatial-pipeline-skill/scripts/08_spatial_trajectory.py` Moran's I | `improved = I_spatial > I_expr`，两侧都是 nan 时 `nan > nan` 为 `False` | 状态声称"**平滑损害了一致性**"，而两个量**一个都没算出来** |

**修法：抽纯函数。** 抽出来才能被标定脚本**直接调**（139 项断言，秒级）；内联在流水线里只能靠跑整条流水线验证（25 分钟一轮），而"跑一轮要 25 分钟"会让人**不去验证**。

```python
summarize_morans_pair(i_expr, i_spatial, ndigits=4) -> {"defined","improved","gain","expression_only","spatially_smoothed","note"}
spatial_z_score(near_mean, perm_means) -> (null_mu, null_sd, z, z_reason)
summarize_reconstruction_error(errors, max_err) -> {"error_defined","median","n_unreliable","frac_unreliable",...}
common.finite_round(x, n)   # 三态：nan/inf/None -> None；有限 -> 四舍五入；0.0 保留
```

`finite_round` 的两个坑：**`0.0` 是合法值、必须保留**（判空要用 `is not None` 而不是 `if not x`）；**`np.float64` 是 `float` 子类、`np.float32` 不是** —— float32 会走到 `json.dump` 的 `default` 回调。

**本仓同族现场（已修，2026-09-26 收口）：** `scripts/05_trajectory.py` 的
`off = [float(cmat.loc[a_, b_]) for i, a_ in enumerate(cv_names) for b_ in cv_names[i + 1:]]` 在 `cv_names` 只有 1 个方法时是**空列表**，于是 `mean_rho` / `min_rho` 记 `nan` —— `mean_rho < 0.3` 那条限制**静默不触发**（限制被跳过而不是被报告）、`round(mean_rho, 4)` 落盘成 `null`（**写盘前是 `nan`**）、日志与限制文案打出 `+nan`。

修法与 spatial 侧同形，也抽成纯函数：

```python
method_correlation_stats(cv_names, cmat) -> (mean_rho, min_rho, state, note)
# state ∈ {"ok", "single_method", "undefined"}
```

**为什么是三个状态而不是"有值 / nan"两个。** `single_method`（分母里只有一个方法，**没有"方法间一致性"这个量**）与 `undefined`（有两个以上方法，但离对角相关**全是非有限值** —— 各方法退化成常数列）是**两种互不相干的处境，排查方向不同**：前者要去看方向参考是不是把方法都剔掉了，后者要去看那些方法为什么退化成常数列。压成一个 `nan` 就等于把排查线索丢掉 —— 这与规则 5「'不确定'与'算不出来'要分开报」、规则 4「失败原因必须是原因码不能是布尔量」是同一条。

配套改动四处，**缺一处这类缺陷就能再犯一次**：

| 位置 | 改动 | 为什么 |
|---|---|---|
| 产出端 | `method_correlation_mean_offdiag` / `_min_offdiag` 走 `common.finite_round` 而不是裸 `round` | `round(nan, 4)` 是 `nan`；靠 `_scrub_nonfinite` 兜底只是**写盘那一刻**变成 `null`，内存里的消费者读到的仍是 `nan` |
| 产出端 | 新增 `method_correlation_state` + `method_correlation_undefined_note`（**正常态与 `no_consensus` 两处都写**） | `null` 本身不说明**为什么**是 `null` |
| 判据 | `if mean_rho is not None and mean_rho < 0.3:` | 加 `is not None` 守卫；**不能用真值判断**（`0.0` 是有效值，而 `0.0 < 0.3` 必须触发那条限制） |
| 消费端 | `main_analysis.py` 的「轨迹方法间一致性已量化」判据由 `_mc_ok` 参与，detail 里带 state 与原因 | **只把状态打进 detail 是装饰** —— 判据本身必须由它参与，否则"算不出来"照样绿 |

**消费端这一处值得单独说。** 旧判据是 `bool(cv) and mean_rho is not None` —— 在"一致性算不出来"时，产出端写的是裸 `NaN`，`json.load` 读回来是 `nan`，于是 `nan is not None` 为真、**这条判据照样是绿的**。也就是说：**这条判据从落地起就对它本该抓的那类处境失效**，而它在干净产物上永远是绿的（正向标定看不出来）。

---

### 30.2 Form B：只写不读的状态字段

字段写进 JSON、而验收层**一个消费者都没有**时，**坏值与"字段不存在"长得一模一样**。全仓扫出 **15 个**这样的字段（`used_counts_layer` / `matrix_source` / `counts_layer` / `full_gene_counts_available` / `hvg_fallback` / `counts_check` / `n_spots_dropped_no_coords` / `dropped_no_coords_examples` / `coord_coverage_note` / `max_genes_cap` / `gene_selection` / `n_genes_dropped` / `n_proportion_values_truncated` / `celltype_spot_alignment` / `spatial_smoothing_improves_coherence`），补了 **12 条注册表探针**：

```python
dict(cid, base="data"|"results", f=<文件>, path=<元组>, kind, severity, opt, fn=<lambda>, good=<文本|lambda>, bad=<文本|lambda>)
```

- `fn` 返回 `None` ⇒ **不适用**（PASS）
- `opt=False` = **无条件写**，缺失 ⇒ 判红「**产生端不再写了**」
- `opt=True` = **条件写**，缺失 ⇒ PASS

12 条：`status:input_is_counts` / `coords_dropped` / `hvg_flavor` / `full_gene_counts` / `svg_gene_subset` / `svg_gene_selection` / `svg_genes_dropped` / `proportions_truncated` / `deconv_matrix_source` / `niche_spot_alignment` / `smoothing_improves` / `morans_I_defined`。

**三个守卫：**

1. `chk(cid, kind, ok, detail, severity="required")` —— **第二个位置参数是 `kind` 不是 `severity`**。第一版把 severity 值塞进 kind 槽，`bad_root` / `not_applicable` / `missing_pca` 全被记成 `required`；两个同名同型的参数相邻，**传错不报错**。
2. 父状态白名单 `_PROBE_PARENT_OK = (None, "ok")` + `_PARENT_NOT_EXECUTED` 映射 12 个非 ok 字面量；**未知状态跳过但打印原始值**（不要静默放行）。
3. `_dig_present(obj, path) -> (found, value)` —— `_dig()` 对"**键不存在**"与"**值为 None**"返回同一个 `None`，而 `hvg_fallback=None` 的意思是**没有发生回退**（好事）。

本仓同族：`manifest_summary` 的 `n_human_review` / `human_review_confirmed`。

---

### 30.3 Form C：恒真判据

`manifest:human_review` 的 `ok` 曾**写死 `True`**、`required: False` —— "**一个都没登记**"被写成"**全部已确认**"。两仓同形。

修法：`_n_hr > 0` 才可能为真；`human_review_confirmed` 只列 `status in ("confirmed","overridden","not_needed")` 的节点；默认 `pending` **不算失败但必须可见**（与规则 16 第 2 条同一条理由 —— 判成 FAIL 会让每个 job 都红，反而没人看）。

**恒真判据比"没有这个判据"更糟**：它占着"已检查"的位置，让人不再去找真正该检查的东西。

---

### 30.4 三个实现坑（第三个坑有两个方向）

1. **消费端标定看不见产生端。** 反向标定里把产生端的 `morans_I_defined` 键名改掉，消费端 12 条探针**全绿** —— 因为消费端读的是标定用的 JSON 夹具，**产生端写什么它根本不知道**。补的判据必须是**源码级**的：`spatial-pipeline-skill/scripts/08_spatial_trajectory.py:404` 必须是 `"morans_I_defined": _mi_pair_defined`。
2. **`_code_only` 不能查字典键名。** `_code_only` 剥 COMMENT+STRING，而**键名本身就是 STRING token**，一起被剥掉 ⇒ 永远找不到。要用只剥 COMMENT 的 `_no_comment`。**同一个文件里两种剥离策略各服务一条判据**（E-68 第三处、E-69 结构组各踩一次）。
3. **tokenize 的 token 是无空格拼接的。** 判据写成 `'"morans_I_defined": _mi_pair_defined'`（带空格）**永远匹配不上**；实测 `with space: False` / `without space: True`，上下文是 `..."morans_I_defined":_mi_pair_defined,...`。

---

### 30.5 收口

- `common.write_json` 加 `allow_nan=False` + `_scrub_nonfinite`（递归 dict / list / tuple / numpy 类型；**dict 键必须能当 `str`**，否则 `json.dump` 的 `default` 回调崩）。此前它会写出**裸 `NaN`** —— 那是**非法 JSON 字面量**，Python 的 `json.load` 能读回去，其它语言的解析器不能。
- `lib/alignment.py` 是唯一绕过 `common.write_json` 的写盘点，已收口。
- `common.read_json` 改成**损坏时抛异常**（原来 `except Exception: return None` 把"文件缺失"与"文件损坏"混成一件事），另加 `read_json_or_none` 给"确实可能没有"的调用点。

---

### 30.6 标定与防复发（本仓这一处：正向 37 项 + 反向 6 类）

**spatial 侧**（`D:\tmp\_q28\`）正向：`calib_e69_probes.py` **139 项 / 0 失败**、`calib_e69.py` 100 项、`calib_e69_consumer.py` 36 项。
反向：`neg_e69_probes.py` **11 类注入 → 符合预期 10 类 / 不符合 0 类**、`exit=0`、末尾 `源码已还原: True`。

**本仓侧**（`D:\tmp\_e69\`）正向：`calib_scrna_traj.py` **37 项 / 0 失败**（三种 state 各取到、`0.0` 不被当成"算不出来"、部分 nan 报出剔除对数、旧实现在同一输入下给出 `nan` 且 `nan < 0.3` 静默为假、旧消费端 `mean_rho is not None` 对 `nan` 判绿）。
反向：`neg_scrna_traj.py` **6 类注入 → 符合预期 6 类 / 不符合 0 类**、`exit=0`、末尾 `源码已还原: True`（注入项：调用点退回裸 `nan` / 去掉 `is not None` 守卫 / 产出端退回裸 `round` / 去掉 `method_correlation_state` / `finite_round` 不再挡 nan / 消费端不再用状态判红）。

**反向标定第一轮有 3 类注回去却全绿** —— 说明那 3 条判据当时**压根不存在**（产出端用没用 `finite_round`、状态字段还在不在、消费端有没有拿状态参与判红，一条都没查）。补上判据后才是 6/6。**这正是反向标定的全部价值：正向 37 项全绿时，你并不知道哪几条判据是空的。**

防复发：

1. **正向标定（干净产物）看不出 Form A/B/C 任何一形** —— 必须反向标定：逐个把原缺陷注回去，确认判据真的会红。
2. **判红必须伴随非空 FAIL 摘要。** 子进程少 `cwd` / `PYTHONIOENCODING` 时按 cp936 读 UTF-8 源码抛 `UnicodeDecodeError`，**崩溃的非零退出会被误读成"判红"**。
3. 断言必须 **None-safe**（`"x" in why` 在 `why is None` 时抛 `TypeError`）。
4. 断言必须**复刻调用点的构造方式**（`np.concatenate` 传空数组不抛、传空列表抛）。
5. 标定要覆盖**整条回退链**，不是只标定其中一环。
6. `nan` 不能参与比较；`or` 不能兜底 nan（`nan` 是 truthy）。
7. 凡"空值有语义"处（顶层路径、根节点、默认分支）判据要显式写，不要靠 truthiness。
8. 新字段要**同时**加产生端和消费端，否则就是下一个 Form B。
9. 检查器查"有没有"之外，还要查"**取到的值是不是平凡值**"（恒 0 / 恒真）。
10. **抽纯函数**是让标定可执行的前提 —— 抽不出来就只能跑整条流水线，而"要跑 25 分钟"会让人不去验证。

台账：`governance/15_ERROR_LEDGER.md` E-69。
