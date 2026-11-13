# 工程规范 — 分析代码规范（Linux 运行面）

> version 1 | 2026-10-13
> 本规范规定 `src/` 下三条线分析代码怎么写。新步骤动工前先读本文件。
> 与 [isolation.md](isolation.md)（三线隔离）、[directories.md](directories.md)（文件放哪）配套；
> 隔离红线优先于本规范任何条目。

## 第一条铁律：代码在 Linux 上跑

**所有分析代码（`src/<line_id>/scripts/` 下的步骤脚本）运行环境是算力节点的
Linux 容器，不是 Windows，也不是本地 PC。**

因此一切按 Linux 运行面规范写，Windows 习惯一律不带进来：

- **shebang + 严格模式**：`#!/usr/bin/env bash` + `set -euo pipefail`，不许缺。
- **换行必须是 LF**，不许出现 CRLF：CRLF 会让 awk 的数值比较退化成字符串比较，
  统计结果全错（坑 #103）。本地编辑后先确认行尾再部署。
- **不写任何 Windows 专有语法**：无 PowerShell/cmd 语义、无 `\` 路径分隔符、
  无盘符、无 `C:/` 路径。路径一律 `/` 分隔的 POSIX 绝对路径。
- **路径与文件名纯英文**：Linux 系统内不出现中文字符（directories.md 命名红线，
  含"Linux 系统内一律英文路径与文件名"条）。中文只允许进注释与报告正文。
- **可执行位**：脚本部署后必须 755（`tools/ops/06_deploy_scripts.ps1` 已带 chmod + `test -x` 复核）。
- **不用 em dash（—）**：注释与报告用逗号、冒号、句号、括号。

ops/ 下的 PowerShell 脚本跑在本地 PC 侧，只做链路运维与部署传输入口，
**不是分析代码**，不受本条约束；但它们下发的远端命令整体用单引号包裹，
避免 `$var` 被 PowerShell 本地展开（坑 #98）。

## 步骤脚本结构（生产模式）

每个步骤一个自包含可执行文件，bash 包装 + python3 heredoc：

```
scripts/<step_name>
├── bash 段：seed 契约检查、参数解析、输入分流
└── python3 - <<'PYSCRIPT' ... PYSCRIPT   # 计算与自校验， heredoc 用引号防变量展开
```

- **seed 契约**：`BIOANALYS_SEED` 必须由引擎注入，缺失即 exit 3，不设默认值。
- **参数解析**：只认 `--input`（可多个）/ `--output` / `--seed`，其余参数 exit 2；
  未知参数立即报错退出，不静默忽略。
- **输入分流**：按 `case` 模式匹配时模式**前后都带 `*`**：引擎传进来的是绝对路径，
  前缀匹配永远不中（坑 #99）。
- **计算段**：环境隔离用 `H5AD_FILE` 等环境变量传路径，不拼命令行。
- **自校验内建**：输入契约（obs/var/uns/layers/obsm 逐项核对）、核心数值口径
  （如库大小归一化后精确等于 target_sum）、产物结构（shape/稀疏性/dtype/对账），
  任一不过 `fail()` 非零退出。**自校验不过不许落盘。**
- **报告**：每步产 `*_report.txt`，口径、阈值、诊断量、计时、产物清单全部如实记录，
  包括"如实报告"类数字（如 mito 入选 HVG 数），不粉饰。

## 壳层书写纪律（写母版与运维脚本时同样适用）

三条来自 2026-10-25 实测的约定：同一难度的活，违反时三次失败，遵守后一次通过。
完整提示词见 `workspace/handoff/archive/2026-10-25_shell_discipline_prompt.md`，可整段贴给新 agent。

1. **多行内容落盘成文件再执行。** 超二十行的代码写进 `.cjs` / `.ps1` / `.md` / `.json`，
   再用 `node 文件` / `pwsh -File 文件` 跑。命令行里出现 `${...}`、`$_`、`$env:`、
   反引号、嵌套花括号任意一种，即改为落盘：PowerShell 双引号先做变量展开，
   这些字符会被提前吃掉，到目标程序时已是残串。
2. **here-string 包裹字面量。** 需要原文义不变地交一段内容给程序时用 `@'` 与 `'@`，
   两者各自单独占一行，闭定界符顶格。单引号定界不展开变量，双引号定界 `@"` 仍会展开。
3. **短字符串默认单引号。** 只有确定要展开变量时才用双引号。

shell 母版写完用 `bash -n <脚本>` 查语法（bash 位于 `D:\Git\bin\bash.exe`，
该目录需在 PATH 中）。母版契约由 `BioanalysDesktop/tests/master_conventions.test.js` 兜底，
本地没有 bash 时靠它间接验证。

## 退出码约定

| 码 | 含义 |
|---|---|
| 0 | 成功，产物与台账齐 |
| 2 | 参数错误 / 输入契约不满足（缺输入、文件不存在、无法识别） |
| 3 | `BIOANALYS_SEED` 未注入（引擎契约破坏） |
| 1 | 自校验失败（计算或产物不合口径） |

包装器专用退出码：`4` = 同步骤锁互斥，`5` = 启动前内存余量不足，任务未执行。

## Python 侧约定

- 栈：anndata / scanpy / numpy / scipy（线3 空间分析按需 squidpy，装前先 dry-run，
  确认不连带降级既有依赖，坑 #107）。**不新增依赖**，除非已写入该线 `envs/manifest.txt`
  并在报告注明版本。
- `ad.settings.allow_write_nullable_strings = False` 每个 python 段开头显式关掉
  （pandas 3 与 2 默认值不同，不关会在写 obs 字符串列时崩，坑 #107）。
- 随机性只来自 `BIOANALYS_SEED`（PCA/neighbors/UMAP/采样的 `random_state` 全部传入），
  保证同 seed 可复现；诊断用采样固定 seed。
- 数值口径写常量不放裸字面量（`TARGET_SUM = 1e4`、`N_TOP_GENES = 2000`），
  跨步骤共享口径先在上游报告里定稿，下游引用不改动。
- 注释只解释"为什么"：口径决策、绕过某库 bug 的写法、与论文/等价工具的对应关系。
  不写"这是什么"。

## 算力与内存红线（硬件受限运行面硬约束）

本分析流水线运行于物理内存 12GB（容器默认上限 7.5GB，单任务最高 8.0GB，由 `BIO_MEM_LIMIT_GB` 控制，纯 ARM CPU）的真机设备上。
为保证长时间批处理稳定运行，对各步骤的算法选型与内存消耗实行刚性熔断限制：

- **空间与单细胞解卷积必须保持轻量 `scipy.optimize.nnls`**：
  单步内存消耗必须控制在 $\le 2.2\text{ GB}$ 以内，运行时间控制在 $6\text{ 分钟}$ 以内。
- **严禁在生产流水线上引入重型模型**：
  严禁引入 `cell2location`、`Stereoscope`、`DestVI`（依赖 PyTorch 变分推断，纯 ARM CPU 上运行需数小时且内存易突破 $10\text{ GB}$）
  或原生 `RCTD`（在内存展开数万单细胞稠密矩阵，容易触发 Android LMK 强杀退出码 137）。
- **零算力增耗的解卷积精度优化路线**：
  若需优化解卷积质量，一律通过数据流优化实现，算力 0 增加：
  1. **输入端 Marker 严格提纯**：每个细胞类型严格仅挑选特异性最高的 Top 20 个 Marker 基因，剔除泛表达管家基因，降低优化求解矩阵维度；
  2. **输出端阈值硬截断**：吸收 RCTD 杂音拒绝理念，对输出丰度低于 $0.05$（$5\%$）的极微弱预测值直接归零，消除背景弥散假阳性。
- **参考代码与生产代码职责分离**：
  `可用代码预筛选/06.补充学习代码/` 下的 RCTD、Seurat 等 R 代码，仅作为离线 PC 知识沉淀与外部科研客户交付模版，严禁并入真机生产执行流。

## 与台账、隔离的接口

- 台账由引擎（BioanalysDesktop server）落 `system_logs/`，脚本只管 stdout/stderr 与退出码；
  但脚本输出必须包含可解析的关键行（`step_executed=<step>` 开头，关键节点 flush 输出），
  供 SSE 日志流与巡检解析。
- 跨线引用（输入/输出路径指向本线目录外）由审计前置拒绝，脚本自身不绕过：
  输入路径出现他线字符串即非零退出并说明，不"帮忙"找文件。

## 扩展名必须与解释器声明一致（2026-11-02 新增，机器可判定）

文件扩展名是**给工具看的执行器约定**，不是备注。首行 shebang 与扩展名不一致，会让三处同时出错：
静态门禁按扩展名把文件送进错误的语法检查器、编辑器 LSP 挂错语言服务、
`bridge_guard.sh` 按扩展名推断执行器。这类错误不会在运行时立刻炸，而是静默跑偏。

实测案例：`src/*/checkpoint/checkpoint_*_bridge.py` 三个文件首行是 `#!/usr/bin/env bash`，
内容是合法 bash（316 行），却被命名为 `.py`。`python -m py_compile` 报 `SyntaxError: unmatched ')'`。
已改名 `.sh`。

| 首行声明 | 允许的扩展名 |
|---|---|
| `#!/usr/bin/env bash` / `sh` / `zsh` | `.sh` |
| `#!/usr/bin/env python` / `python3` | `.py` |
| `#!/usr/bin/env Rscript` | `.R` |
| `#!/usr/bin/env node` | `.js` / `.mjs` / `.cjs` |

判据：`node tools/dev/static_gate.mjs` 的 `shebang` 门禁（读前 128 字节，逐文件比对）。
新增脚本时若扩展名写错，这条门禁会在 CI 与本地同时拦下。
