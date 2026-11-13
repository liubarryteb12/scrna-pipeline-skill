# 工程规范: 目录结构

> version 3.0 | 2026-10-26
> 本规范规定"文件放哪"。新会话接手先读本文件，不要猜。

## 一句话原则

**本地控制面保留代码、治理台账与轻量回搬的项目产物；原始大数据、运行环境和计算都在服务器侧。**
单个项目的数据、图件、调试、交付与台账统一收口在 `workspace/projects/<项目>/`。

## 本地工程目录（本仓库）

```
AI辅助生信分析流程工作台/
├── README.md                # 项目入口：定位、架构、快速开始、算力接口
├── .env.example             # 算力配置模板；真实值只放未跟踪的 .env.local
├── configs/                 # 算力接口示例配置，不保存真实设备信息
├── _private/                # 本机隐私区：设备档案、历史交接件；已 gitignore，不入库
├── SPEC.md / ROUTER.yaml / SCOPE.md / GLOSSARY.md /
│   MIGRATION.md / INDEX.yaml / LICENSES.md / AGENT.md
│                             # 治理文件（2026-11-08 整体重构 W1 上移，契约见 MIGRATION.md）
├── core/ hooks/ schemas/ tools/ workspace/ ledger/
│                             # 参考结构骨架目录（分波迁入，占位 README 见各目录）
├── skills/                   # SOP 规约资产（W7 起迁入，stable 与 in-progress 分级）
│   └── stable/{bio,dev,ops}/ # 已验收技能，每个技能含 skill.yaml 与 prompt.md
├── workspace/               # 数据与制品工作区（客户数据与控制面物理隔离）
│   ├── projects/            # 本地项目唯一根，结构与远端项目语义对应
│   │   ├── INDEX.md         # 全项目机读摘要入口
│   │   └── <项目>/          # 自包含项目目录，见下节
│   └── handoff/             # 跨会话交接、决策、历史派单与调研快照
│       ├── handoff.md
│       ├── SPEC.md
│       ├── SPEC_DECISIONS.md
│       ├── dispatch/
│       ├── topic_survey/
│       └── archive/
├── ledger/                  # 统一知识台账（追加式 JSONL 与运维/坑位台账）
│   ├── actions.jsonl / changes.jsonl / errors.jsonl / fixes.jsonl / pitfalls.jsonl
│   ├── ops_ledger.md / pitfalls_ledger.md / 踩坑台账.md
│   └── sync/                # inbox、outbox、daemon.log
├── outputs/
│   └── system_logs/         # 仅存回搬的机器执行台账，与 JSONL 台账独立
├── references/               # 参考资料唯一存放地，内容只读，结构可整理（约定见下）
├── docs/                    # 工程规范与架构文档
│   ├── directories.md       # 本文件
│   ├── coding-standards.md  # 分析代码规范（Linux 运行面，tools/bio/ 动工前必读）
│   ├── isolation.md         # 三线隔离约束与校验
│   ├── ledger-schema.md     # 任务台账 JSON schema
│   └── runs/                # 历史运行记录与验证证据
├── tools/
│   ├── bio/                 # 三条线分析母版与注册表（W3 已迁入，git 历史保留）
│   │   ├── line1_bulk/
│   │   ├── line2_scrna/
│   │   ├── line3_spatial/
│   │   ├── workflows.json
│   │   └── figure_templates.json
│   ├── dev/                 # 门禁、库与工单编排（W5a 已迁入）
│   ├── reporting/           # 报告、交付与验证链（W5b 已迁入）
│   └── ops/                 # 运维入口与算力面脚本（W5b 已迁入）
└── BioanalysDesktop/        # 管理员工作台应用（W6 评估保留）
```

## 路径双轨

仓库内路径与设备路径语义分开，二者不互相推导：

| 面 | 项目根 | 台账根 | 同步目录 |
|---|---|---|---|
| 仓库 | `workspace/projects/` | `ledger/` | `ledger/sync/` |
| 设备 | `/mnt/ssd/bio_data/projects/` | `/mnt/ssd/bio_data/system_logs/` | 设备侧运行目录 |

仓库侧改目录只影响本地开发与回搬视图。设备侧固定布局由 `docs/isolation.md` 与
`tools/bio/` 注册表约束，不随仓库目录改名。

## 服务器侧目录

```
/mnt/ssd/bio_data/                      # 算力节点数据盘根
├── projects/                           # 项目根：一个项目一套三线，项目之间物理隔离
│   ├── default/                        # 默认项目（历史三线数据迁入处）
│   │   ├── line1_bulk/                 # GEO bulk 线
│   │   ├── line2_scrna/                # 单细胞 + 虚拟敲除线
│   │   ├── line3_spatial/              # 单细胞空间 + 虚拟敲除线
│   │   └── project.json                # 项目元数据（名称 / 描述 / 创建时间）
│   └── <新项目>/...
│       ├── <line>/scripts/             # 【代理自管】本项目适配后的脚本副本（方案A，框架不创建）
│       └── deliverables/               # 【代理自管】客户交付包三大件（框架不创建）
├── bin/<line>/                         # 脚本部署区（代码本体，全局一份，所有项目共用）
├── envs/<line>/                        # 环境声明（跟代码走，全局一份）
├── public_references/                  # 公共参考（跨项目唯一合法共享读入口）
└── system_logs/                        # 全局台账 JSON（条目带 project_id 区分归属）
```

每条线内部数据子目录固定为 `raw_data/`、`processed/`、`results/`、`figures/`；
脚本 `bin/` 与环境声明 `envs/` 在顶层全局一份，不随项目复制（同一线的所有项目共用同一份代码与环境）。

**两个代理自管目录（不在框架 `LINE_SUBDIRS` 内，工作台不创建、不感知）**：

| 目录 | 用途 | 谁创建 |
|---|---|---|
| `projects/<项目>/<line>/scripts/` | 客户数据需适配时，从 `bin/<line>/` 复制母版到此改，跑通后留档；不改母版 | 代理 `mkdir -p` |
| `projects/<项目>/deliverables/` | 客户交付包三大件（适配代码 / 图表与表格 / 解读报告） | 代理 `mkdir -p` |

交付包结构与打包规则见 `docs/delivery_package_spec.md`（v1.0，**实现现状见该文「零」节**）：打包由服务器侧 `tools/reporting/build_client_report.py` 承担，工作台入口是 `POST /api/report/build`（`BioanalysDesktop/server/report.js`，含存储水位 / 进行中任务 / 失败台账三道前置门禁与 `.sha256` 旁挂）；原规划的 `tools/ops/09_package_delivery.ps1` 已作废，不再实现。

适配副本不经工作台任务引擎（它只认全局 `bin/<line>/`），而是直跑运维包装器：
`bash /mnt/ssd/bio_data/ops/run_throttled_task.sh <task_id> <line> <step> <R|Python> <项目内脚本绝对路径>`。

本地回搬入口（`tools/ops/07_pull_outputs.ps1`）按项目写入
`workspace/projects/<项目>/data/<line>/{processed,results}` 与
`workspace/projects/<项目>/figures/<line>/`；机器执行台账独立落 `outputs/system_logs/`。

## 本地项目存档与 sync 运行态

`workspace/projects/` 是本地项目唯一根，只接轻量回传与项目台账，不替代算力面工作目录。
每个项目由 `tools/ops/init_project_scaffold.py` 幂等创建：

```text
workspace/projects/<项目>/
├── project.json
├── handoff.md
├── README.md
├── INDEX.json
├── data/<line>/{processed,results}
├── figures/<line>
├── figures/adapt
├── deliverables/
├── debug/
└── ledger/{actions,changes,errors,fixes,pitfalls}.jsonl
```

建项命令：

```powershell
python tools/ops/init_project_scaffold.py --root . --project <项目> --line <line> [--ticket <id>]
```

`figures/<line>/` 收真实结果图，`figures/adapt/` 收项目内模板派生副本，两者不得混用。
脚手架重复运行只补缺失项，不覆盖项目元数据、README、交接文档或已有台账内容。

全局 `ledger/` 与项目级 `ledger/` 分离：前者只收跨项目结论和全局动作，
后者只收本项目事实。`outputs/system_logs/` 或服务器侧 `system_logs/` 是独立机器台账，
不与 JSONL 台账合并。

`ledger/sync/` 是工作区内的 sync 运行态，只放 `inbox/`、`outbox/`、`daemon.log`。
它是 spool，不是数据仓库。

存量 `default`、`PRJ_20261019_IS_lactylation` 及分散项目资产由 `WP-2.1` 执行零损迁移。
迁移完成前，旧 `outputs/` 项目产物只作兼容存量，新任务一律按 `workspace/projects/<项目>/`
落盘并单独记录迁移前后的文件数与 SHA256。

## 归档分类（写文件前先判断类别）

| 类别 | 落点 |
|---|---|
| 交接上下文 | `workspace/handoff/` |
| 工程规范、架构文档 | `docs/` |
| 运维脚本 | `tools/ops/` |
| 分析代码 | `tools/bio/<line_id>/` |
| 项目专属调试与日志 | `workspace/projects/<项目>/debug/` |
| 跨项目临时探针 | `cn_debug/` |
| 运行记录、验证报告 | `docs/runs/` |

运行记录命名：`YYYY-MM-DD_[类型]_[简短说明]`，类型 ∈ `initial` 基线 / `debug` 调试 / `failed` 失败现场 / `verification` 验证 / `formal` 正式产物。

参考资料是唯一例外：它按内容只读，按结构可整理。

## references/ 的读写规则

`references/` 是参考资料唯一存放地，分两层对待：

内容层只读：已入库的参考文件（实测经验、规范摘录、可用文献、出图代码）正文不改写、
不删除。要改先复制到工作副本，改完作为新文件入库，原文件留痕。

结构层可整理：目录归类、可用性判定、把无用件移到待整理区，这类动作合法，
是维护动作不是破坏。整理前先出盘点清单，整理后留台账。

待整理区固定为 `references/09.待整理/`，存放判定为无用或暂不可用的文件。
移入时保留原始相对路径，便于回溯。待整理区内的文件不参与分析，不引用，定期清。

命名沿用参考区既有风格：`NN.中文主题名`。新增顶类目时编号顺延，不重排已有编号。

### references/ 子目录编号与用途

盘点快照与导航见 `references/README.md`（2026-11-07 盘点，逐文件明细在
`references/09.待整理/inventory_20261107.tsv`）。

| 编号 | 中文名 | 文件数 | 用途 |
|---|---:|---:|---|
| `01` | 工程规范 | 188 | mattpocock-skills 上游工程规范参考 |
| `02` | 分析流水线构建参考与可用参考文献 | 49 | 文献 PDF，机读源 `citations.json`，索引由 `tools/ops/gen_reference_index.py` 生成 |
| `03` | 出图规范 | 2 | `sci.md` 规范正本与单图代码模板包 |
| `04` | 正式使用生信代码 | 0 | 预留空目录，待正式启用代码入驻 |
| `05` | 可参考生信代码 | 3044 | `可用代码预筛选/`（可用与源码）与 `library/` 图录；在役出图证据在 `可用代码预筛选/源码/` |
| `06` | 测试用原始数据 | 22 | 三线回归测试数据集四套（L1 两套、L2/L3 各一套），分类见其 README |
| `07` | Visium_v1_载玻片布局 | 1 | Visium v1 全局 barcode 坐标布局表（git 跟踪） |
| `09` | 待整理 | 2204 | 隔离区，重复件、备份与迁移留痕，保留原始相对路径 |

历史沿革：原顶层 `00a.可参考代码库` 与 `05.可参考出图代码` 已并入
`05.可参考生信代码/可用代码预筛选/源码/`，消费方路径与 gitignore 白名单已同步
（见 `docs/REFERENCE_CLEANUP_LEDGER.md`）；
原 `06.Visium_v1_载玻片布局` 改编号为 `07`。

新增顶类目编号顺延，不重排已有编号。

项目分组及对应目录维护在 `references/项目分组表.tsv`。

## 命名与写法红线

- 标识符、路径用英文；面向用户的说明用中文。
- **Linux 系统内（服务器侧）一律英文路径与文件名**，不出现中文字符：本地既有中文名的
  用户资产（如 `test_use_data/` 下的目录与 `分组情况.txt`）推到服务器前必须换成英文名，
  映射关系集中在引用处（当前在 `tools/ops/05_push_testdata.ps1` 的 `$RemoteNameMap`）。
- 注释只解释"为什么"，不写"这是什么"。
- 注释和文档不用破折号，用逗号、冒号、句号替代。
- 变量命名不用单字母（数学符号和循环计数除外）。
- 不自动重命名、移动或覆盖用户原有文件。
- linux上运行的代码必须符合linux的命名规范。

## 服务器侧执行面

| 面 | 别名 | 用途 | 不该做什么 |
|---|---|---|---|
| 基础用户空间 | 可配置别名 | 文件传输、通道管理、服务看护 | 不跑分析计算 |
| Linux 容器 | 可配置别名 | R/Python 科学计算执行 | 不作为文件长期存放地 |

## 环境隔离

- 容器内 R/Python 环境按线固定：每条线一个独立 conda env 或独立 R 库路径，记录主语言与版本到该线 `envs/`。
- 版本、随机种子写进该线 `envs/manifest.txt`，分析中途不更换。
- 跨线不共用编译产物、不共用 conda env。
