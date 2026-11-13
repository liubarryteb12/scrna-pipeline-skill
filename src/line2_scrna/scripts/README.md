# scripts/ 与 checkpoint/ 说明

## 目录约定

1. **一个文件一个步骤**：文件名就是任务面板里的步骤名，只使用字母、数字与 `. _ -`。参数通过 `--input` / `--output` / `--seed` 传入。
2. **必须带 shebang 且可执行**：bash 包装或直接使用 `Rscript` / `python3`。部署通道会设置 755，本地也保持可执行位。
3. **工作目录是 `<line>/processed/`**：任务引擎进入该目录后执行脚本，脚本用相对路径读写中间产物，或用 `--input` / `--output` 传绝对路径。
4. **只使用本线 line_id**：不得出现另外两条线的字面量，提交前由隔离审计校验。
5. **桥上检查点脚本**：命名 `checkpoint_*.sh`，放 `checkpoint/`，必须校验基因名、矩阵方向、细胞名/样本名、embedding 四项。
6. **注释只解释为什么**：不使用 em dash，变量不用单字母。

## 部署

生产环境由运维层把本目录部署到算力节点的 `bin/line2_scrna/`，任务引擎已把该目录注入 PATH。本地改脚本后必须重新部署才生效。

数据、产物与报告落在 `<BIO_DATA_ROOT>/projects/<project>/line2_scrna/`，不进入本仓库。
