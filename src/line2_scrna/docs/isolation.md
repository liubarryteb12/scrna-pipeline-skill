# 三线隔离约束

> version 2 | 2026-11-08
> 术语定义见 `GLOSSARY.md`，本文件讲"约束是什么 + 怎么执行"。

## 约束内容

三条分析线：`line1_bulk`、`line2_scrna`、`line3_spatial`。

**默认互不交叉**（架构级硬约束，不是工作习惯）。唯一例外是本章末尾定义的受控桥接：

1. 默认不互相取参考。线3 做解卷积时，参考取自 `public_references/` 或线3 自有单细胞分辨率数据。
   需要从线2 导出细胞类型签名时，必须走「受控桥接」：线2 只导出冻结的签名表，
   线3 只读取 `public_references/` 下的冻结表，不读取线2 的 `raw_data/`、`processed/` 或 `h5ad`。
2. 不合并数据集。三条线的 `raw_data/`、`processed/`、`results/` 物理分离。
3. 不跨线验证。线1 的外部验证用独立 GEO 数据集，不拿线2/线3 的结果当验证集。

**虚拟敲除不是第四条线**：它是线2 和线3 各自的末端模块，产出 DR genes。

**R/Python 混用只发生在线内**：同一条线内部可以按工具分工换语言，
例如线2 用 Scanpy 做 QC 聚类、用 scTenifoldKnk 做虚拟敲除。这是工具层混用，不是分析线交叉。

## 目录层怎么保证

物理隔离靠目录边界，不靠人记得：

```
/mnt/ssd/bio_data/
├── projects/
│   ├── default/
│   │   ├── line1_bulk/{raw_data,processed,results,figures}/
│   │   ├── line2_scrna/{raw_data,processed,results,figures}/
│   │   └── line3_spatial/{raw_data,processed,results,figures}/
│   └── <新项目>/...          ← 项目之间同样物理隔离
├── bin/<line>/              ← 脚本部署区（全局一份）
├── envs/<line>/             ← 环境声明（全局一份）
└── public_references/       ← 跨项目唯一合法共享读入口
```

每条线的分析代码工作时，工作目录必须设在该项目该线的 `processed/` 内；
脚本的输入输出路径只允许出现：本项目本线目录、`public_references/`、`system_logs/`。

**两层隔离**：项目之间（不同 `<项目>` 目录互不引用）与线之间（同项目内三条线互不引用）都是硬约束。

## 执行层校验

约束不能只写在文档里，由两个脚本分管执行：

1. **代码静态检查**（PC 侧，`tools/ops/02_check_isolation.ps1`）：
   `src/<line>/` 下的脚本，路径字符串中不得出现其他两条线的目录名或 line_id。
2. **台账检查**（服务器侧，`tools/ops/03_check_logs_remote.sh`，PC 侧一键入口 `tools/ops/04_check_logs.ps1`）：
   `system_logs/` 下任务台账 JSON 的 `line_id` 必须是三个合法值之一，`project_id` 必须合法
   （缺省视为 `default`，兼容迁移前的历史台账）；`io_files` 中的路径必须属于该 line_id 或
   `public_references/` 或 `system_logs/`，且路径中的项目段（`/projects/<项目>/`）必须与
   `project_id` 一致（跨项目引用一律判违规）。
3. **运行时路径断言**（后端 `server/files.js` / `figures.js` / `upload.js`）：
   所有浏览、出图、上传路径都必须落在 `projects/<项目>/<线>/` 或 `public_references/` 之内，
   别的项目、别的线、`..` 穿越一律拒绝。

## 桥接检查点（线内混用时必过）

跨语言传文件（如 Python 出 `.h5ad` 给 R 读）后，必须校验四项，任一不符即停：

- 基因名（symbol 是否一致、有没有带版本号后缀）
- 矩阵方向（基因×细胞 还是 细胞×基因）
- 细胞名是否保留
- embedding 是否保留

校验脚本放 `src/<line>/` 下，命名为 `checkpoint_*.R` 或 `checkpoint_*.py`。

## bulk 线的外部验证冻结

独立 GEO 验证集只做冻结评估：不重新筛基因、不重新调参。
验证集内部可以去批次，但**禁止把训练集和验证集合并校正后再评估**。
多外部数据集、跨平台必须分开报告，关注表达方向一致性。

## 受控桥接：线2 细胞类型签名到线3 解卷积

本例外只适用于「同一项目内，线2 已经完成细胞类型注释，线3 需要用细胞类型表达签名做解卷积」这一种场景。
它不是开放跨线读取，也不是允许线3 复用线2 的原始表达对象。

### 合法产物

线2 侧只允许导出一个冻结签名表，落点固定为：

```
public_references/celltype_signatures/<project_id>/signatures.tsv
public_references/celltype_signatures/<project_id>/bridge_manifest.json
```

`signatures.tsv` 首列为 `gene_symbol`，其余每列一个细胞类型，数值域必须在报告中写明。
`bridge_manifest.json` 至少记录：项目 ID、源线、目标线、导出脚本、原始产物名、sha256、
细胞类型列表、基因数、物种、基因名规范、值域、正式 ortholog 表或 symbol 映射方式、seed、导出时间。

### 执行规则

1. 线2 侧导出脚本只允许写 `public_references/celltype_signatures/<project_id>/`，不允许把线3 路径写进脚本。
2. 线3 侧解卷积脚本只接收冻结签名表的路径作为输入，不允许出现线2 目录名、线2 的 h5ad 名或线2 的 `processed/` 路径。
3. 运行时台账必须把冻结签名表记录为 `public_references/` 输入，并在备注里写桥接 ID。
4. 物种不同源时必须记录正式 ortholog 表；仅大写 symbol 同名匹配时，报告必须保留 `map_mode=symbol_exact_upper` 与匹配率，
   禁止把近似映射表述为正式 ortholog。
5. 解卷积报告必须写细胞类型数、匹配基因数、未匹配数、失败 spot 数、线性混合假设和参考来源哈希。
6. 任何新增桥接产物类型，先更新本文件与 `workflows.json` 的 `allowed_cross_series_bridges`，再写代码。

### 历史基准

`PRJ_DCM_BENCH_20261105` 的适配器历史台账里保留了线2 原始 `signatures.tsv` 路径。该记录用于复核当时执行事实，
不作为后续项目的输入契约。后续项目一律按本章的 `public_references/` 冻结签名表执行。
