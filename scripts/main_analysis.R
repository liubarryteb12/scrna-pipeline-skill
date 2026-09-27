#!/usr/bin/env python3 → R 移植（r_version.md §7 / R-08）
# main_analysis.R — 编排全部步骤 + 验收清单（main_analysis.py 的 R 版镜像）
#
# **验收清单是这一步存在的理由。** 单细胞流水线的失败模式大多是
# "job 绿了但结果不对"：某个可选步骤静默跳过、某张图是空白的、
# 某个中间产物没写出来。这些在 CI 日志里和成功长得一模一样。
#
# 所以最后强制检查（与 Python 版 checks[].id 逐条对齐，r_version.md §7）：
#   - 必需步骤是否都 ok（含 result_status 消费，E-56）
#   - 必需产物文件是否都在（且非空）
#   - 图是否都真的产出（声明扫描，M11；动态槽位，E-70）
#   - 嵌套 status 里有没有 failed（Q-26 / E-48）
# **任何必需项失败 → 退出码 1 → CI 红。** 可选步骤失败只在报告里标出。
#
# 运行方式：必须从**仓库根**运行 `Rscript scripts/main_analysis.R`
# （record/load_config 的产物路径是相对路径 results/<id>、data/<id>；
#   run_00_fetch 的 repo_root 在有 --file= 时自行推导，与 cwd 无关）。
#
# R 版适配点（与 Python 逐条对照，改动处都有注释）：
#   1. INPUT_FILES 中间对象用 .rds（r_version.md §5：中间对象不跨语言，
#      文件名扩展允许不同）；REQUIRED_FILES 的 results 侧文件名两版相同。
#   2. 步骤函数经 sys.source 载入独立环境后按名取出（9 个脚本的尾部守卫
#      `sys.nframe()==0L || SCRNA_STEP_MAIN` 保证 source 安全）；
#      source 前 options(scrna.script_path=…) 让 get_script_path 可解析。
#   3. named_tools 用 R 版登记表（requireNamespace 实测 available）：
#      Python 那张"R 包一个都装不上"的表对 R 版不成立，照抄是写假话；
#      edgeR / Slingshot 的 name_taken 两条（PyPI 同名无关包）在 R 版无意义，
#      登记表里没有 kind="name_taken"，验收层的 squat 判据随之不适用。
#   4. 08 的方法自洽判据按 R 版状态结构适配：R 版没有顶层 method 串，
#      是 method_bits 向量（08_virtual_perturbation.R L544-556）；
#      engine 字段叫 engine_configured；tenifold 失败原因在 method_bits 里
#      （"Tenifold 引擎未产出：…"），没有独立的 tenifold 子状态 dict。
#   5. 空敲除判据消费 R 版状态键 n_empty_knockouts（复数，Python 是
#      n_empty_knockout）；基因名单不在状态里，从产物
#      virtual_perturbation_tenifold.csv 的 empty_knockout 列读。
#   6. part3 契约键名：R 版是 counts_layer（Python 是 layers['counts']）。
#   7. chk() 的**第二个位置参数是 kind 不是 severity**（E-53）——
#      两个参数名分明，severity 有默认值且在末位。

# ---------------------------------------------------------------------------
# 路径与公共层载入（先于一切；必须在任何步骤之前把 common.R 载入 globalenv）
# ---------------------------------------------------------------------------
.scripts_dir <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1L]), mustWork = FALSE))
  else getOption("scrna.main_scripts_dir", "scripts")
})
REPO <- dirname(.scripts_dir)
source(file.path(.scripts_dir, "lib", "common.R"))

# get_script_path 的 source 场景兜底：01/03 用 dirname(dirname(gsp)) 到仓库根、
# 06/07 用 file.path(gsp, "..", "assets")、08 用 file.path(gsp, "lib", …)。
# 给它一个 scripts/ 下的文件路径，三处算法都解析到正确位置。
options(scrna.script_path = normalizePath(
  file.path(.scripts_dir, "main_analysis.R"), mustWork = FALSE))

# ---------------------------------------------------------------------------
# 表结构（与 main_analysis.py L38-132 逐条对齐）
# ---------------------------------------------------------------------------

# (步骤 id, 模块文件, 函数名, 是否必需, 中文名)
STEPS <- list(
  list(id = "fetch", mfile = "00_fetch.R", fn = "run_00_fetch",
       required = TRUE,  label = "取数与校验"),
  list(id = "qc", mfile = "01_qc.R", fn = "run_01_qc",
       required = TRUE,  label = "质量控制"),
  list(id = "integrate", mfile = "02_integrate.R", fn = "run_02_integrate",
       required = TRUE,  label = "标准化与降维"),
  list(id = "cluster_annotate", mfile = "03_cluster_annotate.R",
       fn = "run_03_cluster_annotate", required = TRUE, label = "聚类与注释"),
  list(id = "pseudobulk_de", mfile = "04_pseudobulk_de.R",
       fn = "run_04_pseudobulk_de", required = FALSE, label = "拟bulk 差异表达"),
  list(id = "trajectory", mfile = "05_trajectory.R", fn = "run_05_trajectory",
       required = FALSE, label = "轨迹推断"),
  list(id = "communication", mfile = "06_communication.R",
       fn = "run_06_communication", required = FALSE, label = "细胞通讯"),
  list(id = "grn", mfile = "07_grn.R", fn = "run_07_grn",
       required = FALSE, label = "转录因子调控"),
  list(id = "virtual_perturbation", mfile = "08_virtual_perturbation.R",
       fn = "run_08_virtual_perturbation",
       required = FALSE, label = "虚拟敲除/过表达（§1.7/§1.8 保留框架）")
)

# 每步会写的状态文件。**跑之前先删掉** —— 否则步骤崩溃时旧文件还在，
# 下游"产物存在"检查读的是**上一轮的**结果，会给出虚假的通过。
# 实测踩过：07_grn 因漏 import 崩了，而 grn_status.json 是上一轮的，
# 验收照样 40 项全绿。
STEP_STATUS_FILES <- list(
  "qc"                   = "qc_status.json",
  "integrate"            = "integration_status.json",
  "cluster_annotate"     = "cluster_status.json",
  "pseudobulk_de"        = "pseudobulk_status.json",
  "trajectory"           = "trajectory_status.json",
  "communication"        = "communication_status.json",
  "grn"                  = "grn_status.json",
  "virtual_perturbation" = "virtual_perturbation_status.json"
)

# **M10（R-03 裁决）**：这张表的键必须与 STEPS 的 id 集合一致；fetch 是唯一
# 豁免（它写 data_dir/dataset_info.json，取数步骤的产物是数据本身，
# "状态"没有独立载体）。启动期对齐断言 —— 将来加步骤忘了登记，**启动即报错**
# （否则 unlink 会落到结果目录上，IsADirectoryError 被 tryCatch 接住 ⇒
# 步骤莫名失败）。
.STEPS_WITHOUT_STATUS_FILE <- "fetch"
.MISSING_STATUS_FILES <- setdiff(
  vapply(STEPS, function(s) s$id, ""),
  c(names(STEP_STATUS_FILES), .STEPS_WITHOUT_STATUS_FILE))
if (length(.MISSING_STATUS_FILES)) {
  stop(sprintf(paste0(
    "STEPS 里的步骤既没有登记状态文件、也不在豁免集里: %s",
    " —— 没有登记的话跑前删旧状态文件的逻辑会落到结果目录上（审计 M10）。"),
    paste(.MISSING_STATUS_FILES, collapse = ", ")), call. = FALSE)
}

# 文档 §2「本部分人工复核节点」。**默认 pending，不是 confirmed** ——
# 自动化流水线不能替人签字。验收里作为**可见但不阻断**的项列出。
HUMAN_REVIEW_NODES <- list(
  list(node = "celltype_labels", label = "细胞类型注释最终标签", required = TRUE),
  list(node = "cluster_resolution", label = "聚类分辨率选择依据", required = TRUE),
  list(node = "pseudobulk_design", label = "拟bulk 差异分析设计", required = FALSE),
  list(node = "trajectory_direction", label = "拟时序轨迹方向确认（marker 验证）", required = TRUE),
  list(node = "trajectory_branches", label = "拟时序分支点的生物学解释", required = TRUE),
  list(node = "virtual_perturbation_targets",
       label = "虚拟扰动靶基因的生物学合理性", required = FALSE)
)

# 需要登记哈希的输入（相对 data_dir）。**R 版适配：中间对象是 .rds**
# （r_version.md §5——中间对象不跨语言，语义标签一致、扩展允许不同）。
INPUT_FILES <- list(
  list(name = "raw.rds",          desc = "原始计数矩阵", required = TRUE),
  list(name = "dataset_info.json", desc = "数据集元信息", required = TRUE),
  list(name = "qc_filtered.rds",  desc = "QC 后矩阵",    required = TRUE),
  list(name = "integrated.rds",   desc = "整合后矩阵",   required = TRUE),
  list(name = "clustered.rds",    desc = "聚类后矩阵",   required = TRUE)
)

# 必需产物（相对 results_dir）。这些文件名两版相同（状态 JSON schema 同、
# CSV 同名，r_version.md §5）。
REQUIRED_FILES <- list(
  list(name = "state.json",                        desc = "步骤状态",        required = TRUE),
  list(name = "qc_status.json",                    desc = "QC 结论",         required = TRUE),
  list(name = "integration_status.json",           desc = "降维与整合结论",   required = TRUE),
  list(name = "cluster_status.json",               desc = "聚类与注释结论",   required = TRUE),
  list(name = "markers_all.csv",                   desc = "marker 基因表",    required = TRUE),
  list(name = "cluster_resolution_scan.csv",       desc = "分辨率扫描",       required = TRUE),
  list(name = "celltype_annotation.csv",           desc = "细胞类型注释",     required = TRUE),
  list(name = "celltype_scores.csv",               desc = "细胞类型打分矩阵",  required = TRUE),
  list(name = "qc_cells.csv",                      desc = "每细胞 QC 指标",   required = TRUE),
  list(name = "pseudobulk_status.json",            desc = "拟bulk 状态",      required = FALSE),
  list(name = "trajectory_status.json",            desc = "轨迹状态",         required = FALSE),
  list(name = "communication_status.json",         desc = "通讯状态",         required = FALSE),
  list(name = "grn_status.json",                   desc = "GRN 状态",         required = FALSE),
  list(name = "virtual_perturbation_status.json",  desc = "虚拟扰动状态",     required = FALSE)
)

# 必需的图（相对 figures_dir，不含扩展名）。
#
# **M11（R-03 裁决）：这张表不再单独承担判据** —— 真正的判据是
# declared_figures()（从源码扫声明）。这张表退化成**说明文字**；
# 下面的 figures:documented 检查专门盯"扫出来的图有没有说明"。
REQUIRED_FIGURES <- list(
  c("02-01-01-unit1-genes-detected",     "过滤前 QC: genes detected"),
  c("02-01-01-unit2-total-counts",       "过滤前 QC: total counts"),
  c("02-01-01-unit3-mito-fraction",      "过滤前 QC: mito fraction"),
  c("02-01-01-unit4-ribo-fraction",      "过滤前 QC: ribo fraction"),
  c("02-01-01-unit5-hb-fraction",        "过滤前 QC: hb fraction"),
  c("02-01-02-unit1-qc-scatter-thresholds", "QC 阈值散点"),
  c("02-02-01-unit1-hvg-selection",      "高变基因选择"),
  c("02-02-02-unit1-pca-variance-ratio", "PCA 方差解释"),
  c("02-02-03-unit1-batch-mixing",       "批次混合前后对比"),
  c("02-03-01-unit1-cluster-resolution-scan", "分辨率扫描曲线"),
  c("02-03-02-unit1-umap-clusters",      "UMAP 聚类图"),
  c("02-03-03-unit1-markers-dotplot",    "marker 点图"),
  c("02-03-04-unit1-celltype-scores-heatmap", "细胞类型打分热图"),
  c("02-05-01-unit1-paga-graph",         "PAGA 连接图（M11 补）"),
  c("02-05-02-unit1-trajectory-method-correlation", "轨迹方法一致性"),
  c("02-05-03-unit1-trajectory-modules-heatmap", "轨迹模块热图"),
  c("02-05-03-unit2-trajectory-module-profiles", "轨迹模块轮廓"),
  c("02-05-04-unit1-pseudotime-consensus", "共识拟时序 UMAP"),
  c("02-05-04-unit2-pseudotime-dpt",     "DPT 拟时序 UMAP"),
  c("02-05-04-unit3-celltype-on-umap",   "细胞类型 UMAP（M6，条件产出）"),
  c("02-05-04-unit4-pseudotime-principal-path", "拟时序主路径+root"),
  c("02-05-05-unit1-pseudotime-by-cluster", "各簇拟时序箱线图"),
  c("02-05-05-unit2-pseudotime-ridgeline", "拟时序山脊图"),
  c("02-06-01-unit1-communication-heatmap", "细胞通讯热图（M11 补）"),
  c("02-07-01-unit1-tf-activity-vs-pseudotime", "TF 活性沿拟时序"),
  c("02-07-02-unit1-tf-activity-heatmap", "TF 活性热图"),
  c("02-07-03-unit1-tf-specificity-scatter", "TF 特异性散点"),
  c("02-08-01-unit1-virtual-perturbation-effect", "虚拟扰动效应")
)
.FIG_DESC <- vapply(REQUIRED_FIGURES, function(x) x[2L], "")
names(.FIG_DESC) <- vapply(REQUIRED_FIGURES, function(x) x[1L], "")

PART <- "02"

# ---------------------------------------------------------------------------
# 源码扫描辅助（declared_figures / dynamic_fig_bases，M11）
# ---------------------------------------------------------------------------

# 逐行剥注释，**引号内不剥** —— 与 tools/check_fig_names.mjs 的 stripComments
# 同义。口径必须一致：门禁层用 JS 那份扫"声明了哪些图"，验收层用这份扫，
# 两边算法不同就会对同一份源码给出不同的图名集合，而**没有任何东西能发现
# 它们不一致**（E-62 同源：纯文本扫描必须先抹掉非代码区域）。
.strip_comments_r <- function(src) {
  lines <- strsplit(src, "\n", fixed = TRUE)[[1L]]
  out <- character(length(lines))
  for (i in seq_along(lines)) {
    line <- lines[i]
    chs <- strsplit(line, "", fixed = TRUE)[[1L]]
    q <- NA_character_; cut <- NA_integer_
    for (j in seq_along(chs)) {
      ch <- chs[j]
      if (!is.na(q)) { if (identical(ch, q)) q <- NA_character_ }
      else if (ch %in% c('"', "'")) q <- ch
      else if (ch == "#") { cut <- j; break }
    }
    out[i] <- if (is.na(cut)) line else substr(line, 1L, cut - 1L)
  }
  paste(out, collapse = "\n")
}

.script_sources <- function() {
  files <- sort(Sys.glob(file.path(REPO, "scripts", "[0-9][0-9]_*.R")))
  # main_analysis.R 不匹配 [0-9][0-9]_*.R 前缀 —— 与 Python 版同语义：
  # 本文件的验收层会引用图名，那是"检查对象"不是"出图声明"。
  lapply(files, function(p) {
    list(file = p,
         src = .strip_comments_r(paste(readLines(p, warn = FALSE), collapse = "\n")))
  })
}

.FIG_NAME_RE <- "^02-\\d{2}-\\d{2}-unit\\d+-[a-z0-9]+(?:-[a-z0-9]+)*$"

# 扫出**声明要出**的静态图名（字符串字面量）。M11。
# 含 sprintf 占位符（%）的模板串跳过 —— R 版的运行时拼名对应物是 sprintf，
# 由 DYNAMIC_FIG_BASES 声明豁免（Python 版跳的是含 { } 的 f-string）。
declared_figures <- function() {
  out <- character(0)
  pat <- '"02-[^"]*"'
  for (e in .script_sources()) {
    for (m in regmatches(e$src, gregexpr(pat, e$src, perl = TRUE))[[1L]]) {
      nm <- substr(m, 2L, nchar(m) - 1L)
      nm <- sub("\\.(pdf|png)$", "", nm)
      if (grepl("%", nm, fixed = TRUE)) next
      if (grepl(.FIG_NAME_RE, nm, perl = TRUE)) out <- c(out, nm)
    }
  }
  sort(unique(out))
}

# 扫出 `DYNAMIC_FIG_BASES <- list("<图号>" = <张数>)` 声明（R 写法，L67/L25）。
# 键补全成完整前缀 `02-<模块号>-<图号>`。这是**槽位上限**，不是精确值：
# 少出合法（只有 top 8 个 marker 基因时 8 个槽位里出不满），
# 但**一张都没有说明那段循环整段没跑**。
dynamic_fig_bases <- function() {
  out <- numeric(0)
  pat <- 'DYNAMIC_FIG_BASES\\s*<-\\s*list\\(([^)]*)\\)'
  for (e in .script_sources()) {
    m <- regmatches(e$src, regexpr(pat, e$src, perl = TRUE))
    if (!length(m)) next
    body <- sub(pat, "\\1", m, perl = TRUE)
    pairs <- regmatches(body, gregexpr('"\\d{2}"\\s*=\\s*\\d+', body, perl = TRUE))[[1L]]
    for (p in pairs) {
      key <- gsub('"', "", regmatches(p, regexpr('"\\d{2}"', p, perl = TRUE))[[1L]])
      val <- suppressWarnings(as.integer(
        gsub("[^0-9]", "", sub('^"\\d{2}"\\s*=\\s*', "", p, perl = TRUE))))
      if (is.na(val)) next
      mod <- substr(basename(e$file), 1L, 2L)
      out[[sprintf("%s-%s-%s", PART, mod, key)]] <- val
    }
  }
  out
}

# 条件产出的图：**图名 -> (状态文件, 判据路径, 能力就位时的取值, 说明)**。
# **M6（R-03 裁决）**：`05_trajectory.R` 的 unit3 在 `celltype` 列缺失时
# 回退成按 `leiden` 着色，然后**跳过保存**。这里把它变成可判定的：
# `celltype` 列存在时这张图**变成必需**，不存在时允许缺失但**必须可见**。
# 语义与 Python 版一致：**不是"已知缺陷白名单"**，且是**自愈**的。
# （R 适配：R 版 cluster_status.json 的 annotation.status 是两态
#   ok/not_possible —— 判据读的正是这个键，语义不变。）
CONDITIONAL_FIGURES <- list(
  "02-05-04-unit3-celltype-on-umap" = list(
    status_file = "cluster_status.json", path = c("annotation", "status"),
    ready = "ok",
    why = paste0("03 步骤的细胞类型注释没成功时没有 `celltype` 列，",
                 "此时 unit3 会退回 leiden 且按约定不落盘")),
  "02-02-03-unit1-batch-mixing" = list(
    status_file = "integration_status.json", path = "has_batch_key",
    ready = TRUE,
    why = "单样本数据没有批次，没有「前后对比」可画 —— 这张图本就不该存在")
)

# 按路径取值；任一层缺失返回 NULL（**不抛异常**）。
.dig <- function(obj, path) {
  cur <- obj
  for (k in path) {
    if (!is.list(cur)) return(NULL)
    if (is.character(k)) {
      cur <- if (!is.null(names(cur)) && k %in% names(cur)) cur[[k]] else NULL
    } else if (is.numeric(k)) {
      cur <- if (k >= 1L && k <= length(cur)) cur[[k]] else NULL
    } else return(NULL)
    if (is.null(cur)) return(NULL)
  }
  cur
}

# 构造一条验收项（Q-27 范式，与 Python 侧同签名）。
#
# **第二个位置参数是 `kind` 不是 `severity`**（E-53）—— 两者同名同型、
# 位置相邻，传错不会报错，只会把 severity 值写进 kind 槽。所以 kind
# 有名、severity 有默认值且在末位。
#
# `severity` 三档：required（失败即红，默认）/ content（内容正确性）/
# info（只可见，永不判红）。
chk <- function(cid, kind, ok, detail, severity = "required") {
  list(id = cid, kind = kind, item = sprintf("[%s] %s", kind, cid),
       ok = isTRUE(ok), required = severity != "info",
       severity = severity, detail = detail)
}

has_file <- function(p) {
  file.exists(p) && !dir.exists(p) && file.info(p)$size > 0
}

# 嵌套 status 里，哪些取值算"崩了"。其余（not_applied / not_done /
# not_configured / skipped …）都是**设计如此地没做**，只可见、不阻断。
#
# **注意这个集合只管嵌套字段，不要和 common.R 的 STEP_ABORT_VALUES 合并**：
# 顶层返回值里的 bad_root / insufficient_methods 是"这一步没做成"，
# 而嵌套字段里同名或近义的值往往只是"某个可选能力没成"——
# （08 的 tenifold.status 可以是 timeout 或 no_candidates，顶层仍是 ok）。
# 判红会让每个 job 都红，反而没人看（台账元规则 ④）。
NESTED_FAILED_VALUES <- c("failed", "error", "fail")

# 递归产出所有名为 `status` 的字段：list(path, value, 同级 reason)。
# Q-26 / E-48：grn_status.json 的**顶层** status 是 "ok"，而里面
# regulon_vs_pseudotime.status 是 "failed" —— 顶层全绿、里面已经崩了。
.iter_nested_status <- function(obj, path = character(0)) {
  out <- list()
  if (!is.list(obj)) return(out)
  nms <- names(obj)
  if (!is.null(nms) && "status" %in% nms) {
    reason <- obj$reason %||% obj$message %||% ""
    out[[length(out) + 1L]] <-
      list(path = path, value = obj$status, reason = reason)
  }
  if (is.null(nms)) {
    # JSON 数组 → 无名 list：递归子元素，路径用下标（Python 用 str(i)）
    for (i in seq_along(obj)) {
      if (is.list(obj[[i]])) {
        out <- c(out, .iter_nested_status(obj[[i]], c(path, as.character(i))))
      }
    }
  } else {
    for (nm in nms) {
      if (is.list(obj[[nm]])) {
        out <- c(out, .iter_nested_status(obj[[nm]], c(path, nm)))
      }
    }
  }
  out
}

# 状态里的引擎清单与方法自述必须自洽（K-01b 语义，R 版适配）。
#
# 原来这条查的是"方法串里必须出现 'scTenifoldKnk'"来证明**没跑**它；
# 真引擎接进来之后那条判据会**把正确的结果判红**。判据换成自洽性：
#   * 记了 engines_used 就必须逐个在方法自述里被点名；
#   * 没跑 tenifold 时，方法自述必须写明**为什么没跑**（而不是假装跑了）；
#   * PerturbNet 始终没跑，方法自述必须仍然否掉它。
#
# **R 版适配**（08_virtual_perturbation.R L544-575 的状态结构）：
#   - 没有顶层 method 串，是 method_bits **向量** —— 判据对拼接串做；
#   - engine 字段叫 engine_configured；
#   - tenifold 失败原因没有独立子 dict，写在 method_bits 的
#     "Tenifold 引擎未产出：…" 一条里 —— 它就是 td.reason 的 R 版载体。
.vp_engine_consistent <- function(vp) {
  if (!length(vp)) return(FALSE)
  used <- as.character(vp$engines_used %||% character(0))
  bits <- as.character(vp$method_bits %||% character(0))
  method_str <- paste(bits, collapse = "；")
  if (!nzchar(method_str)) return(FALSE)
  # R 版 method_bits 的点名措辞："一阶引擎：…"（L546）/"Tenifold 引擎（…）"（L551）
  names_map <- c(first_order = "一阶", tenifold = "Tenifold 引擎")
  for (e in used) {
    token <- if (!is.na(names_map[e]) && length(names_map[e]) &&
                 !is.null(names_map[[e]])) names_map[[e]] else e
    if (!grepl(token, method_str, fixed = TRUE)) return(FALSE)
  }
  if (!("tenifold" %in% used)) {
    eng <- vp$engine_configured %||% "both"
    td_reasoned <- any(grepl("Tenifold 引擎未产出", bits, fixed = TRUE))
    # 没跑 tenifold：要么配置里就没选它，要么它失败并留了原因（method_bits）
    if (eng %in% c("tenifold", "both") && !td_reasoned) return(FALSE)
  }
  grepl("PerturbNet", method_str, fixed = TRUE)
}

# ---------------------------------------------------------------------------
# §2 点名工具的缺口登记 —— **R 版登记表**（R-08 裁决，b456/b464）
#
# Python 版 NAMED_TOOLS 是"R 包一个都装不上"的表，对 R 版**不成立**：
# R 版就是 R，SoupX / SCTransform / scran / DESeq2 / Monocle3 / CellChat
# 都装得上 —— 照抄就是写假话（common.R named_tools_note 的同一裁决）。
# 所以这里：
#   * r_package 类用 requireNamespace **实测** available（每加一个探测式判断
#     都要能回答"目标存在时它返回什么"—— requireNamespace 对已装包返回 TRUE）；
#   * 去掉 edgeR / Slingshot 两条（kind=name_taken 是 PyPI 同名无关包，
#     R 版没有这个处境）；剩 11 个工具；
#   * reason 保留"这条工具在本流水线的实际处境"，**每个不可用都有原因**。
.NAMED_TOOLS_R <- list(
  list(name = "SoupX", section = "§2.1", kind = "r_package", pkg = "SoupX",
       reason = paste0("ambient RNA 校正（§2.1）。R 版可装（requireNamespace 实测）；",
                       "本流水线 01 走启发式、ambient_rna=heuristic_only，未接 SoupX ",
                       "—— 登记的是能力，不是「已用」")),
  list(name = "CellBender", section = "§2.1", kind = "deps",
       reason = "GPU 导向的去噪工具：GitHub 托管 runner 无 GPU，两版都跑不了"),
  list(name = "SCTransform", section = "§2.2", kind = "r_package", pkg = "Seurat",
       reason = paste0("SCTransform 归一化（§2.2）。包随 Seurat 可装；",
                       "R 版 02 按配置 flavor 走 LogNormalize 路线，未走 SCTransform")),
  list(name = "scran", section = "§2.2", kind = "r_package", pkg = "scran",
       reason = "scran 尺寸因子（§2.2）。R 版 02 未用 scran；包可装"),
  list(name = "SingleR", section = "§2.4", kind = "needs_reference",
       reason = paste0("SingleR 注释需要参考数据集（§2.4）。R 版 03 走 marker 打分，",
                       "避免引入参考集依赖；包可装但缺口是参考集本身")),
  list(name = "DESeq2", section = "§2.5", kind = "r_package", pkg = "DESeq2",
       reason = paste0("拟bulk 差异表达（§2.5）。R 版 04 **已接入** DESeq2",
                       "（pydeseq2 的 R 原版，方法学等级 A）；包可装")),
  list(name = "Monocle3", section = "§2.6", kind = "r_package", pkg = "monocle3",
       reason = paste0("轨迹（§2.6）。monocle3 不在 CRAN（GitHub 装包）；",
                       "R 版 05 用 destiny(DPT) + slingshot，未接 monocle3")),
  list(name = "CytoTRACE2", section = "§2.6", kind = "not_on_pypi",
       reason = paste0("CytoTRACE2 不在 PyPI/CRAN。R 版 05 用**自实现**的 ",
                       "CytoTRACE 风格 GCS（与 Python 版逐行同式），",
                       "并如实标注非官方包")),
  list(name = "scVelo", section = "§2.6", kind = "needs_layers",
       reason = paste0("RNA 速率需要 spliced/unspliced 两套计数（§2.6）；",
                       "输入是 10x filtered 矩阵只有一套计数，两版都记 not_done")),
  list(name = "CellChat", section = "§2.7", kind = "r_package", pkg = "CellChat",
       reason = paste0("细胞通讯（§2.7）。R 版可装（Bioconductor/GitHub）；",
                       "06 按 §3 裁决走 自建 + liana 双路，未接 CellChat（R-17 评估）")),
  list(name = "pySCENIC", section = "§2.8", kind = "needs_resources",
       reason = paste0("GRN 的 motif 剪枝需要 GB 级 cisTarget 排名库（§2.8）；",
                       "07 共表达推断**两版都刻意不做** motif 剪枝，方法保持一致"))
)

probe_named_tools_r <- function() {
  out <- list()
  for (t in .NAMED_TOOLS_R) {
    avail <- if (identical(t$kind, "r_package")) {
      tryCatch(requireNamespace(t$pkg, quietly = TRUE), error = function(e) FALSE)
    } else FALSE
    out[[t$name]] <- list(available = avail, kind = t$kind,
                          section = t$section, reason = t$reason)
  }
  out
}

# ---------------------------------------------------------------------------
# 步骤函数载入：source 到独立环境，按名取出 run_* 函数
# ---------------------------------------------------------------------------
.step_fns <- new.env(parent = globalenv())
for (.s in STEPS) {
  .e <- new.env(parent = globalenv())
  sys.source(file.path(REPO, "scripts", .s$mfile), envir = .e)
  assign(.s$fn, get(.s$fn, envir = .e, inherits = FALSE), envir = .step_fns)
}
rm(.s, .e)

# ---------------------------------------------------------------------------
# run_all —— 编排 + 验收（main_analysis.py L400-1263 的 R 镜像）
# ---------------------------------------------------------------------------
run_all <- function(cfg, only = NULL) {
  .set_orchestrated(TRUE)
  res_dir <- cfg$output$results_dir
  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

  log_info(strrep("=", 68L))
  log_info(sprintf("单细胞流水线 | 数据集 %s", cfg$dataset_id))
  log_info(sprintf("结果目录 %s", res_dir))
  log_info(strrep("=", 68L))

  # ---- 模块零：建立本轮运行清单（§0.3 / §0.4）-----------------------------
  # **必须在任何步骤之前建，且先清掉上一轮** —— 清单描述的是本轮。
  init_manifest(cfg)
  capture_versions(cfg)
  record_params(cfg, list(
    seed = cfg$analysis$seed,
    qc = cfg$qc %||% list(),
    cluster = cfg$cluster %||% list(),
    integration = cfg$integration %||% list(),
    trajectory = cfg$trajectory %||% list(),
    # **这三段原先漏了。** 漏掉的可选步骤参数意味着那几步的结果无法被复现。
    communication = cfg$communication %||% list(),
    grn = cfg$grn %||% list(),
    perturbation = cfg$perturbation %||% list()
  ))
  for (h in HUMAN_REVIEW_NODES) {
    record_human_review(cfg, h$node, required = h$required, status = "pending",
                        note = sprintf("%s —— 需人工确认，本轮自动化未确认", h$label))
  }

  # ---- §2 点名工具的缺口登记 ----------------------------------------------
  # **为什么放在清单里，而不是各步骤的 status 文件里：** "R 包用不上"是
  # 整轮运行的属性，不是某一步的属性。它必须存在，因为产物**看不出来**：
  # 每一步都有东西产出，所以"点名的方法哪些没用上"得自己说出来。
  named <- probe_named_tools_r()
  n_avail <- sum(vapply(named, function(x) isTRUE(x$available), TRUE))
  record_decision(
    cfg, "named_tools",
    "文档 §2.1–§2.8 点名的工具，哪些真的用上了？",
    sprintf("%d/%d 个可用；R 版直接使用原生 R 包/内置替代实现", n_avail, length(named)),
    evidence = named_tools_note())
  record_params(cfg, list(named_tools = named))
  log_info("")
  log_info("§2 点名工具的使用情况（缺口要自己说出来，这是缺口不是已覆盖）：")
  for (nm in names(named)) {
    info <- named[[nm]]
    mark <- if (isTRUE(info$available)) "可用" else "未使用"
    log_info(sprintf("  [%s] %s %s（%s）—— %s", mark, info$section, nm,
                     info$kind, substr(info$reason, 1L, 66L)))
  }
  log_info(sprintf("运行清单：%s", manifest_path(cfg)))

  failed_required <- list()
  # 溢出记录同样要**跑前清空**（与状态文件同理）：否则这一轮修好了，
  # 上一轮留下的记录还在，验收会把已经修好的图判红。
  .ovf_stale <- file.path(res_dir, "figure_overflow.json")
  if (file.exists(.ovf_stale) && !dir.exists(.ovf_stale)) unlink(.ovf_stale)

  for (s in STEPS) {
    sid <- s$id
    if (!is.null(only) && !sid %in% only) {
      log_info(sprintf("--- 跳过 %s (%s)：不在 --steps 里", s$label, sid))
      next
    }
    log_info("")
    log_info(sprintf("--- %s (%s)%s", s$label, sid,
                     if (isTRUE(s$required)) "" else "  [可选]"))
    # 先删本步的状态文件：崩溃时不留旧文件冒充本轮结果。
    # M10：**显式取键再判空**，不用默认空串 —— 空串会让 file.path(res_dir, "")
    # 等于结果目录本身。上面的启动断言已保证键存在，这里仍然显式判空，
    # 两层各管一件事（断言管"配置一致"，判空管"这一行不会删到目录"）。
    # **R 陷阱**：fetch 不在这张表里（豁免），`[[` 对不存在的名字抛
    # subscript out of bounds（Python 的 .get(sid) 返回 None）——
    # 而这一行不在 tryCatch 里，会崩掉整个 run_all。先 %in% 再取。
    fname <- if (sid %in% names(STEP_STATUS_FILES)) STEP_STATUS_FILES[[sid]] else NULL
    if (!is.null(fname)) {
      stale <- file.path(res_dir, fname)
      if (file.exists(stale) && !dir.exists(stale)) unlink(stale)
    }
    t0 <- Sys.time()
    err <- NULL; res <- NULL
    tryCatch(
      res <- .step_fns[[s$fn]](cfg),
      error = function(e) err <<- sprintf("%s: %s", class(e)[1L], conditionMessage(e)))
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (is.null(err)) {
      # **返回值必须接住**（E-56）。步骤函数有一条"跑完了、但结果是
      # 『没做成』"的返回路径（bad_root / insufficient_methods /
      # no_candidates / no_tfs_in_data / no_usable_pseudobulk …），
      # 它们是 write_json + log_warn + return，**不抛异常**。
      record_step(cfg, sid, "ok", seconds = secs, required = s$required,
                  result_status = result_status_of(res))
    } else {
      record_step(cfg, sid, "failed", seconds = secs, message = err,
                  required = s$required, result_status = "failed")
      if (isTRUE(s$required)) {
        log_error(sprintf("%s 失败（必需）: %s", s$label, err))
        failed_required[[length(failed_required) + 1L]] <-
          list(id = sid, label = s$label, error = err)
      } else {
        # **可选步骤失败不让 job 变红，但必须显眼。**
        log_warn(sprintf("%s 失败（可选，不影响 job 结论）: %s", s$label, err))
      }
    }
  }

  # ---- 模块零：登记输入哈希（§0.4）----------------------------------------
  # 放在所有步骤之后 —— 可选步骤的产物这轮有没有，跑完才知道。
  data_dir <- cfg$output$data_dir
  for (f in INPUT_FILES) {
    e <- record_input(cfg, file.path(data_dir, f$name),
                      label = sprintf("%s (%s)", f$desc, f$name),
                      required = f$required)
    if (identical(e$status, "missing") && isTRUE(f$required)) {
      log_warn(sprintf("输入缺失：%s (%s)", f$desc, f$name))
    }
  }
  msum <- manifest_summary(cfg)
  log_info(sprintf("输入登记 %d 项%s", msum$n_inputs,
                   if (length(msum$inputs_missing))
                     sprintf("，其中缺失 %d 项", length(msum$inputs_missing))
                   else "，全部就位"))

  # =========================================================================
  # 验收清单
  # =========================================================================
  log_info("")
  log_info(strrep("=", 68L))
  log_info("验收清单")
  log_info(strrep("=", 68L))

  fig_dir <- cfg$output$figures_dir
  checks <- list()

  # ---- 步骤 status（failed 不等于"正确地跳过"——AGENTS 规则 14）-----------
  steps_list <- read_state(cfg)$steps %||% list()
  smap <- list()
  for (s in steps_list) smap[[s$id]] <- s
  for (s in STEPS) {
    sid <- s$id
    if (!is.null(only) && !sid %in% only) next
    st <- if (sid %in% names(smap)) smap[[sid]] else list()
    # R 陷阱：list 的 [[ 对不存在的名字抛 subscript out of bounds
    # （Python 的 .get(sid, {}) 返回 {}）—— %||% 兜不了 [[ 的抛错。
    status <- st$status %||% "not_run"
    ok <- identical(status, "ok")
    # not_configured/not_applicable/disabled/not_run 是**设计如此**，
    # 可选步骤这样算通过；但 failed 是**崩了**。
    if (!isTRUE(s$required) &&
        status %in% c("not_configured", "not_applicable", "disabled", "not_run")) {
      ok <- TRUE
    }
    checks[[length(checks) + 1L]] <- list(
      item = sprintf("步骤 %s", s$label), ok = ok,
      required = isTRUE(s$required), detail = as.character(status))

    # ---- 步骤函数自己返回的 status 也要看（E-56）--------------------------
    # 光记不读等于没记 —— E-48 的形态（异常/自述失败被降级成没人读的字段）。
    rs <- st$result_status
    if (is.null(rs) || identical(tolower(as.character(rs)), "ok")) next
    verdict <- classify_step_result(rs)
    abort <- identical(verdict, "abort")
    checks[[length(checks) + 1L]] <- list(
      item = sprintf("步骤 %s 自述状态 = %s", s$label, as.character(rs)),
      ok = !abort,
      required = abort,
      detail = if (abort)
        sprintf("**步骤跑完了但结果是失败**：%s 返回 status=%s", sid, as.character(rs))
      else
        sprintf("设计如此地没做（%s 返回 status=%s），可见不阻断", sid, as.character(rs)))
  }

  for (f in REQUIRED_FILES) {
    ok <- has_file(file.path(res_dir, f$name))
    checks[[length(checks) + 1L]] <- list(
      item = sprintf("产物 %s (%s)", f$desc, f$name), ok = ok,
      required = isTRUE(f$required),
      detail = if (ok) "存在" else "**缺失或为空**")
  }

  # 图属于哪一步，由图名里的模块号决定（`02-07-01-…` 的 `07` → `07_grn.R`
  # → `grn`）。**图不能无条件要求产出** —— 可选步骤没配时它本就不该有图。
  mod2sid <- list()
  for (s in STEPS) mod2sid[[substr(s$mfile, 1L, 2L)]] <- s$id
  # **先给全部步骤填默认值 "not_run"**：`[[` 对不存在的名字会直接报
  # subscript out of bounds（Python 的 `.get(sid, {})` 有默认值，R 没有）——
  # 某个可选步骤从未运行时 state.json 里没有它的条目，验收层查
  # step_status[[owner]] 会把整个验收层崩掉，而不是判它"未运行"。
  step_status <- setNames(rep("not_run", length(STEPS)),
                          vapply(STEPS, function(s) s$id, ""))
  for (s in steps_list) step_status[[s$id]] <- s$status %||% "not_run"
  step_required <- vapply(STEPS, function(s) isTRUE(s$required), TRUE)
  names(step_required) <- vapply(STEPS, function(s) s$id, "")

  fig_names <- sort(sub("\\.png$", "", list.files(fig_dir, pattern = "\\.png$")))
  declared <- declared_figures()
  dyn <- dynamic_fig_bases()
  declared_set <- declared

  owner_of <- function(fname) {
    mod <- strsplit(fname, "-", fixed = TRUE)[[1L]]
    mod <- if (length(mod) >= 2L) mod[2L] else ""
    # [[ 陷阱：mod 不在表里时 mod2sid[[mod]] 抛错，先 %in% 再取
    owner <- if (nzchar(mod) && mod %in% names(mod2sid)) mod2sid[[mod]] else NULL
    if (is.null(owner)) return(list(owner = NA_character_, status = "not_run"))
    list(owner = owner, status = step_status[[owner]] %||% "not_run")
  }

  # ---- 静态图：声明过的每一张都要在（可选步骤没跑则豁免但可见）------------
  missing <- waived <- skipped_optional <- character(0)
  for (fname in declared) {
    ow <- owner_of(fname)
    if (has_file(file.path(fig_dir, paste0(fname, ".png")))) next
    # [[ 陷阱：绝大多数图名不在条件表里，CONDITIONAL_FIGURES[[fname]] 会抛
    # subscript out of bounds（Python 的 dict.get 返回 None）—— 先 %in% 再取
    cond <- if (fname %in% names(CONDITIONAL_FIGURES)) CONDITIONAL_FIGURES[[fname]] else NULL
    if (!is.null(cond)) {
      sp <- file.path(res_dir, cond$status_file)
      got <- NULL
      if (file.exists(sp)) got <- .dig(read_json(sp), cond$path)
      hit <- if (isTRUE(cond$ready)) identical(got, TRUE) else identical(got, cond$ready)
      if (!hit) {
        waived <- c(waived, sprintf("%s（%s；%s %s=%s）", fname, cond$why,
                                    cond$status_file,
                                    paste(cond$path, collapse = "."),
                                    if (is.null(got)) "NULL" else as.character(got)))
        next
      }
    }
    if (is.na(ow$owner) || (!isTRUE(step_required[[ow$owner]]) &&
                            !identical(ow$status, "ok"))) {
      # 可选步骤没跑（或配置关闭）→ 不要求这张图，但**必须可见**；
      # 模块号反查不到步骤（ow$owner 为 NA）同样不判 missing ——
      # step_required[[NA]] 是 subscript out of bounds，不能让验收层崩。
      skipped_optional <- c(skipped_optional,
                            sprintf("%s（步骤 %s = %s）", fname, ow$owner, ow$status))
      next
    }
    missing <- c(missing, fname)
  }

  notes <- character(0)
  if (length(waived))
    notes <- c(notes, sprintf("%d 张条件图本轮不适用：%s", length(waived),
                              paste(waived, collapse = ", ")))
  if (length(skipped_optional))
    notes <- c(notes, sprintf("%d 张属于未运行的可选步骤：%s",
                              length(skipped_optional),
                              paste(skipped_optional, collapse = ", ")))
  tail_notes <- if (length(notes)) paste(notes, collapse = "；") else ""
  checks[[length(checks) + 1L]] <- chk(
    "figures:declared", "required", length(missing) == 0L,
    if (length(missing) == 0L)
      paste0(sprintf("源码声明 %d 张静态图，全部产出", length(declared)),
             if (nzchar(tail_notes)) paste0("；", tail_notes) else "")
    else paste0("**声明了但没产出** ", paste(missing, collapse = ", "),
                if (nzchar(tail_notes)) paste0("（", tail_notes, "）") else "",
                sprintf(" —— 实际产出 %d 张", length(fig_names))))

  # ---- 动态图名：每组前缀至少 1 张 ----------------------------------------
  # **每组的 产出/槽位 都要报出来（E-70）。** 槽位是上限、少出合法，
  # 于是"3 个槽位只出 1 张"与"3 个槽位出 3 张"在这条判据眼里完全一样。
  # 判据本身不改，但把比值**摆出来**（E-69 Form B：有人读才算）。
  dyn_missing <- character(0); dyn_report <- character(0)
  for (base in names(dyn)) {
    n_slots <- dyn[[base]]
    got <- fig_names[startsWith(fig_names, paste0(base, "-")) &
                       !fig_names %in% declared_set]  # Q-27 排静态图
    dyn_report <- c(dyn_report, sprintf("%s %d/%d", base, length(got), n_slots))
    if (!length(got)) {
      dyn_missing <- c(dyn_missing,
                       sprintf("%s（声明 %d 个槽位，实际 0 张）", base, n_slots))
    }
  }
  checks[[length(checks) + 1L]] <- chk(
    "figures:dynamic", "required", length(dyn_missing) == 0L,
    if (length(dyn_missing) == 0L)
      sprintf("%d 组动态图名共 %d 个槽位，各自至少产出 1 张（产出/槽位：%s）",
              length(dyn), sum(unlist(dyn)), paste(dyn_report, collapse = ", "))
    else paste0("**动态图名整组没产出** ", paste(dyn_missing, collapse = ", "),
                " —— 槽位是上限不是精确值，少出合法，",
                "**但一张都没有说明那段循环整段没跑**",
                sprintf("（产出/槽位：%s）", paste(dyn_report, collapse = ", "))))

  # 保留计数作为**下限兜底**：声明扫描本身失效时这条还能拦住"一张图都没有"。
  checks[[length(checks) + 1L]] <- chk(
    "figures:count", "required", length(fig_names) >= 8L,
    sprintf("%d 张图（要求 >=8）", length(fig_names)))

  # 扫出来的图名必须有说明 —— 否则报告里会出现一个只有文件名、没人知道它
  # 想表达什么的条目。**这条不判红**：判红会让"加了新图"变成 CI 失败。
  no_desc <- setdiff(declared, names(.FIG_DESC))
  checks[[length(checks) + 1L]] <- chk(
    "figures:documented", "required", TRUE,
    if (!length(no_desc)) sprintf("%d 张声明图都有中文说明", length(declared))
    else sprintf("%d 张声明图缺中文说明（不影响正确性，但报告里只有文件名）：%s",
                 length(no_desc), paste(no_desc, collapse = ", ")),
    severity = "info")

  # ---- 模块零：运行清单（§0.3 / §0.4）-------------------------------------
  # 清单缺项不是"分析错了"，而是"这轮跑出来的东西没法追溯"。
  # **人工复核未确认不算失败**（默认就是 pending，设计如此）但必须可见。
  checks[[length(checks) + 1L]] <- list(
    item = "运行清单存在（run_manifest.json）",
    ok = isTRUE(msum$present), required = TRUE,
    detail = if (isTRUE(msum$present))
      sprintf("%d 个包版本、%d 项输入、%d 条决策",
              msum$n_versions %||% 0L, msum$n_inputs %||% 0L,
              msum$n_decisions %||% 0L)
    else "**缺失**")
  if (isTRUE(msum$present)) {
    # **三态，不是二态。** versions_enumeration：ok / failed / unknown。
    # failed 是本轮唯一该判红的；unknown = 清单由旧版本代码写出，
    # 读历史 artifact 不是缺陷 —— 与 E-64 同一条教训：
    # **"没跑"和"跑了没问题"必须长得不一样**，所以 unknown 只可见。
    enum <- msum$versions_enumeration %||% "unknown"
    enum_err <- msum$versions_enumeration_error
    checks[[length(checks) + 1L]] <- list(
      item = "版本枚举成功（枚举器未抛异常）",
      ok = !identical(enum, "failed"),
      required = identical(enum, "failed"),
      detail = if (identical(enum, "ok"))
        sprintf("%d 个已安装包", msum$n_versions %||% 0L)
      else if (identical(enum, "failed"))
        sprintf("**枚举失败**：%s —— versions 只有 %d 项，不足以复现本轮",
                if (is.null(enum_err)) "<NA>" else as.character(enum_err),
                msum$n_versions %||% 0L)
      else
        sprintf(paste0("**无法判断**：清单里没有 versions_enumeration 字段",
                       "（该清单由旧版本代码写出）；versions 有 %d 项。",
                       "本轮代码写出的清单会带这个字段。"),
                msum$n_versions %||% 0L))
    # 关键工具解析率**可见但不阻断**：判「至少一个解析出来」，
    # 未解析的名单完整列出。
    nkr <- msum$n_key_resolved %||% 0L
    nkt <- msum$n_key_total %||% 0L
    unres <- msum$key_unresolved %||% character(0)
    checks[[length(checks) + 1L]] <- list(
      item = "关键工具版本有记录（§0.3）",
      ok = nkt > 0L && nkr > 0L, required = TRUE,
      detail = sprintf("%d/%d 个关键工具解析出版本%s", nkr, nkt,
                       if (length(unres))
                         paste0("；未解析（查过了，未安装）: ",
                                paste(unres, collapse = ", "))
                       else ""))
    # **只看 required 的缺失。** 可选项本来就可以不存在。
    checks[[length(checks) + 1L]] <- list(
      item = "输入哈希已登记且必需项无缺失",
      ok = (msum$n_inputs %||% 0L) >= length(INPUT_FILES) &&
        length(msum$inputs_missing_required %||% character(0)) == 0L,
      required = TRUE,
      detail = paste0(
        sprintf("%d 项", msum$n_inputs %||% 0L),
        if (length(msum$inputs_missing_required %||% character(0)))
          paste0("，必需缺失 ",
                 paste(msum$inputs_missing_required, collapse = ","))
        else "，必需项齐全",
        if (length(msum$inputs_missing %||% character(0)))
          paste0("（可选缺失 ", paste(msum$inputs_missing, collapse = ","), "）")
        else ""))
    # **不能直接拿 pending 判红。** record_human_review() 默认记 pending
    # （自动化流水线不能替人签字），任何一轮跑完 pending 都非空 ——
    # 拿它判红 = 每个 job 都红 = E-29「永远红的门禁等于没有门禁」。
    # 要判的是 n_human_review == 0：说明登记循环压根没跑到（那是缺陷）。
    pend <- msum$human_review_pending %||% character(0)
    n_hr <- msum$n_human_review %||% 0L
    n_conf <- length(msum$human_review_confirmed %||% character(0))
    checks[[length(checks) + 1L]] <- list(
      item = sprintf("人工复核节点已登记（%d/%d，其中 %d 个待确认，不阻断 job）",
                     n_hr, length(HUMAN_REVIEW_NODES), length(pend)),
      ok = n_hr > 0L, required = FALSE,
      detail = if (n_hr > 0L)
        paste0(sprintf("登记 %d 个", n_hr),
               if (n_conf > 0L) sprintf("，已确认 %d 个", n_conf) else "",
               if (length(pend))
                 sprintf("，待确认 %d 个：%s", length(pend), paste(pend, collapse = ", "))
               else "（全部已确认，或本轮无人签字但节点已登记）")
      else
        # pending 为空时千万别让读者以为「人签过字了」。
        paste0("**清单里一个人工复核节点都没有** —— HUMAN_REVIEW_NODES（",
               length(HUMAN_REVIEW_NODES),
               " 个）那个登记循环没跑到（不是「全部已确认」：已确认会体现在",
               " human_review_confirmed 里）"))
    # ---- §0.2 跨语言转换：**空数组必须被解释** --------------------------
    # 空数组和"忘了记"长得一模一样。判据不是"必须有记录"，
    # 而是"空的话必须有解释"。用 read_manifest 而不是拼文件名 ——
    # 文件名只该有一处（MANIFEST_NAME 常量）。
    mfull <- read_manifest(cfg)
    cl <- msum$n_cross_language %||% 0L
    dec_nodes <- vapply(mfull$decisions %||% list(),
                        function(d) as.character(d$node %||% ""), "")
    checks[[length(checks) + 1L]] <- list(
      item = "§0.2 跨语言交接已登记或有解释",
      ok = cl > 0L || "cross_language" %in% dec_nodes,
      required = FALSE,
      detail = if (cl > 0L)
        sprintf("%d 条跨语言转换记录（交接表 + 丢失字段）", cl)
      else if ("cross_language" %in% dec_nodes)
        "0 条，**但决策链里已说明本轮没有交接**"
      else
        "**0 条且没有任何解释** —— 读者无法区分『本轮没配』和『忘了记』")
  }

  # 可选步骤的"没做"要在报告里可见 —— 但"崩了"和"没做"必须分开（M12）。
  # ok 反映真实取值（classify_step_result 只有 abort 算崩了），
  # severity="info" 保证永不判红（判红会让每个 job 都红，元规则 ④）。
  for (pair in list(
      c("pseudobulk_de", "pseudobulk_status.json"),
      c("trajectory", "trajectory_status.json"),
      c("communication", "communication_status.json"),
      c("grn", "grn_status.json"),
      c("virtual_perturbation", "virtual_perturbation_status.json"))) {
    sid <- pair[1L]; fname <- pair[2L]
    d <- read_json(file.path(res_dir, fname))
    if (is.null(d) || !is.list(d)) next
    st <- d$status
    note <- substr(as.character(d$reason %||% d$method %||% ""), 1L, 150L)
    ok <- !identical(classify_step_result(st), "abort")
    checks[[length(checks) + 1L]] <- chk(
      sprintf("step_status:%s", sid), "info", ok,
      sprintf("可选步骤状态 %s = %s；%s", sid, as.character(st), note),
      severity = "info")
  }

  # ---- 内嵌 status 扫描（Q-26 / E-48 / E-56）------------------------------
  # 判红只留给 failed/error/fail；not_applied/not_done 是设计如此地没做。
  # **E-56：顶层路径不能滤掉** —— path 为空（length 0）时 key 是 "status"，
  # 05/04/08 的失败恰恰写在顶层。R 里 length(path)==0 是显式判空，
  # 不是 truthiness（Python 的 `if p` 恒假坑照抄就复刻缺陷）。
  status_files <- sort(list.files(res_dir, pattern = "status\\.json$",
                                  full.names = TRUE))
  for (stf in status_files) {
    d <- read_json(stf)
    if (is.null(d) || !is.list(d)) next
    for (it in .iter_nested_status(d)) {
      v <- tolower(as.character(it$value))
      if (!v %in% NESTED_FAILED_VALUES) next
      key <- if (length(it$path)) paste0(paste(it$path, collapse = "."), ".status")
             else "status"
      why <- substr(as.character(it$reason %||% ""), 1L, 200L)
      checks[[length(checks) + 1L]] <- list(
        item = sprintf("%s → %s = failed", basename(stf), key),
        ok = FALSE, required = TRUE,
        detail = sprintf("**内嵌失败**：%s",
                         if (nzchar(why)) why else "(无 reason)"))
    }
  }

  # ---- 内容超出画布（Q-26 / E-49）------------------------------------------
  # save_fig 的溢出检测升级成**产物级**：写进 figure_overflow.json，
  # 这里读它并判红。（R 版 save_fig 目前不做溢出落盘 —— 本轮天然 PASS，
  # 判据保留：R 侧接入溢出检测后自动生效。）
  ovf <- read_json(file.path(res_dir, "figure_overflow.json")) %||% list()
  ovf_figs <- ovf$figures %||% list()
  if (length(ovf_figs)) {
    for (nm in names(ovf_figs)) {
      info <- ovf_figs[[nm]]
      info_txt <- if (is.list(info))
        as.character(jsonlite::toJSON(info, auto_unbox = TRUE))
      else as.character(info)
      checks[[length(checks) + 1L]] <- list(
        item = sprintf("图 %s 内容超出画布", nm),
        ok = FALSE, required = TRUE,
        detail = paste0("**会被静默裁掉**：", info_txt))
    }
  } else {
    checks[[length(checks) + 1L]] <- list(
      item = "没有图的内容超出画布", ok = TRUE, required = FALSE,
      detail = "figure_overflow.json 无记录（或无图被检测出溢出）")
  }

  # ---- 内容级检查：状态文件在 ≠ 结果是对的 --------------------------------
  tj <- read_json(file.path(res_dir, "trajectory_status.json")) %||% list()
  if (identical(tj$status, "ok")) {
    n_m <- as.integer(tj$n_methods %||% 0L)
    checks[[length(checks) + 1L]] <- list(
      item = "轨迹用了 >=2 种方法交叉验证",
      ok = n_m >= 2L, required = TRUE,
      detail = if (n_m >= 2L)
        sprintf("%d 种：%s", n_m, paste(names(tj$methods_ok), collapse = ","))
      else
        paste0(sprintf("**只有 %d 种**", n_m),
               " —— 单一方法的拟时序是某个算法的一次输出，不是数据里的结构"))
    cv <- tj$cross_validated_methods %||% list()
    # **E-69 Form A 的消费端。** 产出端用 finite_round 落 null（R 里是 NULL），
    # 这里读到 NULL —— **「算不出来」不是「一致性不低」**。
    # mc_state（single_method / undefined）与原因文案必须一起看。
    mean_rho <- tj$method_correlation_mean_offdiag
    mc_state <- tj$method_correlation_state
    mc_note <- tj$method_correlation_undefined_note
    mc_ok <- length(cv) > 0L && !is.null(mean_rho)
    checks[[length(checks) + 1L]] <- list(
      item = "轨迹方法间一致性已量化（且排除方向参考）",
      ok = mc_ok, required = TRUE,
      detail = if (mc_ok)
        paste0(sprintf("交叉验证 %d 种，平均 rho=%+.4f，参考方法 %s",
                       length(cv), as.numeric(mean_rho),
                       as.character(tj$direction_reference_method %||% "")),
               if (!is.null(mc_note)) sprintf("（%s）", as.character(mc_note)) else "")
      else
        paste0(sprintf("**方法间一致性算不出来**（state=%s）：%s",
                       as.character(mc_state %||% ""),
                       if (!is.null(mc_note)) as.character(mc_note)
                       else "缺 cross_validated_methods 或一致性数值"),
               " —— 「算不出来」不是「一致性不低」"))
    ds <- tj$direction_source
    checks[[length(checks) + 1L]] <- list(
      item = "轨迹方向来源已写明",
      ok = length(ds) > 0L && nzchar(as.character(ds)), required = TRUE,
      detail = if (length(ds) > 0L && nzchar(as.character(ds)))
        as.character(ds) else "**缺失**")
    sv <- tj$scvelo %||% list()
    checks[[length(checks) + 1L]] <- list(
      item = "scVelo 的不可得已如实记录",
      ok = length(sv$status %||% NULL) > 0L, required = FALSE,
      detail = as.character(sv$status %||% "**缺失**"))
    lim_n <- length(tj$limitations %||% character(0))
    checks[[length(checks) + 1L]] <- list(
      item = "轨迹方法学限定已写明（>=4 条）",
      ok = lim_n >= 4L, required = TRUE,
      detail = sprintf("%d 条", lim_n))
    # **02-05-03 已按单图原则拆成两张**：unit1 = 热图、unit2 = 模块曲线。
    # 改名时**必须同步这里** —— 否则验收会查一个不存在的文件而静默变 false。
    for (fg in list(
        c("02-05-02-unit1-trajectory-method-correlation", "方法一致性矩阵"),
        c("02-05-03-unit1-trajectory-modules-heatmap", "沿轨迹基因模块热图"),
        c("02-05-03-unit2-trajectory-module-profiles", "各模块拟时序曲线"))) {
      ok <- has_file(file.path(fig_dir, paste0(fg[1L], ".png")))
      checks[[length(checks) + 1L]] <- list(
        item = sprintf("图 %s (%s.png)", fg[2L], fg[1L]), ok = ok,
        required = TRUE, detail = if (ok) "存在" else "**缺失**")
    }
  }

  # ---- 文档指定工具的落地情况（§2.4 CellTypist / §2.7 LIANA）-------------
  # **不阻断 job，但没跑成就要红字显示。** 规则 4：可选步骤的"没做"
  # 必须和"做了没问题"长得不一样。
  cj <- read_json(file.path(res_dir, "cluster_status.json")) %||% list()
  ann <- cj$annotation %||% list()
  ct <- ann$celltypist %||% list()
  ct_cmp <- ann$celltypist_vs_marker %||% list()

  # ---- §0.2 交接的**产出侧**：Part 3 拿 clustered 参考当输入 --------------
  # 交接的两端都要能被核对 —— 只检查消费侧的话，产出侧悄悄丢掉 counts 层
  # 不会被任何人发现。（R 适配：契约键名是 counts_layer，
  # Python 版是 layers['counts']。）
  p3 <- cj$part3_reference %||% list()
  p3c <- p3$contract %||% list()
  p3_missing <- character(0)
  if (!isTRUE(p3c$counts_layer)) p3_missing <- c(p3_missing, "counts 层")
  if (is.null(p3c$celltype_column)) p3_missing <- c(p3_missing, "obs 细胞类型列")
  checks[[length(checks) + 1L]] <- list(
    item = "§0.2 Part 3 参考导出的契约（counts 层 + 细胞类型列）",
    ok = length(p3) > 0L && length(p3_missing) == 0L, required = FALSE,
    detail = if (length(p3) > 0L && length(p3_missing) == 0L)
      sprintf("契约完整：counts 层 + `%s` 列，%s 细胞 x %s 基因",
              as.character(p3c$celltype_column),
              as.character(p3$n_cells %||% ""), as.character(p3$n_genes %||% ""))
    else if (length(p3) > 0L)
      sprintf("**契约缺**：%s —— Part 3 拿这份文件当参考时会静默算错或 KeyError",
              paste(p3_missing, collapse = ", "))
    else
      "**没有 part3_reference 记录** —— cluster_status.json 是上一轮的旧文件，或这一步没跑完")
  checks[[length(checks) + 1L]] <- list(
    item = "CellTypist 自动注释（§2.4）",
    ok = identical(ct$status, "ok"), required = FALSE,
    detail = if (identical(ct$status, "ok"))
      sprintf("ok：%s 种标签，模型 %s，%s 基因输入",
              as.character(ct$n_labels %||% ""), as.character(ct$model %||% ""),
              as.character(ct$n_genes_input %||% ""))
    else
      sprintf("**未跑成**（%s：%s）—— marker 打分仍在，但少了第二条独立证据",
              as.character(ct$status %||% ""),
              substr(as.character(ct$reason %||% ""), 1L, 120L)))
  if (identical(ct$status, "ok")) {
    # **"恒为 0 的一致率"必须显示成红的**（审计 S3 / E-58）：
    # 有对比但一个簇都没映射上 → 红；有映射 → 报真的一致率。
    cmp_ok <- isTRUE(ct_cmp$compared)
    n_map <- ct_cmp$n_mapped
    no_map <- cmp_ok && identical(as.integer(n_map %||% 0L), 0L) &&
      (as.integer(ct_cmp$n_clusters %||% 0L) > 0L)
    checks[[length(checks) + 1L]] <- list(
      item = "CellTypist 与 marker 注释的一致性已量化（§2.4）",
      ok = cmp_ok && !no_map, required = FALSE,
      detail = if (cmp_ok && !no_map)
        sprintf("簇层面一致 %s/%s（%s）；%s 个簇的词表无映射，不计入分母",
                as.character(ct_cmp$n_agree_mapped %||% ""),
                as.character(n_map %||% ""),
                as.character(ct_cmp$agreement_frac_mapped %||% ""),
                as.character(ct_cmp$n_unmapped %||% ""))
      else if (no_map)
        paste0(sprintf("**%s 个簇一个都没映射上** —— `%s` 缺条目或词表已变；",
                       as.character(ct_cmp$n_clusters %||% ""),
                       as.character(ct_cmp$mapping_source %||% "")),
               "此时字符串全等恒为 0，**不能读成『两条路不一致』**")
      else
        sprintf("**未对比**：%s", as.character(ct_cmp$reason %||% "")))
  }

  cm <- read_json(file.path(res_dir, "communication_status.json")) %||% list()
  li <- cm$liana %||% list()
  li_cmp <- cm$liana_vs_builtin %||% list()
  checks[[length(checks) + 1L]] <- list(
    item = "LIANA rank_aggregate（§2.7 指定的主工具）",
    ok = identical(li$status, "ok"), required = FALSE,
    detail = if (identical(li$status, "ok"))
      sprintf("ok：%s 行，v%s", as.character(li$n_rows %||% ""),
              as.character(li$version %||% ""))
    else
      sprintf("**未跑成**（%s：%s）—— 自建共表达打分仍在，但它不是 consensus rank aggregate",
              as.character(li$status %||% ""),
              substr(as.character(li$reason %||% ""), 1L, 120L)))
  if (identical(li$status, "ok")) {
    checks[[length(checks) + 1L]] <- list(
      item = "LIANA 与自建打分的差异已量化（§2.7）",
      ok = isTRUE(li_cmp$compared), required = FALSE,
      detail = if (isTRUE(li_cmp$compared))
        sprintf("共同组合 %s，Spearman rho=%s，top25 重叠 %s/25",
                as.character(li_cmp$n_common_combinations %||% ""),
                as.character(li_cmp$spearman_rho %||% ""),
                as.character(li_cmp$top25_overlap %||% ""))
      else
        sprintf("**未对比**：%s", as.character(li_cmp$reason %||% "")))
  }

  # ---- §1.7 / §1.8 虚拟扰动（保留框架）------------------------------------
  # 这一节的**重点是"没做什么"**。**R 版适配**：R 版状态里没有 tools 段
  # （Python 的 vp_tools 探针表）—— R 版探针（tools_probe）是局部的，
  # tenifold 脚本缺失的原因走 method_bits，所以本段没有 tools 一条；
  # 引擎自洽与失败原因见 _vp_engine_consistent。
  vp <- read_json(file.path(res_dir, "virtual_perturbation_status.json")) %||% list()
  if (identical(vp$status, "ok")) {
    ts <- vp$target_source
    checks[[length(checks) + 1L]] <- list(
      item = "虚拟扰动的候选来源已写明（Part 1 交接 vs 内部回退）",
      ok = length(ts) > 0L && nzchar(as.character(ts)), required = FALSE,
      detail = paste0(
        sprintf("%s，%s 个目标（celltype_key=%s）", as.character(ts %||% ""),
                as.character(vp$n_targets %||% ""),
                as.character(vp$celltype_key %||% "")),
        # 内部回退时 signature_alignment 为空是预期的
        if (length(ts) > 0L && startsWith(as.character(ts), "internal"))
          "（**内部回退：没有 Part 1 签名，signature_alignment 为空是预期的**）"
        else ""))
    lim_n <- length(vp$limitations %||% character(0))
    checks[[length(checks) + 1L]] <- list(
      item = "虚拟扰动的方法学限定已写明（>=5 条）",
      ok = lim_n >= 5L, required = FALSE, detail = sprintf("%d 条", lim_n))
    checks[[length(checks) + 1L]] <- list(
      item = "虚拟扰动声明的方法与状态里的引擎一致",
      # 防止"报了个数就被当成因果预测"：状态里记了哪些引擎，
      # 方法自述就必须点名哪些，且必须仍然否掉没跑的 PerturbNet（K-01b）。
      ok = .vp_engine_consistent(vp), required = FALSE,
      detail = sprintf("engines_used=%s；%s",
                       paste(as.character(vp$engines_used %||% ""), collapse = "+"),
                       substr(paste(as.character(vp$method_bits %||% ""),
                                    collapse = "；"), 1L, 100L)))
    # Tenifold 的产出或失败原因已记录（失败也必须留下原因 ——
    # 否则"没跑"和"跑了没结果"长得一样）。R 版载体是 tenifold_ran +
    # method_bits 里的 "Tenifold 引擎未产出：…"。
    td_ran <- isTRUE(vp$tenifold_ran)
    td_reasoned <- any(grepl("Tenifold 引擎未产出",
                             as.character(vp$method_bits %||% character(0)),
                             fixed = TRUE))
    checks[[length(checks) + 1L]] <- list(
      item = "scTenifoldKnk（R 引擎）的产出或失败原因已记录",
      ok = td_ran || td_reasoned, required = FALSE,
      detail = if (td_ran) sprintf("tenifold_ran=TRUE，图 %s",
                                   if (isTRUE(vp$tenifold_figure_drawn)) "已画" else "退回一阶")
      else if (td_reasoned) "tenifold 未产出（原因见 method_bits）"
      else "**两者都没有** —— tenifold_ran=FALSE 且 method_bits 无失败原因")
    # 空敲除（出度为 0 的候选基因）必须被点名。包的敲除方式是"把网络里
    # 该基因那一行清零"；出度为 0 时清一行全 0 等于没敲，"距离"只剩浮点
    # 噪声（实测 ~1e-16）—— 不报错、不给 NA，读表人只会得出"敲除无影响"。
    # **R 版适配**：状态键是 n_empty_knockouts（Python 是 n_empty_knockout）；
    # 基因名单不在状态里，从产物 virtual_perturbation_tenifold.csv 读
    # （empty_knockout 列是结构证据，不从整行 NA 反推）。
    if (td_ran) {
      n_empty <- as.integer(vp$n_empty_knockouts %||% 0L)
      empty_genes <- character(0)
      if (n_empty > 0L) {
        td_csv <- file.path(res_dir, "virtual_perturbation_tenifold.csv")
        if (file.exists(td_csv)) {
          tdd <- tryCatch(utils::read.csv(td_csv, stringsAsFactors = FALSE),
                          error = function(e) NULL)
          if (!is.null(tdd) && all(c("gene", "empty_knockout") %in% names(tdd))) {
            empty_genes <- as.character(tdd$gene[tdd$empty_knockout %in% TRUE])
          }
        }
      }
      checks[[length(checks) + 1L]] <- list(
        item = "scTenifoldKnk 的空敲除（网络里出度为 0 的候选）已点名",
        ok = (n_empty == 0L) || (length(empty_genes) == n_empty),
        required = FALSE,
        detail = if (n_empty == 0L) "没有空敲除"
        else sprintf("%d 个：%s（**不是效应为 0，是该网络表达不了该扰动**）",
                     n_empty, paste(empty_genes, collapse = "、")))
    }
    vp_fig <- "02-08-01-unit1-virtual-perturbation-effect.png"
    ok_fig <- has_file(file.path(fig_dir, vp_fig))
    checks[[length(checks) + 1L]] <- list(
      item = sprintf("图 虚拟敲除效应 (%s)", vp_fig), ok = ok_fig,
      required = FALSE, detail = if (ok_fig) "存在" else "**缺失**")
  }

  # ---- §2 点名工具的缺口登记 -------------------------------------------
  # **判据是"理由写了没有"，不是"工具跑了没有"。** 将来某个工具能装了，
  # 这条应该依然 PASS —— 如实记录不该被惩罚成失败。
  # **变量名是 mfull 不是 msum**：两个结构键完全不同，同名会让下一次改动
  # 时 .get 静默通过（Python 侧 L1186-1194 的教训照抄）。
  mfull <- read_manifest(cfg)
  named <- (mfull$params %||% list())$named_tools
  if (!is.list(named) || length(named) == 0L) {
    checks[[length(checks) + 1L]] <- list(
      item = "§2 点名工具的缺口已登记（清单 named_tools）",
      ok = FALSE, required = FALSE,
      detail = "清单里没有 named_tools —— 读者会以为 §2 点名的方法都用上了")
  } else {
    no_reason <- names(named)[vapply(named, function(i)
      !nzchar(as.character((i %||% list())$reason %||% "")), TRUE)]
    no_kind <- names(named)[vapply(named, function(i)
      !nzchar(as.character((i %||% list())$kind %||% "")), TRUE)]
    checks[[length(checks) + 1L]] <- list(
      item = "§2 点名工具的缺口已登记（清单 named_tools）",
      ok = length(no_reason) == 0L && length(no_kind) == 0L,
      required = FALSE,
      detail = if (length(no_reason) == 0L && length(no_kind) == 0L)
        sprintf("%d 个点名工具已登记，%d 个当前可用", length(named),
                sum(vapply(named, function(i) isTRUE(i$available), TRUE)))
      else sprintf("无理由 %s；无分类 %s",
                   paste(no_reason, collapse = ", "),
                   paste(no_kind, collapse = ", ")))
    # R 版没有 name_taken kind（PyPI 同名无关包在 R 版不适用，见登记表注释），
    # Python 版的 squat 判据在此**不适用而非省略** —— r_version.md §7
    # 的 R_CHK_DIFF 需要标记这条差异。
    # 决策链也要留痕（§0.4）
    dec <- mfull$decisions %||% list()
    dec_has <- any(vapply(dec, function(d) identical(d$node, "named_tools"), TRUE))
    checks[[length(checks) + 1L]] <- list(
      item = "点名工具的使用情况进了决策链（§0.4）",
      ok = dec_has, required = FALSE,
      detail = sprintf("decisions 里 %d 条，named_tools %s", length(dec),
                       if (dec_has) "在" else "**不在**"))
  }

  # ---- 汇总与三态日志 ------------------------------------------------------
  n_fail <- 0L
  for (c in checks) {
    mark <- if (isTRUE(c$ok)) "PASS" else if (isTRUE(c$required)) "FAIL" else "INFO"
    if (!isTRUE(c$ok) && isTRUE(c$required)) n_fail <- n_fail + 1L
    if (isTRUE(c$ok) && isTRUE(c$required)) {
      log_info(sprintf("  [%s] %s", mark, c$item))
    } else if (!isTRUE(c$ok)) {
      log_error(sprintf("  [%s] %s  %s", mark, c$item, c$detail))
    } else {
      log_info(sprintf("  [%s] %s  %s", mark, c$item, c$detail))
    }
  }

  summary <- list(
    dataset_id = cfg$dataset_id,
    n_checks = length(checks),
    n_failed_required = n_fail,
    n_info = sum(vapply(checks, function(c) identical(c$severity, "info"), TRUE)),
    checks = checks,
    steps_failed_required = failed_required,
    manifest = mfull,
    verdict = if (n_fail == 0L) "ok" else "failed"
  )
  write_json(file.path(res_dir, "acceptance.json"), summary)

  log_info("")
  if (n_fail == 0L) {
    log_info(sprintf("验收通过：%d 项检查，0 项必需失败", length(checks)))
  } else {
    log_error(sprintf("验收失败：%d 项必需检查未通过", n_fail))
  }
  if (n_fail == 0L) 0L else 1L
}

# ---------------------------------------------------------------------------
# 入口（--config / --steps 逗号拆分；守卫与 9 个步脚本同款）
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  args <- parse_args()
  cfg <- load_config(args$config)
  only <- if (is.null(args$steps) || !nzchar(args$steps)) NULL
          else trimws(strsplit(args$steps, ",", fixed = TRUE)[[1L]])
  quit(save = "no", status = run_all(cfg, only))
}
