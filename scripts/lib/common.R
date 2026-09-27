# lib/common.R — 单细胞流水线公共库（R 版，R-06，2026-09-26）
#
# 与 scripts/lib/common.py 同一套心智模型，函数一一对应（对照表见下）。
# 立项文档：references/r_version.md（R-05）。
#
#   - `parse_args()` **故意没有默认配置** —— 多数据集下静默默认到其中某一个，
#     正是"跑错数据集"的来源。
#   - 日志走 log_info/log_warn/log_error，不用裸 cat/printf。
#   - 路径一律从 cfg 派生（dataset_id -> results/<id>、data/<id>），不硬编码。
#   - 步骤失败必须 stop() 让编排器记录，不要 tryCatch 后静默继续。
#   - 可选步骤的失败要在自己的状态 JSON 里留 reason。
#
# ## 与 Python 版的三处**语言级**差异（立项文档 §5 / §11）
#
# 1. **"没有这个值"在 R 里是 `NA`，不是 `NULL`。** jsonlite 把 `NULL`
#    序列化成 `{}`（空对象），把 `NA` 在 `na="null"` 下序列化成 `null` ——
#    Python 版写的是 `null`，所以 R 版一律用 `NA` 表示"算不出来/没查到"，
#    **用 `NULL` 表示"这个键不存在"**（会被序列化丢掉）。两者混用会产出
#    结构不同但都合法的 JSON，下游 `is not None` / `is.na` 判据就全错了。
# 2. **中间对象是 `.rds`，不是 `.h5ad`**（立项 §5：不走 basilisk/conda，
#    跨版本对照只走 CSV）。
# 3. **没有 matplotlib 的 bbox_inches 概念。** ggsave 写出的就是给定尺寸，
#    Python 版 `save_fig(tight=)` 的语义在 R 版不存在；"内容被裁掉"的兜底
#    （`figure_overflow.json`）也没有对应物 —— ggplot 的图如果放不下会
#    **挤在一起**而不是裁掉，这类问题靠 `check_figures.mjs` 的墨迹/边距判据
#    兜（不是靠写盘记录）。
#
# ## R1 纪律
# 本文件是**分析公共库**，本地不执行（R1：分析只在 GitHub Actions 上跑）。
# 本地允许的操作只有静态检查与 `--selftest`（CI 已验证的模式，见
# `scripts/lib/tenifold_knk.R`）—— selftest 只做 JSON/配置/状态文件的往返，
# 不做任何分析计算。

# ============================================================================
# 确定性 —— **必须在载入任何重型包之前执行**
# ============================================================================
# 与 common.py 的 env 段同一组变量：浮点末位分叉有两个独立来源 ——
# 线程调度（多线程归约的求和顺序）与 OpenBLAS 运行期按 CPU 型号选内核。
# R 走的 OpenBLAS 同样读这些环境变量；R 自己的 RNG 用 set_seed() 钉。
# 规则 12 的口径：可声称的是**结构与量级可复现**，不是逐字节一致。
for (v in c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
            "NUMBA_NUM_THREADS")) {
  if (is.na(Sys.getenv(v, unset = NA))) Sys.setenv(v = "1")
}
# 规则 12 的第二个来源：OpenBLAS 运行期按 CPU 型号选 SIMD 内核。
# Python 版 common.py:53 同款；R 在 Windows 上链接的也是 OpenBLAS（R 4.2+）。
if (is.na(Sys.getenv("OPENBLAS_CORETYPE", unset = NA))) {
  Sys.setenv(OPENBLAS_CORETYPE = "Haswell")
}
if (is.na(Sys.getenv("PYTHONHASHSEED", unset = NA))) Sys.setenv(PYTHONHASHSEED = "0")

# R 自身的 RNG 也要钉：默认 Mersenne-Twister + Inversion 是稳定的，
# 但 `sample()` 的特殊分支（sample(1:n) vs sample(x)）随版本变过 ——
# 这里把 kind 显式写死，不赌 R 版本默认值。
# **第三个参数（sample_kind）只接受 "Rejection"/"Rounding"** ——
# 曾写过的 "Reproducible" 不是合法 choice，R 4.6（CI 实跑 run 36308546495）
# 直接 stop("'NA' is not a valid choice")：R 把非法值替换成 NA 再校验，
# 报错文本里看到的是 NA 不是原值，grep 原值是搜不到的。9 个步骤每步都调
# set_seed → 全部以同一错误失败。
options(stringsAsFactors = FALSE)
ScrnaRNGKind <- function() {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
}

# ============================================================================
# 日志
# ============================================================================
.scrna_orchestrated <- function() isTRUE(getOption("scrna.orchestrated", FALSE))

.set_orchestrated <- function(flag = TRUE) options(scrna.orchestrated = isTRUE(flag))

.log <- function(level, msg) {
  ts <- format(Sys.time(), "%H:%M:%S")
  cat(sprintf("[%s] %-7s %s\n", ts, level, msg), file = if (level == "ERROR") stderr() else stdout())
}
log_info <- function(...) .log("INFO", paste0(...))
log_warn <- function(...) .log("WARN", paste0(...))
log_error <- function(...) .log("ERROR", paste0(...))

# ============================================================================
# 小工具
# ============================================================================
# Python 的 `a or b`（缺省替换）。**只在"键不存在/值为 NULL"时给缺省**，
# 不要拿它替换 NA —— NA 是"有这个键，值算不出来"，替换它会掩盖信息。
`%||%` <- function(a, b) if (is.null(a)) b else a

# 硬依赖版（缺了就 stop）。geo 仓同名 `require_pkg` 是**可选包**版
# （缺了 log_warn + 返回 FALSE）—— 两个语义都存在，R 版两处分开起名：
# 硬依赖用 `need_pkg`（本文件），可选包将来接 geo 的 `require_pkg` 语义。
# 好奇它为什么分开：把"缺硬依赖"降级成告警，步骤会在没有地基的情况下
# 继续跑，产出一个"看起来正常"的空壳 —— 比 crash 糟得多（E-48 同族）。
need_pkg <- function(pkg, why = "") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("缺少 R 包 %s%s —— 这是 R 版流水线的硬依赖，装法见 scripts/r_deps.R",
                 pkg, if (nzchar(why)) paste0("（", why, "）") else ""), call. = FALSE)
  }
}

# ============================================================================
# 参数与配置
# ============================================================================
# 用法: Rscript scripts/main_analysis.R --config assets/config.pbmc3k.yml
#       [--steps fetch,qc,...]（调试用，只跑指定步骤）
# **没有默认配置** —— 与 common.parse_args 同一理由。
parse_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  get_val <- function(flag) {
    i <- match(flag, argv)
    if (is.na(i)) return(NULL)
    if (i >= length(argv)) stop(sprintf("%s 缺少取值", flag), call. = FALSE)
    argv[i + 1]
  }
  unknown <- setdiff(grep("^--", argv, value = TRUE),
                     c("--config", "--steps", "--selftest", "--help", "-h"))
  if (length(unknown)) {
    stop(sprintf("未知参数: %s（认识: --config --steps --selftest --help）",
                 paste(unknown, collapse = ", ")), call. = FALSE)
  }
  list(config = get_val("--config"), steps = get_val("--steps"),
       selftest = "--selftest" %in% argv)
}

load_config <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    stop(sprintf("配置文件不存在: %s", path %||% "<NULL>"), call. = FALSE)
  }
  need_pkg("yaml", "读 assets/config.*.yml")
  cfg <- yaml::read_yaml(path)
  if (!is.list(cfg) || is.null(cfg$dataset_id) || !nzchar(as.character(cfg$dataset_id))) {
    stop(sprintf("配置缺少 dataset_id: %s", path), call. = FALSE)
  }
  did <- as.character(cfg$dataset_id)
  # 目录由 dataset_id 派生 —— 一个数据集一个产物目录，互不覆盖。
  cfg$output$results_dir <- cfg$output$results_dir %||% file.path("results", did)
  cfg$output$data_dir    <- cfg$output$data_dir    %||% file.path("data", did)
  cfg$output$figures_dir <- cfg$output$figures_dir %||% file.path("results", did, "figures")
  # 300 dpi 是投稿图的底线（150 dpi 下 183 mm 宽只有 1080 px，印刷发虚）。
  # setdefault 只在键缺失时生效 —— config.pbmc3k.yml 里显式写了 figure_dpi，
  # 所以这里只是兜底；两处都要 300（common.py:127-140 的实测教训）。
  cfg$analysis$seed       <- cfg$analysis$seed %||% 20260919
  cfg$analysis$figure_dpi <- cfg$analysis$figure_dpi %||% 300
  cfg
}

ensure_dirs <- function(cfg) {
  for (k in c("results_dir", "data_dir", "figures_dir")) {
    dir.create(cfg$output[[k]], recursive = TRUE, showWarnings = FALSE)
  }
  invisible(NULL)
}

set_seed <- function(cfg) {
  seed <- as.integer(cfg$analysis$seed)
  ScrnaRNGKind()
  set.seed(seed)
  seed
}

# ============================================================================
# 非有限值 —— E-69 Form A 的修复必须先于一切数值状态字段
# ============================================================================
# Python 版 `finite_round`（common.py:161-180）同名同语义。**裸 NaN 参与
# 任何比较都是 FALSE 且不报错** —— `if (mean_rho < 0.3)` 对"算不出来"
# 静默判假，"算不出来"被当成"算出来不低"。所以「算不出来」应当在
# **产出它的地方**就标成 NA（写盘时 jsonlite 再换成 null 是最后一道兜底，
# 内存里的消费者读到的已经是 NA）。调用点判空一律 `is.na()` / `!is.finite()`，
# **不能用 `is.null()` 或真值判断**（`0` 是有效值，`FALSE` 也是）。
finite_round <- function(x, ndigits = 4L) {
  if (is.null(x) || length(x) == 0L || is.na(x)) return(NA_real_)
  x <- suppressWarnings(as.numeric(x))
  if (length(x) != 1L || is.na(x) || !is.finite(x)) return(NA_real_)
  round(x, ndigits)
}

# ============================================================================
# JSON —— 合法 JSON 是硬约束（M9）
# ============================================================================
# Python 版 `_scrub_nonfinite`：把 NaN/Infinity 递归换成 None 并返回命中数。
# R 版换成 NA（jsonlite 在 na="null" 下把 NA 写成 null）。**不删键** ——
# "这个量算不出来"本身是信息，删掉它读者只会以为没这个字段。
.scrub_nonfinite <- function(obj) {
  hits <- 0L
  walk <- function(x) {
    if (is.list(x)) {
      out <- lapply(x, walk)
      names(out) <- names(x)
      out
    } else if (is.numeric(x)) {
      bad <- !is.finite(x) & !is.na(x)
      if (any(bad)) hits <<- hits + sum(bad)
      x[bad] <- NA_real_
      x
    } else x
  }
  out <- walk(obj)
  list(obj = out, hits = hits)
}

# Python 的 `default=` 等价物：把 R 的类型归一成 jsonlite 能写的东西。
# data.frame -> records（逐行 list），Matrix/矩阵 -> 不支持（先落 CSV）。
.json_default <- function(obj) {
  if (inherits(obj, "data.frame")) {
    if (nrow(obj) == 0L) return(list())
    return(lapply(seq_len(nrow(obj)), function(i) as.list(obj[i, , drop = FALSE])))
  }
  as.character(obj)
}

write_json <- function(path, obj) {
  need_pkg("jsonlite")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  scrubbed <- .scrub_nonfinite(obj)
  if (scrubbed$hits > 0L) {
    log_warn(sprintf("write_json: %s 有 %d 个 NaN/Infinity，已写成 null（非法 JSON 会让 node/jq/R 读者解析失败）",
                     path, scrubbed$hits))
  }
  # **原子写**（geo 仓 common.R:143-148 同款）：先写 .tmp 再改名。
  # 半截的 state.json / *_status.json 不是"少一条数据"，而是让对应验收项
  # **静默消失**（read_json 对坏内容抛错、验收 `if not d: continue` 跳过）
  # —— M8 的教训，写的时候防比读的时候炸便宜。
  tmp <- paste0(path, ".tmp")
  jsonlite::write_json(scrubbed$obj, tmp, auto_unbox = TRUE, pretty = TRUE,
                       na = "null", digits = NA, default = .json_default)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) {
    stop(sprintf("write_json: 原子改名失败 %s -> %s", tmp, path), call. = FALSE)
  }
  invisible(path)
}

read_json <- function(path) {
  # 不存在 -> NULL；**内容坏了 -> stop()**（M8：坏文件让验收项静默消失，
  # 比报错糟得多）。判据产物（*_status.json）一律走这条；需要"坏也当没有"
  # 的调用点显式用 read_json_or_none()，让意图在调用处可见。
  need_pkg("jsonlite")
  if (is.null(path) || !file.exists(path)) return(NULL)
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

read_json_or_none <- function(path) {
  tryCatch(read_json(path),
           error = function(e) {
             log_warn(sprintf("read_json_or_none: %s 存在但解析失败（%s）—— 当作不存在处理",
                              path, conditionMessage(e)))
             NULL
           })
}

sha256_file <- function(path) {
  need_pkg("digest", "输入哈希")
  digest::digest(file = path, algo = "sha256")
}

# ============================================================================
# 步骤状态（state.json）
# ============================================================================
state_path <- function(cfg) file.path(cfg$output$results_dir, "state.json")
read_state <- function(cfg) read_json(state_path(cfg)) %||% list(steps = list())

# 步骤函数自述状态的翻译（E-56；与 common.py:308-345 同名同值）。
# 这两个集合**只管步骤函数顶层返回值**；嵌套字段里的同一批词不能照搬
# （tenifold 的 status=timeout 是可选能力没成，不是这一步崩了）。
STEP_ABORT_VALUES <- c(
  # 通用崩溃
  "failed", "error", "fail",
  # 05 轨迹：根簇不在簇列表 / 没产出 leiden / 成功方法不足 2 种 / 方向全不可判
  "bad_root", "missing_clusters", "insufficient_methods", "no_consensus",
  # 08 虚拟敲除
  "no_candidates",
  # 07 GRN
  "no_tfs_in_data", "no_regulons",
  # 06 通讯
  "no_pairs_in_data", "no_signal",
  # 04 拟bulk（missing_counts / mixed_group 是**契约被破坏**，要判红）
  "no_usable_pseudobulk", "no_celltype_testable", "column_missing",
  "missing_counts", "mixed_group"
)
# 设计如此地没做（配置关了 / 输入不支持 / 环境缺包 / **部分成功**）—— 可见不阻断。
STEP_SKIP_VALUES <- c(
  "disabled", "not_configured", "not_applicable", "not_run", "not_applied",
  "not_done", "not_available", "not_possible", "heuristic_only",
  "package_missing", "unavailable", "needs_reference", "partial"
)

result_status_of <- function(res) {
  # 步骤函数**返回值**里自述的 status（E-56）。不是 dict（list）或没有
  # status 键 => "ok"（正常跑完、没有自述失败）。
  if (is.list(res) && !is.null(res$status)) return(as.character(res$status))
  "ok"
}

classify_step_result <- function(result_status) {
  # ok / abort / skip。**未登记的取值归 skip 不判红** —— 一个没见过的词
  # 有可能是新加的"设计如此"，判红会让 job 变红、让人开始忽略告警。
  v <- tolower(trimws(as.character(result_status %||% "")))
  if (!nzchar(v) || v == "ok") return("ok")
  if (v %in% STEP_ABORT_VALUES) return("abort")
  "skip"
}

record_step <- function(cfg, id, status, seconds = NA_real_, message = "",
                        required = TRUE, result_status = NULL) {
  st <- read_state(cfg)
  st$steps <- Filter(function(s) !identical(s$id, id), st$steps %||% list())
  entry <- list(id = id, status = status, required = isTRUE(required),
                message = message)
  if (!is.null(result_status)) {
    # 与 `status`（编排器记的"有没有崩"）是**两件事**：status="ok" +
    # result_status="bad_root" 完全可能，那正是 E-56 要让它显形的组合。
    entry$result_status <- as.character(result_status)
  }
  if (!is.na(seconds)) entry$seconds <- round(as.numeric(seconds), 1)
  st$steps <- c(st$steps, list(entry))
  st$dataset_id <- cfg$dataset_id
  st$updated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  write_json(state_path(cfg), st)
  invisible(NULL)
}

# ============================================================================
# 运行清单（run_manifest.json）
# ============================================================================
# state.json 记"这一步跑没跑成"，每步重写；manifest 记"本轮是在什么条件下
# 跑出来的"，是证据，写入后不该再变。混在一起会让后者被前者覆盖。
MANIFEST_NAME <- "run_manifest.json"

# R 版的关键包清单（Python 版 KEY_PACKAGES 是它自己的依赖清单 ——
# 两版各自记录**自己实际用的**环境，键名不要求一致，但
# n_key_resolved / n_key_total / versions_enumeration 这组结构键一致）。
# **2026-09-26 更正**：`scFates` 从 CRAN/Bioconductor 均不可得（R-07 实测
# CRAN 404，Python scFates 从未发布 R 版）——R 版 05 用 slingshot 替代
#（r_version.md §11 未决 1 的裁决）。清单里删掉 scFates、补 slingshot。
KEY_PACKAGES <- c(
  # 本仓 R 版按立项文档 §3 选定的方法栈
  "Seurat", "harmony", "sva", "scDblFinder", "scater", "SingleR", "celldex",
  "DESeq2", "destiny", "slingshot", "bluster", "cowplot", "liana", "scTenifoldKnk",
  # 基础设施
  "yaml", "jsonlite", "digest", "ggplot2", "Matrix", "patchwork"
)

manifest_path <- function(cfg) file.path(cfg$output$results_dir, MANIFEST_NAME)
read_manifest <- function(cfg) read_json(manifest_path(cfg)) %||% list()

init_manifest <- function(cfg, language = "r") {
  # 建立本轮清单骨架。**会清掉上一轮的内容** —— 清单描述的是本轮。
  m <- list(
    dataset_id = cfg$dataset_id, language = language,
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    seed = cfg$analysis$seed,
    versions = list(), key_versions = list(),
    inputs = list(), params = list(),
    decisions = list(), human_review = list(), cross_language = list()
  )
  write_json(manifest_path(cfg), m)
  m
}

.manifest_append <- function(cfg, key, entry) {
  m <- read_manifest(cfg)
  m[[key]] <- c(m[[key]] %||% list(), list(entry))
  m$dataset_id <- cfg$dataset_id
  m$updated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  write_json(manifest_path(cfg), m)
}

record_input <- function(cfg, path, label = "", required = TRUE) {
  # 文件不存在时**记 missing 而不是抛异常** —— 可选步骤的产物这轮未必有。
  entry <- list(label = if (nzchar(label)) label else basename(path),
                path = as.character(path), required = isTRUE(required))
  if (file.exists(path) && !dir.exists(path)) {
    entry$sha256 <- sha256_file(path)
    entry$bytes <- as.integer(file.info(path)$size)
    entry$status <- "present"
  } else entry$status <- "missing"
  .manifest_append(cfg, "inputs", entry)
  entry
}

record_params <- function(cfg, params) {
  m <- read_manifest(cfg)
  m$params <- c(m$params %||% list(), params %||% list())
  m$seed <- cfg$analysis$seed
  m$dataset_id <- cfg$dataset_id
  m$updated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  write_json(manifest_path(cfg), m)
  invisible(NULL)
}

record_decision <- function(cfg, node, question, answer, evidence = "") {
  # evidence 要写**支持这个选择的实际数字**，不是"因为这是通行做法"。
  .manifest_append(cfg, "decisions", list(
    node = node, question = question, answer = answer, evidence = evidence,
    at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")))
}

record_human_review <- function(cfg, node, required = TRUE, status = "pending", note = "") {
  # **默认 pending 而不是 confirmed。** 自动化流水线不能替人签字。
  .manifest_append(cfg, "human_review", list(
    node = node, required = isTRUE(required), status = status, note = note,
    at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")))
}

record_cross_language <- function(cfg, src, dst, format_, before = list(),
                                  after = list(), lost = character(0),
                                  tool = "", note = "") {
  # 只走 CSV；桥接工具限 zellkonverter / anndata2ri，**禁止 sceasy**。
  .manifest_append(cfg, "cross_language", list(
    src = src, dst = dst, format = format_, tool = tool,
    before = before, after = after, lost_fields = sort(unique(lost)),
    note = note, at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")))
}

capture_versions <- function(cfg, key_packages = KEY_PACKAGES) {
  # §0.3 版本记录。R 版枚举 installed.packages()（等价于 importlib 全量枚举，
  # 不起子进程）；枚举失败**必须落盘**成 versions_enumeration="failed"，
  # 不能只打 WARN（R-03 裁决 L12：日志没有人读）。
  full <- list()
  enum_error <- NA_character_
  tryCatch({
    ip <- installed.packages()
    for (i in seq_len(nrow(ip))) full[[ip[i, "Package"]]] <- ip[i, "Version"]
  }, error = function(e) {
    enum_error <<- sprintf("%s: %s", class(e)[1], conditionMessage(e))
    log_warn(sprintf("枚举已安装包失败（%s）—— versions 会不完整", enum_error))
  })

  key <- list()
  for (p in key_packages) key[[p]] <- full[[p]] %||% NA_character_

  m <- read_manifest(cfg)
  m$versions <- full
  m$key_versions <- key
  m$n_packages <- length(full)
  m$versions_enumeration <- if (is.na(enum_error)) "ok" else "failed"
  if (!is.na(enum_error)) m$versions_enumeration_error <- enum_error
  m$n_key_resolved <- sum(!vapply(key, function(v) is.na(v) || is.null(v), logical(1)))
  m$n_key_total <- length(key)
  m$r_version <- as.character(getRversion())
  m$platform <- R.version$platform
  # `r` 块：Python 版经 probe_r_packages 写同一结构（rscript/r_version/
  # packages/reason）。R 版就是 R，rscript 恒可用 —— 结构保持一致，
  # 让消费端不用分语言判读。
  m$r <- list(rscript = TRUE, r_version = as.character(getRversion()),
              packages = key, reason = NA)
  m$dataset_id <- cfg$dataset_id
  m$updated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  write_json(manifest_path(cfg), m)

  missing <- names(key)[vapply(key, function(v) is.na(v) || is.null(v), logical(1))]
  if (length(missing)) {
    log_warn(sprintf("关键 R 包未安装（%d，R %s）：%s", length(missing),
                     as.character(getRversion()), paste(missing, collapse = ", ")))
  } else {
    log_info(sprintf("关键 R 包全部就位（%d 个，R %s）", length(key),
                     as.character(getRversion())))
  }
  invisible(key)
}

manifest_summary <- function(cfg) {
  # 给验收用的一行摘要。**必需项缺失和可选项缺失分开报**（"设计如此"
  # 不当"出错"），「登记了几个」与「确认了几个」分开给（E-29）。
  m <- read_manifest(cfg)
  if (!length(m)) return(list(present = FALSE))
  ins <- m$inputs %||% list()
  miss <- Filter(function(i) identical(i$status, "missing"), ins)
  key <- m$key_versions %||% list()
  res <- vapply(key, function(v) !(is.na(v) || is.null(v)), logical(1))
  hr <- m$human_review %||% list()
  list(
    present = TRUE,
    n_versions = length(m$versions %||% list()),
    versions_enumeration = m$versions_enumeration %||% "unknown",
    versions_enumeration_error = m$versions_enumeration_error %||% NA,
    n_key_resolved = m$n_key_resolved %||% sum(res),
    n_key_total = m$n_key_total %||% length(key),
    key_unresolved = sort(names(key)[!res]),
    n_inputs = length(ins),
    inputs_missing = sort(vapply(miss, function(i) i$label %||% "?", character(1))),
    inputs_missing_required = sort(vapply(
      Filter(function(i) isTRUE(i$required), miss),
      function(i) i$label %||% "?", character(1))),
    n_decisions = length(m$decisions %||% list()),
    human_review_pending = sort(vapply(
      Filter(function(h) identical(h$status, "pending"), hr),
      function(h) h$node, character(1))),
    n_human_review = length(hr),
    human_review_confirmed = sort(vapply(
      Filter(function(h) h$status %in% c("confirmed", "overridden", "not_needed"), hr),
      function(h) h$node, character(1))),
    n_cross_language = length(m$cross_language %||% list())
  )
}

# Python 版 NAMED_TOOLS 表（§2 点名但 Python 版用不了的工具）在 R 版的对应物：
# R 版**就是** R，SoupX / SCTransform / scran / DESeq2 / edgeR / Monocle3 /
# Slingshot / CellChat 都装得上 —— 那张"用不了"的表对 R 版不成立，照抄就是
# 写假话。仍然用不了的只剩 Python 侧那三个与语言无关的缺口。
named_tools_note <- function() {
  paste0("Python 版 §2.1–§2.8 点名而用不了的 R 包（SoupX / SCTransform / scran / ",
         "DESeq2 / edgeR / Monocle3 / Slingshot / CellChat），R 版按立项文档 §3 ",
         "接上了其中与本流水线相关的那部分（scDblFinder / DESeq2 / destiny / ",
         "slingshot / liana）。仍然没有的：CellBender（GPU 导向，GitHub runner 无 GPU）、",
         "scVelo（需要 spliced/unspliced 层）、pySCENIC 的 motif 剪枝（GB 级排名库）",
         " —— **07 共表达推断在两版都刻意不做 motif 剪枝**，方法保持一致。",
         "Palantir 无 R 实现，05 轨迹的第四个方法槽位在 R 版丢弃（r_version.md §3）。",
         "这是缺口，不是「已覆盖」。")
}

# ============================================================================
# 出图样式与调色板
# ============================================================================
# 规范来源：scientific-agent-skills/skills/scientific-visualization
# （assets/publication.mplstyle 与 assets/color_palettes.py，K-Dense，MIT）。
# PAL 的十六进制值与 Python 版**逐位相同**（check_palette.mjs 是同一个文件），
# 挑选理由（删 yellow/sky_blue、10 色判据、grey 排最后）见 common.py:1011-1079。
PAL <- list(
  orange = "#E69F00", green = "#009E73", blue = "#0072B2",
  vermillion = "#D55E00", purple = "#CC79A7", black = "#000000",
  primary = "#0072B2", highlight = "#D55E00", muted = "#999999",
  cyan = "#17BECF", maroon = "#A50F15", indigo = "#5B4FCF", grey = "#A8A8A8"
)
# 分类循环色必须 >= 类别数（10 簇实测），顺序与 Python 版一致 ——
# 前 5 个不变让 <=5 类别的图外观不变。
PAL_CYCLE <- c(PAL$blue, PAL$vermillion, PAL$green, PAL$purple, PAL$black,
               PAL$orange, PAL$cyan, PAL$maroon, PAL$indigo, PAL$grey)

# ggplot 的"出版样式"：对应 assets/publication.mplstyle 的关键项。
# **幂等**（theme 对象不携带状态）。每张图都要加上它 —— save_fig 会自动补。
theme_paper <- function(base_size = 8) {
  need_pkg("ggplot2")
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(size = base_size),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 1),
      plot.title = ggplot2::element_text(size = base_size),
      plot.subtitle = ggplot2::element_text(size = base_size - 1),
      legend.position = "right",            # 规则 31：框外右侧
      legend.title = ggplot2::element_text(size = base_size - 1),
      legend.text = ggplot2::element_text(size = base_size - 1),
      strip.background = ggplot2::element_blank()
    )
}

scale_colour_pal <- function(name = NULL, ...) {
  ggplot2::scale_colour_manual(values = PAL_CYCLE, name = name, ...)
}
scale_fill_pal <- function(name = NULL, ...) {
  ggplot2::scale_fill_manual(values = PAL_CYCLE, name = name, ...)
}

# 毫米 -> 英寸。投稿图的尺寸单位是毫米，不是英寸（89 / 136 / 183 三档）。
# **别写 `c(..., recycle0 = FALSE)`**：传给 `c()` 的命名参数会变成
# **带名字的元素**而不是选项 —— `mm(89)` 会返回 `c(89=89, recycle0=0)`
# 两个值，几何参数悄悄错掉且不报错。逐个 `/` 即可。
mm <- function(...) {
  v <- c(...)
  v / 25.4
}

W_SINGLE    <- mm(89)    # 单栏
W_ONE_HALF  <- mm(136)   # 单栏半
W_DOUBLE    <- mm(183)   # 双栏（= 满版宽）

grid_size <- function(n_cols, n_rows, col_mm = 11.4, row_mm = 8.6,
                      min_w_mm = 136.0, max_w_mm = 183.0,
                      min_h_mm = 60.0, max_h_mm = 200.0) {
  # 按"列数 x 行数"算热图尺寸，全程毫米，**高度也封顶**（common.py:1140-1162：
  # 单位混用不会报错，只会把"装不装得进一页"变成没人检查的猜测）。
  w_mm <- min(max_w_mm, max(min_w_mm, col_mm * max(1, n_cols) + 20.0))
  h_mm <- min(max_h_mm, max(min_h_mm, row_mm * max(1, n_rows) + 14.0))
  list(width = mm(w_mm), height = mm(h_mm))
}

save_fig <- function(cfg, name, plot = NULL, width = NULL, height = NULL) {
  # 保存一张图为 PNG + PDF，返回写出的路径。
  # **PDF 与 PNG 都要出**（PDF 矢量可再编辑；PNG 给 check_figures 量像素）。
  # PDF 用 cairo_pdf（字体内嵌，对应 mpl 的 Type 42）。
  # R 版没有 tight/bbox_inches —— ggsave 写出的就是给定尺寸（见文件头注释 3）。
  need_pkg("ggplot2")
  figdir <- cfg$output$figures_dir
  dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
  p <- plot
  if (is.null(p)) stop("save_fig: plot 为空（R 版不读 last_plot，显式传图）", call. = FALSE)
  if (inherits(p, "ggplot") || inherits(p, "patchwork")) p <- p + theme_paper()
  w <- width %||% W_ONE_HALF
  h <- height %||% mm(76)
  dpi <- as.integer(cfg$analysis$figure_dpi %||% 300)
  out <- c(file.path(figdir, paste0(name, ".png")),
           file.path(figdir, paste0(name, ".pdf")))
  ggplot2::ggsave(out[1], p, width = w, height = h, units = "in", dpi = dpi,
                  device = grDevices::png)
  ggplot2::ggsave(out[2], p, width = w, height = h, units = "in",
                  device = grDevices::cairo_pdf)
  out
}

# ============================================================================
# 小工具
# ============================================================================
df_to_records <- function(df) {
  # data.frame -> JSON 安全的 records（NA 变 null）。
  need_pkg("jsonlite")
  if (is.null(df) || nrow(df) == 0L) return(list())
  jsonlite::fromJSON(jsonlite::toJSON(df, na = "null"), simplifyVector = FALSE)
}

# ============================================================================
# marker dotplot（common.py:1197-1358 同语义；plot_marker_dotplot +
# build_marker_dotplot_figure 合成一个 ggplot 构建器）
#
# **为什么抽在 common 而不是内联在脚本里（common.py docstring 的教训）**：
# 内联会导致"脚本一份、验证脚本照抄一份"，两份镜像不同步 —— 改一处
# 另一处不知道，"全 PASS"是假的。脚本与验证调用同一份代码。
#
# 语义逐条对齐 Python 版：
#   - 颜色 = **真 z-score**（按基因跨簇标准化），不是 scanpy 旧版
#     standard_scale="var" 的 min-max 0-1（那是标度矛盾：副标题写 z-score
#     但 0-1 不可能有负值）。RdBu_r 对称色标，0 = 簇间平均水平。
#   - 点大小 = 表达细胞比例（frac > 0 计表达）。
#   - vmax = max(1.5, p95(|z|))——个别极端基因不把色标撑得两极分化；
#     色标刻度按实际 vmax 生成（L11），不写死 [-1,0,1]。
#   - 图例放**主图下方横带**（Seurat do_DotPlot 范式，用户第七/八轮反馈）：
#     左半点大小四档（25/50/75/100%），右半横向色标。竖着塞右侧窄列
#     两块图例加标题没有不受挤的排法（实测调四轮仍相撞）。
#   - 列数上限 24 由调用方负责（双栏 183mm / 每基因约 7mm）。
#
# :param frac_df: data.frame，rowname=分组，列=基因，值 0..1（表达比例）
# :param z_df:    data.frame，同形状，值 = 按基因 z-score 后的表达
# :param group_label: y 轴标签（"Leiden cluster"）
# :param title / subtitle: 总标题 / 主图第二行说明
# :return: ggplot 对象（save_fig 自动补 theme_paper）
# ============================================================================
build_marker_dotplot_figure <- function(frac_df, z_df, group_label, title,
                                        subtitle, fig_width = NULL,
                                        fig_height_mm = 126) {
  need_pkg("ggplot2")
  w <- fig_width %||% W_DOUBLE
  genes <- colnames(z_df)
  groups <- rownames(z_df)
  ng <- length(genes); ngrp <- length(groups)

  # 长表：geom_point 的 size/colour 双标度映射
  long <- data.frame(
    gene = factor(rep(genes, times = ngrp), levels = genes),
    group = factor(rep(groups, each = ng), levels = rev(groups)),
    frac = as.numeric(t(frac_df[groups, genes, drop = FALSE])),
    z = as.numeric(t(z_df[groups, genes, drop = FALSE])))
  # 对称归一上限（Python common.py:1243 同式）
  vmax <- max(1.5, stats::quantile(abs(long$z), 0.95, na.rm = TRUE))
  vmax <- round(as.numeric(vmax), 1)
  long$size_pt <- 3 + long$frac * 5.5   # 面积观感对应 size_min=10..size_max=170

  p <- ggplot2::ggplot(long, ggplot2::aes(x = gene, y = group,
                                          size = frac, colour = z)) +
    ggplot2::geom_point(shape = 21, stroke = 0.3,
                        fill = "transparent") +
    ggplot2::scale_colour_gradient2(
      low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
      midpoint = 0, limits = c(-vmax, vmax),
      name = "Mean Expression",
      breaks = c(-vmax, 0, vmax), labels = c("Low", "Mid", "High")) +
    ggplot2::scale_size(range = c(0.5, 6), limits = c(0, 1),
                        breaks = c(0.25, 0.5, 0.75, 1),
                        labels = c("25", "50", "75", "100"),
                        name = "Percent Expressed (%)") +
    ggplot2::scale_x_discrete(labels = sort_unique_labels(genes)) +
    ggplot2::labs(x = "gene", y = group_label, title = title,
                  subtitle = subtitle) +
    # theme_paper 是本文件自定义主题，不是 ggplot2 导出 —— 不能加 ggplot2::
    # 前缀（"not an exported object"，CI run7 实跑抓到）。
    theme_paper(base_size = 8) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 90, hjust = 1, vjust = 0.4, size = 7),
      legend.box = "horizontal",
      legend.position = "bottom") # LEGEND_DIFF: 底部横带 = do_DotPlot 范式（R-06 定版，双图例横排）
  p
}
# 基因列原顺序展示（ggplot scale_x_discrete 的 labels 需与 levels 对齐，
# 这里直接用因子 levels 保持传入顺序，标签原样）
sort_unique_labels <- function(genes) genes

# ============================================================================
# 自检（--selftest）
# ============================================================================
# 与 tenifold_knk.R --selftest 同一模式：**只做 JSON/配置/状态文件往返**，
# 不做任何分析计算，本地与 CI 都能跑。判据都是往返一致或显式失败，
# 一个静默的空通过都不要（E-63：分支永不执行时错永不显形）。
selftest_common <- function() {
  cases <- list(); ok <- 0L; failed <- 0L
  check <- function(label, cond, detail = "") {
    if (isTRUE(cond)) { ok <<- ok + 1L; cat(sprintf("  [OK]   %s\n", label))
    } else { failed <<- failed + 1L
             cat(sprintf("  [FAIL] %s  <-  %s\n", label, detail)) }
  }
  tmp <- file.path(tempdir(), paste0("scrna_common_selftest_",
                                     format(Sys.time(), "%H%M%OS")))
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)

  # 1) load_config：dataset_id 派生目录 + 兜底值
  yml <- file.path(tmp, "c.yml")
  writeLines(c("dataset_id: selftest", "analysis:", "  seed: 7"), yml)
  cfg <- load_config(yml)
  check("load_config 派生三个目录",
        identical(cfg$output$results_dir, "results/selftest") &&
          identical(cfg$output$data_dir, "data/selftest") &&
          identical(cfg$output$figures_dir, "results/selftest/figures"),
        paste(capture.output(str(cfg$output)), collapse = " "))
  check("load_config 兜底 seed/dpi",
        identical(as.integer(cfg$analysis$seed), 7L) &&
          identical(as.integer(cfg$analysis$figure_dpi), 300L))
  check("load_config 拒绝缺 dataset_id", {
    bad <- file.path(tmp, "bad.yml"); writeLines("analysis:", bad)
    inherits(tryCatch(load_config(bad), error = function(e) e), "error")
  })
  check("load_config 拒绝不存在的文件",
        inherits(tryCatch(load_config(file.path(tmp, "nope.yml")),
                          error = function(e) e), "error"))

  # 2) write_json / read_json：非有限值 -> null，且**结构可往返**
  p <- file.path(tmp, "a.json")
  write_json(p, list(a = 1.5, b = c(1, NaN, Inf), c = list(d = NA_real_)))
  txt <- paste(readLines(p, warn = FALSE), collapse = "")
  check("write_json 把 NaN/Inf 写成 null",
        !grepl("NaN|Infinity|-Inf", txt), txt)
  back <- read_json(p)
  check("read_json 往返结构一致",
        identical(as.numeric(back$a), 1.5) && length(back$b) == 3L &&
          is.null(back$b[[2]]) && is.null(back$b[[3]]) && is.null(back$c$d),
        paste(capture.output(str(back)), collapse = " "))

  # 3) read_json：不存在 -> NULL；坏内容 -> stop（M8）
  check("read_json 不存在返回 NULL", is.null(read_json(file.path(tmp, "nope.json"))))
  badp <- file.path(tmp, "bad.json"); writeLines("{oops", badp)
  check("read_json 坏内容抛错",
        inherits(tryCatch(read_json(badp), error = function(e) e), "error"))
  check("read_json_or_none 坏内容返回 NULL",
        is.null(read_json_or_none(badp)))

  # 4) 步骤状态：result_status 与编排器 status 是两件事（E-56）
  cfg$output$results_dir <- file.path(tmp, "results/selftest")
  record_step(cfg, "qc", "ok", seconds = 1.2, result_status = "bad_root")
  st <- read_state(cfg)
  e1 <- st$steps[[1]]
  check("record_step 记下两件事",
        identical(e1$status, "ok") && identical(e1$result_status, "bad_root"),
        paste(capture.output(str(e1)), collapse = " "))
  check("classify：abort 值判 abort / 未知值判 skip / 空判 ok",
        classify_step_result("bad_root") == "abort" &&
          classify_step_result("weird_new_word") == "skip" &&
          classify_step_result("") == "ok")
  record_step(cfg, "qc", "ok", seconds = 2.0)   # 同 id 覆盖，不追加
  check("record_step 同 id 覆盖", length(read_state(cfg)$steps) == 1L)

  # 5) manifest：登记/摘要/计数
  init_manifest(cfg, language = "r")
  record_input(cfg, yml, label = "配置", required = TRUE)
  record_input(cfg, file.path(tmp, "missing.csv"), label = "可选项", required = FALSE)
  record_decision(cfg, "selftest", "q", "a", evidence = "没有量化依据")
  record_human_review(cfg, "celltype_labels")
  record_params(cfg, list(k = "v"))
  s <- manifest_summary(cfg)
  check("manifest_summary 分开报必需/可选缺失",
        isTRUE(s$present) && identical(s$inputs_missing_required, "配置") &&
          identical(s$inputs_missing, "可选项"),
        paste(capture.output(str(s)), collapse = " "))
  check("manifest_summary 分开报登记数与确认数",
        identical(s$n_human_review, 1L) && identical(s$human_review_pending, "celltype_labels") &&
          length(s$human_review_confirmed) == 0L)
  check("manifest_summary 计数键齐全",
        all(c("n_versions", "n_key_resolved", "n_key_total", "versions_enumeration",
              "n_decisions", "n_cross_language") %in% names(s)))
  capture_versions(cfg, key_packages = c("jsonlite", "不存在的包名xyz"))
  s2 <- manifest_summary(cfg)
  check("capture_versions：key_versions 记 NA 而不是丢键",
        s2$n_key_total == 2L && "不存在的包名xyz" %in% s2$key_unresolved,
        paste(capture.output(str(s2[c("n_key_resolved","n_key_total","key_unresolved")])), collapse = " "))

  # 6) mm / grid_size
  check("mm() 换算（含 E-61 形态回归：c() 命名参数陷阱）",
        length(mm(25.4)) == 1L && abs(mm(25.4) - 1) < 1e-12 &&
          abs(W_DOUBLE - 183 / 25.4) < 1e-12 && length(W_SINGLE) == 1L,
        sprintf("mm(25.4) 长度 %d", length(mm(25.4))))
  check("finite_round：非有限值落 NA / 0 保留 / 位数正确",
        is.na(finite_round(NaN)) && is.na(finite_round(Inf)) &&
          is.na(finite_round(NULL)) && identical(finite_round(0), 0) &&
          identical(finite_round(0.123456), 0.1235) && is.na(finite_round("abc")))
  gs <- grid_size(n_cols = 24, n_rows = 24)
  check("grid_size 高度封顶 200mm",
        gs$height <= mm(200) + 1e-9 && gs$width >= mm(136) - 1e-9,
        sprintf("w=%.2f h=%.2f in", gs$width, gs$height))

  # 7) df_to_records：NA -> null
  df <- data.frame(a = c(1, NA), b = c("x", "y"))
  r <- df_to_records(df)
  check("df_to_records NA 变 null", length(r) == 2L && is.null(r[[2]]$a))

  n <- ok + failed
  cat(sprintf("自检%s（%d 个用例，通过 %d，失败 %d）\n",
              if (failed == 0L) "通过" else "失败", n, ok, failed))
  failed == 0L
}

if (identical(Sys.getenv("SCRNA_COMMON_SELFTEST"), "1") ||
    ("--selftest" %in% commandArgs(trailingOnly = TRUE))) {
  good <- selftest_common()
  quit(save = "no", status = if (good) 0L else 1L)
}
