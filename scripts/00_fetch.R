#!/usr/bin/env Rscript
# =============================================================================
# 00_fetch.R — 取数、读入、硬校验（R 版，与 00_fetch.py 同号同语义）
#
# 方法学等级（references/r_version.md §3）：**B 换实现同思路**。
# 差异点与对齐义务：
#   - 读入：`Seurat::Read10X()`（DenseMatrix/DGEMatrix 均可）替换
#     `sc.read_10x_mtx`。**坐标约定硬约束：rownames = 基因符号、
#     colnames = 细胞条码** —— 与 Python 版一致，否则后续全部错位。
#     重复基因符号照 Python 版 `var_names_make_unique()` 的语义
#     （首个保留原名，其余加 `-1` 后缀）。
#   - 下载：`download.file(mode="wb")` + `options(HTTPUserAgent=...)`。
#     Python 侧的教训（00_fetch.py:49-91）在 R 侧**要重新踩一遍**：
#     R 默认 UA 是 `R/x.y.z`，10x 的 CDN 同样会 403 —— 所以 UA 必须伪装。
#   - 计数校验 `validate_counts`：与 Python 版同构（两批抽查、结论一致才算数），
#     抽样行号逐位对齐（等距主判据 + 中段连续交叉核对，审计 S10）。
#
# 硬门禁与落盘字段与 Python 版逐项相同：
#   MIN_CELLS / MIN_GENES / MAX_CELLS_WARN、extract_filter（M14）、
#   max_cells_warn_gate / n_cells_over_warn_gate（L18）。
#
# 中间对象走 .rds（references/r_version.md §5）：raw.rds 而非 raw.h5ad。
# 两版产物目录同名但中间对象格式不同 —— 对照请用 CSV 产物。
# =============================================================================

# 独立 Rscript 运行时由下方 .load_common_for_local() 把 lib/common.R 载入
# globalenv；被 main_analysis.R source() 时 common.R 已先载入，此步是空操作。

# ---- 硬门禁常量（00_fetch.py:31-37 逐值相同）--------------------------------
MIN_CELLS <- 50L
MIN_GENES <- 50L
MAX_CELLS_WARN <- 200000L
# 为什么是 20 万（L18）：本地实测 ~2 万细胞用 ~1.5 GB、CI runner 7 GB 内存，
# 20 万约 15 GB —— 已经**超过 runner 内存**，只告警不拒绝是因为有些分析
# （只做 QC + 聚类）确实跑得完。这个依据原先只存在于作者脑子里，
# 现在写下来，否则没人知道这个数字能不能改。

.load_common_for_local <- function() {
  # 独立 Rscript 运行时：把 lib/common.R source 进本环境。
  argv0 <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(argv0)) {
    self <- normalizePath(sub("^--file=", argv0[1]), mustWork = FALSE)
    lib <- file.path(dirname(self), "lib", "common.R")
    if (file.exists(lib) && !exists("record_step", envir = globalenv())) {
      sys.source(lib, envir = globalenv())
    }
  }
  invisible(NULL)
}
.load_common_for_local()

# ---------------------------------------------------------------------------
# 登记表（00_fetch.py:40-46 同语义）
# ---------------------------------------------------------------------------
load_registry <- function(repo_root) {
  p <- file.path(repo_root, "assets", "datasets.yml")
  if (!file.exists(p)) return(list())
  need_pkg("yaml", "读 assets/datasets.yml")
  y <- yaml::read_yaml(p)
  y$datasets %||% list()
}

# ---------------------------------------------------------------------------
# 带重试的下载（00_fetch.py:49-91 同语义）
#
# **必须带 User-Agent。** 10x 的 CDN 对默认 UA（Python-urllib / R-x.y）直接
# 403，而浏览器请求同 URL 是 200 —— 实测踩过，报错只有 "HTTP 403"，
# 很容易误判成"链接失效了"。R 侧用 options(HTTPUserAgent=...) 伪装。
# 已存在且非空则跳过（CI 里 data/ 不持久化，本地会命中）。
# ---------------------------------------------------------------------------
download <- function(url, dest, retries = 3L) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest) && file.info(dest)$size > 0) {
    log_info(sprintf("已缓存，跳过下载: %s (%.1f MB)",
                     basename(dest), file.info(dest)$size / 1e6))
    return(invisible(dest))
  }
  old_ua <- options(HTTPUserAgent = sprintf(
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0 Safari/537.36 R-fetch"))
  on.exit(options(old_ua), add = TRUE)
  last <- NULL
  for (attempt in seq_len(retries)) {
    ok <- tryCatch({
      log_info(sprintf("下载 (%d/%d): %s", attempt, retries, url))
      tmp <- paste0(dest, ".part")
      status <- download.file(url, tmp, mode = "wb", quiet = TRUE,
                              method = "libcurl", cacheOK = FALSE)
      if (status != 0L) stop(sprintf("download.file 返回状态 %d", status), call. = FALSE)
      if (!file.exists(tmp) || file.info(tmp)$size == 0) stop("下载到 0 字节", call. = FALSE)
      if (file.exists(dest)) unlink(dest)
      if (!file.rename(tmp, dest)) stop("改名失败", call. = FALSE)
      log_info(sprintf("下载完成: %s (%.1f MB)",
                       basename(dest), file.info(dest)$size / 1e6))
      TRUE
    }, error = function(e) {
      last <<- sprintf("%s: %s", class(e)[1], conditionMessage(e))
      log_warn(sprintf("下载失败 (%d/%d): %s", attempt, retries, last))
      FALSE
    })
    if (isTRUE(ok)) return(invisible(dest))
  }
  stop(sprintf("下载 %s 失败（重试 %d 次）: %s", url, retries, last), call. = FALSE)
}

# ---------------------------------------------------------------------------
# 解压 10x tar.gz 并读入（00_fetch.py:94-145 同语义）
#
# **10x 的 tar 里目录名带版本后缀**（filtered_gene_bc_matrices/hg19/），
# 不能写死路径 —— 找含 matrix.mtx 的**最浅**目录。
#
# M14（R-03 裁决）：解压方式必须留痕。R 的 untar() 没有路径穿越过滤参数，
# 等价于 Python 侧的 `none_fallback` 降级路径 —— **这是事实，不是缺陷**，
# 但必须写进产物，不能让"过滤过"与"没过滤"长得一样。10x 官方 tar 包风险低。
# ---------------------------------------------------------------------------
read_10x_tar <- function(tar_path, cache_dir) {
  extract_to <- file.path(cache_dir, "10x_extracted")
  marker <- file.path(extract_to, ".done")
  extract_filter <- "untar_no_filter"
  if (!file.exists(marker)) {
    dir.create(extract_to, recursive = TRUE, showWarnings = FALSE)
    log_info(sprintf("解压 %s -> %s", basename(tar_path), extract_to))
    untar(tar_path, exdir = extract_to)
    writeLines("ok", marker)
  } else {
    log_info(sprintf("已解压，跳过: %s", extract_to))
  }

  # 找含 matrix.mtx(.gz) 的目录，取路径最浅的
  hits <- list.files(extract_to, pattern = "^matrix\\.mtx(\\.gz)?$",
                     recursive = TRUE, full.names = TRUE)
  if (!length(hits)) {
    stop(sprintf("%s 里找不到 matrix.mtx —— 不是 10x 计数矩阵格式",
                 basename(tar_path)), call. = FALSE)
  }
  # 深度按路径分隔段数算（"/" 计数比 normalizePath 的平台差异可靠）
  depth <- vapply(hits, function(p) length(strsplit(gsub("\\\\", "/", p),
                                                    "/", fixed = TRUE)[[1]]), integer(1))
  mtx_file <- hits[[which.min(depth)]]
  mtx_dir <- dirname(mtx_file)
  log_info(sprintf("10x 矩阵目录: %s", mtx_dir))

  need_pkg("Seurat", "读 10x 三件套")
  m <- Seurat::Read10X(data.dir = mtx_dir, var.names.repair = "unique",
                       strip.suffix = FALSE)
  # Read10X 返回单个矩阵或按 feature 类型分层的 list（Gene Expression 在前）。
  if (is.list(m)) {
    # 取第一个（"Gene Expression"），与 scanpy read_10x_mtx 的 gene 层一致。
    m <- m[[1]]
  }
  # scanpy 语义：细胞 x 基因（obs x var）。Seurat 给 基因 x 细胞 —— 转置。
  # **稀疏矩阵不物化成 dense**：dgCMatrix 转置在 Matrix 内部完成，内存安全。
  counts <- Matrix::t(m)
  list(counts = counts, extract_filter = extract_filter)
}

read_h5ad_counts <- function(path) {
  # kind=h5ad 的 R 侧读入。**不引 zellkonverter/basilisk**（立项 §5 的裁决：
  # conda 运行时又慢又脆）。h5ad 的 X 若是 CSR 稀疏矩阵，HDF5 布局是
  # /X/data + /X/indices + /X/indptr —— 用 rhdf5 手工拼回 dgCMatrix。
  # 已 normalize 的 h5ad 会被 validate_counts 拦下（与 Python 版同一道门）。
  need_pkg("rhdf5", "读 h5ad（kind=h5ad 时）")
  rhdf5::h5disableLibLoc()
  h5 <- tryCatch(rhdf5::H5Fopen(path), error = function(e) {
    stop(sprintf("h5ad 打不开（%s）。R 版读 h5ad 需要 rhdf5；若数据来自 scanpy ",
                 conditionMessage(e)), call. = FALSE)
  })
  on.exit(rhdf5::H5Fclose(h5), add = TRUE)
  rd <- rhdf5::h5read(h5, "obs")          # 细胞元数据（data.frame）
  vd <- rhdf5::h5read(h5, "var")          # 基因元数据
  obs_names <- if (!is.null(rd$`_index`)) as.character(rd$`_index`) else NULL
  var_names <- if (!is.null(vd$`_index`)) as.character(vd$`_index`) else NULL
  if (is.null(obs_names) || is.null(var_names)) {
    stop("h5ad 缺 obs/_index 或 var/_index —— 不是合法的 AnnData 文件", call. = FALSE)
  }
  if (rhdf5::h5exists(h5, "X/data")) {
    vals <- as.numeric(rhdf5::h5read(h5, "X/data"))
    idxs <- as.integer(rhdf5::h5read(h5, "X/indices"))
    ptrs <- as.integer(rhdf5::h5read(h5, "X/indptr"))
    shape <- rhdf5::h5read(h5, "X/shape")
    # scanpy 的 CSR：行 = 细胞。Matrix::sparseMatrix(i,j,p=) 直接按 CSR 建。
    counts <- Matrix::sparseMatrix(p = ptrs, j = idxs, x = vals,
                                   dims = as.integer(shape), index1 = FALSE)
  } else {
    xd <- rhdf5::h5read(h5, "X")
    counts <- Matrix::Matrix(as.numeric(xd), nrow = length(obs_names),
                             ncol = length(var_names), sparse = TRUE)
  }
  dimnames(counts) <- list(obs_names, var_names)
  counts
}

read_10x_h5_counts <- function(path) {
  need_pkg("Seurat", "读 10x h5")
  m <- Seurat::Read10X_h5(filename = path, use.names = TRUE,
                          unique.names = TRUE)
  if (is.list(m)) m <- m[[1]]
  Matrix::t(m)   # 基因 x 细胞 -> 细胞 x 基因，与 Python 版坐标约定一致
}

# ---------------------------------------------------------------------------
# 计数性质校验（00_fetch.py:158-243 同语义，两批抽查逐位对齐）
#
# **这是本步骤存在的理由。** 把已 normalize 的矩阵当计数喂进去，QC 指标失真、
# HVG 前提被破坏、pseudobulk 的负二项模型不成立，而所有图和数字看起来都正常。
#
# **抽查的是"跨全表的非零元素"，不是"前 200 个细胞"**（审计 S10）。
# 两批抽法不同才有意义：等距抽样覆盖全表，中段连续块才能暴露"数据按样本
# 拼接、各样本处理方式不同"。占比差 >1 个百分点即判不一致。
# ---------------------------------------------------------------------------
validate_counts <- function(counts, cfg) {
  n_obs <- nrow(counts)

  .vals_of <- function(rows) {
    if (!length(rows)) return(numeric(0))
    sub <- counts[rows, , drop = FALSE]
    # 只取非零元素：稀疏矩阵零是隐式的，`which(..., arr.ind)` 不会物化
    # 全矩阵；dense 情形 Matrix::as.vector 也安全。
    if (methods::is(sub, "sparseMatrix")) {
      v <- sub@x
    } else {
      v <- as.numeric(as.vector(sub))
    }
    v[is.finite(v)]
  }
  .strided <- function(k = 200L) {
    # 跨全表等距抽 k 行 —— **主判据**
    if (n_obs == 0L) return(list(vals = numeric(0), n = 0L))
    step <- max(1L, n_obs %/% k)
    rows <- seq.int(1L, n_obs, by = step)
    if (length(rows) > k) rows <- rows[seq_len(k)]
    list(vals = .vals_of(rows), n = length(rows))
  }
  .contiguous <- function(start, k = 200L) {
    # 从 `start` 起抽**连续** k 行 —— **交叉核对**
    if (n_obs == 0L) return(list(vals = numeric(0), n = 0L))
    rows <- seq.int(start, min(start + k - 1L, n_obs))
    list(vals = .vals_of(rows), n = length(rows))
  }

  b1 <- .strided()
  b2 <- .contiguous(n_obs %/% 2L)
  if (!length(b1$vals)) {
    return(list(is_counts = FALSE, reason = "矩阵里没有有限值",
                n_cells_checked = 0L))
  }

  close_to_int <- function(v) abs(v - round(v)) < 1e-8
  is_int <- all(close_to_int(b1$vals))
  min_v <- min(b1$vals)
  max_v <- max(b1$vals)
  frac_int <- mean(close_to_int(b1$vals))

  frac_int2 <- if (length(b2$vals)) mean(close_to_int(b2$vals)) else NULL
  consistent <- is.null(frac_int2) || abs(frac_int - frac_int2) < 0.01

  ok <- is_int && min_v >= 0 && consistent
  if (!consistent) {
    reason <- sprintf(paste0("**两批抽查结论不一致**（等距抽样整数值占比 %.1f%% vs ",
                             "中段连续抽样 %.1f%%）—— 数据可能按样本拼接且各样本",
                             "处理方式不同，此时『是/不是计数』的二值判定不可靠"),
                      100 * frac_int, 100 * frac_int2)
  } else if (ok) {
    reason <- sprintf(paste0("非负整数 -> 判定为原始计数（跨 %d 个细胞等距抽 %d 行 + ",
                             "中段连续抽 %d 行，两批一致）"),
                      n_obs, b1$n, b2$n)
  } else {
    reason <- sprintf(paste0("**不是整数计数**（最小 %.4g，最大 %.4g，",
                             "整数值占比 %.1f%%）—— 可能已经 normalize/log 过"),
                      min_v, max_v, 100 * frac_int)
  }

  list(is_counts = ok, reason = reason, min = finite_round(min_v),
       max = finite_round(max_v),
       frac_integer = round(frac_int, 6),
       frac_integer_second_batch = if (is.null(frac_int2)) NULL else round(frac_int2, 6),
       sampling_consistent = consistent,
       sampling_note = paste0("第一批 = 跨全表等距抽样；第二批 = 从中段起连续抽样。",
                              "两批抽法不同，占比若差 >1 个百分点即判为不一致"),
       n_cells_total = n_obs,
       n_cells_checked = b1$n + b2$n,
       sampled_values = length(b1$vals) + length(b2$vals))
}

# ---------------------------------------------------------------------------
# 主流程（00_fetch.py:246-339 同语义）
# ---------------------------------------------------------------------------
run_00_fetch <- function(cfg) {
  ensure_dirs(cfg)
  set_seed(cfg)
  argv0 <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  repo_root <- if (length(argv0)) {
    dirname(dirname(normalizePath(sub("^--file=", argv0[1]))))
  } else {
    normalizePath(".")
  }
  data_dir <- cfg$output$data_dir
  cache_dir <- file.path(data_dir, "cache")

  src <- cfg$source %||% list()
  if (!is.null(src$dataset)) {
    reg <- load_registry(repo_root)
    key <- as.character(src$dataset)
    if (is.null(reg[[key]])) {
      stop(sprintf("source.dataset='%s' 不在 assets/datasets.yml 里；可选: %s",
                   key, paste(sort(names(reg)), collapse = ", ")), call. = FALSE)
    }
    entry <- reg[[key]]
    extra <- src[names(src) != "dataset"]
    for (nm in names(extra)) entry[[nm]] <- extra[[nm]]
    log_info(sprintf("数据集 '%s': %s", key, entry$note %||% ""))
  } else if (!is.null(src$url)) {
    entry <- as.list(src)
  } else {
    stop("配置缺少 source.dataset 或 source.url", call. = FALSE)
  }

  kind <- entry$kind
  url <- entry$url
  if (is.null(kind) || is.null(url)) {
    stop("数据源条目缺少 kind/url", call. = FALSE)
  }

  # ---- 取数 --------------------------------------------------------------
  extract_filter <- NULL
  counts <- NULL
  if (identical(kind, "10x_tar")) {
    tar_path <- download(url, file.path(cache_dir, basename(url)))
    rd <- read_10x_tar(tar_path, cache_dir)
    counts <- rd$counts
    extract_filter <- rd$extract_filter
  } else if (identical(kind, "10x_h5")) {
    h5_path <- download(url, file.path(cache_dir, basename(url)))
    counts <- read_10x_h5_counts(h5_path)
  } else if (identical(kind, "h5ad")) {
    h5_path <- download(url, file.path(cache_dir, basename(url)))
    counts <- read_h5ad_counts(h5_path)
  } else {
    stop(sprintf("不支持的 kind: %s（可选 10x_tar / 10x_h5 / h5ad）", kind),
         call. = FALSE)
  }

  # var_names_make_unique 同语义（Read10X 已用 unique.names 时是幂等的）
  if (anyDuplicated(colnames(counts))) {
    colnames(counts) <- make.unique(colnames(counts), sep = "-1")
  }
  if (nrow(counts) == 0L || ncol(counts) == 0L) {
    stop(sprintf("读入后为空: %d 细胞 x %d 基因", nrow(counts), ncol(counts)),
         call. = FALSE)
  }

  # ---- 硬校验 ------------------------------------------------------------
  counts_check <- validate_counts(counts, cfg)
  log_info(counts_check$reason)
  if (!isTRUE(counts_check$is_counts)) {
    stop(paste0("输入矩阵不是原始整数计数，拒绝继续。\n",
                "  依据: ", counts_check$reason, "\n",
                "  为什么必须停: QC 指标会失真、HVG 前提被破坏、\n",
                "  pseudobulk 的负二项模型不成立 —— 而所有图和数字看起来都正常。\n",
                "  修法: 用原始计数矩阵（10x 的 matrix.mtx / filtered_feature_bc_matrix）。"),
         call. = FALSE)
  }

  if (nrow(counts) < MIN_CELLS) {
    stop(sprintf("细胞数 %d < %d，单细胞分析无从谈起", nrow(counts), MIN_CELLS),
         call. = FALSE)
  }
  if (ncol(counts) < MIN_GENES) {
    stop(sprintf("基因数 %d < %d", ncol(counts), MIN_GENES), call. = FALSE)
  }
  if (nrow(counts) > MAX_CELLS_WARN) {
    log_warn(sprintf("细胞数 %d 较大，CI 上可能超出内存/时限", nrow(counts)))
  }

  # ---- 落盘 --------------------------------------------------------------
  raw_path <- file.path(data_dir, "raw.rds")
  saveRDS(counts, raw_path)
  log_info(sprintf("已写出 %s（%d 细胞 x %d 基因）",
                   raw_path, nrow(counts), ncol(counts)))

  info <- list(
    dataset_id = cfg$dataset_id,
    source_kind = kind,
    source_url = url,
    source_note = entry$note %||% "",
    organism = entry$organism,
    tissue = entry$tissue,
    condition = entry$condition,
    n_cells_raw = nrow(counts),
    n_genes_raw = ncol(counts),
    counts_check = counts_check,
    min_cells_gate = MIN_CELLS,
    min_genes_gate = MIN_GENES,
    # **M14**：解压时有没有做路径穿越过滤。R 的 untar() 没有过滤参数，
    # `untar_no_filter` 是降级路径的如实记录。
    extract_filter = extract_filter,
    # **L18（R-03 裁决）**：阈值必须落盘，否则事后想解释 CI 为什么慢
    # 没有任何依据。
    max_cells_warn_gate = MAX_CELLS_WARN,
    n_cells_over_warn_gate = nrow(counts) > MAX_CELLS_WARN
  )
  write_json(file.path(data_dir, "dataset_info.json"), info)
  info
}

# ---------------------------------------------------------------------------
# 入口（00_fetch.py:342-356 同语义；E-56：result_status 与 status 分开记）
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L ||
    identical(Sys.getenv("SCRNA_STEP_MAIN"), "00_fetch")) {
  args <- parse_args()
  cfg <- load_config(args$config)
  t0 <- Sys.time()
  tryCatch({
    out <- run_00_fetch(cfg)
    record_step(cfg, "fetch", "ok", as.numeric(difftime(Sys.time(), t0, units = "secs")),
                result_status = result_status_of(out))
    out
  }, error = function(e) {
    record_step(cfg, "fetch", "failed",
                as.numeric(difftime(Sys.time(), t0, units = "secs")),
                message = conditionMessage(e), result_status = "failed")
    stop(e)
  })
}
