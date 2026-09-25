#!/usr/bin/env Rscript
# =============================================================================
# tenifold_knk.R — 真 · scTenifoldKnk 虚拟敲除（R 引擎）
#
# ## 为什么有这个脚本
#
# `08_virtual_perturbation.py` 原来只有一条路径：在共表达 GRN 上做**一阶单跳
# 线性传播**。那是近似，不是 scTenifoldKnk。本脚本把规范 §1.7 点名的
# scTenifoldKnk 真正跑起来，与一阶近似**并列落盘**，让两法可比。
#
# ## 两处"读源码才发现"的关键事实
#
# K-01 立项时把这两条都登记错了（已修正，见 02_TASKLIST.md 的 K-01b）：
#
# 1. **矩阵方向是 genes × cells，不是 cells × genes。**
#    源码 `R/scTenifoldKnk.R` L10 与 `R/scQC.R` L9 逐字写着
#    "Raw counts matrix with cells as columns and genes (symbols) as rows."
#    写反了**不会报错** —— `makeNetworks` 照跑，只是把细胞当基因、基因当
#    细胞，算出一个垃圾网络，而所有状态字段看起来都正常。所以本脚本在入口
#    做显式断言（`assert_orientation`），并且自测里专门验"转置能被拦住"。
#
# 2. **应该用 `transcriptomeWide = TRUE`。** 源码 L183–L221：该模式下
#    **只建一次 WT 网络、批量扰动多个基因**；默认的单基因模式每扰动一个
#    基因就重建一次网络，慢 N 倍。
#
# ## 为什么 `qc = FALSE`
#
# 包默认 `qc = TRUE`，会用 `qc_minLibSize = 1000` / `qc_minPCT = 0.05` /
# `qc_maxMTratio = 0.1` 再过滤一遍。本仓库的 `01_qc.py` 已经做过 QC，
# 这里再滤一次会**改变细胞集** —— 那么 tenifold 与一阶近似用的就不是同一批
# 细胞，两法的数不可比。所以显式传 `qc = FALSE`，并在 meta 里记下来。
#
# ## 确定性
#
# 包在源码里**硬编码了四处 `set.seed(1)`**（L163/L172/L208/L231），这是包内
# 的设计。所以给定同一个输入矩阵与同一个包版本，结果是逐位确定的 ——
# 这一点与仓库其它"设了种子也没用"的随机源不同，但仍要在 meta 里记版本。
#
# ## 输入格式：裸文本矩阵，不是 CSV
#
#   <input>   第 1 行 = tab 分隔的基因名（n_genes 个）
#             第 2..n_genes+1 行 = 一个基因一行，tab 分隔的细胞计数
#   <meta>    JSON，含 n_genes / n_cells / genes（有序）/ gene_sums（每基因总和）
#
# **为什么不用 CSV**：pbmc3k 规模下这是 10^6–10^7 个整数，`read.csv` 要几十秒；
# 而它唯一的额外好处（自带列名）我们不需要，因为基因顺序本来就记在 meta 里。
# 更要紧的是 `read.csv` 会把细胞条码里的 `-`/`.` 改写成 `X` 前缀（要靠
# `check.names = FALSE` 兜），裸矩阵没有这个坑。
#
# ## 方向见证：`gene_sums`
#
# 光靠"行数 == n_genes 且列数 == n_cells"**拦不住方阵转置** —— 转置后方阵
# 行列数不变。所以 meta 里额外记了每个基因的总计数（`gene_sums`），读进来后
# 用 `rowSums` 复核：如果值被转置了，行和就变成了"每个细胞的文库大小"，
# 与基因总计数完全不同，立刻暴露。自测里专门造了这个退化情形。
#
# ## 用法
#
#   Rscript tenifold_knk.R --input in_matrix.txt --meta in_meta.json \
#                          --targets targets.txt --out_dir /tmp/out
#   Rscript tenifold_knk.R --selftest          # 合成数据自测，不需要任何输入
#
# ## 输出
#
#   <out_dir>/tenifold_perturbation_distances.csv  扰动基因 × 网络基因 的距离
#   <out_dir>/tenifold_meta.json                   版本 / 维度 / 参数 / 耗时
# =============================================================================

# ---- 参数解析 ---------------------------------------------------------------
# 不引 optparse：多一个依赖就多一个装不上的理由。参数表很小，手写够用，
# 而且**未知参数直接报错**（静默忽略拼错的参数会让人以为改动生效了）。
DEFAULTS <- list(
  input = NULL, meta = NULL, targets = NULL, out_dir = NULL,
  n_net = "10", n_cells = "500", n_comp = "3",
  td_k = "3", ma_n_dim = "2", n_cores = "1",
  selftest = FALSE
)

parse_args <- function(argv) {
  out <- DEFAULTS
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (identical(a, "--selftest")) {
      out$selftest <- TRUE
      i <- i + 1L
      next
    }
    if (!grepl("^--", a)) {
      stop("无法识别的参数（不是 -- 开头）: ", a, call. = FALSE)
    }
    key <- sub("^--", "", a)
    if (!key %in% names(DEFAULTS)) {
      stop("未知参数: ", a, call. = FALSE)
    }
    if (i + 1L > length(argv)) {
      stop("参数 ", a, " 缺少取值", call. = FALSE)
    }
    out[[key]] <- argv[[i + 1L]]
    i <- i + 2L
  }
  out
}

as_int <- function(x, name) {
  v <- suppressWarnings(as.integer(x))
  if (is.na(v) || v < 1L) {
    stop("参数 ", name, " 必须是正整数，收到: ", x, call. = FALSE)
  }
  v
}

# ---- 方向断言（本脚本存在的核心理由之一）------------------------------------
# 判据是**四条互相独立的事实**都要对上：
#   1. 行数 == Python 记的 n_genes
#   2. 列数 == Python 记的 n_cells
#   3. 行名**逐位**等于 Python 记的基因顺序
#   4. `rowSums(mat)` 等于 Python 记的 `gene_sums`
#
# 1)+2) 拦不住方阵转置（行列数不变）；3) 能拦住，但前提是基因名恰好也是
# 细胞名才会漏 —— 而 4) 是**数值证据**：转置后行和变成文库大小，数量级都
# 不一样。四条一起上，转置不可能静默通过。
assert_orientation <- function(mat, meta, label = "输入矩阵") {
  if (is.null(meta)) {
    stop("缺少输入 meta —— 没有它就无法断言矩阵方向（genes x cells）。",
         call. = FALSE)
  }
  if (!all(c("n_genes", "n_cells", "genes") %in% names(meta))) {
    stop("输入 meta 缺少 n_genes / n_cells / genes 字段", call. = FALSE)
  }
  n_genes <- as.integer(meta$n_genes)
  n_cells <- as.integer(meta$n_cells)
  if (nrow(mat) != n_genes || ncol(mat) != n_cells) {
    stop(sprintf(
      "%s 维度不符：收到 %d x %d，meta 记的是 %d genes x %d cells。%s",
      label, nrow(mat), ncol(mat), n_genes, n_cells,
      "**scTenifoldKnk 要求 genes x cells（基因行、细胞列）** —— 传反了不会报错，只会算出一个垃圾网络。"
    ), call. = FALSE)
  }
  expect <- as.character(meta$genes)
  if (!identical(rownames(mat), expect)) {
    bad <- which(rownames(mat) != expect)[1]
    stop(sprintf(
      "%s 的行名与 meta 记的基因顺序不一致（第 %s 位：收到 %s，期望 %s）。%s",
      label, if (is.na(bad)) "?" else as.character(bad),
      if (is.na(bad)) "?" else rownames(mat)[bad],
      if (is.na(bad)) "?" else expect[bad],
      "行名挂错说明矩阵可能被转置过。"
    ), call. = FALSE)
  }
  # 数值证据：每个基因的总计数。转置后这里会变成每个细胞的文库大小。
  if (!is.null(meta$gene_sums)) {
    gs <- as.numeric(meta$gene_sums)
    if (length(gs) != n_genes) {
      stop(sprintf("meta$gene_sums 长度 %d != n_genes %d", length(gs), n_genes),
           call. = FALSE)
    }
    got <- as.numeric(Matrix::rowSums(mat))
    bad <- which(abs(got - gs) > 1e-6)[1]
    if (!is.na(bad)) {
      stop(sprintf(
        paste0("%s 的每基因总计数对不上（第 %d 位基因 %s：收到 %.6g，期望 %.6g）。%s"),
        label, bad, expect[bad], got[bad], gs[bad],
        "这是**矩阵被转置**的数值证据 —— 行和变成了细胞的文库大小。"
      ), call. = FALSE)
    }
  }
  invisible(TRUE)
}

# ---- 读输入 -----------------------------------------------------------------
read_counts <- function(path, meta_path) {
  if (!file.exists(path)) stop("找不到输入矩阵: ", path, call. = FALSE)
  if (!file.exists(meta_path)) stop("找不到输入 meta: ", meta_path, call. = FALSE)
  meta <- jsonlite::fromJSON(meta_path)

  # 裸文本矩阵：第 1 行是基因名，其后每行一个基因。
  # `scan` 比 `read.csv` 快一个数量级，且不做任何列名改写。
  n_genes <- as.integer(meta$n_genes)
  n_cells <- as.integer(meta$n_cells)
  if (is.na(n_genes) || is.na(n_cells)) {
    stop("meta 里的 n_genes / n_cells 不是整数", call. = FALSE)
  }
  con <- file(path, "r")
  on.exit(close(con), add = TRUE)
  header <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1]]
  vals <- scan(con, what = double(), sep = "\t", quiet = TRUE,
               comment.char = "", na.strings = "NA")
  if (length(vals) != n_genes * n_cells) {
    stop(sprintf("矩阵元素个数 %d != n_genes * n_cells = %d",
                 length(vals), n_genes * n_cells), call. = FALSE)
  }
  # scan 按行优先填 —— 与写出来的顺序一致，所以 byrow = TRUE
  mat <- matrix(vals, nrow = n_genes, ncol = n_cells, byrow = TRUE)
  rownames(mat) <- header
  colnames(mat) <- sprintf("cell_%d", seq_len(n_cells))
  assert_orientation(mat, meta, "输入矩阵")
  if (anyNA(mat)) stop("输入矩阵里有 NA", call. = FALSE)
  if (any(mat < 0)) stop("输入矩阵里有负值 —— CPM 归一化会算出负的库大小", call. = FALSE)
  storage.mode(mat) <- "double"
  list(mat = Matrix::Matrix(mat, sparse = TRUE), meta = meta)
}

# ---- 主流程 -----------------------------------------------------------------
run_tenifold <- function(input, meta_path, targets_path, out_dir,
                         n_net, n_cells, n_comp, td_k, ma_n_dim, n_cores) {
  t0 <- Sys.time()

  for (pkg in c("scTenifoldKnk", "scTenifoldNet", "Matrix", "jsonlite")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("R 包缺失: ", pkg, " —— 请检查 CI 的 R 依赖安装步骤", call. = FALSE)
    }
  }
  if (!file.exists(targets_path)) stop("找不到候选基因表: ", targets_path, call. = FALSE)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  targets <- readLines(targets_path, warn = FALSE)
  targets <- trimws(targets)
  targets <- targets[nzchar(targets)]
  targets <- unique(targets)
  if (length(targets) == 0L) stop("候选基因表是空的", call. = FALSE)

  inp <- read_counts(input, meta_path)
  mat <- inp$mat
  cat(sprintf("[tenifold] 输入: %d genes x %d cells\n", nrow(mat), ncol(mat)))
  cat(sprintf("[tenifold] 候选基因 %d 个: %s\n", length(targets),
              paste(utils::head(targets, 20), collapse = ", ")))

  # 包在源码里会自己 stop，但在这里先报一次，报错信息里能带上"是谁不在"
  missing <- setdiff(targets, rownames(mat))
  if (length(missing) > 0L) {
    stop("以下候选基因不在输入矩阵里: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  cat(sprintf(
    "[tenifold] 参数: transcriptomeWide=TRUE qc=FALSE nNet=%d nCells=%d nComp=%d td_K=%d ma_nDim=%d nCores=%d\n",
    n_net, n_cells, n_comp, td_k, ma_n_dim, n_cores
  ))

  res <- scTenifoldKnk::scTenifoldKnk(
    countMatrix = mat,
    gKO = targets,
    transcriptomeWide = TRUE,
    qc = FALSE,                      # 见文件头"为什么 qc = FALSE"
    nc_nNet = n_net,
    nc_nCells = n_cells,
    nc_nComp = n_comp,
    td_K = td_k,
    ma_nDim = ma_n_dim,
    nCores = n_cores
  )

  pd <- res$perturbationDistances
  if (is.null(pd)) stop("scTenifoldKnk 没有返回 perturbationDistances", call. = FALSE)
  # 自检返回结构 —— 这是本脚本与包之间的契约，破了要立刻知道
  if (!identical(rownames(pd), targets)) {
    stop("perturbationDistances 的行名与候选基因顺序不一致", call. = FALSE)
  }
  if (!identical(colnames(pd), rownames(mat))) {
    stop("perturbationDistances 的列名与输入基因顺序不一致", call. = FALSE)
  }

  dist_csv <- file.path(out_dir, "tenifold_perturbation_distances.csv")
  utils::write.csv(pd, dist_csv, row.names = TRUE)

  # 网络本身也留一份"形状证据"：WT 网络的维度与非零边数。
  # **不落盘完整网络** —— 那是 nGenes x nGenes 的稀疏矩阵，落盘只是占地方；
  # 真正被下游用的是距离矩阵。
  wt <- res$tensorNetworks$WT
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  meta <- list(
    engine = "scTenifoldKnk",
    engine_version = as.character(utils::packageVersion("scTenifoldKnk")),
    scTenifoldNet_version = as.character(utils::packageVersion("scTenifoldNet")),
    r_version = R.version.string,
    transcriptome_wide = TRUE,
    qc = FALSE,
    qc_note = "01_qc.py 已做过 QC；这里不再滤，保证与一阶近似用同一批细胞",
    n_genes = nrow(mat),
    n_cells = ncol(mat),
    n_targets = length(targets),
    target_genes = targets,
    wt_network_dim = as.integer(dim(wt)),
    wt_network_nonzero = as.integer(sum(wt != 0)),
    distances_non_na = as.integer(sum(!is.na(pd))),
    distances_total = as.integer(length(pd)),
    params = list(n_net = n_net, n_cells = n_cells, n_comp = n_comp,
                  td_k = td_k, ma_n_dim = ma_n_dim, n_cores = n_cores),
    internal_seed = "包内硬编码 set.seed(1)（scTenifoldKnk.R L163/L172/L208/L231）",
    elapsed_sec = round(elapsed, 1),
    distances_csv = basename(dist_csv)
  )
  meta_path_out <- file.path(out_dir, "tenifold_meta.json")
  writeLines(jsonlite::toJSON(meta, auto_unbox = TRUE, pretty = TRUE,
                              null = "null"), meta_path_out)
  cat(sprintf("[tenifold] 完成: %d 个候选 x %d 个网络基因，耗时 %.1f s\n",
              length(targets), nrow(mat), elapsed))
  cat(sprintf("TENIFOLD_OK %s\n", dist_csv))
  invisible(meta)
}

# ---- 自测 -------------------------------------------------------------------
# 本地没有 R（R1：本地零执行），所以这段只能在 CI 里跑。它验三件事：
#   1. 包真的能装、能跑（不是"import 成功"就算）
#   2. 返回结构与我们的契约一致（行名 = 候选基因顺序，列名 = 输入基因顺序）
#   3. **方向断言真的能拦住转置矩阵** —— 这条最重要，因为写反了包本身不报错
run_selftest <- function() {
  cat("[selftest] 造合成计数矩阵（负二项，含 mt- 基因）...\n")
  set.seed(42)
  n_genes <- 60L
  n_cells <- 150L
  cnt <- matrix(stats::rnbinom(n_genes * n_cells, size = 20, prob = 0.7),
                nrow = n_genes, ncol = n_cells)
  rownames(cnt) <- c(sprintf("G%03d", seq_len(n_genes - 5L)),
                     sprintf("MT-%d", seq_len(5L)))
  colnames(cnt) <- sprintf("C%04d", seq_len(n_cells))
  storage.mode(cnt) <- "double"

  meta <- list(n_genes = n_genes, n_cells = n_cells, genes = rownames(cnt),
               gene_sums = as.numeric(Matrix::rowSums(cnt)))

  # --- 断言 1：正确方向必须通过 ---
  assert_orientation(cnt, meta, "selftest 正例")
  cat("[selftest] 方向断言：正例通过\n")

  # --- 断言 2：转置必须被拦住（这是本脚本存在的理由）---
  caught <- tryCatch({
    assert_orientation(t(cnt), meta, "selftest 转置")
    FALSE
  }, error = function(e) TRUE)
  if (!caught) {
    stop("方向断言没有拦住转置矩阵 —— 保护失效，写反了会静默算出垃圾网络",
         call. = FALSE)
  }
  cat("[selftest] 方向断言：转置被拦住\n")

  # --- 断言 3：方阵转置也要被拦住（行列数相同的退化情形）---
  # 只查维度的实现会在这里漏掉。这里把**行名与数值证据都抹平**，逼出最难的一种：
  # 方阵 + 基因名与细胞名同集合，只有 gene_sums 能发现。
  sq <- cnt[seq_len(40L), seq_len(40L)]
  sq_meta <- list(n_genes = 40L, n_cells = 40L, genes = rownames(sq),
                  gene_sums = as.numeric(Matrix::rowSums(sq)))
  caught_sq <- tryCatch({
    assert_orientation(t(sq), sq_meta, "selftest 方阵转置")
    FALSE
  }, error = function(e) TRUE)
  if (!caught_sq) {
    stop("方向断言没拦住**方阵**转置 —— 只查维度的实现在这里会漏", call. = FALSE)
  }
  cat("[selftest] 方向断言：方阵转置被拦住\n")

  # --- 断言 3b：方阵 + 行列名相同，只有数值证据能拦 ---
  # 把 sq 的列名改成与行名相同（模拟"基因名恰好也是细胞名"的退化输入），
  # 此时前三条判据全过，唯一能发现转置的就是 gene_sums。
  sq2 <- sq
  colnames(sq2) <- rownames(sq2)
  sq2_meta <- list(n_genes = 40L, n_cells = 40L, genes = rownames(sq2),
                   gene_sums = as.numeric(Matrix::rowSums(sq2)))
  caught_num <- tryCatch({
    assert_orientation(t(sq2), sq2_meta, "selftest 数值证据")
    FALSE
  }, error = function(e) TRUE)
  if (!caught_num) {
    stop("方向断言没拦住「行列名相同的方阵转置」—— 说明 gene_sums 数值复核没生效",
         call. = FALSE)
  }
  cat("[selftest] 方向断言：行列名相同时靠 gene_sums 拦住\n")

  # --- 断言 3c：gene_sums 被篡改也必须被拦住（证明它真在比对）---
  bad_meta <- meta
  bad_meta$gene_sums <- meta$gene_sums * 1.5
  caught_gs <- tryCatch({
    assert_orientation(cnt, bad_meta, "selftest 篡改 gene_sums")
    FALSE
  }, error = function(e) TRUE)
  if (!caught_gs) {
    stop("gene_sums 对不上却没报错 —— 数值复核形同虚设", call. = FALSE)
  }
  cat("[selftest] 方向断言：gene_sums 不符被拦住\n")

  # --- 断言 4：真的跑一遍 scTenifoldKnk ---
  targets <- c("G001", "G002", "G003")
  cat("[selftest] 跑 scTenifoldKnk（3 网络 / 80 细胞 / 60 基因）...\n")
  res <- scTenifoldKnk::scTenifoldKnk(
    countMatrix = Matrix::Matrix(cnt, sparse = TRUE),
    gKO = targets,
    transcriptomeWide = TRUE,
    qc = FALSE,
    nc_nNet = 3L, nc_nCells = 80L, nc_nComp = 3L,
    td_K = 3L, ma_nDim = 2L, nCores = 1L
  )
  pd <- res$perturbationDistances
  if (is.null(pd)) stop("selftest: 没有 perturbationDistances", call. = FALSE)
  if (!identical(dim(pd), c(3L, n_genes))) {
    stop("selftest: 距离矩阵维度是 ",
         paste(dim(pd), collapse = "x"), "，期望 3x", n_genes, call. = FALSE)
  }
  if (!identical(rownames(pd), targets)) {
    stop("selftest: 行名不是候选基因的原顺序", call. = FALSE)
  }
  if (!identical(colnames(pd), rownames(cnt))) {
    stop("selftest: 列名不是输入基因的原顺序", call. = FALSE)
  }
  if (anyNA(pd)) {
    stop("selftest: 距离矩阵里有 NA —— dRegulation 没有覆盖全部基因",
         call. = FALSE)
  }
  if (all(pd == 0)) {
    stop("selftest: 距离全为 0 —— 敲除没有产生任何扰动，链路可能断了",
         call. = FALSE)
  }
  cat(sprintf("[selftest] 距离矩阵 %d x %d，非零 %d 个，范围 [%.4f, %.4f]\n",
              nrow(pd), ncol(pd), sum(pd != 0), min(pd), max(pd)))

  # --- 断言 5：读写契约往返 ---
  # 这一段验的是**本脚本与 Python 之间的文件格式契约**，不是 scTenifoldKnk。
  # Python 侧写矩阵的格式必须与 read_counts 的读法一致 —— 格式对不上时
  # 最典型的失败是"行列数对但元素错位"，而它不会报错，只会算出错误的数。
  cat("[selftest] 验读写契约（裸矩阵往返）...\n")
  tmpdir <- tempfile("tenifold_selftest_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
  mat_path <- file.path(tmpdir, "matrix.txt")
  meta_path <- file.path(tmpdir, "meta.json")
  # 按 Python 的写法逐行写：第 1 行基因名，其后每行一个基因的细胞计数。
  # **必须逐行显式拼** —— `cat(matrix, sep="\n")` 走的是列优先展开，
  # 写出来会变成"每个细胞一行"，行列数恰好还对，是最难发现的那种错位。
  writeLines(c(
    paste(rownames(cnt), collapse = "\t"),
    vapply(seq_len(n_genes), function(i) {
      paste(format(cnt[i, ], scientific = FALSE, trim = TRUE), collapse = "\t")
    }, character(1))
  ), mat_path)
  writeLines(jsonlite::toJSON(
    list(n_genes = n_genes, n_cells = n_cells, genes = rownames(cnt),
         gene_sums = as.numeric(Matrix::rowSums(cnt))),
    auto_unbox = TRUE), meta_path)
  back <- read_counts(mat_path, meta_path)
  if (!identical(dim(back$mat), c(n_genes, n_cells))) {
    stop("selftest: 往返后维度不对: ", paste(dim(back$mat), collapse = "x"),
         call. = FALSE)
  }
  if (!identical(rownames(back$mat), rownames(cnt))) {
    stop("selftest: 往返后基因顺序变了", call. = FALSE)
  }
  if (max(abs(as.matrix(back$mat) - cnt)) > 1e-9) {
    stop("selftest: 往返后数值不一致 —— 元素错位（byrow 或分隔符有问题）",
         call. = FALSE)
  }
  cat("[selftest] 读写契约：往返逐元素一致\n")

  cat(sprintf("[selftest] 版本: scTenifoldKnk %s / scTenifoldNet %s / %s\n",
              as.character(utils::packageVersion("scTenifoldKnk")),
              as.character(utils::packageVersion("scTenifoldNet")),
              R.version.string))
  cat("SELFTEST OK\n")
  invisible(TRUE)
}

# ---- 入口 -------------------------------------------------------------------
main <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  args <- parse_args(argv)
  if (isTRUE(args$selftest)) {
    run_selftest()
    return(invisible(NULL))
  }
  for (k in c("input", "meta", "targets", "out_dir")) {
    if (is.null(args[[k]])) {
      stop("缺少必需参数 --", k, "（或加 --selftest 跑自测）", call. = FALSE)
    }
  }
  run_tenifold(
    input = args$input, meta_path = args$meta, targets_path = args$targets,
    out_dir = args$out_dir,
    n_net = as_int(args$n_net, "--n_net"),
    n_cells = as_int(args$n_cells, "--n_cells"),
    n_comp = as_int(args$n_comp, "--n_comp"),
    td_k = as_int(args$td_k, "--td_k"),
    ma_n_dim = as_int(args$ma_n_dim, "--ma_n_dim"),
    n_cores = as_int(args$n_cores, "--n_cores")
  )
  invisible(NULL)
}

main()
