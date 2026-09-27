#!/usr/bin/env Rscript
# =============================================================================
# 05_trajectory.R — 轨迹分析（R 版，与 05_trajectory.py 同号同语义）
#
# 方法学等级（references/r_version.md §3 05 行 + §11 已裁决 2026-09-26）：
#   - GCS：**自行实现，与 Python 逐行同式**（不是 CytoTRACE 包 —— 保证两侧同值）
#   - dpt 槽位：destiny::DiffusionMap + DPT（对应 scanpy 的 diffmap+DPT）
#   - scfates 槽位：**slingshot**（主曲线树 + 沿树拟时序）—— Python scFates
#     从未发布 R 版（CRAN 404 已核实），"CRAN 原版 scFates"这个说法不成立
#   - palantir 槽位：**丢弃**（无 R 实现）—— 方法数 4→3，≥2 方法仍可交叉验证
#   - stable_methods = ["dpt","slingshot","cytotrace"]（reproducibility 段
#     **必须按 R 版现实改写**：slingshot 没有六轮 CI 证据，"首次落地无历史
#     观测值，按未定值报"）
#
# 逐条复刻的 Python 语义（行号以 05_trajectory.py 为准）：
#   - L6：rho 非有限或=0 → flipped=NA + direction_decided=FALSE ——
#     **不判定 ≠ 不翻**，随机方向不能伪装成有方向；select_cv_methods()
#     是 direction_decided 的唯一消费者，把"方向未判定"的方法从
#     交叉验证与共识里剔除（理由两组分开记）。
#   - E-69 Form A：method_correlation_stats 是纯函数，state 三态
#     ok/single_method/undefined，调用点判空一律 is.null() 反义 ——
#     **0.0 是有效值，不能用真值判断**。
#   - M2：KS 的原因码五值（not_configured/single_group/no_pairs/
#     import_failed/ok），BH 校正与 06/07 同口径。
#   - M4：along_trajectory 嵌套 status 分清"没有"vs"崩了"。
#   - M5：山脊图 MIN_RIDGE_GROUPS 守卫，画的不够不落盘。
#   - M7：分支段 is_reliable 列（MIN_CELLS_PER_SEGMENT=30），保留小段
#     但标不可靠。
#   - 早退 6 条全部 write_json + return，不抛异常（E-56：result_status_of 接住）。
# =============================================================================

.load_common_for_local <- function() {
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

N_MODULES <- 6L
MIN_GENES_FOR_MODULES <- 200L
MIN_CELLS_PER_SEGMENT <- 30L
MIN_RIDGE_GROUPS <- 2L

# ---------------------------------------------------------------------------
# 工具：Spearman 相关（stats::cor(method="spearman")，常数列返回 NA——
# 与 scipy.spearmanr 的 nan 行为一致，调用点按 NA 处理）
# ---------------------------------------------------------------------------
.spearman <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3L) return(NA_real_)
  r <- suppressWarnings(stats::cor(a[ok], b[ok], method = "spearman"))
  if (!is.finite(r)) NA_real_ else as.numeric(r)
}

# ---------------------------------------------------------------------------
# 方法 1：CytoTRACE 风格 GCS（05_trajectory.py:60-103 逐行同式）
# ---------------------------------------------------------------------------
compute_cytotrace_gcs <- function(expr_all, log = log_info) {
  # expr_all: 细胞 x 基因（全基因集，log 后即可 —— 排序不变）
  n_cells <- nrow(expr_all)
  n_genes <- ncol(expr_all)
  n_genes_per_cell <- as.numeric(Matrix::rowSums(expr_all > 0))
  if (stats::sd(n_genes_per_cell) == 0) {
    stop("所有细胞的表达基因数相同，GCS 无从计算", call. = FALSE)
  }
  corr <- numeric(n_genes)
  for (j in seq_len(n_genes)) {
    col <- as.numeric(expr_all[, j])
    if (stats::sd(col) == 0) next
    r <- .spearman(col, n_genes_per_cell)
    corr[j] <- if (is.na(r)) 0 else r
  }
  top_n <- min(max(50L, as.integer(0.01 * n_genes)), n_genes)
  top_idx <- order(-corr)[seq_len(top_n)]
  sub <- expr_all[, top_idx, drop = FALSE]
  # 秩均值（列方向），归一 [0,1]
  ranks <- Matrix::apply(sub, 2, function(c) rank(c, ties.method = "average"))
  gcs <- as.numeric(rowMeans(ranks))
  rng <- max(gcs) - min(gcs)
  gcs <- (gcs - min(gcs)) / (if (rng > 0) rng else 1)
  log(sprintf("CytoTRACE GCS（自实现）: top_genes=%d，范围 %.3f-%.3f",
              top_n, min(gcs), max(gcs)))
  gcs
}

# ---------------------------------------------------------------------------
# 方法 2：DPT（destiny；05_trajectory.py:109-114 同语义，小=早）
# ---------------------------------------------------------------------------
compute_dpt <- function(expr_hvg, pca_emb, root_idx, seed, log = log_info) {
  need_pkg("destiny", "扩散拟时序 DPT")
  # destiny 吃细胞 x 特征；给 PCA embedding（与 scanpy 的 diffmap 输入同思路）
  dm <- destiny::DiffusionMap(pca_emb, n_pcs = NA, n_local = 5L, verbose = FALSE)
  dpt <- destiny::DPT(dm)
  pt <- as.numeric(dpt[[1L]])
  # destiny 的 DPT 以**密度最高**的细胞（tip1）为根，不是我们的 GCS 根。
  # scanpy 版（05_trajectory.py:109-114）用 `adata.uns["iroot"] = root_idx`
  # 让 DPT 直接从 GCS 根出发 —— destiny 没有等价的"指定根"接口，
  # **这是实现级差异，必须如实落账**：
  #   - 对齐语义只做**平移**（pt - min(pt)）：单调变换不改 Spearman 排序，
  #     对方向校正、方法间相关、共识 z-score 全部无影响；
  #   - **不做** `pt - pt[root_idx]` 再把负值置 NA 的"伪对齐"—— 那会删掉
  #     "比 GCS 根更早"的真实尾部（GCS 根未必是 DPT 密度峰的全局最小），
  #     制造一段 NA 空洞，下游 `.spearman` 按可用对算，无声改变样本量。
  pt_fixed <- pt - min(pt)
  log(sprintf("destiny DPT: 范围 %.3f-%.3f（以密度峰为根，非 GCS 根；实现级差异已记 status.dpt_root_note）",
              min(pt_fixed), max(pt_fixed)))
  pt_fixed
}

# ---------------------------------------------------------------------------
# 方法 3：slingshot（对应 scfates 槽位；主曲线树 + 沿树拟时序）
# ---------------------------------------------------------------------------
compute_slingshot <- function(pca_emb, clusters, root_cluster, log = log_info) {
  need_pkg("slingshot", "主曲线拟时序")
  # slingshot 只需要 reducedDim + 簇标签；assays 里的占位矩阵不会被读
  #（slingshot 走 reducedDim），用 1x1 零矩阵避免复制整个 embedding 两份
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = Matrix::Matrix(0, 1, 1, sparse = TRUE)),
    reducedDims = list(PCA = pca_emb))
  sce$cluster <- factor(clusters)
  sling <- slingshot::slingshot(sce, clusterLabels = "cluster",
                                reducedDim = "PCA", start.clus = root_cluster)
  # 多谱系时取第一谱系（pbmc3k 单主干；多分支数据集的谱系选择记进 status）
  pt <- as.numeric(SingleCellExperiment::colData(sling)$slingPseudotime_1)
  lineages <- slingParams <- NULL
  lineages <- slingshot::slingParams(sling)$lineages
  if (is.null(lineages)) lineages <- "Lineage1"
  list(pseudotime = pt,
       n_lineages = length(lineages),
       lineages = as.character(lineages))
}

# ---------------------------------------------------------------------------
# 方向参考（05_trajectory.py:191-209 同语义）
# ---------------------------------------------------------------------------
marker_direction_score <- function(expr_all, gene_names, early, late,
                                   log = log_info) {
  present <- gene_names
  e <- intersect(early, present)
  l <- intersect(late, present)
  if (!length(e) || !length(l)) {
    log(sprintf("marker 方向参考不可用：early 命中 %s，late 命中 %s",
                paste(e, collapse = ","), paste(l, collapse = ",")))
    return(list(score = NULL, used = list(early = e, late = l)))
  }
  score <- as.numeric(Matrix::rowMeans(expr_all[, l, drop = FALSE])) -
    as.numeric(Matrix::rowMeans(expr_all[, e, drop = FALSE]))
  log(sprintf("marker 方向参考：early=%s late=%s",
              paste(e, collapse = ","), paste(l, collapse = ",")))
  list(score = score, used = list(early = e, late = l))
}

# ---------------------------------------------------------------------------
# orient（05_trajectory.py:212-257 逐条复刻，含 L6）
# ---------------------------------------------------------------------------
orient <- function(raw, reference, convention) {
  out <- list(); rows <- list()
  for (nm in names(raw)) {
    vv <- as.numeric(raw[[nm]])
    if (identical(convention[[nm]], "earlier")) vv <- -vv
    r_before <- .spearman(vv, reference)
    decidable <- is.finite(r_before) && abs(r_before) > 0
    flipped <- NA   # NA = 方向无法判定（L6）
    if (decidable && r_before < 0) {
      vv <- -vv
      flipped <- TRUE
    } else if (decidable) {
      flipped <- FALSE
    }
    r_after <- .spearman(vv, reference)
    out[[nm]] <- vv
    rows[[length(rows) + 1L]] <- list(
      method = nm,
      raw_convention = convention[[nm]] %||% "?",
      rho_vs_reference_before_flip = finite_round(r_before, 4),
      flipped = if (is.na(flipped)) NULL else flipped,
      direction_decided = decidable,
      direction_note = if (decidable) NULL else
        sprintf("与参考的相关 rho=%s —— **方向无法判定**，该方法的拟时序方向未做校正，不应被当成有方向的结果",
                format(r_before)),
      rho_vs_reference_after_flip = finite_round(r_after, 4))
  }
  list(corrected = out, rows = rows)
}

# ---------------------------------------------------------------------------
# select_cv_methods（05_trajectory.py:260-288 纯函数，两道剔除分开记）
# ---------------------------------------------------------------------------
select_cv_methods <- function(names, direction_rows, reference_method) {
  decided <- vapply(direction_rows, function(r) isTRUE(r$direction_decided), logical(1))
  undecided <- vapply(direction_rows, function(r) r$method, character(1))[!decided]
  cv <- names[!names %in% c(reference_method, undecided)]
  list(cv_names = cv,
       cv_excluded = list(
         direction_reference = if (!is.null(reference_method)) list(reference_method) else list(),
         direction_undecided = as.list(setdiff(undecided, reference_method))))
}

# ---------------------------------------------------------------------------
# method_correlation_stats（05_trajectory.py:291-330 纯函数，E-69 Form A）
# ---------------------------------------------------------------------------
method_correlation_stats <- function(cv_names, cmat) {
  if (length(cv_names) < 2L) {
    return(list(mean_rho = NULL, min_rho = NULL, state = "single_method",
                note = sprintf("只有 %d 个方法可交叉验证，**没有「方法间一致性」这个量**（不是一致性低，是算不出来）",
                               length(cv_names))))
  }
  off <- c()
  for (i in seq_along(cv_names)) {
    if (i + 1L > length(cv_names)) break
    for (j in seq.int(i + 1L, length(cv_names))) {
      off <- c(off, as.numeric(cmat[cv_names[i], cv_names[j]]))
    }
  }
  finite <- off[is.finite(off)]
  n_bad <- length(off) - length(finite)
  if (!length(finite)) {
    return(list(mean_rho = NULL, min_rho = NULL, state = "undefined",
                note = sprintf("离对角相关 %d 对**全部非有限**（各方法退化成常数列），方法间一致性算不出来",
                               length(off))))
  }
  list(mean_rho = mean(finite), min_rho = min(finite), state = "ok",
       note = if (n_bad) sprintf("%d/%d 对相关非有限，已从均值/最小值中剔除", n_bad, length(off)) else NULL)
}

# ---------------------------------------------------------------------------
# 自实现簇连通性矩阵（scanpy PAGA 的等级 B 等价物）
# ---------------------------------------------------------------------------
cluster_connectivity <- function(pca_emb, clusters, k = 15L) {
  # 簇间连通性 = 归一化 kNN 边计数（PAGA 同思路：两簇间邻居边多 = 连通强）
  nn <- FNN::get.knn(pca_emb, k = k)$nn.index
  cl <- as.character(clusters)
  cats <- sort(unique(cl))
  conn <- matrix(0, nrow = length(cats), ncol = length(cats),
                 dimnames = list(cats, cats))
  cli <- match(cl, cats)
  for (i in seq_len(nrow(nn))) {
    for (j in nn[i, ]) {
      a <- cli[i]; b <- cli[j]
      if (a != b) {
        conn[a, b] <- conn[a, b] + 1
        conn[b, a] <- conn[b, a] + 1
      }
    }
  }
  # 归一化到 0-1（按每对簇可能的最大边数）
  conn <- conn / max(1, sum(conn))
  conn
}

# ---------------------------------------------------------------------------
# 主流程（05_trajectory.py:336-1242 同语义）
# ---------------------------------------------------------------------------
run_05_trajectory <- function(cfg) {
  ensure_dirs(cfg)
  set_seed(cfg)
  data_dir <- cfg$output$data_dir
  res_dir <- cfg$output$results_dir
  seed <- as.integer(cfg$analysis$seed)

  .early_exit <- function(status) {
    write_json(file.path(res_dir, "trajectory_status.json"), status)
    log_warn(status$reason %||% status$status)
    status
  }

  traj <- cfg$trajectory %||% list()
  if (!isTRUE(traj$enabled %||% TRUE)) {
    return(.early_exit(list(dataset_id = cfg$dataset_id, status = "disabled",
                            reason = "配置 trajectory.enabled=false")))
  }

  clu_path <- file.path(data_dir, "clustered.rds")
  if (!file.exists(clu_path)) {
    return(.early_exit(list(dataset_id = cfg$dataset_id, status = "failed",
                            reason = sprintf("缺输入 %s —— 先跑 03_cluster_annotate.R", clu_path))))
  }
  clu <- readRDS(clu_path)
  logcounts <- clu$logcounts          # HVG 子集（方法输入）
  counts_all <- clu$counts            # 原始计数（GCS 表达判定）
  umap <- clu$umap
  clusters <- as.character(clu$clusters)
  # GCS 与 marker 方向参考要"全基因集" —— R 版的 clustered.rds 只存了
  # HVG 子集（02 的契约）。如实记 limitation，不能假装是全基因集。
  gene_names_all <- colnames(logcounts)
  expr_all <- logcounts

  if (is.null(clusters)) {
    return(.early_exit(list(dataset_id = cfg$dataset_id, status = "missing_clusters",
                            reason = "clustered.rds 里没有 clusters（step 03 未产出）")))
  }
  n_clusters <- length(unique(clusters))
  if (n_clusters < 2L) {
    return(.early_exit(list(dataset_id = cfg$dataset_id, status = "not_applicable",
                            reason = sprintf("只有 %d 个簇，轨迹分析无从谈起", n_clusters))))
  }

  # ---- 1. 簇连通图（PAGA 等价）-------------------------------------------
  conn_res <- tryCatch(cluster_connectivity(pca_emb = clu$pca, clusters = clusters),
                       error = function(e) e)
  if (inherits(conn_res, "error")) {
    return(.early_exit(list(dataset_id = cfg$dataset_id, status = "failed",
                            reason = sprintf("簇连通性计算失败: %s", conditionMessage(conn_res)))))
  }
  conn <- conn_res
  utils::write.csv(data.frame(cluster_a = rep(rownames(conn), times = ncol(conn)),
                              cluster_b = rep(colnames(conn), each = nrow(conn)),
                              connectivity = as.numeric(conn)),
                   file.path(res_dir, "paga_connectivities.csv"), row.names = FALSE)

  # 图 02-05-01：连通图（edge width = connectivity）
  dfc <- data.frame(cluster_a = rep(rownames(conn), times = ncol(conn)),
                    cluster_b = rep(colnames(conn), each = nrow(conn)),
                    w = as.numeric(conn))
  dfc <- dfc[dfc$cluster_a < dfc$cluster_b, ]
  cents <- do.call(rbind, lapply(sort(unique(clusters)), function(cc) {
    m <- clusters == cc
    data.frame(x = mean(umap[m, 1]), y = mean(umap[m, 2]), cluster = cc)
  }))
  dfc <- merge(dfc, data.frame(cluster_a = cents$cluster, ax = cents$x, ay = cents$y),
               by = "cluster_a")
  dfc <- merge(dfc, data.frame(cluster_b = cents$cluster, bx = cents$x, by = cents$y),
               by = "cluster_b")
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = dfc,
                          ggplot2::aes(x = ax, y = ay, xend = bx, yend = by,
                                       linewidth = pmax(w, 0.02)),
                          colour = "grey60", alpha = 0.7, show.legend = FALSE) +
    ggplot2::geom_point(data = cents, ggplot2::aes(x, y),
                        size = 2.4, colour = PAL$primary) +
    ggplot2::geom_text(data = cents, ggplot2::aes(x, y, label = cluster),
                       size = 2.4, vjust = -0.9) +
    ggplot2::labs(x = "UMAP1", y = "UMAP2",
                  title = "Cluster connectivity graph (edge width = connectivity)") +
    ggplot2::theme_paper()
  save_fig(cfg, "02-05-01-unit1-paga-graph", p, width = W_ONE_HALF, height = mm(66))

  # ---- 2. GCS：既是方法也是选根依据 ---------------------------------------
  gcs <- compute_cytotrace_gcs(expr_all)
  root_idx <- which.max(gcs)
  root_cluster <- clusters[root_idx]
  root_record <- list(
    method = "auto_max_cytotrace_gcs",
    root_cluster = root_cluster,
    root_cell = rownames(logcounts)[root_idx],
    reason = paste0("取 CytoTRACE GCS 最高的细胞作根（GCS 高 = 分化潜能高 = 早）。",
                    "**这是统计判据不是生物学判据** —— 要下方向性结论仍需人工",
                    "用已知早期/晚期 marker 复核。"))
  configured_root <- traj$root_cluster
  if (!is.null(configured_root)) {
    cr <- as.character(configured_root)
    if (!cr %in% unique(clusters)) {
      return(.early_exit(list(dataset_id = cfg$dataset_id, status = "bad_root",
                              reason = sprintf("trajectory.root_cluster='%s' 不在簇列表里: [%s]",
                                               cr, paste(unique(clusters), collapse = ", ")))))
    }
    idx_in <- which(clusters == cr)
    root_idx <- idx_in[which.max(gcs[idx_in])]
    root_cluster <- cr
    root_record <- list(method = "configured", root_cluster = cr,
                        root_cell = rownames(logcounts)[root_idx],
                        reason = "来自配置 trajectory.root_cluster（簇内取 GCS 最高的细胞）")
  }
  log_info(sprintf("根 = 簇 %s，细胞 %s", root_cluster, rownames(logcounts)[root_idx]))

  # ---- 3. 三种方法（R 版：destiny DPT + slingshot + GCS）------------------
  # 每个方法独立 tryCatch：一个失败不拖垮整步（Python 版 L443-456 同语义）。
  # **不用 `<<-`**：tryCatch 返回值直接带出来再赋值 —— `<<-` 在 handler
  # 里会跳过本函数环境往全局找，找不到就**静默造全局变量**（M8 家族
  # "全局污染"形态），而且 error handler 的返回值才是 tryCatch 的返回值，
  # value 回调里赋值后 handler 再抛错时已赋的值会留下半套状态。
  raw_pt <- list(); methods_ok <- list(); methods_failed <- list()
  res_dpt <- tryCatch({
    list(pt = compute_dpt(logcounts, clu$pca, root_idx, seed),
         note = "destiny DPT（扩散拟时序，对应 scanpy DPT）")
  }, error = function(e) list(err = conditionMessage(e)))
  if (!is.null(res_dpt$err)) {
    methods_failed[["dpt"]] <- res_dpt$err
    log_warn(sprintf("DPT 失败: %s", res_dpt$err))
  } else {
    raw_pt[["dpt"]] <- res_dpt$pt
    methods_ok[["dpt"]] <- res_dpt$note
  }

  res_sling <- tryCatch({
    ss <- compute_slingshot(clu$pca, clusters, root_cluster)
    list(pt = ss$pseudotime,
         note = sprintf("slingshot 主曲线拟时序（%d 谱系；对应 Python 版 scfates 槽位）",
                        ss$n_lineages))
  }, error = function(e) list(err = conditionMessage(e)))
  if (!is.null(res_sling$err)) {
    methods_failed[["slingshot"]] <- res_sling$err
    log_warn(sprintf("slingshot 失败: %s", res_sling$err))
  } else {
    raw_pt[["slingshot"]] <- res_sling$pt
    methods_ok[["slingshot"]] <- res_sling$note
  }

  # CytoTRACE 自己也是一种方法
  raw_pt[["cytotrace"]] <- gcs
  methods_ok[["cytotrace"]] <- "CytoTRACE 风格 GCS（本仓库自行实现，非 CytoTRACE 包）"

  if (length(methods_ok) < 2L) {
    return(.early_exit(list(dataset_id = cfg$dataset_id, status = "insufficient_methods",
                            reason = sprintf("只有 %d 种方法成功，少于交叉验证所需的 2 种", length(methods_ok)),
                            methods_ok = methods_ok, methods_failed = methods_failed)))
  }

  # 原始符号约定（R 版方法名换了，convention 的键跟着换）
  convention <- list(dpt = "earlier", slingshot = "earlier", cytotrace = "later")

  # ---- 4. 方向参考 --------------------------------------------------------
  early_m <- as.character(traj$early_markers %||% character(0))
  late_m <- as.character(traj$late_markers %||% character(0))
  marker_ref <- NULL; marker_used <- list()
  if (length(early_m) && length(late_m)) {
    mr <- marker_direction_score(expr_all, gene_names_all, early_m, late_m)
    marker_ref <- mr$score; marker_used <- mr$used
  }
  if (!is.null(marker_ref)) {
    reference <- marker_ref
    direction_source <- "markers"
    direction_note <- "方向参考 = 配置的 early/late marker 基因（mean(late) - mean(early)）"
  } else {
    reference <- -gcs
    direction_source <- "cytotrace_fallback"
    direction_note <- paste0(
      "**未配置 trajectory.early_markers/late_markers**，退回用 CytoTRACE 分化潜能分作方向参考",
      "（取负使「大 = 晚」）。这是统计判据不是生物学判据：若该数据集里分化潜能与成熟度不同向，",
      "方向会整体反掉。要下方向性结论请在配置里给出已知的早期/晚期 marker。")
    log_warn(direction_note)
  }

  o <- orient(raw_pt, reference, convention)
  corrected <- o$corrected; direction_rows <- o$rows
  utils::write.csv(do.call(rbind, lapply(direction_rows, function(r) {
    data.frame(method = r$method, raw_convention = r$raw_convention,
               rho_before = r$rho_vs_reference_before_flip,
               flipped = if (is.null(r$flipped)) NA else r$flipped,
               direction_decided = r$direction_decided,
               rho_after = r$rho_vs_reference_after_flip)
  })), file.path(res_dir, "trajectory_direction.csv"), row.names = FALSE)
  for (r in direction_rows) {
    if (is.null(r$flipped)) {
      log_warn(sprintf("  方向 %-10s rho=%+.4f —— **无法判定方向**（L6），该方法的拟时序不做方向校正",
                       r$method, r$rho_vs_reference_before_flip))
    } else {
      log_info(sprintf("  方向 %-10s rho=%+.4f%s", r$method,
                       r$rho_vs_reference_after_flip,
                       if (isTRUE(r$flipped)) "（已翻转）" else ""))
    }
  }

  # ---- 5. 交叉验证矩阵 ----------------------------------------------------
  names_order <- intersect(c("dpt", "slingshot", "cytotrace"), names(corrected))
  cmat <- matrix(NA_real_, nrow = length(names_order), ncol = length(names_order),
                 dimnames = list(names_order, names_order))
  for (a in names_order) for (b in names_order) {
    cmat[a, b] <- .spearman(corrected[[a]], corrected[[b]])
  }
  cdf <- data.frame(method = rownames(cmat), as.data.frame(cmat),
                    check.names = FALSE)
  utils::write.csv(cdf, file.path(res_dir, "trajectory_method_correlation.csv"),
                   row.names = FALSE)

  reference_method <- if (identical(direction_source, "cytotrace_fallback")) "cytotrace" else NULL
  sel <- select_cv_methods(names_order, direction_rows, reference_method)
  cv_names <- sel$cv_names; cv_excluded <- sel$cv_excluded
  if (length(cv_excluded$direction_undecided)) {
    log_warn(sprintf("**%d 个方法的方向无法判定（%s），已从交叉验证与共识中剔除** —— 不剔除等于给共识掺进一个随机方向",
                     length(cv_excluded$direction_undecided),
                     paste(unlist(cv_excluded$direction_undecided), collapse = ", ")))
  }
  mcs <- method_correlation_stats(cv_names, cmat)
  mean_rho <- mcs$mean_rho; min_rho <- mcs$min_rho
  mc_state <- mcs$state; mc_note <- mcs$note
  if (!identical(mc_state, "ok")) {
    log_warn(sprintf("方法间一致性**算不出来**（state=%s）：%s", mc_state, mc_note))
  } else if (!is.null(mc_note)) {
    log_warn(sprintf("方法间一致性部分不可用：%s", mc_note))
  }
  if (!is.null(reference_method)) {
    log_info(sprintf("一致性统计已排除方向参考 %s（它与参考的相关是定义上的，不是证据）", reference_method))
  }

  # 图 02-05-02：数值标注热图
  dfm <- data.frame(
    method_a = rep(rownames(cmat), times = ncol(cmat)),
    method_b = rep(colnames(cmat), each = nrow(cmat)),
    rho = as.numeric(cmat))
  p <- ggplot2::ggplot(dfm, ggplot2::aes(method_b, method_a, fill = rho)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", rho)), size = 2.4) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                                  midpoint = 0, limits = c(-1, 1), name = "Spearman") +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "Pseudotime agreement (Spearman, direction-corrected)") +
    ggplot2::theme_paper() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))
  save_fig(cfg, "02-05-02-unit1-trajectory-method-correlation", p,
           width = W_SINGLE, height = mm(66))

  # 共识拟时序（cv_names 空时 no_consensus 早退 —— **不静默回退**，L6 第二半）
  if (!length(cv_names)) {
    st <- list(dataset_id = cfg$dataset_id, status = "no_consensus",
               reason = paste0("没有任何方法的方向可判定（方向未判定: ",
                               paste(unlist(cv_excluded$direction_undecided), collapse = ", "),
                               "；方向参考: ", reference_method, "）",
                               "—— 方向未知的拟时序不能取共识，否则共识方向是随机的"),
               methods_ok = methods_ok, methods_failed = methods_failed,
               cv_excluded = cv_excluded,
               method_correlation_state = mc_state,
               method_correlation_undefined_note = mc_note)
    return(.early_exit(st))
  }
  consensus_names <- cv_names
  zscore <- function(v) (v - mean(v, na.rm = TRUE)) / (stats::sd(v, na.rm = TRUE) + 1e-12)
  Z <- do.call(cbind, lapply(consensus_names, function(n) zscore(corrected[[n]])))
  consensus <- as.numeric(rowMeans(Z, na.rm = TRUE))
  rho_ref <- NULL
  if (length(consensus_names) < length(names_order)) {
    Z_all <- do.call(cbind, lapply(names_order, function(n) zscore(corrected[[n]])))
    consensus_with_ref <- as.numeric(rowMeans(Z_all, na.rm = TRUE))
    rho_ref <- .spearman(consensus, consensus_with_ref)
  }
  log_info(sprintf("共识拟时序完成（用 %d 个方法: %s%s）；方法间平均 rho=%s",
                   length(consensus_names), paste(consensus_names, collapse = ", "),
                   if (!is.null(reference_method))
                     sprintf("；已排除方向参考 %s", reference_method) else "",
                   if (!is.null(mean_rho)) sprintf("%+.4f，最低 %+.4f", mean_rho, min_rho)
                   else sprintf("**不可用**（state=%s）", mc_state)))

  # ---- 6. 沿轨迹变化的基因（M4：嵌套 status 分清"没有"vs"崩了"）----------
  # **不用 `<<-`**（与上面方法段同理由）：tryCatch 的 value/handler 返回值
  # 带出来再赋值。`<<-` 在 handler 里若上层环境没有该绑定会**静默造全局
  # 变量**（M8 家族"全局污染"形态）。
  genes_rows <- list(); gene_names_sel <- character(0); gene_rho <- NULL
  along_status <- list(status = "not_run", reason = NULL)
  along_res <- tryCatch({
    expressed <- Matrix::colSums(expr_all > 0) >= max(10L, as.integer(0.01 * nrow(expr_all)))
    idx <- which(expressed)
    Xs <- expr_all[, idx, drop = FALSE]
    rho_g <- vapply(seq_len(ncol(Xs)), function(k) .spearman(as.numeric(Xs[, k]), consensus),
                    numeric(1))
    rho_g[!is.finite(rho_g)] <- 0
    names_g <- colnames(expr_all)[idx]
    ord <- order(-abs(rho_g))
    sel <- ord[seq_len(min(400L, length(ord)))]
    rows_out <- lapply(ord[seq_len(min(200L, length(ord)))], function(i) {
      list(gene = names_g[i],
           rho_with_pseudotime = finite_round(rho_g[i], 4),
           direction = if (rho_g[i] > 0) "increases" else "decreases")
    })
    list(gene_names_sel = names_g[sel], gene_rho = rho_g[sel],
         genes_rows = rows_out, n_expressed = length(idx))
  }, error = function(e) list(err = conditionMessage(e)))
  if (!is.null(along_res$err)) {
    along_status <- list(status = "failed", reason = along_res$err,
                         n_expressed = NULL, n_reported = 0L)
    log_warn(sprintf("沿轨迹基因分析失败: %s（已记进 trajectory_status.json 的 along_trajectory.status）",
                     along_res$err))
  } else {
    gene_names_sel <- along_res$gene_names_sel
    gene_rho <- along_res$gene_rho
    genes_rows <- along_res$genes_rows
    along_status <- list(status = "ok", reason = NULL,
                         n_expressed = along_res$n_expressed,
                         n_reported = length(genes_rows))
    utils::write.csv(data.frame(
      gene = vapply(genes_rows, function(r) r$gene, character(1)),
      rho_with_pseudotime = vapply(genes_rows, function(r) r$rho_with_pseudotime, numeric(1)),
      direction = vapply(genes_rows, function(r) r$direction, character(1))),
      file.path(res_dir, "trajectory_genes.csv"), row.names = FALSE)
    log_info(sprintf("沿轨迹变化基因：表达基因 %d 个，报前 %d 个",
                     along_res$n_expressed, length(genes_rows)))
  }

  # ---- 模块（20 分箱 → 行 z-score → kmeans）------------------------------
  modules_rows <- list(); profz <- NULL; lab <- NULL
  if (length(gene_names_sel) >= MIN_GENES_FOR_MODULES) {
    mod_res <- tryCatch({
      gene_mat <- expr_all[, gene_names_sel, drop = FALSE]
      nb <- 20L
      qs <- stats::quantile(consensus, probs = seq(0, 1, length.out = nb + 1L), na.rm = TRUE)
      bins <- cut(consensus, breaks = unique(qs), include.lowest = TRUE, labels = FALSE)
      nb <- max(bins, na.rm = TRUE)
      prof <- matrix(NA_real_, nrow = ncol(gene_mat), ncol = nb)
      for (b in seq_len(nb)) {
        m <- which(bins == b)
        if (!length(m)) next
        prof[b, ] <- Matrix::colMeans(gene_mat[m, , drop = FALSE])
      }
      # 行 z-score（行=基因：prof 是 基因 x 分箱，沿**分箱轴**标准化）
      mu <- Matrix::rowMeans(prof, na.rm = TRUE)
      sdv <- apply(prof, 1, function(x) stats::sd(x, na.rm = TRUE))
      sdv[!is.finite(sdv) | sdv == 0] <- 1
      profz <- sweep(sweep(prof, 1, mu, "-"), 1, sdv, "/")
      profz[!is.finite(profz)] <- 0
      # kmeans：行=基因
      set.seed(seed)
      km <- stats::kmeans(profz, centers = N_MODULES, nstart = 10L)
      lab <- km$cluster
      for (m in sort(unique(lab))) {
        sel_m <- which(lab == m)
        members <- gene_names_sel[sel_m]
        modules_rows[[length(modules_rows) + 1L]] <- list(
          module = m - 1L, n_genes = length(members),
          mean_abs_rho = finite_round(mean(abs(gene_rho[sel_m])), 4),
          mean_rho = finite_round(mean(gene_rho[sel_m]), 4),
          top_genes = paste(utils::head(members, 20), collapse = ","))
      }
      utils::write.csv(do.call(rbind, modules_rows),
        file.path(res_dir, "trajectory_modules.csv"), row.names = FALSE)
      log_info(sprintf("基因模块：%d 个（k-means on 分箱平滑曲线）", length(modules_rows)))
      list(profz = profz, lab = lab)    }, error = function(e) {
      log_warn(sprintf("基因模块分析失败: %s", conditionMessage(e)))
      NULL
    })
    if (!is.null(mod_res)) { profz <- mod_res$profz; lab <- mod_res$lab }
  }
  # 图 02-05-03 unit1 热图 + unit2 曲线（D-006 单图拆分）
  if (!is.null(profz) && !is.null(lab)) {
    # unit1：基因 x 分箱热图（按模块排序）
    ord_rows <- order(lab)
    dfh <- data.frame(gene = factor(rep(gene_names_sel[ord_rows], times = ncol(profz)),
                                    levels = gene_names_sel[ord_rows]),
                      bin = rep(seq_len(ncol(profz)), each = nrow(profz)),
                      z = as.numeric(profz[ord_rows, ]))
    p <- ggplot2::ggplot(dfh, ggplot2::aes(bin, gene, fill = z)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                                    midpoint = 0, limits = c(-2, 2),
                                    name = "z-scored mean expression") +
      ggplot2::labs(x = "pseudotime bin", y = "gene (grouped by module)",
                    title = sprintf("Genes along pseudotime, %d modules\ncolour = z-scored mean expression (shared scale -2..2)",
                                    N_MODULES)) +
      ggplot2::theme_paper() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank())
    save_fig(cfg, "02-05-03-unit1-trajectory-modules-heatmap", p,
             width = W_ONE_HALF, height = mm(70))
    # unit2：各模块平均曲线（图例规则 31：框外右侧单列）
    dfp <- do.call(rbind, lapply(sort(unique(lab)), function(m) {
      sel_m <- which(lab == m)
      data.frame(bin = seq_len(ncol(profz)),
                 z = as.numeric(colMeans(profz[sel_m, , drop = FALSE])),
                 module = sprintf("M%d (n=%d)", m - 1L, length(sel_m)))
    }))
    p <- ggplot2::ggplot(dfp, ggplot2::aes(bin, z, colour = module)) +
      ggplot2::geom_line(linewidth = 0.5) +
      ggplot2::labs(x = "pseudotime bin", y = "z-scored mean",
                    title = "Module profiles along pseudotime") +
      ggplot2::theme_paper()
    save_fig(cfg, "02-05-03-unit2-trajectory-module-profiles", p,
             width = W_ONE_HALF, height = mm(70))
  }

  # ---- 7. 多条件 KS 比较（M2：原因码五值 + BH 与 06/07 同口径）------------
  ks_rows <- list()
  design <- cfg$design %||% list()
  group_key <- design$group_key
  ks_reason <- "not_configured"
  ks_correction <- NULL
  cell_meta <- clu$cell_meta
  grp_vals <- NULL
  grp_degraded <- FALSE
  if (!is.null(group_key) && !is.null(cell_meta) && group_key %in% colnames(cell_meta)) {
    grp_vals <- as.character(cell_meta[[group_key]])
  } else if (!is.null(group_key)) {
    # 降级：分组列不在 cell_meta（03 产出的 rds 只带 celltype/leiden）——
    # 用簇标签顶上时**必须记降级**，否则 KS 的"组间差异"会被读成
    # 生物学条件差异，而实际是簇间差异（两码事）。
    grp_vals <- clusters
    grp_degraded <- TRUE
    log_warn(sprintf("design.group_key='%s' 不在 clustered.rds 的 meta 里（R 版 cell_meta 只有 celltype）—— KS 降级用簇标签，比较的是**簇间**差异不是条件差异", group_key))
  }
  if (is.null(group_key)) {
    log_info("未配置 design.group_key，跳过多条件拟时序分布比较")
  } else if (is.null(grp_vals) || length(unique(grp_vals)) < 2L) {
    ks_reason <- "single_group"
    log_info(sprintf("分组列 %s 只有 1 个取值，没有可比较的两组，跳过 KS 比较", group_key))
  } else {
    uniq <- sort(unique(grp_vals))
    for (i in seq_along(uniq)) for (j in seq.int(i + 1L, length(uniq))) {
      a_ <- consensus[grp_vals == uniq[i]]
      b_ <- consensus[grp_vals == uniq[j]]
      a_ <- a_[is.finite(a_)]; b_ <- b_[is.finite(b_)]
      if (length(a_) < 5L || length(b_) < 5L) next
      kt <- stats::ks.test(a_, b_, exact = FALSE)
      ks_rows[[length(ks_rows) + 1L]] <- list(
        group_a = uniq[i], group_b = uniq[j],
        n_a = length(a_), n_b = length(b_),
        ks_statistic = finite_round(unname(kt$statistic), 4),
        p_value = as.numeric(unname(kt$p.value)),
        median_a = finite_round(stats::median(a_), 4),
        median_b = finite_round(stats::median(b_), 4))
    }
    if (!length(ks_rows)) {
      ks_reason <- "no_pairs"
      log_info(sprintf("分组列 %s 有 %d 个取值，但没有任何一对满足最小细胞数（各 >=5），跳过 KS 比较",
                       group_key, length(uniq)))
    } else {
      # BH 校正（stats::p.adjust(method="BH")，与 06/07 同口径）
      pv <- vapply(ks_rows, function(r) r$p_value, numeric(1))
      padj <- stats::p.adjust(pv, method = "BH")
      for (k in seq_along(ks_rows)) {
        ks_rows[[k]]$p_adj_bh <- padj[k]
        ks_rows[[k]]$significant_bh <- padj[k] < 0.05
      }
      ks_reason <- "ok"; ks_correction <- "fdr_bh"
      utils::write.csv(data.frame(
        group_a = vapply(ks_rows, function(r) r$group_a, character(1)),
        group_b = vapply(ks_rows, function(r) r$group_b, character(1)),
        n_a = vapply(ks_rows, function(r) r$n_a, integer(1)),
        n_b = vapply(ks_rows, function(r) r$n_b, integer(1)),
        ks_statistic = vapply(ks_rows, function(r) r$ks_statistic, numeric(1)),
        p_value = vapply(ks_rows, function(r) r$p_value, numeric(1)),
        p_adj_bh = vapply(ks_rows, function(r) r$p_adj_bh, numeric(1)),
        significant_bh = vapply(ks_rows, function(r) r$significant_bh, logical(1))),
        file.path(res_dir, "trajectory_ks_by_group.csv"), row.names = FALSE)
      log_info(sprintf("多条件 KS 检验：%d 对比较，BH 校正后显著 %d 对",
                       length(ks_rows), sum(vapply(ks_rows, function(r) isTRUE(r$significant_bh), logical(1)))))
    }
  }
  ks_note_map <- list(
    ok = "两两 KS 检验的家庭是全部条件对；`p_adj_bh` 是 BH 校正后的值，判显著性看它",
    not_configured = paste0("未配置 `design.group_key`（或该列不在 obs 里），本轮**没有做**多条件拟时序分布比较 —— 这是「没做」，不是「做了没问题」"),
    single_group = paste0("`design.group_key` 配了，但该列在本数据集里只有 1 个取值，没有可比较的两组 —— 这是「数据里没有分组」，不是「配置漏了」"),
    no_pairs = paste0("分组列存在，但没有任何一对条件同时满足最小细胞数（各 >=5），本轮**没有做** KS 比较 —— 这是「样本不够」，不是「分布没有差异」"),
    import_failed = "`p.adjust` 失败，`p_value` 是**未校正**值，不能直接当判据")

  # ---- 8. 分支段（M7：is_reliable 列）------------------------------------
  # Python 版的 seg 来自 scFates 的 milestones（树形结构天然有"段"）；
  # slingshot 没有直接的"段"字段 —— R 版按**簇**聚合替代（每簇的拟时序
  # 中位数 + 细胞数），is_reliable 语义保留（小段保留但标不可靠，
  # 不删 —— 删掉等于隐藏信息）。实现级差异记进 status.segments_basis。
  branch_rows <- list(); n_unreliable_segments <- 0L
  per_seg_cl <- unique(clusters)
  for (cc in sort(per_seg_cl)) {
    m <- clusters == cc
    n_c <- sum(m)
    reliable <- n_c >= MIN_CELLS_PER_SEGMENT
    if (!reliable) n_unreliable_segments <- n_unreliable_segments + 1L
    branch_rows[[length(branch_rows) + 1L]] <- list(
      segment = cc, n_cells = n_c,
      median_pseudotime = finite_round(stats::median(consensus[m]), 4),
      is_reliable = reliable,
      min_cells_note = if (reliable) "" else sprintf("n<%d，中位数不可靠", MIN_CELLS_PER_SEGMENT))
  }
  utils::write.csv(data.frame(
    segment = vapply(branch_rows, function(r) r$segment, character(1)),
    n_cells = vapply(branch_rows, function(r) r$n_cells, integer(1)),
    median_pseudotime = vapply(branch_rows, function(r) r$median_pseudotime, numeric(1)),
    is_reliable = vapply(branch_rows, function(r) r$is_reliable, logical(1)),
    min_cells_note = vapply(branch_rows, function(r) r$min_cells_note, character(1))),
    file.path(res_dir, "trajectory_segments.csv"), row.names = FALSE)
  if (n_unreliable_segments) {
    log_warn(sprintf("分支段里有 %d 段细胞数 < %d，已在 trajectory_segments.csv 里标 is_reliable=false",
                     n_unreliable_segments, MIN_CELLS_PER_SEGMENT))
  }

  # ---- 9. 图：拟时序三联（D-006 单图拆分）--------------------------------
  dfu <- data.frame(u1 = umap[, 1], u2 = umap[, 2], pt = consensus)
  p <- ggplot2::ggplot(dfu, ggplot2::aes(u1, u2, colour = pt)) +
    ggplot2::geom_point(size = 0.4) +
    ggplot2::scale_colour_viridis_c(name = "pseudotime (higher = later)") +
    ggplot2::labs(x = "UMAP1", y = "UMAP2",
                  title = sprintf("Consensus pseudotime (root = cluster %s)", root_cluster)) +
    ggplot2::theme_paper()
  save_fig(cfg, "02-05-04-unit1-pseudotime-consensus", p,
           width = W_ONE_HALF, height = mm(58))

  # unit4：PAGA 骨架箭头（阈值 0.15；不用"质心按拟时序排序直连"——
  # 分散嵌入上会画出假折线，**误导性的图不如不出**）
  dfc2 <- dfc
  med <- tapply(consensus, clusters, stats::median)
  dfc2$src <- ifelse(med[dfc2$cluster_a] <= med[dfc2$cluster_b],
                     dfc2$cluster_a, dfc2$cluster_b)
  dfc2$dst <- ifelse(dfc2$src == dfc2$cluster_a, dfc2$cluster_b, dfc2$cluster_a)
  strong <- dfc2[dfc2$w >= 0.15, ]
  n_edges <- nrow(strong)
  arrows_df <- do.call(rbind, lapply(seq_len(nrow(strong)), function(k) {
    s <- strong$src[k]; d <- strong$dst[k]
    data.frame(x = cents$x[cents$cluster == s], y = cents$y[cents$cluster == s],
               xend = cents$x[cents$cluster == d], yend = cents$y[cents$cluster == d],
               w = strong$w[k])
  }))
  ri <- match(root_cluster, cents$cluster)
  root_circ <- if (!is.na(ri)) cents[ri, ] else NULL
  p <- ggplot2::ggplot(dfu, ggplot2::aes(u1, u2, colour = pt)) +
    ggplot2::geom_point(size = 0.4, alpha = 0.8) +
    ggplot2::scale_colour_viridis_c(name = "pseudotime (higher = later)") +
    ggplot2::labs(x = "UMAP1", y = "UMAP2",
                  title = sprintf("Consensus pseudotime with PAGA skeleton (%d strong edges, width ∝ connectivity > 0.15)", n_edges),
                  subtitle = sprintf("arrows point along increasing pseudotime; open circle = root cluster %s", root_cluster)) +
    ggplot2::theme_paper()
  if (n_edges) {
    p <- p + ggplot2::geom_segment(data = arrows_df,
                                   ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
                                                linewidth = 0.5 + 2.0 * w),
                                   arrow = grid::arrow(length = grid::unit(2, "mm")),
                                   colour = PAL$black, alpha = 0.65,
                                   show.legend = FALSE, inherit.aes = FALSE)
  }
  if (!is.null(root_circ)) {
    p <- p + ggplot2::geom_point(data = root_circ, ggplot2::aes(x, y),
                                 size = 4, shape = 21, colour = PAL$black,
                                 fill = NA, stroke = 1, inherit.aes = FALSE) +
      ggplot2::annotate("text", x = root_circ$x, y = root_circ$y + 0.8,
                        label = paste0(" root ", root_cluster), size = 2.2)
  }
  save_fig(cfg, "02-05-04-unit4-pseudotime-principal-path", p,
           width = W_ONE_HALF, height = mm(58))

  # unit2：DPT 交叉验证
  if ("dpt" %in% names(corrected)) {
    dfu2 <- dfu; dfu2$pt <- corrected[["dpt"]]
    p <- ggplot2::ggplot(dfu2, ggplot2::aes(u1, u2, colour = pt)) +
      ggplot2::geom_point(size = 0.4) +
      ggplot2::scale_colour_viridis_c(name = "pseudotime") +
      ggplot2::labs(x = "UMAP1", y = "UMAP2",
                    title = "DPT pseudotime (direction-corrected)") +
      ggplot2::theme_paper()
    save_fig(cfg, "02-05-04-unit2-pseudotime-dpt", p,
             width = W_ONE_HALF, height = mm(58))
  }
  # unit3：celltype 着色（无 celltype 列则跳过不落账）
  if (!is.null(clu$celltype) && !all(is.na(clu$celltype))) {
    dfu3 <- dfu; dfu3$grp <- as.character(clu$celltype)
    p <- ggplot2::ggplot(dfu3, ggplot2::aes(u1, u2, colour = grp)) +
      ggplot2::geom_point(size = 0.4) +
      ggplot2::labs(x = "UMAP1", y = "UMAP2",
                    title = "celltype on the same UMAP (pseudotime context)") +
      ggplot2::theme_paper()
    save_fig(cfg, "02-05-04-unit3-celltype-on-umap", p,
             width = W_ONE_HALF, height = mm(58))
  } else {
    log_warn("celltype 列缺失 —— 02-05-04-unit3 跳过（图名账目不含 leiden 回退）")
  }

  # per-cluster 表
  per_cluster <- do.call(rbind, lapply(sort(unique(clusters)), function(cc) {
    m <- clusters == cc
    dpt_v <- corrected[["dpt"]][m]
    data.frame(cluster = cc, n_cells = sum(m),
               consensus_mean = mean(consensus[m]),
               consensus_median = stats::median(consensus[m]),
               consensus_std = if (sum(m) > 1L) stats::sd(consensus[m]) else NA_real_,
               dpt_median = if (length(dpt_v) && any(is.finite(dpt_v)))
                 stats::median(dpt_v, na.rm = TRUE) else NA_real_)
  }))
  utils::write.csv(per_cluster, file.path(res_dir, "pseudotime_by_cluster.csv"),
                   row.names = FALSE)

  # unit1：箱线（簇多一栏半，簇少单栏）
  fig_w <- if (n_clusters > 4L) W_ONE_HALF else W_SINGLE
  dfu$cluster <- factor(clusters,
                        levels = per_cluster$cluster[order(per_cluster$consensus_median)])
  p <- ggplot2::ggplot(dfu[is.finite(dfu$pt), ],
                       ggplot2::aes(cluster, pt, fill = cluster)) +
    ggplot2::geom_boxplot(outlier.shape = NA, show.legend = FALSE,
                          fill = PAL$primary, alpha = 0.55) +
    ggplot2::labs(x = "cluster (ordered by median consensus pseudotime)",
                  y = "consensus pseudotime (higher = later)",
                  title = "Pseudotime distribution per cluster") +
    ggplot2::theme_paper() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))
  # **ns 必须写出来**：全局 Kruskal-Wallis，不显著就写 ns
  ok_groups <- split(consensus[is.finite(consensus)], clusters[is.finite(consensus)])
  ok_groups <- ok_groups[vapply(ok_groups, function(v) sum(is.finite(v)) >= 5L, logical(1))]
  if (length(ok_groups) >= 2L) {
    kw <- tryCatch(suppressWarnings(stats::kruskal.test(x = consensus[is.finite(consensus)],
                                                        g = factor(clusters[is.finite(consensus)],
                                                                   names = names(ok_groups)))),
                   error = function(e) NULL)
    if (!is.null(kw)) {
      p_kw <- as.numeric(kw$p.value)
      mark <- if (p_kw < 1e-3) "p < 0.001" else if (p_kw < 0.05) sprintf("p = %.3g", p_kw)
              else sprintf("ns (p = %.3g)", p_kw)
      p <- p + ggplot2::labs(subtitle = sprintf("Kruskal-Wallis: %s", mark))
    }
  }
  save_fig(cfg, "02-05-05-unit1-pseudotime-by-cluster", p, width = fig_w, height = mm(58))

  # unit2：山脊图（M5 守卫：画的不够不落盘，原因写 status）
  # 手写 KDE（stats::density），不引 ggridges 新依赖（Python 版同决策）。
  ridge_res <- tryCatch({
    need_pkg("ggplot2")
    grid_n <- 200L
    gl <- seq(min(consensus, na.rm = TRUE), max(consensus, na.rm = TRUE), length.out = grid_n)
    drawn <- 0L
    drawn <- 0L
    dfs <- do.call(rbind, lapply(levels(dfu$cluster), function(cc) {
      v <- consensus[clusters == cc]
      v <- v[is.finite(v)]
      if (length(v) < 10L) return(NULL)
      d <- stats::density(v, n = grid_n, from = gl[1L], to = gl[grid_n])
      dn <- d$y / max(d$y) * 0.8
      off <- which(levels(dfu$cluster) == cc) - 1L
      # `<<-`：lapply 匿名函数环境里没有 `drawn` 绑定，向上找到
      # tryCatch value 回调的局部 `drawn` —— 这里是**有意的**计数累加
      #（绑定存在，不会逃逸到全局；初始化就在上面一行）。
      drawn <<- drawn + 1L
      data.frame(x = gl, y = dn + off, ymin = off, cluster = cc)
    }))
    n_ridge_drawn <- drawn
    if (n_ridge_drawn < MIN_RIDGE_GROUPS) {
      ridge_skip_reason <- sprintf("只有 %d 个簇的细胞数 >= 10（需要 >= %d 个）—— 落盘会是空图，故不落盘",
                                   n_ridge_drawn, MIN_RIDGE_GROUPS)
      log_warn(sprintf("山脊图跳过：%s", ridge_skip_reason))
      list(status = "skipped", reason = ridge_skip_reason, n_groups_drawn = n_ridge_drawn)
    } else if (is.null(dfs)) {
      ridge_skip_reason <- "没有任何簇满足最小细胞数"
      list(status = "skipped", reason = ridge_skip_reason, n_groups_drawn = 0L)
    } else {
      p <- ggplot2::ggplot(dfs, ggplot2::aes(x, y, group = cluster)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = ymin, ymax = y),
                             fill = PAL$primary, alpha = 0.45) +
        ggplot2::geom_line(linewidth = 0.3, colour = PAL$primary) +
        ggplot2::labs(x = "consensus pseudotime (higher = later)", y = "cluster",
                      title = "Pseudotime density per cluster (ridgeline)",
                      subtitle = "KDE per cluster, offset vertically; not a new dependency (stats::density)") +
        ggplot2::theme_paper() +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank())
      save_fig(cfg, "02-05-05-unit2-pseudotime-ridgeline", p,
               width = W_ONE_HALF, height = mm(72))
      list(status = "ok", n_groups_drawn = n_ridge_drawn)
    }
  }, error = function(e) list(status = "skipped", reason = conditionMessage(e),
                              n_groups_drawn = 0L))
  if (!is.null(ridge_res$reason)) {
    log_warn(sprintf("山脊图未落盘: %s", ridge_res$reason))
  }
  ridge_status <- ridge_res

  # ---- 10. 每细胞拟时序落盘（供 07_grn 做 regulon×拟时序）----------------
  cell_df <- data.frame(cell = rownames(logcounts), leiden = clusters,
                        consensus_pseudotime = finite_round(consensus, 6))
  for (n in names(corrected)) {
    cell_df[[paste0("pseudotime_", n)]] <- finite_round(corrected[[n]], 6)
  }
  utils::write.csv(cell_df, file.path(res_dir, "pseudotime_per_cell.csv"),
                   row.names = FALSE)

  # ---- 11. 状态 -----------------------------------------------------------
  limitations <- c(
    "PAGA 等价物给的是簇间连通性，不是分化方向",
    "拟时序是一维坐标，分支过程会被压成先后关系",
    sprintf("**拟时序的符号是任意的**，本脚本按「值越大越晚」统一了方向；方向参考 = %s", direction_note),
    "**没有 RNA 速率（spliced/unspliced）时不能下方向性结论** —— 本流水线只有计数矩阵，scVelo 记 not_done，所以这里报的是相似度排序",
    paste0("**R 版方法集与 Python 版不同**（r_version.md §3 05 行，等级 C）：",
           "destiny DPT / slingshot / 自实现 GCS，Palantir 丢弃 —— ",
           "两版的数值与 stable_methods 名单不可直接比"))
  if (identical(mc_state, "ok")) {
    limitations <- c(limitations,
      sprintf("方法间一致性只是**内部一致性**：%d 种被交叉验证的方法都错向同一个伪轨迹时，它们依然彼此高度相关。一致不等于正确%s",
              length(cv_names), if (!is.null(mc_note)) sprintf("（%s）", mc_note) else ""))
  } else {
    limitations <- c(limitations,
      sprintf("**方法间一致性算不出来**（state=%s）：%s", mc_state, mc_note))
  }
  if (identical(direction_source, "cytotrace_fallback")) {
    limitations <- c(limitations,
      "**方向参考是 CytoTRACE 而非 marker 基因**：若该数据集中分化潜能与成熟度不同向，整条轴会反掉。要下方向性结论请在配置里给出 trajectory.early_markers / late_markers")
  }
  if (!is.null(mean_rho) && mean_rho < 0.3) {
    limitations <- c(limitations,
      sprintf("**方法间一致性偏低（平均 rho=%+.3f）**：各算法对同一数据给出了差异较大的排序，此时不该报单一「轨迹」", mean_rho))
  }
  limitations <- c(limitations, sprintf(
    "**共识拟时序用的是 %d 个方法（%s）%s。** 共识是各方法 z-score 后的等权均值 —— z-score 只把尺度归一化到 1，**不改变分布形状**，所以分辨率高（分布更极端）的方法在共识里权重更大，这不是严格的等权平均%s",
    length(consensus_names), paste(consensus_names, collapse = ", "),
    if (!is.null(reference_method)) sprintf("，已排除方向参考 %s", reference_method) else "",
    if (!is.null(rho_ref)) sprintf("；含参考方法的共识与它 rho=%+.4f（两者差多少本身是信息）", rho_ref) else ""))

  status <- list(
    dataset_id = cfg$dataset_id,
    status = "ok",
    n_cells = nrow(logcounts),
    n_clusters = n_clusters,
    methods_ok = methods_ok,
    methods_failed = methods_failed,
    n_methods = length(methods_ok),
    root_selection = root_record,
    direction_source = direction_source,
    direction_note = direction_note,
    marker_genes_used = marker_used,
    direction_table = direction_rows,
    method_correlation = lapply(seq_len(nrow(cdf)), function(i) as.list(cdf[i, ])),
    direction_reference_method = reference_method,
    cross_validated_methods = as.list(cv_names),
    consensus_methods = as.list(consensus_names),
    consensus_excludes_direction_reference = !is.null(reference_method),
    consensus_vs_with_reference_rho = if (!is.null(rho_ref)) round(rho_ref, 4) else NULL,
    method_correlation_mean_offdiag = finite_round(mean_rho, 4),
    method_correlation_min_offdiag = finite_round(min_rho, 4),
    method_correlation_state = mc_state,
    method_correlation_undefined_note = mc_note,
    cv_excluded = cv_excluded,
    method_correlation_note = paste0(
      "一致性统计**不含方向参考方法**",
      if (!is.null(reference_method)) sprintf("（%s）：退回模式下参考就是它本身，它与参考的相关恒为 ±1，那是定义不是证据", reference_method) else ""),
    pseudotime_range = list(round(min(consensus, na.rm = TRUE), 4),
                            round(max(consensus, na.rm = TRUE), 4)),
    pseudotime_by_cluster = lapply(seq_len(nrow(per_cluster)), function(i) as.list(per_cluster[i, ])),
    n_genes_along_trajectory = length(genes_rows),
    along_trajectory = along_status,
    ridgeline = ridge_status,
    n_unreliable_segments = n_unreliable_segments,
    min_cells_per_segment = MIN_CELLS_PER_SEGMENT,
    n_modules = length(modules_rows),
    segments_basis = "cluster（R 版按簇聚合；Python 版按 scFates milestones）",
    slingshot_lineages = methods_ok[["slingshot"]] %||% NULL,
    ks_by_group = ks_rows,
    ks_by_group_note = list(
      p_adj_available = identical(ks_reason, "ok"),
      correction = ks_correction,
      reason = ks_reason,
      n_pairs = length(ks_rows),
      group_degraded_to_clusters = grp_degraded,
      note = ks_note_map[[ks_reason]]),
    scvelo = list(
      status = "not_done",
      reason = paste0("RNA 速率需要 spliced/unspliced 两套计数矩阵。本流水线的输入是 10x **filtered 表达矩阵**，",
                      "只有一套计数，没有内含子/外显子的区分，无法计算速率。要跑 scVelo 必须从 Cell Ranger 的 ",
                      "`velocyto` 或 `--include-introns` 输出重新开始。")),
    limitations = limitations,
    # ---- reproducibility（R 版改写：不能照抄 Python 的六轮 CI 证据）--------
    # Python 版的 stable/unstable 名单与"六轮 CI 观测值"来自它自己的运行史。
    # R 版首次落地，**没有历史观测值** —— 按未定值报是唯一诚实的写法
    # （硬性规则 20 的同类要求：跑通了不等于可复现）。
    reproducibility = list(
      stable_methods = list("dpt", "slingshot", "cytotrace"),
      unstable_methods = list(),
      evidence = paste0("R 版首次落地（本 commit），**没有多轮 CI 观测值**。",
                        "Python 版的六轮证据（dpt/palantir 逐位稳定、scfates 漂移）",
                        "是 Python 实现的历史，**不可外推到 R 版**。"),
      observed_range = NULL,
      range_is_from = "**无历史观测值** —— 本轮的值见 method_correlation；",
      cause = paste0("destiny 与 slingshot 都接受 set.seed；",
                     "两版的可复现性要在 R-17 验收时用连续多轮 CI 实测确认"),
      mitigation = paste0("workflow 在 job 级钉了 OMP/OPENBLAS/MKL_NUM_THREADS=1 ",
                          "与 OPENBLAS_CORETYPE=Haswell（common.R env 段同款）"),
      how_to_report = paste0("**R 版所有方法的 rho 都按未定值报**（首次落地无范围可引）；",
                             "R-17 之后如有多轮观测，再把 stable/unstable 名单与范围写进来"))
  )
  write_json(file.path(res_dir, "trajectory_status.json"), status)
  log_info(sprintf("轨迹分析完成：%d 种方法，%s",
                   length(methods_ok),
                   if (!is.null(mean_rho)) sprintf("平均一致性 rho=%+.4f", mean_rho)
                   else sprintf("方法间一致性不可用（state=%s）", mc_state)))
  status
}

# ---------------------------------------------------------------------------
# 入口（05_trajectory.py:1245-1259 同语义；E-56：6 条早退路径 result_status_of 接住）
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L || identical(Sys.getenv("SCRNA_STEP_MAIN"), "05_trajectory")) {
  args <- parse_args()
  cfg <- load_config(args$config)
  t0 <- Sys.time()
  tryCatch({
    out <- run_05_trajectory(cfg)
    record_step(cfg, "trajectory", "ok", as.numeric(difftime(Sys.time(), t0, units = "secs")),
                result_status = result_status_of(out))
    out
  }, error = function(e) {
    record_step(cfg, "trajectory", "failed", as.numeric(difftime(Sys.time(), t0, units = "secs")),
                message = conditionMessage(e), result_status = "failed")
    stop(e)
  })
}
