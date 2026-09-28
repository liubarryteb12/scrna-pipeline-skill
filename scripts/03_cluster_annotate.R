#!/usr/bin/env Rscript
# =============================================================================
# 03_cluster_annotate.R — 邻居图、UMAP、Leiden 聚类、marker 基因、细胞类型打分
# （R 版，与 03_cluster_annotate.py 同号同语义）
#
# **注释是打分提示，不是结论。** 见 assets/celltype_markers.yml 开头的说明。
# 每个簇的 assignment 都带 `score_margin`（第一名与第二名的差）与
# `margin_state`（这个差算不算得出来）——
# margin 小的 assignment 不该被当成结论，而这一点只有把 margin 写出来
# 才看得出来。只报一个类型名等于把不确定性藏起来。
#
# `margin_state` 有四个取值，**它们的处理建议完全不同**（E-66 第二版同族）：
#   - ok               margin >= 0.05，assignment 可当结论
#   - low_margin       有第二名但差 < 0.05，两条路真的分不开
#   - single_celltype  **没有第二名可比较**（该簇只落进一个候选类型）
#   - margin_undefined 有第二名但分数含 nan/inf，差算不出来
# 后两者 `assignment_confident` 是 `null`（**不是 `false`**）：它们不是
# "不确定"，是"这个指标在这里不适用"—— 要补签名基因，不是去比对两条路。
#
# **两条独立的注释路径**（文档 §2.4）：
#   1. marker 签名打分（score_genes scanpy 语义复刻 + 簇均值）→ celltype_annotation.csv
#   2. CellTypist 预训练模型 → **R 侧无对应包，如实记 not_run**（r_version.md §3）
#
# 方法学等级（references/r_version.md §3）：聚类 A/B（bluster 图聚类等价
# leiden igraph flavor）、marker A（wilcoxon 同法）、打分 C（score_genes
# scanpy 语义复刻：AddModuleScore 在全基因集上 bins 必塌，run21/23R/25 实锤，
# 改 rank 整除分桶逐条复刻 scanpy/_score_genes.py —— 见 score_genes_scanpy_style）。
#
# **签名基因必须在全基因集里找**（03_cluster_annotate.py:96-103 的教训）：
# 实测用 HVG 找会把 CD3D 判缺失、两簇分数逐位相同 margin 恰为 0 ——
# 签名被削到无法区分。R 版从 raw.rds（全基因集 counts）取签名基因，
# 在 log 化矩阵上打分。
# =============================================================================

.load_common_for_local <- function() {
  argv0 <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(argv0)) {
    self <- normalizePath(sub("^--file=", "", argv0[1]), mustWork = FALSE)
    lib <- file.path(dirname(self), "lib", "common.R")
    if (file.exists(lib) && !exists("record_step", envir = globalenv())) {
      sys.source(lib, envir = globalenv())
    }
  }
  invisible(NULL)
}
.load_common_for_local()

# ---------------------------------------------------------------------------
# 簇标签排序键（03_cluster_annotate.py:50-69 同语义；M18）
# 标签是纯数字时按数值排（数字在前），混合标签按字符串排 —— 对任何输入
# 都返回可比较的键，不抛异常。
# ---------------------------------------------------------------------------
cluster_sort_key <- function(x) {
  s <- as.character(x)
  ifelse(grepl("^[0-9]+$", s),
         list(as.integer(s)), list(s))   # 逐元素比较时 R 的 order 数字在前
}
# order() 用的键：数字簇 -> (0, 数值)，其它 -> (1, 字符串)
cluster_order <- function(labels) {
  s <- as.character(labels)
  is_num <- grepl("^[0-9]+$", s)
  num_val <- ifelse(is_num, as.integer(s), NA_integer_)
  str_val <- ifelse(is_num, "", s)
  order(ifelse(is_num, 0L, 1L), num_val, str_val)
}

# 灰底 marker UMAP：DYNAMIC_FIG_BASES 声明与 Python 版同款（8 张运行时命名）
DYNAMIC_FIG_BASES <- list("03" = 8L)

# 脚本定位：独立 Rscript 用 --file=；被 main_analysis.R source 时由入口
# 先设置 options(scrna.script_path=...)（与 01_qc.R:43-47 同款；01/06/07/08
# 都有本定义，03 之前漏了 —— load_signatures 引用未定义名字，单独跑
# Rscript 03 时会在 L70 NameError）。
get_script_path <- function() {
  argv0 <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(argv0)) return(normalizePath(sub("^--file=", "", argv0[1]), mustWork = FALSE))
  getOption("scrna.script_path", "")
}

load_signatures <- function(cfg) {
  p <- file.path(dirname(dirname(get_script_path())), "assets", "celltype_markers.yml")
  need_pkg("yaml")
  doc <- yaml::read_yaml(p)
  name <- (cfg$analysis %||% list())$celltype_markers %||% "default"
  sigs <- doc$signatures %||% list()
  if (is.null(sigs[[name]])) {
    stop(sprintf("celltype_markers='%s' 不在 %s 里；可选: %s",
                 name, basename(p), paste(sort(names(sigs)), collapse = ", ")),
         call. = FALSE)
  }
  sigs[[name]]
}

# ---------------------------------------------------------------------------
# score_genes scanpy 语义复刻（根修 run21-25 实锤的 bins 塌缩问题）
#
# **为什么不再用 AddModuleScore**：Seurat 的 bin 是 cut(mean.expr, breaks=ctrl+1)
# 分位数切桶，全基因集上重复分位值多 → `Insufficient data values to produce
# 24 bins` 任何 ctrl 档位都塌（run21/23R/25 三轮实锤，降级链 25→1 无效）。
# scanpy score_genes 的分桶是 **rank(method="min") 整除 n_items**（scanpy 源码
# scanpy/tools/_score_genes.py:274-275）—— 重复表达值落同桶是正常路径、
# 空桶只是 WARN 跳过，**结构上不会塌**。R 侧逐条复刻该语义：
#   1. obs_avg = gene_pool 每基因平均表达（稀疏列均值）
#   2. n_items = round(len(obs_avg) / (n_bins-1))；obs_cut = rank(min) 整除 n_items
#   3. 对签名基因命中的每个桶：取同桶**非签名**基因，空桶 WARN 跳过，
#      多于 ctrl_size 个时随机抽 ctrl_size 个（seed 固定可复现）
#   4. score = mean(签名基因) - mean(对照基因)，按细胞
# 与 scanpy 数值不可逐位比（稀疏 nanmean 与 Matrix colMeans 的边界差异），
# 但语义同源：都是"签名基因平均表达 - 表达匹配对照平均表达"。
# ---------------------------------------------------------------------------
score_genes_scanpy_style <- function(log_counts_all, genes, n_bins = 25L,
                                     ctrl_size = 50L, seed = 42L) {
  stopifnot(length(genes) > 0L)
  gene_pool <- rownames(log_counts_all)
  genes <- intersect(genes, gene_pool)
  if (!length(genes)) {
    stop("score_genes: 签名基因一个都不在数据里", call. = FALSE)
  }
  obs_avg <- Matrix::rowMeans(log_counts_all)
  obs_avg <- obs_avg[is.finite(obs_avg)]
  gene_pool <- names(obs_avg)
  genes <- intersect(genes, gene_pool)
  if (!length(genes)) {
    stop("score_genes: 签名基因平均表达全部非有限", call. = FALSE)
  }
  n_items <- max(1L, round(length(obs_avg) / (n_bins - 1L)))
  obs_cut <- floor(rank(obs_avg, ties.method = "min") / n_items)
  names(obs_cut) <- gene_pool
  is_sig <- gene_pool %in% genes
  set.seed(seed)
  ctrl_sel <- character(0)
  for (cut_v in unique(obs_cut[genes])) {
    r_genes <- gene_pool[obs_cut == cut_v & !is_sig]
    if (!length(r_genes)) {
      log_warn(sprintf("score_genes: 桶 %d 无对照基因（跳过，scanpy 同为 WARN）", cut_v))
      next
    }
    if (ctrl_size < length(r_genes)) {
      r_genes <- sample(r_genes, ctrl_size)
    }
    ctrl_sel <- c(ctrl_sel, r_genes)
  }
  if (!length(ctrl_sel)) {
    stop("score_genes: 所有桶都没有对照基因（gene_pool 过小）", call. = FALSE)
  }
  sig_idx <- match(genes, gene_pool)
  ctrl_idx <- match(ctrl_sel, gene_pool)
  as.numeric(Matrix::rowMeans(log_counts_all[sig_idx, , drop = FALSE])) -
    as.numeric(Matrix::rowMeans(log_counts_all[ctrl_idx, , drop = FALSE]))
}

# ---------------------------------------------------------------------------
# 细胞类型打分（03_cluster_annotate.py:84-189 同语义；margin 三态逐条复刻）
# 返回 list(per_cell=矩阵, assign=data.frame, diag=list)；打分不可行时
# assign=NULL 且 diag 带 reason。
# ---------------------------------------------------------------------------
score_celltypes <- function(log_counts_all, clusters, sig, seed = 42L) {
  celltypes <- sig$celltypes %||% list()
  all_genes <- rownames(log_counts_all)
  present <- list(); missing_in_data <- list(); not_in_hvg <- list()
  for (ct in names(celltypes)) {
    marks <- as.character(celltypes[[ct]]$markers %||% character(0))
    genes <- intersect(marks, all_genes)
    if (length(genes)) present[[ct]] <- genes
    miss <- setdiff(marks, all_genes)
    if (length(miss)) missing_in_data[[ct]] <- miss
    # R 版直接在全基因集上打分，not_in_hvg 不适用 —— 保留键、记 empty，
    # schema 与 Python 版对齐（审计 S2：两个字段刻意分开）。
  }
  if (!length(present)) {
    return(list(per_cell = NULL, assign = NULL,
                diag = list(reason = "签名里的基因一个都不在数据里")))
  }
  # 诊断（run21 实锤沿用的输入自检：签名基因必须真的在矩阵行名里）
  log_info(sprintf("打分诊断: present %d 类型 / %d 基因; all_genes 头3=%s; CD3D 在 all_genes=%s",
                   length(present), length(unlist(present)),
                   paste(utils::head(all_genes, 3L), collapse = ","),
                   "CD3D" %in% all_genes))

  # **打分实现 = score_genes_scanpy_style（本文件上方，scanpy 语义复刻）**
  # run21/23R/25 三轮实锤 AddModuleScore 在全基因集上 bins 必塌
  # （cut() 分位数切桶对重复值敏感），rank 整除分桶（scanpy 原语义）结构上
  # 不会塌。放弃 Seurat 路线 = 03 的方法学等级从 A/B 降 C（r_version.md §3
  # 03 行待同步）；每细胞打分逐签名独立算（与 sc.tl.score_genes 循环同构）。
  # vapply 按列堆叠 → sc_mat = 签名x细胞（lg_all 列序 == clusters 序，调用点 539
  # 已保证：lg_all colnames = rownames(raw_sub) = common_cells = names(clusters)）。
  # 显式命名细胞后按 names(clusters) 重排，防静默错位（与 06 colLabels 同纪律）。
  sc_mat <- vapply(seq_along(present), function(j) {
    score_genes_scanpy_style(log_counts_all, present[[j]], n_bins = 25L,
                             ctrl_size = 50L, seed = seed)
  }, numeric(ncol(log_counts_all)))
  rownames(sc_mat) <- colnames(log_counts_all)
  colnames(sc_mat) <- names(present)
  if (!is.null(names(clusters))) {
    miss_cells <- sum(!names(clusters) %in% rownames(sc_mat))
    if (miss_cells > 0L) {
      stop(sprintf("打分矩阵缺 %d 个细胞的行（细胞名对齐失败）", miss_cells),
           call. = FALSE)
    }
    sc_mat <- sc_mat[names(clusters), , drop = FALSE]
  }

  # 按簇均值
  ucl <- unique(clusters)
  ucl <- ucl[cluster_order(ucl)]
  per_cell <- t(vapply(ucl, function(cl) {
    idx <- which(clusters == cl)
    colMeans(sc_mat[idx, , drop = FALSE])
  }, numeric(length(present))))
  rownames(per_cell) <- ucl
  colnames(per_cell) <- names(present)

  rows <- vector("list", length(ucl))
  for (i in seq_along(ucl)) {
    cl <- ucl[i]
    s <- sort(per_cell[i, ], decreasing = TRUE)
    top <- names(s)[1]
    if (length(s) <= 1L) {
      margin_state <- "single_celltype"; margin_out <- NULL; confident <- NULL
    } else {
      margin <- as.numeric(s[1] - s[2])
      if (!is.finite(margin)) {
        margin_state <- "margin_undefined"; margin_out <- NULL; confident <- NULL
      } else if (margin < 0.05) {
        margin_state <- "low_margin"; margin_out <- round(margin, 4); confident <- FALSE
      } else {
        margin_state <- "ok"; margin_out <- round(margin, 4); confident <- TRUE
      }
    }
    rows[[i]] <- list(
      cluster = as.character(cl),
      assigned = top,
      top_score = round(as.numeric(s[1]), 4),
      runner_up = if (length(s) > 1L) names(s)[2] else NULL,
      runner_up_score = if (length(s) > 1L) round(as.numeric(s[2]), 4) else NULL,
      score_margin = margin_out,
      margin_state = margin_state,
      assignment_confident = confident)
  }
  assign <- do.call(rbind, lapply(rows, function(r) {
    data.frame(cluster = r$cluster, assigned = r$assigned,
               top_score = r$top_score,
               runner_up = r$runner_up %||% NA_character_,
               runner_up_score = r$runner_up_score %||% NA_real_,
               score_margin = r$score_margin,
               margin_state = r$margin_state,
               assignment_confident = r$assignment_confident,
               stringsAsFactors = FALSE)
  }))
  diag <- list(
    n_celltypes_scored = length(present),
    celltypes = sort(names(present)),
    missing_markers = missing_in_data,
    not_in_hvg = list(),
    # **同样刻意分开**（E-66 第二版）：只有 low_margin 才叫"不确定"，
    # single_celltype / margin_undefined 单列 —— 它们不是"不确定"，
    # 是"这个指标算不出来"。
    n_clusters_low_margin = sum(assign$margin_state == "low_margin"),
    n_clusters_margin_undefined = sum(assign$margin_state %in%
                                        c("single_celltype", "margin_undefined")),
    margin_state_counts = as.list(table(assign$margin_state)))
  list(per_cell = per_cell, assign = assign, diag = diag)
}

# ---------------------------------------------------------------------------
# CellTypist（R 侧无对应包 —— r_version.md §3 裁决：如实记 not_run）
# ---------------------------------------------------------------------------
try_celltypist_r <- function(cfg) {
  list(attempted = TRUE, status = "not_run",
       reason = paste0("CellTypist 是 Python 包，R 版无对应实现（r_version.md §3 ",
                       "裁决）；本路径 R 版记 not_run，不退回 marker 打分、不假装"),
       model = (cfg$analysis %||% list())$celltypist_model %||% "Immune_All_Low.pkl")
}

# ---------------------------------------------------------------------------
# marker 点图的 top3 挑选（03_cluster_annotate.py:454-482 同语义；M24）
# ---------------------------------------------------------------------------
top_markers_per_cluster <- function(markers, top_n = 3L) {
  # markers: data.frame(cluster, gene, score)
  ordered <- character(0); kept <- list(); dropped <- list()
  cls <- unique(as.character(markers$cluster))
  cls <- cls[cluster_order(cls)]
  for (c in cls) {
    sub <- markers[as.character(markers$cluster) == c, , drop = FALSE]
    sub <- sub[order(-sub$score), , drop = FALSE]
    sub <- utils::head(sub, top_n)
    kept[[c]] <- character(0)
    for (g in sub$gene) {
      if (g %in% ordered) {
        owner <- names(kept)[vapply(kept, function(v) g %in% v, logical(1))][1]
        dropped[[length(dropped) + 1L]] <-
          list(cluster = c, gene = g, kept_by_cluster = owner)
        next
      }
      ordered <- c(ordered, g)
      kept[[c]] <- c(kept[[c]], g)
    }
  }
  list(ordered = ordered,
       diag = list(
         per_cluster_kept = as.list(vapply(kept, length, integer(1))),
         n_columns = length(ordered),
         n_clusters = length(unique(as.character(markers$cluster))),
         top_n_per_cluster_requested = top_n,
         shared_dropped = dropped,
         why = paste0("同一基因可能是多个簇的前 3，而点图的列不能重复 —— ",
                      "所以列数少于 `top_n × n_clusters`。被去掉的那些记在 ",
                      "`shared_dropped` 里，不是静默丢弃")))
}

# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
run_03_cluster_annotate <- function(cfg) {
  ensure_dirs(cfg)
  set_seed(cfg)
  data_dir <- cfg$output$data_dir
  res_dir <- cfg$output$results_dir

  need_pkg("Seurat")
  ig_path <- file.path(data_dir, "integrated.rds")
  if (!file.exists(ig_path)) {
    stop(sprintf("缺输入 %s —— 先跑 02_integrate.R", ig_path), call. = FALSE)
  }
  ig <- readRDS(ig_path)
  counts <- ig$counts          # HVG 子集 counts（细胞 x 基因）
  logcounts <- ig$logcounts
  use_rep <- ig$use_rep %||% "pca"
  pca <- ig$pca
  if (identical(use_rep, "harmony")) {
    # harmony 的 embedding 存在 Seurat 对象里；02_integrate.R 只落了 pca ——
    # 这里用 pca 兜底并如实告警（下游一致性影响由验收层 honesty 记）。
    log_warn("use_rep=harmony 但 R 版中间对象只落了 pca embedding —— 用 pca 继续")
    use_rep <- "pca"
  }
  log_info(sprintf("读入 %d 细胞 x %d 基因（use_rep=%s）",
                   nrow(counts), ncol(counts), use_rep))
  # 02 落盘的中间对象是 细胞x基因（Python adata 语义）；Seurat API 要 **基因x细胞**。
  # CI run11（36340148115）实锤：直接 CreateSeuratObject(counts=细胞x基因) 会把基因
  # 当细胞建对象，随后 obj[["pca"]] <- pca_dr 因细胞名不匹配报
  # `Cannot add new cells with [[<-`。统一在这里转一次，下游 Seurat 调用全用这个。
  counts_gxc <- Matrix::t(counts)   # 基因x细胞（Seurat API 方向）

  rd <- cfg$reduce
  n_pcs <- as.integer(rd$n_pcs)
  n_pcs_use <- rd$n_pcs_use %||% n_pcs

  # ---- 1. 邻居 + UMAP -----------------------------------------------------
  set.seed(as.integer(cfg$analysis$seed))
  k <- as.integer(rd$n_neighbors)
  # Seurat 链（可参考代码 99 验证过的路）：CreateSeuratObject 带上 PCA
  # embedding 作 DimReduc -> RunUMAP(reduction="pca", dims=1:n_pcs_use)。
  # **不是** RunUMAP(reduction.model=...) —— 那是"用已训好的 UMAP 模型投影
  # 新数据"的接口，给它传 PCA 矩阵语义就错了。
  obj <- Seurat::CreateSeuratObject(counts = counts_gxc)   # 基因x细胞
  emb <- as.matrix(pca[, seq_len(min(n_pcs_use, ncol(pca))), drop = FALSE])
  colnames(emb) <- sprintf("PC_%d", seq_len(ncol(emb)))
  rownames(emb) <- rownames(pca)
  pca_dr <- Seurat::CreateDimReducObject(
    embeddings = emb, key = "PC_", assay = "RNA",
    global = TRUE)
  obj[["pca"]] <- pca_dr
  obj <- Seurat::RunUMAP(obj, reduction = "pca",
                         dims = seq_len(ncol(emb)), verbose = FALSE)
  umap <- Seurat::Embeddings(obj, "umap")
  log_info(sprintf("UMAP 完成（n_neighbors=%d, n_pcs=%d）", k, n_pcs_use))

  # ---- 2. 聚类（bluster 图聚类等价 leiden；多分辨率扫描）-------------------
  # 让"选 1.0"变成有依据的决定，而不是默认值。
  # bluster 的 NNGraphParam 本身没有 resolution 概念（leiden 在图上直接
  # 优化模块度）；分辨率语义用**边的 jaccard 权重 + k** 近似 —— 扫描靠
  # `cluster.args$partition_args` 透传 leiden 的 resolution。
  need_pkg("bluster")
  .clus_at <- function(res) {
    # clusterRows 聚**行** = 观测：emb 是 细胞xPC，行已是细胞 —— 直接传，
    # **不能 t()**（t 后行=PC，会去聚 40 个主成分，返回长度 40 与 2574 细胞静默错配）。
    # jaccard 权重经 makeSNNGraph 的 **type="jaccard"**（run13 实锤：weights.type
    # 是虚构参数名 → unused argument）。leiden 的 objective_function/resolution
    # 经 cluster.args 透传 igraph::cluster_leiden（拼写已对 igraph 源码核实）。
    bluster::clusterRows(emb, bluster::NNGraphParam(
      k = k, type = "jaccard", cluster.fun = "leiden",
      cluster.args = list(objective_function = "modularity",
                          resolution = res)))
  }
  scan_res <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0)
  clus_saved <- vector("list", length(scan_res))
  names(clus_saved) <- as.character(scan_res)
  scan_rows <- lapply(scan_res, function(r) {
    clus <- clus_saved[[as.character(r)]] %||% {
      clus_saved[[as.character(r)]] <<- .clus_at(r)
      clus_saved[[as.character(r)]]
    }
    list(resolution = r, n_clusters = length(unique(clus)))
  })
  scan <- do.call(rbind, lapply(scan_rows, as.data.frame))
  utils::write.csv(scan, file.path(res_dir, "cluster_resolution_scan.csv"),
                   row.names = FALSE)
  log_info(paste("分辨率扫描:", paste(sprintf("%s->%d", scan$resolution,
                                               scan$n_clusters), collapse = ", ")))

  # 扫描图（02-03-01-unit1，W_SINGLE×58mm；axvline = 配置使用值）
  df <- data.frame(res = scan$resolution, n = scan$n_clusters)
  p <- ggplot2::ggplot(df, ggplot2::aes(res, n)) +
    ggplot2::geom_line(colour = PAL$primary) +
    ggplot2::geom_point(colour = PAL$primary, size = 1.5) +
    ggplot2::geom_vline(xintercept = as.numeric(rd$resolution),
                        colour = PAL$highlight, linetype = "dashed") +
    ggplot2::labs(x = "Leiden resolution", y = "number of clusters",
                  title = "Cluster count vs resolution") +
    theme_paper()
  save_fig(cfg, "02-03-01-unit1-cluster-resolution-scan", p,
           width = W_SINGLE, height = mm(58))

  # ---- 3. 用配置的分辨率定稿 ----------------------------------------------
  # 复用扫描时存下的聚类结果（clus_saved），**不再重算一遍** ——
  # 重算既慢，又有与扫描表不一致的风险（leiden 有随机性）。
  res_used <- as.numeric(rd$resolution)
  idx <- which.min(abs(scan_res - res_used))
  clusters <- as.character(clus_saved[[idx]])
  names(clusters) <- rownames(counts)   # 行名 = 细胞名（counts 是 细胞x基因）
  n_clusters <- length(unique(clusters))
  log_info(sprintf("最终聚类: %d 个簇（resolution=%.1f -> 扫描档 %.1f）",
                   n_clusters, res_used, scan_res[idx]))

  # UMAP 簇图（02-03-02-unit1，W_ONE_HALF×84mm；簇号白底标注 = 第二线索）
  ord <- cluster_order(unique(clusters))
  ucl <- unique(clusters)[ord]
  dfc <- data.frame(u1 = umap[, 1], u2 = umap[, 2], leiden = factor(clusters, levels = ucl))
  # 簇中心标注
  cent <- do.call(rbind, lapply(ucl, function(cl) {
    m <- clusters == cl
    data.frame(x = mean(umap[m, 1]), y = mean(umap[m, 2]), lab = cl)
  }))
  # 图例规则 31：框外右侧单列
  leg <- cowplot::get_legend(
    ggplot2::ggplot(dfc, ggplot2::aes(u1, u2, colour = leiden)) +
      ggplot2::geom_point(size = 0.4) + theme_paper())
  p <- ggplot2::ggplot(dfc, ggplot2::aes(u1, u2, colour = leiden)) +
    ggplot2::geom_point(size = 0.4, alpha = 0.75) +
    ggplot2::scale_colour_manual(values = rep_len(PAL_CYCLE, max(length(ucl), length(PAL_CYCLE)))) +
    ggplot2::geom_point(data = cent, ggplot2::aes(x, y), colour = "black",
                        alpha = 0, size = 0.001) +   # 保持图层占位
    ggplot2::geom_text(data = cent, ggplot2::aes(x, y, label = lab),
                       inherit.aes = FALSE, fontface = "bold", size = 2.4,
                       colour = "black") +
    ggplot2::labs(x = "UMAP1", y = "UMAP2",
                  title = sprintf("Leiden clusters (n=%d, resolution=%s)",
                                  n_clusters, rd$resolution)) +
    theme_paper() + ggplot2::theme(legend.position = "none")
  p <- patchwork::wrap_plots(p, leg, ncol = 2, widths = c(3, 0.5))
  save_fig(cfg, "02-03-02-unit1-umap-clusters", p,
           width = W_ONE_HALF, height = mm(84))

  # ---- 4. Marker 基因（Seurat FindAllMarkers，wilcoxon 同法）--------------
  obj <- Seurat::CreateSeuratObject(counts = counts_gxc)   # 基因x细胞（Seurat API 方向）
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj@meta.data$leiden <- factor(clusters, levels = ucl)
  Seurat::Idents(obj) <- "leiden"
  wm <- Seurat::FindAllMarkers(obj, only.pos = TRUE, test.use = "wilcox",
                               logfc.threshold = 0, min.pct = 0,
                               verbose = FALSE)
  markers <- data.frame(
    cluster = as.character(wm$cluster), gene = wm$gene, score = wm$avg_log2FC)
  # Python 版按 scores 列取 top；Seurat 无 scores 列，用 avg_log2FC 排序
  markers <- markers[order(match(markers$cluster, ucl), -markers$score), , drop = FALSE]
  utils::write.csv(markers, file.path(res_dir, "markers_all.csv"), row.names = FALSE)
  top_n <- as.integer((cfg$analysis %||% list())$top_markers %||% 25L)
  log_info(sprintf("marker 表: %d 行（每簇前 %d）", nrow(markers), top_n))

  # ---- 5. marker 点图 + 灰底 marker UMAP ----------------------------------
  top3res <- top_markers_per_cluster(markers, 3L)
  top3 <- top3res$ordered
  top3_diag <- top3res$diag
  all_genes <- colnames(logcounts)   # logcounts 是 细胞x基因，基因名在**列**
  top3 <- head(intersect(top3, all_genes), 24L)
  top3_diag$n_columns_after_cap <- length(top3)
  top3_diag$width_cap_note <- paste0(
    "列数上限 24 是按双栏 183mm / 每基因约 7mm 算的；",
    "被上限截掉的不影响每簇至少一个代表基因")
  log_info(sprintf("marker 点图: %d 列（%d 个簇，每簇请求 3 个）；跨簇重复去掉 %d 个",
                   top3_diag$n_columns, top3_diag$n_clusters,
                   length(top3_diag$shared_dropped)))

  if (length(top3)) {
    # **真 z-score**（按基因跨簇标准化），不是 min-max 0-1 —— 详见
    # common.R build_marker_dotplot_figure docstring ① 的标度矛盾。
    groups <- clusters
    frac <- matrix(0, nrow = length(ucl), ncol = length(top3),
                   dimnames = list(ucl, top3))
    mean_expr <- frac
    for (ri in seq_along(ucl)) {
      idx_g <- which(groups == ucl[ri])
      blk <- as.matrix(logcounts[idx_g, top3, drop = FALSE])
      frac[ri, ] <- colMeans(blk > 0)
      mean_expr[ri, ] <- colMeans(blk)
    }
    sd_v <- apply(mean_expr, 2, stats::sd)
    sd_v[sd_v == 0] <- 1
    zmat <- sweep(sweep(mean_expr, 2, colMeans(mean_expr), "-"), 2, sd_v, "/")

    fig <- build_marker_dotplot_figure(
      frac_df = as.data.frame(frac), z_df = as.data.frame(zmat),
      group_label = "Leiden cluster", title = "Top markers per cluster",
      subtitle = paste0("rows = Leiden clusters (identities: celltype_annotation.csv)\n",
                        "dot size = fraction of cells expressing the gene"))
    save_fig(cfg, "02-03-03-unit1-markers-dotplot", fig,
             width = W_DOUBLE, height = mm(126))

    # **灰底 marker UMAP 网格**（差距清单 #19，文献范式）：
    # dotplot 给"哪个簇高表达"，灰底 UMAP 给"高表达在哪块区域"。
    # 单图原则：每个基因一张（02-03-03 unit2..unit9）。
    for (gi in seq_along(top3)) {
      if (gi > 8L) break
      g <- top3[gi]
      expr <- as.numeric(logcounts[, g])
      dfm <- data.frame(u1 = umap[, 1], u2 = umap[, 2],
                        expr = ifelse(expr > 0, expr, NA_real_))
      pm <- ggplot2::ggplot(dfm, ggplot2::aes(u1, u2)) +
        ggplot2::geom_point(size = 0.4, alpha = 0.25, colour = "grey80") +
        ggplot2::geom_point(data = dfm[!is.na(dfm$expr), ],
                            ggplot2::aes(u1, u2, colour = expr), size = 0.5) +
        ggplot2::scale_colour_viridis_c(name = "expression") +
        ggplot2::labs(title = sprintf("%s on UMAP (grey = not detected)", g),
                      x = "UMAP1", y = "UMAP2") +
        theme_paper()
      save_fig(cfg, sprintf("02-03-03-unit%d-marker-%s", gi + 1L, tolower(g)),
               pm, width = W_SINGLE, height = mm(58))
    }
  }

  # ---- 6. 细胞类型打分 ----------------------------------------------------
  sig <- load_signatures(cfg)
  # **签名基因在全基因集（raw.rds）里找**，不是 HVG —— 见文件头教训。
  raw_path <- file.path(data_dir, "raw.rds")
  raw_counts <- if (file.exists(raw_path)) readRDS(raw_path) else counts
  # 打分在 log 化全基因集矩阵上做（签名基因交集 + 补全其它基因做对照 bin）
  # raw_counts（raw.rds）是 **细胞x基因**（00_fetch 落盘）→ 细胞名在 **rownames**；
  # logcounts 是 细胞x基因 → 细胞名在 rownames。run20 实锤：写成 colnames(raw_counts)
  # 得到基因名 ∩ 细胞名 = 空 → 签名基因"一个都不在数据里" → 打分 not_possible。
  common_cells <- intersect(rownames(raw_counts), rownames(logcounts))
  raw_sub <- raw_counts[common_cells, , drop = FALSE]
  # log 化（与 02_integrate 的 target_sum 一致）。产出保持 **基因x细胞**：
  # score_celltypes 的 CreateSeuratObject 要 Seurat API 方向，且
  # all_genes <- rownames(log_counts_all) 取基因名（run20 实锤：外面多套一层 t()
  # 让 lg_all 变 细胞x基因 → rownames 全是细胞名 → 同样全 miss）。
  target_sum <- (cfg$norm %||% list())$target_sum %||% 1e4
  lg_all <- log1p(Matrix::t(raw_sub) /
                    pmax(Matrix::colSums(raw_sub), 1) * target_sum)
  st <- score_celltypes(lg_all, clusters[common_cells], sig, seed = 42L)
  per_cell <- st$per_cell; assign <- st$assign; diag <- st$diag
  annot_record <- list(
    signature = (cfg$analysis %||% list())$celltype_markers %||% "default",
    signature_description = sig$description,
    signature_note = sig$note)
  if (is.null(assign)) {
    annot_record <- c(annot_record, list(status = "not_possible"), diag)
    log_warn(sprintf("细胞类型打分不可行: %s", diag$reason %||% ""))
  } else {
    utils::write.csv(data.frame(cluster = rownames(per_cell), per_cell,
                                check.names = FALSE),
                     file.path(res_dir, "celltype_scores.csv"), row.names = FALSE)
    utils::write.csv(assign, file.path(res_dir, "celltype_annotation.csv"),
                     row.names = FALSE)
    annot_record <- c(annot_record, diag, list(
      status = "ok",
      assignments = df_to_records(assign)))
    log_info(sprintf(paste0("细胞类型打分: %d 种类型；%d/%d 个簇 margin<0.05",
                            "（assignment 不确定）；margin 状态分布 %s"),
                     diag$n_celltypes_scored, diag$n_clusters_low_margin,
                     n_clusters,
                     paste(names(diag$margin_state_counts),
                           diag$margin_state_counts, collapse = ", ")))
    if (diag$n_clusters_low_margin > 0) {
      log_warn(sprintf(paste0("margin<0.05 的 %d 个簇其 assignment 不应被当成结论",
                              " —— 见 celltype_annotation.csv 的 ",
                              "margin_state='low_margin'"),
                       diag$n_clusters_low_margin))
    }
    if (diag$n_clusters_margin_undefined > 0) {
      log_warn(sprintf(paste0("另有 %d 个簇的 margin **算不出来**",
                              "（margin_state='single_celltype' / ",
                              "'margin_undefined'）—— 这不是「不确定」，",
                              "是「这个指标在这里不适用」：候选细胞类型只有一个，",
                              "没有第二名可比。处理方式与 low_margin 不同：",
                              "要补签名基因，不是去比对两条路"),
                       diag$n_clusters_margin_undefined))
    }

    # 打分热图（02-03-04-unit1）：grid_size 尺寸 + viridis + x 标签旋转
    gs <- grid_size(length(colnames(per_cell)), length(rownames(per_cell)))
    dfh <- data.frame(ct = factor(rep(colnames(per_cell), each = nrow(per_cell)),
                                  levels = colnames(per_cell)),
                      cl = factor(rep(rownames(per_cell), times = ncol(per_cell)),
                                  levels = rev(rownames(per_cell))),
                      score = as.numeric(per_cell))
    ph <- ggplot2::ggplot(dfh, ggplot2::aes(ct, cl, fill = score)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c(name = "score") +
      ggplot2::labs(x = "cell type signature", y = "cluster",
                    title = "Mean signature score per cluster") +
      theme_paper() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_fig(cfg, "02-03-04-unit1-celltype-scores-heatmap", ph,
             width = mm(gs$width), height = mm(gs$height))
  }

  # ---- 7. CellTypist（R 版记 not_run）-------------------------------------
  ct_info <- try_celltypist_r(cfg)
  ct_cmp <- list(compared = FALSE,
                 reason = "CellTypist 在 R 版 not_run，对比无对象")
  annot_record$celltypist <- ct_info
  annot_record$celltypist_vs_marker <- ct_cmp

  # ---- 8. 落盘（R 走 .rds）------------------------------------------------
  # **logcounts_all 透传**（02_integrate.R 存的全基因集 log 表达，
  # Python 版 adata.raw 的对应物）：06_communication 的配体-受体基因
  # 几乎不进 HVG，缺了它 06 会得到大量假阴性。03 不用它，但必须原样
  # 带到 clustered.rds，否则 06 只能拿 HVG 子集冒充全基因集。
  out <- file.path(data_dir, "clustered.rds")
  saveRDS(list(counts = counts, logcounts = logcounts,
               logcounts_all = ig$logcounts_all %||% NULL,
               clusters = clusters,
               umap = umap, hvg = colnames(counts),
               celltype = if (!is.null(assign))
                 setNames(assign$assigned[match(clusters, assign$cluster)], names(clusters))
               else NULL,
               pca = pca, use_rep = use_rep, n_comps = n_pcs), out)
  log_info(sprintf("已写出 %s（logcounts_all %s）", out,
                   if (!is.null(ig$logcounts_all))
                     sprintf("%d 基因（全基因集）", ncol(ig$logcounts_all))
                   else "缺失——06 通讯分析将退化为 HVG 子集（假阴性风险）"))

  # §0.2 跨部分交接的产出侧：这一份 rds 就是 Part 3 的参考（R 版边界）。
  # 契约三条（counts / celltype 列 / 基因集），语义与 Python 版相同；
  # R 版中间对象不走 .h5ad（r_version.md §5），Part 3 R 版消费本文件。
  has_celltype <- !is.null(assign)
  n_hvg <- ncol(counts)
  part3_ref <- list(
    path = out,
    contract = list(
      counts_layer = TRUE,   # counts 原始计数在 clustered.rds 的 counts 字段
      celltype_column = if (has_celltype) "celltype" else NULL),
    n_cells = nrow(counts),
    n_genes = n_hvg,
    counts_layer_source = "02_integrate 存的原始计数（HVG 子集）",
    how_part3_uses_it = paste0("Part 3 R 版的 05_deconvolution.R 读 ",
                               "clustered.rds 的 counts + celltype，按类型求均值",
                               "得到参考谱 S，再用 NNLS 解每个 spot 的组成"),
    limitations = c(
      if (n_hvg < 10000) paste0("**基因集是 HVG 子集** —— 本文件来自 ",
                                "integrated.rds（", n_hvg, " 个高变基因），不是全基因集。",
                                "Part 3 的参考谱因此只覆盖两边的共同基因，",
                                "共同基因太少时 Part 3 会直接报错") else NULL,
      if (!has_celltype) "**obs 里没有细胞类型列** —— Part 3 配 celltype_key 时会报错" else NULL))
  if (!has_celltype) {
    log_warn(sprintf(paste0("Part 3 参考契约不完整：细胞类型列=%s",
                            " —— 见 cluster_status.json 的 part3_reference"),
                     "null"))
  }

  status <- list(
    dataset_id = cfg$dataset_id,
    n_cells = nrow(counts),
    n_clusters = n_clusters,
    resolution = res_used,
    resolution_scan = df_to_records(scan),
    n_neighbors = k,
    n_pcs_used = as.integer(n_pcs_use),
    use_rep = use_rep,
    n_markers_rows = nrow(markers),
    # **M24**：点图的列数为什么少于 3 × 簇数。
    marker_dotplot_columns = top3_diag,
    cluster_sizes = as.list(table(clusters)[ucl]),
    annotation = annot_record,
    part3_reference = part3_ref,
    status = "ok")
  write_json(file.path(res_dir, "cluster_status.json"), status)
  status
}

# ---------------------------------------------------------------------------
# 入口（03_cluster_annotate.py:803-815 同语义；E-56）
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L || identical(Sys.getenv("SCRNA_STEP_MAIN"), "03_cluster_annotate")) {
  args <- parse_args()
  cfg <- load_config(args$config)
  t0 <- Sys.time()
  tryCatch({
    out <- run_03_cluster_annotate(cfg)
    record_step(cfg, "cluster_annotate", "ok",
                as.numeric(difftime(Sys.time(), t0, units = "secs")),
                result_status = result_status_of(out))
    out
  }, error = function(e) {
    record_step(cfg, "cluster_annotate", "failed",
                as.numeric(difftime(Sys.time(), t0, units = "secs")),
                message = conditionMessage(e), result_status = "failed")
    stop(e)
  })
}
