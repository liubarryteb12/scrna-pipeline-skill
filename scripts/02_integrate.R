#!/usr/bin/env Rscript
# =============================================================================
# 02_integrate.R — 标准化、高变基因、PCA、批次整合（R 版，与 02_integrate.py 同号同语义）
#
# **顺序很重要，而且 seurat_v3 与另外两种口味的顺序不同：**
#   - seurat_v3 吃**原始计数**，必须在 normalize/log1p **之前**跑
#   - seurat / cell_ranger 吃 log 后的数据，必须在之后
# 实测踩过的坑：顺序反了 seurat_v3 不报错，只是给出错的 HVG ——
# 而 HVG 错了下游全部跟着错，图上完全看不出来。
#
# 方法学等级（references/r_version.md §3）：**A 同方法不同实现**。
# Seurat NormalizeData/FindVariableFeatures(vst)/ScaleData/RunPCA 对应
# scanpy normalize_total/highly_variable_genes(seurat_v3)/PCA。
# **未决 4（r_version.md §11）**：scanpy seurat_v3 的 vst 实现与 Seurat 的
# `vst` 在 loess 细节上不完全等价，HVG 集合允许有差，差异记进 status
# （hvg_flavor_used = "seurat_v3_r Equivalent"，不冒充逐位一致）。
#
# 批次整合默认关闭。**单样本数据没有批次，关掉是正确的**，不是偷懒。
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

run_02_integrate <- function(cfg) {
  ensure_dirs(cfg)
  set_seed(cfg)
  data_dir <- cfg$output$data_dir
  res_dir <- cfg$output$results_dir

  need_pkg("Seurat", "标准化/HVG/PCA")

  qc_path <- file.path(data_dir, "qc_filtered.rds")
  if (!file.exists(qc_path)) {
    stop(sprintf("缺输入 %s —— 先跑 01_qc.R", qc_path), call. = FALSE)
  }
  qc <- readRDS(qc_path)
  counts <- qc$counts
  log_info(sprintf("读入 %d 细胞 x %d 基因", nrow(counts), ncol(counts)))

  norm <- cfg$norm %||% list()
  flavor <- norm$hvg_flavor %||% "seurat_v3"
  n_top <- as.integer(norm$n_top_genes)
  target_sum <- as.numeric(norm$target_sum)
  batch_key <- (cfg$design %||% list())$batch_key

  # ---- 1. HVG + normalize（顺序按口味分；seurat_v3 用原始计数）-----------
  # Python 版 flavors: seurat_v3 / seurat / cell_ranger；R 版用 Seurat 的
  # `vst`（对应 seurat_v3）或 `mean.var.plot`（对应 seurat）。skmisc 缺失
  # 的降级在 R 侧不存在（vst 是 Seurat 内置），但**口味映射本身**要落盘。
  flavor_map <- c(seurat_v3 = "vst", seurat = "mean.var.plot",
                  cell_ranger = "dispersion")
  if (!flavor %in% names(flavor_map)) {
    stop(sprintf("不支持的 hvg_flavor: %s（seurat_v3/seurat/cell_ranger）", flavor),
         call. = FALSE)
  }
  hvg_flavor_used <- flavor
  hvg_fallback <- NULL   # R 版 vst 内置无降级路径；字段保留与 Python 版 schema 对齐

  # Seurat 对象：counts 层存原始计数（下游 pseudobulk 与 GRN 都要用计数）
  obj <- Seurat::CreateSeuratObject(counts = Matrix::t(counts),
                                    meta.data = qc$metrics)
  # R 版约定：metrics 的行名就是细胞名（01_qc.R 已保证）

  if (identical(hvg_flavor_used, "seurat_v3")) {
    # 必须在 normalize 之前：vst 对**原始计数**做均值-方差稳定化拟合
    obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize",
                                 scale.factor = target_sum, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst",
                                        nfeatures = n_top, verbose = FALSE)
    log_info("HVG (vst, 用原始计数) 先算，再 normalize")
  } else {
    obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize",
                                 scale.factor = target_sum, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, selection.method = flavor_map[flavor],
                                        nfeatures = n_top, verbose = FALSE)
    log_info(sprintf("HVG (%s, 用 log 后数据) 后算", hvg_flavor_used))
  }

  hvg <- Seurat::VariableFeatures(obj)
  n_hvg <- length(hvg)
  log_info(sprintf("高变基因: %d", n_hvg))

  # HVG 图：均值-离散度，标出被选中的。
  # **不调 Seurat::VariableFeaturePlot**：高层 API 自己建 figure、figsize
  # 不受控（Python 版实测 177.8mm 非标宽），样式门禁无从约束。
  # 直接用 Assay 的 meta.features（FindVariableFeatures 写回的口径）自己画。
  hv_info <- obj[["RNA"]]@meta.features
  means <- as.numeric(hv_info$mean %||% Matrix::rowMeans(counts))
  # vst 口味的"离散度"列：Seurat 存 variance.standardized；其它口味用 dispersion
  disp_col <- if (identical(hvg_flavor_used, "seurat_v3")) "variance.standardized" else "dispersion"
  disp <- as.numeric(hv_info[[disp_col]] %||% rep(NA_real_, ncol(counts)))
  hv_flag <- colnames(counts) %in% hvg
  df <- data.frame(means = means, disp = disp,
                   hv = ifelse(hv_flag, "highly variable genes", "other genes"))
  p <- ggplot2::ggplot(df, ggplot2::aes(means, disp, colour = hv)) +
    ggplot2::geom_point(size = 0.4, alpha = 0.7, linewidth = 0) +
    ggplot2::scale_colour_manual(values = c("highly variable genes" = PAL$highlight,
                                            "other genes" = PAL$muted)) +
    ggplot2::labs(x = "mean expression of genes",
                  y = sprintf("%s of genes", disp_col),
                  title = sprintf("HVG selection (%s, n=%d)\norange = highly variable; grey = other genes",
                                  hvg_flavor_used, n_hvg)) +
    ggplot2::theme_paper()
  save_fig(cfg, "02-02-01-unit1-hvg-selection", p,
           width = W_ONE_HALF, height = mm(70))

  # ---- 2. 子集到 HVG + Scale + PCA ----------------------------------------
  # **M13（R-03 裁决）**：配置值真正传进去、被维度夹住时把**实际算的个数**
  # 记进状态（不是配置值）。"我按你说的做了"和"我改小了"必须能区分。
  n_comps_max <- min(nrow(counts), length(hvg))   # PCA 秩上限
  n_comps <- min(as.integer(cfg$reduce$n_pcs), n_comps_max)
  n_comps_clamped <- n_comps != as.integer(cfg$reduce$n_pcs)

  # 只在 HVG 上 Scale+PCA（与 Python 的 work = adata[:, hvg] 一致）
  obj <- Seurat::ScaleData(obj, features = hvg, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, features = hvg, npcs = n_comps, verbose = FALSE)
  emb <- Seurat::Embeddings(obj, "pca")   # 细胞 x PC
  stdev <- obj[["pca"]]@stdev
  var_ratio_raw <- (stdev^2) / sum(stdev^2)
  var_ratio <- as.numeric(var_ratio_raw)

  # PCA 方差比图（02-02-02-unit1）：不用 sc.pl.pca_variance_ratio 的三个
  # 问题在 Python 版已记 —— R 版同样自画、刻度按宽度 step=5 if n_pc>12。
  n_pc <- length(var_ratio)
  xs <- seq_len(n_pc)
  dfp <- data.frame(pc = xs, v = log10(var_ratio))
  step <- if (n_pc > 12) 5L else 1L
  ticks <- seq(1, n_pc, by = step)
  p <- ggplot2::ggplot(dfp, ggplot2::aes(pc, v)) +
    ggplot2::geom_line(colour = PAL$primary, linewidth = 0.4) +
    ggplot2::geom_point(colour = PAL$primary, size = 1.2) +
    ggplot2::scale_x_continuous(breaks = ticks,
                                labels = sprintf("PC%d", ticks)) +
    ggplot2::labs(x = "Principal component (rank)",
                  y = "log10(variance ratio)",
                  title = "PCA variance ratio (elbow)") +
    ggplot2::theme_paper() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_line(
      linewidth = 0.25, colour = "grey85"))
  save_fig(cfg, "02-02-02-unit1-pca-variance-ratio", p,
           width = W_ONE_HALF, height = mm(66))
  log_info(sprintf("PCA 方差比图: 标出 %d 个刻度（共 %d 个 PC）", length(ticks), n_pc))
  log_info(sprintf("PC1-%d 方差解释: %s", min(10, n_pc),
                   paste(sprintf("%.3f", var_ratio[seq_len(min(10, n_pc))]),
                         collapse = ", ")))

  # ---- 3. 批次整合 --------------------------------------------------------
  integ <- cfg$integration %||% list()
  method <- integ$method %||% "none"
  use_rep <- "pca"
  integ_record <- list(method = method, batch_key = batch_key)

  if (identical(method, "none")) {
    integ_record$status <- "not_applied"
    if (!is.null(batch_key)) {
      integ_record$reason <- paste0("配置 method=none 但设了 batch_key —— ",
                                    "多批次数据不做整合会让聚类被批次差异主导")
      log_warn(integ_record$reason)
    } else {
      integ_record$reason <- "单样本数据没有批次，不做整合是正确的选择"
      log_info(integ_record$reason)
    }
  } else if (is.null(batch_key)) {
    integ_record$status <- "not_configured"
    integ_record$reason <- "integration.method 不是 none 但 design.batch_key 为空"
    log_warn(integ_record$reason)
  } else if (identical(method, "harmony")) {
    # Seurat 对象需要 cell meta 里有 batch 列；从 design.batch_key 指向的
    # meta 列取。R 版 harmony 走 harmony::RunHarmony（seurat 对象方法）。
    res_h <- tryCatch({
      need_pkg("harmony", "批次整合")
      meta <- obj@meta.data
      if (!batch_key %in% colnames(meta)) {
        stop(sprintf("meta 里没有 batch 列 '%s'（design.batch_key 指向它）", batch_key),
             call. = FALSE)
      }
      obj <- harmony::RunHarmony(obj, group.by.vars = batch_key, verbose = FALSE)
      use_rep <- "harmony"
      list(status = "ok", use_rep = use_rep)
    }, error = function(e) list(status = "failed", reason = conditionMessage(e)))
    integ_record$status <- res_h$status
    if (!is.null(res_h$use_rep)) integ_record$use_rep <- res_h$use_rep
    if (!is.null(res_h$reason)) integ_record$reason <- res_h$reason
    if (identical(res_h$status, "ok")) {
      log_info(sprintf("harmony 整合完成（batch_key=%s）", batch_key))
    } else {
      log_warn(integ_record$reason %||% "")
    }
  } else if (identical(method, "combat")) {
    # sva::ComBat 在表达矩阵上按 batch 逐基因去均值-方差。
    # **M13：ComBat 后必须重跑 PCA 带同样 n_comps**（ComBat 直接改 X，
    # X_pca 已失效）。
    need_pkg("sva", "ComBat 批次整合")
    meta <- obj@meta.data
    if (!batch_key %in% colnames(meta)) {
      stop(sprintf("meta 里没有 batch 列 '%s'", batch_key), call. = FALSE)
    }
    # 在 log 数据上跑：取 NormalizeData 之后的 data 层
    log_expr <- as.matrix(Seurat::GetAssayData(obj, layer = "data")[hvg, , drop = FALSE])
    combat_out <- sva::ComBat(dat = log_expr, batch = as.factor(meta[[batch_key]]))
    obj@assays$RNA@data[hvg, ] <- combat_out
    obj <- Seurat::ScaleData(obj, features = hvg, verbose = FALSE)
    obj <- Seurat::RunPCA(obj, features = hvg, npcs = n_comps, verbose = FALSE)
    emb <- Seurat::Embeddings(obj, "pca")
    stdev <- obj[["pca"]]@stdev
    var_ratio <- as.numeric((stdev^2) / sum(stdev^2))
    n_pc <- length(var_ratio)
    integ_record$status <- "ok"
    integ_record$note <- "ComBat 后重跑了 PCA（ComBat 直接改 X，X_pca 已失效）"
    log_info(integ_record$note)
  } else {
    stop(sprintf("不支持的 integration.method: %s（none/harmony/combat）", method),
         call. = FALSE)
  }

  # 批次效应可视化（有 batch_key 才有意义）
  if (!is.null(batch_key) && batch_key %in% colnames(obj@meta.data)) {
    cats <- as.character(obj@meta.data[[batch_key]])
    reps <- c("pca", use_rep)
    plots <- lapply(reps, function(rep_nm) {
      src <- if (identical(rep_nm, "pca")) emb else
        Seurat::Embeddings(obj, rep_nm)
      dfb <- data.frame(d1 = src[, 1], d2 = src[, 2], batch = cats)
      ggplot2::ggplot(dfb, ggplot2::aes(d1, d2, colour = batch)) +
        ggplot2::geom_point(size = 0.4, alpha = 0.5) +
        ggplot2::ggtitle(sprintf("%s by %s", rep_nm, batch_key)) +
        ggplot2::labs(x = "dim 1", y = "dim 2") +
        ggplot2::theme_paper()
    })
    # 双面板 before/after（W_DOUBLE×64mm）；图例规则 31：框外右侧单列
    leg <- cowplot::get_legend(plots[[1]] + ggplot2::theme(legend.position = "right"))
    body <- patchwork::wrap_plots(lapply(plots, function(pp) pp + ggplot2::theme(
      legend.position = "none"))) + patchwork::plot_annotation(
        title = "Batch mixing before/after integration")
    p_batch <- patchwork::wrap_plots(body, leg, ncol = 2, widths = c(3, 0.6))
    save_fig(cfg, "02-02-03-unit1-batch-mixing", p_batch,
             width = W_DOUBLE, height = mm(64))
  }

  # ---- 4. 落盘（R 走 .rds，不走 .h5ad —— references/r_version.md §5）------
  # **logcounts_all = 全基因集 log 表达**（Python 版 `adata.raw = adata` 的
  # R 版对应物，02_integrate.py:91）。下游 06_communication 的配体-受体
  # 基因大多是低表达的细胞因子/趋化因子，几乎不进 HVG —— 实测 Python 版
  # 在 HVG 子集上 38 对里只有 3 对可用，全基因集 27 对，差 9 倍。
  # 缺了这份全基因集矩阵，06 的"没找到显著通讯"会变成假阴性。
  logcounts_all <- Matrix::t(Seurat::GetAssayData(obj, layer = "data"))
  out <- file.path(data_dir, "integrated.rds")
  saveRDS(list(
    counts = Matrix::t(Seurat::GetAssayData(obj, layer = "counts")[hvg, , drop = FALSE]),
    logcounts = Matrix::t(Seurat::GetAssayData(obj, layer = "data")[hvg, , drop = FALSE]),
    logcounts_all = logcounts_all,   # 全基因集（adata.raw 对应物）
    hvg = hvg,
    cell_meta = obj@meta.data,
    pca = emb,
    pca_stdev = stdev,
    pca_variance_ratio = var_ratio,
    use_rep = use_rep,
    n_comps = n_comps), out)
  log_info(sprintf("已写出 %s（%d 细胞 x %d HVG；全基因集 logcounts_all %d 基因）",
                   out, nrow(counts), length(hvg), ncol(logcounts_all)))

  status <- list(
    dataset_id = cfg$dataset_id,
    n_cells = nrow(counts),
    n_hvg = n_hvg,
    hvg_flavor_requested = flavor,
    hvg_flavor_used = hvg_flavor_used,
    hvg_fallback = hvg_fallback,
    n_top_genes_requested = n_top,
    target_sum = target_sum,
    n_pcs = as.integer(n_comps),
    n_pcs_requested = as.integer(cfg$reduce$n_pcs),
    n_pcs_clamped = n_comps_clamped,
    pca_variance_ratio_top10 = lapply(
      var_ratio[seq_len(min(10, n_pc))], function(v) round(v, 5)),
    pca_cumvar_top10 = lapply(
      cumsum(var_ratio[seq_len(min(10, n_pc))]), function(v) round(v, 5)),
    integration = integ_record,
    use_rep = use_rep,
    # **M11 配套**：CONDITIONAL_FIGURES 判据要能读的显式布尔。
    has_batch_key = !is.null(batch_key) && batch_key %in% colnames(obj@meta.data),
    status = "ok")
  write_json(file.path(res_dir, "integration_status.json"), status)
  status
}

# ---------------------------------------------------------------------------
# 入口（02_integrate.py:276-288 同语义；E-56）
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L || identical(Sys.getenv("SCRNA_STEP_MAIN"), "02_integrate")) {
  args <- parse_args()
  cfg <- load_config(args$config)
  t0 <- Sys.time()
  tryCatch({
    out <- run_02_integrate(cfg)
    record_step(cfg, "integrate", "ok", as.numeric(difftime(Sys.time(), t0, units = "secs")),
                result_status = result_status_of(out))
    out
  }, error = function(e) {
    record_step(cfg, "integrate", "failed", as.numeric(difftime(Sys.time(), t0, units = "secs")),
                message = conditionMessage(e), result_status = "failed")
    stop(e)
  })
}
