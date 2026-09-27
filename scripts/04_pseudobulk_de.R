#!/usr/bin/env Rscript
# =============================================================================
# 04_pseudobulk_de.R — 拟bulk 差异表达（DESeq2，R 版，与 04_pseudobulk_de.py 同号同语义）
#
# **为什么必须做拟bulk而不是直接在细胞层面跑检验。**
# 把每个细胞当一个独立样本，是把"细胞数"当成了"样本量"。
# 10000 个细胞来自 3 个供体，自由度不是 9997 而是 2 —— 而 Wilcoxon
# 会给出 p = 1e-300 这样的数字。这是单细胞差异分析里最严重的系统性错误，
# **而且它给出的假阳性看起来极其显著**。
#
# 方法学等级（references/r_version.md §3）：**A 同方法不同实现**。
# pydeseq2 本来就是 DESeq2 的 Python 重实现 —— R 版用原版 DESeq2，
# 方向相反的两版在这里汇合（数值上原版才是"真"参考）。
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

# 一个 (样本, 细胞类型) 组合至少要有这么多细胞才拿去做拟bulk
MIN_CELLS_PER_PSEUDOBULK <- 10L
# 一个细胞类型至少要在这么多个样本里出现，才能做组间比较
MIN_SAMPLES_PER_CELLTYPE <- 2L
# 每组至少这么多个生物学重复 —— 1 个重复算不出组内方差
MIN_REPLICATES_PER_GROUP <- 2L

# ---------------------------------------------------------------------------
# 聚合（04_pseudobulk_de.py:47-95 同语义）
#
# **必须用 counts，不能用 log 后的表达。** DESeq2 的负二项模型要求原始
# 计数 —— 把 log 值喂进去不会报错，只会给出错的离散度估计。
# **没有 counts 时不再静默退回 log 数据**（审计 S6）：返回
# counts_source="missing"，由调用方判成失败写进 status。
# ---------------------------------------------------------------------------
aggregate_pseudobulk <- function(counts, cell_meta, sample_key, celltype_key) {
  # counts: 细胞 x 基因（稀疏）；cell_meta: data.frame 带样本/分组列
  if (is.null(counts)) {
    return(list(mat = NULL, meta = NULL, agg = list(
      reason = paste0("输入没有原始计数 —— DESeq2 的负二项模型要求原始计数，",
                      "喂 log 值不会报错但离散度估计是错的。**拒绝用 log 数据",
                      "冒充计数**；请确认上游 01_qc / 02_integrate 保留了 counts"),
      counts_source = "missing")))
  }
  obs <- data.frame(
    sample = as.character(cell_meta[[sample_key]]),
    celltype = as.character(cell_meta[[celltype_key]]),
    stringsAsFactors = FALSE)
  obs$row <- seq_len(nrow(obs))

  mats <- list(); meta_rows <- list(); k <- 0L
  for (key in split(seq_len(nrow(obs)), paste(obs$sample, obs$celltype, sep = "|"))) {
    g <- obs[key[1L], ]
    s <- g$sample; ct <- g$celltype
    if (length(key) < MIN_CELLS_PER_PSEUDOBULK) next
    k <- k + 1L
    mats[[k]] <- Matrix::colSums(counts[key, , drop = FALSE])
    meta_rows[[k]] <- data.frame(sample = s, celltype = ct,
                                 n_cells = length(key),
                                 row.names = paste(s, ct, sep = "|"))
  }
  if (k == 0L) {
    return(list(mat = NULL, meta = NULL, agg = list(
      reason = sprintf("没有任何 (样本, 细胞类型) 组合达到 %d 个细胞",
                       MIN_CELLS_PER_PSEUDOBULK),
      counts_source = "raw counts")))
  }
  mat <- do.call(rbind, lapply(mats, as.numeric))
  meta <- do.call(rbind, meta_rows)
  list(mat = mat, meta = meta,
       agg = list(n_pseudobulk_samples = k, counts_source = "raw counts"))
}

# ---------------------------------------------------------------------------
# 主流程（04_pseudobulk_de.py:98-304 同语义）
# ---------------------------------------------------------------------------
run_04_pseudobulk_de <- function(cfg) {
  ensure_dirs(cfg)
  set_seed(cfg)
  data_dir <- cfg$output$data_dir
  res_dir <- cfg$output$results_dir

  design <- cfg$design %||% list()
  sample_key <- design$sample_key
  group_key <- design$group_key
  ref_group <- design$reference_group

  status <- list(dataset_id = cfg$dataset_id, status = "not_configured")
  write_status <- function(s) {
    write_json(file.path(res_dir, "pseudobulk_status.json"), s)
  }

  if (is.null(sample_key) || is.null(group_key)) {
    status$reason <- paste0(
      "design.sample_key / design.group_key 为空 —— 拟bulk DE 需要知道",
      "「哪些细胞来自同一个生物学重复」和「怎么分组」。",
      "单细胞数据的样本量单位是供体/病人，不是细胞；",
      "没有这两列就无法构造拟bulk，也不应退化成细胞层面的检验",
      "（那会把细胞数当成样本量，给出 p=1e-300 的假阳性）")
    log_warn(sprintf("跳过拟bulk DE: %s", status$reason))
    write_status(status)
    return(status)
  }

  pkg_ok <- tryCatch({ need_pkg("DESeq2", "拟bulk 差异表达"); TRUE },
                     error = function(e) FALSE)
  if (!pkg_ok) {
    status$status <- "package_missing"
    status$reason <- "DESeq2 未安装（拟bulk 差异表达的负二项检验）"
    log_warn(status$reason)
    write_status(status)
    return(status)
  }

  clu_path <- file.path(data_dir, "clustered.rds")
  if (!file.exists(clu_path)) {
    stop(sprintf("缺输入 %s —— 先跑 03_cluster_annotate.R", clu_path), call. = FALSE)
  }
  clu <- readRDS(clu_path)
  cell_meta <- clu$cell_meta
  for (kv in list(list(sample_key, "sample_key"), list(group_key, "group_key"))) {
    if (!kv[[1]] %in% colnames(cell_meta)) {
      status$status <- "column_missing"
      status$reason <- sprintf("design.%s='%s' 不在 meta 里；实际列: %s",
                               kv[[2]], kv[[1]],
                               paste(utils::head(sort(colnames(cell_meta)), 40),
                                     collapse = ", "))
      log_warn(status$reason)
      write_status(status)
      return(status)
    }
  }

  celltype_key <- if ("celltype" %in% colnames(cell_meta)) "celltype" else "leiden"
  log_info(sprintf("拟bulk: sample_key=%s, group_key=%s, celltype_key=%s",
                   sample_key, group_key, celltype_key))

  agg_res <- aggregate_pseudobulk(clu$counts, cell_meta, sample_key, celltype_key)
  mat <- agg_res$mat; meta <- agg_res$meta; agg <- agg_res$agg
  if (is.null(mat)) {
    # counts_source=="missing" 是**上游契约被破坏**，不是"这批数据不适合
    # 做拟bulk" —— 两者必须长得不一样（审计 S6）。
    if (identical(agg$counts_source, "missing")) {
      status <- utils::modifyList(status, c(list(status = "missing_counts"), agg))
    } else {
      status <- utils::modifyList(status, c(list(status = "no_usable_pseudobulk"), agg))
    }
    log_warn(sprintf("拟bulk 不可行: %s", agg$reason))
    write_status(status)
    return(status)
  }

  # ---- 样本 -> 分组 的映射（审计 M1）---------------------------------------
  # 逐样本枚举真实取值集合；只有恰好 1 个取值才接受，否则记进
  # mixed_group_samples 并让整步以 mixed_group 收尾。
  # **错的结果看起来和真的完全一样 —— 这比抛错危险得多。**
  meta$sample <- as.character(meta$sample)
  sample_groups <- c()
  mixed_group_samples <- list()
  g_by_s <- split(as.character(cell_meta[[group_key]]),
                  as.character(cell_meta[[sample_key]]))
  for (s in unique(meta$sample)) {
    vals <- sort(unique(g_by_s[[s]]))
    if (length(vals) == 1L) {
      sample_groups[s] <- vals
    } else {
      mixed_group_samples[[length(mixed_group_samples) + 1L]] <-
        list(sample = s, groups = as.list(vals))
    }
  }
  if (length(mixed_group_samples)) {
    status$status <- "mixed_group"
    status$reason <- sprintf("样本与分组不是一对一：%d 个样本跨多个分组",
                             length(mixed_group_samples))
    status$mixed_group_samples <- mixed_group_samples
    status$n_samples_checked <- length(unique(meta$sample))
    log_warn(sprintf("拟bulk 不可行：%d 个样本跨多个分组（样本-分组契约被破坏，拒绝按第一个猜分组）",
                     length(mixed_group_samples)))
    write_status(status)
    return(status)
  }
  meta$group <- unname(sample_groups[meta$sample])
  log_info(sprintf("拟bulk 矩阵: %d 个 (样本 x 细胞类型) x %d 基因",
                   nrow(mat), ncol(mat)))

  # ---- 每个细胞类型单独做 --------------------------------------------------
  need_pkg("DESeq2", "拟bulk 差异表达")
  results <- list(); skipped <- list()
  contrasts_used <- list()
  for (ct in sort(unique(meta$celltype))) {
    sub <- meta[meta$celltype == ct, , drop = FALSE]
    if (length(unique(sub$group)) < 2L) {
      skipped[[length(skipped) + 1L]] <- list(celltype = ct, reason = "只有一个分组取值")
      next
    }
    tab <- table(sub$group)
    thin <- tab[tab < MIN_REPLICATES_PER_GROUP]
    if (length(thin) > 0L) {
      skipped[[length(skipped) + 1L]] <- list(
        celltype = ct,
        reason = sprintf("组内生物学重复不足 %d: %s", MIN_REPLICATES_PER_GROUP,
                         paste(names(thin), as.character(thin), sep = "=", collapse = ", ")))
      next
    }
    if (nrow(sub) < MIN_SAMPLES_PER_CELLTYPE * 2L) {
      skipped[[length(skipped) + 1L]] <- list(
        celltype = ct, reason = sprintf("拟bulk 样本只有 %d 个", nrow(sub)))
      next
    }

    idx <- match(rownames(sub), rownames(mat))
    cts <- round(mat[idx, , drop = FALSE])
    colnames(cts) <- colnames(mat)
    rownames(cts) <- rownames(sub)
    # DESeq2 不接受全零基因
    cts <- cts[, colSums(cts) > 0, drop = FALSE]

    # ---- contrast 的分子/分母必须**显式定死**，不能靠行序（审计 S5）----
    # 分组水平按**排序后的字典序**定死（可复现、与行序无关），
    # 再把 numerator / reference 都记进 status 供事后核对。
    levels_all <- sort(unique(as.character(sub$group)))
    if (!is.null(ref_group) && !as.character(ref_group) %in% levels_all) {
      # 指定的参考组在这个细胞类型里不存在 —— DESeq2 会抛错；
      # 显式记成"跳过 + 原因"，与"算失败"区分开（审计 M20）。
      skipped[[length(skipped) + 1L]] <- list(
        celltype = ct,
        reason = sprintf("指定的 reference_group='%s' 在该细胞类型里不存在（实际水平: %s）",
                         ref_group, paste(levels_all, collapse = ", ")))
      log_warn(sprintf("  %s 跳过：reference_group='%s' 不在 [%s]",
                       ct, ref_group, paste(levels_all, collapse = ", ")))
      next
    }
    if (is.null(ref_group)) {
      numerator <- levels_all[length(levels_all)]
      reference <- levels_all[1L]
    } else {
      numerator <- levels_all[levels_all != as.character(ref_group)][1L]
      reference <- as.character(ref_group)
    }
    coldata <- data.frame(group = factor(sub$group, levels = levels_all),
                          row.names = rownames(sub))

    res_ct <- tryCatch({
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = cts, colData = coldata, design = ~group)
      dds <- DESeq2::DESeq(dds, quiet = TRUE)
      res <- DESeq2::results(dds, contrast = c("group", numerator, reference))
      df <- data.frame(gene = rownames(res), as.data.frame(res),
                       stringsAsFactors = FALSE, row.names = NULL)
      df$celltype <- ct   # 第一列语义；构造顺序不敏感（写 CSV 前重排）
      df <- df[, c("celltype", names(df)[names(df) != "celltype"])]
      n_sig <- if ("padj" %in% names(df)) sum(df$padj < 0.05, na.rm = TRUE) else 0L
      log_info(sprintf("  %s: %d 个拟bulk 样本, contrast=%s vs %s, %d 个 padj<0.05",
                       ct, nrow(sub), numerator, reference, n_sig))
      list(df = df, contrast = list(numerator = numerator,
                                    reference = reference,
                                    levels = as.list(levels_all)))
    }, error = function(e) {
      skipped[[length(skipped) + 1L]] <<- list(
        celltype = ct, reason = sprintf("%s: %s", class(e)[1L], conditionMessage(e)))
      log_warn(sprintf("  %s 的 DESeq2 失败: %s", ct, conditionMessage(e)))
      NULL
    })
    if (!is.null(res_ct)) {
      results[[length(results) + 1L]] <- res_ct$df
      contrasts_used[[ct]] <- res_ct$contrast
    }
  }

  if (length(results)) {
    allres <- do.call(rbind, results)
    utils::write.csv(allres, file.path(res_dir, "pseudobulk_de.csv"), row.names = FALSE)
    utils::write.csv(meta, file.path(res_dir, "pseudobulk_samples.csv"), row.names = TRUE)
    status$n_genes_tested <- nrow(allres)
    status$n_celltypes_tested <- length(results)
    # **M20**：有跳过就记 `partial`，让"结果只覆盖了一部分"在顶层可见
    # （不判红 —— 跳过本身是设计如此，判红会让每个 job 都红）。
    status$status <- if (length(skipped)) "partial" else "ok"
    if (length(skipped)) {
      status$partial_reason <- sprintf(
        "%d/%d 个细胞类型被跳过，差异表只覆盖 %d 个",
        length(skipped), length(skipped) + length(results), length(results))
    }
  } else {
    status$status <- "no_celltype_testable"
  }

  status$sample_key <- sample_key
  status$group_key <- group_key
  status$celltype_key <- celltype_key
  status$reference_group <- ref_group
  # 实际用到的对比方向（审计 S5）：没有它就无法事后核对
  # log2FoldChange 的符号指的是哪个方向。
  status$contrast_used <- contrasts_used
  status$contrast_rule <- paste0(
    "分组水平按**字典序**定死；`reference_group` 为空时",
    "numerator=最大的水平、reference=最小的水平。**不看 obs 行序**")
  status$counts_source <- agg$counts_source
  status$min_cells_per_pseudobulk <- MIN_CELLS_PER_PSEUDOBULK
  status$min_replicates_per_group <- MIN_REPLICATES_PER_GROUP
  status$n_pseudobulk_samples <- nrow(mat)
  status$skipped_celltypes <- skipped
  status$design_note <- paste0(
    "样本量单位是生物学重复（供体/病人），不是细胞 —— ",
    "把细胞当样本会给出 p=1e-300 的假阳性")
  write_status(status)
  status
}

# ---------------------------------------------------------------------------
# 入口（04_pseudobulk_de.py:307-321 同语义；E-56：本步骤 5 条早退路径全不抛异常）
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L || identical(Sys.getenv("SCRNA_STEP_MAIN"), "04_pseudobulk_de")) {
  args <- parse_args()
  cfg <- load_config(args$config)
  t0 <- Sys.time()
  tryCatch({
    out <- run_04_pseudobulk_de(cfg)
    record_step(cfg, "pseudobulk_de", "ok",
                as.numeric(difftime(Sys.time(), t0, units = "secs")),
                result_status = result_status_of(out))
    out
  }, error = function(e) {
    record_step(cfg, "pseudobulk_de", "failed",
                as.numeric(difftime(Sys.time(), t0, units = "secs")),
                message = conditionMessage(e), result_status = "failed")
    stop(e)
  })
}
