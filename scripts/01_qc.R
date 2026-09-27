#!/usr/bin/env Rscript
# =============================================================================
# 01_qc.R — 质量控制与过滤（R 版，与 01_qc.py 同号同语义）
#
# 做四件事，每一件都留记录：
#   1. QC 指标（n_genes / total_counts / pct_mt / pct_ribo / pct_hb）
#   2. 硬阈值过滤（可配）
#   3. 双细胞检测（scDblFinder）—— **失败不等于跳过，要写明为什么**
#   4. ambient RNA —— 默认不做，**并写明不做**
#
# 方法学等级（references/r_version.md §3）：**C 真方法差异**。
# scrublet 是"模拟双体 + 分类器"，scDblFinder 是"人工双体 + 建模"；
# 检出数量与阈值口径**会不同** —— qc_status.json 里 doublet_detection
# 的字段两版不可直接比（R-17 验收时按字段语义对照，不按数值对照）。
#
# 逐条复刻的 Python 语义（行号以 01_qc.py 为准）：
#   - M15：血红蛋白用 assets/qc_gene_sets.yml 的**精确基因名**（前缀匹配
#     会误收 HBEGF/RPSA，指标偏了而量级不变看不出异常）；文件缺失 raise。
#   - S9：期望双细胞率按 10x 经验式 0.008*n/1000 反推、夹 [0.01,0.10]，
#     旧公式值一并落盘对照。
#   - D-006 单图原则：五联小提琴拆 5 张独立单图（图名逐字相同）。
#   - M23：过滤链两个维度各自计数、名字与语义一致。
#   - ambient：not_done / heuristic_only 如实写明，不冒充校正。
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
# 脚本定位：独立 Rscript 用 --file=；被 main_analysis.R source 时由入口
# 先设置 options(scrna.script_path=...)（R 没有 __file__，这是最小替代）。
# ---------------------------------------------------------------------------
get_script_path <- function() {
  argv0 <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(argv0)) return(normalizePath(sub("^--file=", "", argv0[1]), mustWork = FALSE))
  getOption("scrna.script_path", "")
}

QC_GENE_SETS_PATH <- file.path(dirname(dirname(get_script_path())), "assets",
                               "qc_gene_sets.yml")

load_qc_gene_sets <- function(organism = "Homo sapiens") {
  if (!file.exists(QC_GENE_SETS_PATH)) {
    stop(sprintf(paste0("缺少 QC 基因集定义 %s —— 它定义 pct_counts_hb / ",
                        "pct_counts_ribo 的语义，不能缺省。（旧版本用前缀匹配，",
                        "会误收 HBEGF / RPSA，见审计 M15）"), QC_GENE_SETS_PATH),
         call. = FALSE)
  }
  need_pkg("yaml", "读 assets/qc_gene_sets.yml")
  doc <- yaml::read_yaml(QC_GENE_SETS_PATH) %||% list()
  is_mouse <- grepl("^(mus|mouse)", tolower(organism))
  pre <- if (is_mouse) "mouse" else "human"
  list(hb = as.character(doc[[paste0(pre, "_hb_genes")]] %||% character(0)),
       ribo = as.character(doc[[paste0(pre, "_ribo_genes")]] %||% character(0)),
       source = "assets/qc_gene_sets.yml")
}

# 旧前缀规则（**只用于记录差异，不再参与判据**）
LEGACY_HB_PREFIXES <- c("HBA", "HBB", "HBD", "HBE", "HBG", "HBM", "HBQ", "HBZ")
LEGACY_RIBO_PREFIXES <- c("RPS", "RPL")

# ---------------------------------------------------------------------------
# QC 指标（01_qc.py:79-131 同语义；Seurat 实现替代 scanpy）
# ---------------------------------------------------------------------------
add_qc_metrics <- function(counts, organism = "Homo sapiens") {
  genes <- colnames(counts)   # R 版矩阵是 细胞 x 基因
  is_mouse <- grepl("^(mus|mouse)", tolower(organism))
  mt_pref <- if (is_mouse) "^mt-" else "^MT-"

  sets <- load_qc_gene_sets(organism)
  hb_genes <- intersect(sets$hb, genes)
  # 核糖体仍走前缀（家族 80+ 成员命名规整），但**误收的 RPSA/Rpsa 排除**
  ribo_hit <- grepl("^(RPS|RPL)", genes) & !(genes %in% c("RPSA", "Rpsa"))
  hb_hit <- genes %in% hb_genes
  mt_hit <- grepl(mt_pref, genes)

  # **counts 布局 细胞 x 基因：每细胞指标 = 行聚合（rowSums），不是 colSums**
  # —— 之前误写 colSums 得到"每基因"长度 32738 的向量，data.frame 长度恰好
  # 匹配 row.names=colnames(counts) 构造成功但语义全反，下游
  # counts[keep_cell, ]（keep_cell 长 32738、行数 2700）报
  # "logical subscript too long"（CI run8 实跑抓到）。
  # run6 曾把 row.names 从 rownames(counts)（细胞名 2700）改成
  # colnames(counts)（基因名 32738）——方向修反：metrics 每行=每细胞，
  # 行名就该是 rownames(counts)=细胞名。两处一起修正。
  lib_size <- Matrix::rowSums(counts)
  n_genes <- Matrix::rowSums(counts > 0)
  n_detected <- function(hit) {
    if (!any(hit)) return(rep(0, nrow(counts)))
    Matrix::rowSums(counts[, hit, drop = FALSE])
  }
  pct <- function(num) ifelse(lib_size > 0, 100 * num / pmax(lib_size, 1), 0)

  list(metrics = data.frame(
         n_genes_by_counts = as.numeric(n_genes),
         total_counts = as.numeric(lib_size),
         pct_counts_mt = as.numeric(pct(n_detected(mt_hit))),
         pct_counts_ribo = as.numeric(pct(n_detected(ribo_hit))),
         pct_counts_hb = as.numeric(pct(n_detected(hb_hit))),
         # row.names 是**细胞名**：counts 布局 细胞x基因，metrics 每行=每细胞。
         row.names = rownames(counts)),
       qc_vars = c("mt", "ribo", "hb")[c(any(mt_hit), any(ribo_hit), any(hb_hit))],
       gene_sets = list(
         source = sets$source,
         n_hb_genes_in_data = length(hb_genes),
         n_ribo_genes_in_data = sum(ribo_hit),
         hb_matching = "exact gene names",
         ribo_matching = "prefix RPS*/RPL* minus RPSA",
         legacy_prefix_extra_hb = sort(setdiff(
           genes[grepl(paste0("^(", paste(LEGACY_HB_PREFIXES, collapse = "|"), ")"),
                       genes)],
           hb_genes)),
         legacy_prefix_extra_ribo = sort(setdiff(genes[grepl("^(RPS|RPL)", genes)],
                                                 genes[ribo_hit])),
         why = paste0("前缀匹配会把 HBEGF（肝素结合 EGF 样生长因子）算进",
                      "血红蛋白、把 RPSA（67 kDa 层粘连蛋白受体前体）算进",
                      "核糖体 —— 指标偏了而量级不变，看不出异常")))
}

# ---------------------------------------------------------------------------
# 双细胞检测（01_qc.py:134-181 同语义；scDblFinder 替代 scrublet，等级 C）
# ---------------------------------------------------------------------------
run_doublet <- function(counts, metrics, cfg) {
  if (!isTRUE(cfg$qc$scrublet %||% TRUE)) {
    return(list(status = "disabled", reason = "配置 qc.scrublet=false"))
  }
  # scDblFinder 依赖 BiocParallel；显式单线程（workflow 的 env 已钉）。
  res <- tryCatch({
    need_pkg("scDblFinder", "双细胞检测")
    n <- nrow(counts)
    # 10x 经验式（审计 S9，两个公式都落盘对照）
    rate <- min(0.10, max(0.01, 0.008 * n / 1000.0))
    rate_old <- min(0.10, max(0.05, 5000 / max(n, 1) * 0.01))
    set.seed(as.integer(cfg$analysis$seed))
    # scDblFinder 要求**基因 x 细胞**（rows=features）；counts 是 细胞x基因，
    # 必须转置（colData 行=细胞，与 metrics 行序一致）。
    out <- scDblFinder::scDblFinder(Matrix::t(counts), dbr = rate, samples = NULL,
                                    BPPARAM = BiocParallel::SerialParam())
    # out 是 SingleCellExperiment；scDblFinder 列在 colData
    cd <- SummarizedExperiment::colData(out)
    is_db <- as.logical(cd$scDblFinder.class %in% c("doublet"))
    scores <- as.numeric(cd$scDblFinder.score)
    n_db <- sum(is_db)
    list(status = "ok",
         expected_doublet_rate = round(rate, 4),
         expected_doublet_rate_10x_rule = round(0.008 * n / 1000.0, 5),
         expected_doublet_rate_old_rule = round(rate_old, 4),
         rate_rule_note = paste0("`expected_doublet_rate` 按 10x 经验式 ",
                                 "0.008 x n/1000 反推，夹在 [0.01, 0.10]；",
                                 "`expected_doublet_rate_old_rule` 是旧公式",
                                 "（n>1000 时恒为下限 0.05）的值，留作对照；",
                                 "scDblFinder 的 dbr 参数接这个值"),
         n_predicted_doublets = n_db,
         frac_predicted_doublets = round(n_db / max(n, 1), 5),
         is_doublet = is_db, doublet_score = scores,
         tool = "scDblFinder",
         method_note = paste0("R 版用 scDblFinder（人工双体 + 建模），",
                              "Python 版用 scrublet（模拟双体 + 分类器）",
                              "—— 方法学等级 C，检出数量两版不可直接比"))
  }, error = function(e) {
    msg <- conditionMessage(e)
    if (grepl("need_pkg|BiocManager|not available|不存在包", msg, ignore.case = TRUE)) {
      list(status = "package_missing", reason = msg,
           fix = "BiocManager::install('scDblFinder')")
    } else {
      list(status = "failed", reason = msg)
    }
  })
  res
}

# ---------------------------------------------------------------------------
# ambient RNA（01_qc.py:184-216 同语义）
# ---------------------------------------------------------------------------
assess_ambient_rna <- function(metrics, cfg) {
  mode <- cfg$qc$ambient_rna %||% "none"
  if (identical(mode, "none")) {
    return(list(status = "not_done",
                reason = paste0("配置 qc.ambient_rna=none。正规做法需要 SoupX（R，需空液滴）",
                                "或 CellBender（需 GPU）；GitHub 托管 runner 无 GPU，",
                                "故未执行，也未用代理指标冒充校正结果")))
  }
  if (identical(mode, "simple")) {
    out <- list(status = "heuristic_only",
                note = paste0("**这不是 ambient RNA 校正**，只是背景水平提示；",
                              "要真正校正需 SoupX/CellBender"))
    v <- metrics$pct_counts_hb
    if (!is.null(v)) {
      out$pct_counts_hb_median <- round(median(v), 4)
      out$pct_counts_hb_p95 <- round(as.numeric(quantile(v, 0.95, na.rm = TRUE)), 4)
    }
    return(out)
  }
  list(status = "unknown_mode", reason = sprintf("qc.ambient_rna='%s' 不是 none/simple", mode))
}

# ---------------------------------------------------------------------------
# 过滤前 QC 图（01_qc.py:250-311 同语义；六个图名逐字相同）
# ---------------------------------------------------------------------------
QC_LABELS <- c(
  n_genes_by_counts = "Genes detected per cell",
  total_counts = "Total counts per cell",
  pct_counts_mt = "Mitochondrial fraction (%)",
  pct_counts_ribo = "Ribosomal fraction (%)",
  pct_counts_hb = "Haemoglobin fraction (%)"
)
# 图名来自查表（字面量，过 check_fig_names 命名门禁；与 01_qc.py:272-278 逐字相同）
QC_UNIT_NAMES <- list(
  n_genes_by_counts = "02-01-01-unit1-genes-detected",
  total_counts = "02-01-01-unit2-total-counts",
  pct_counts_mt = "02-01-01-unit3-mito-fraction",
  pct_counts_ribo = "02-01-01-unit4-ribo-fraction",
  pct_counts_hb = "02-01-01-unit5-hb-fraction"
)

plot_qc_figures <- function(cfg, metrics) {
  need_pkg("ggplot2")
  keys <- intersect(c("n_genes_by_counts", "total_counts", "pct_counts_mt",
                      "pct_counts_ribo", "pct_counts_hb"), colnames(metrics))
  # **单图原则拆分（D-006）**：每个指标一张独立单图（W_SINGLE×58mm）。
  # IQR=0 退化成 strip 散点并注明 no variance（与 Python 版一致）。
  for (k in keys) {
    v <- as.numeric(metrics[[k]])
    label <- unname(QC_LABELS[k]) %||% k
    iqr0 <- (as.numeric(quantile(v, 0.75)) - as.numeric(quantile(v, 0.25))) == 0
    if (iqr0) {
      df <- data.frame(x = runif(length(v), -0.12, 0.12), y = v)
      p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
        ggplot2::geom_point(size = 0.3, alpha = 0.25, colour = PAL$primary) +
        ggplot2::ggtitle(sprintf("%s\nno variance: %d/%d cells",
                                 label, sum(v != 0), length(v)))
    } else {
      df <- data.frame(x = k, y = v)
      p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
        ggplot2::geom_violin(fill = PAL$primary, alpha = 0.6,
                             draw_quantiles = c(0.25, 0.5, 0.75)) +
        ggplot2::ggtitle(label)
    }
    # theme_paper 是 common.R 自定义主题 —— **不是 ggplot2 导出**，
    # 加 ggplot2:: 前缀会报 "not an exported object"（CI run7 实跑抓到）。
    p <- p + theme_paper() +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank())
    save_fig(cfg, QC_UNIT_NAMES[[k]], p, width = W_SINGLE, height = mm(58))
  }

  # 阈值线：total_counts vs n_genes，色 = pct_mt（02-01-02-unit1）
  q <- cfg$qc
  df <- data.frame(total = metrics$total_counts, genes = metrics$n_genes_by_counts,
                   mt = metrics$pct_counts_mt)
  p <- ggplot2::ggplot(df, ggplot2::aes(total, genes, colour = mt)) +
    ggplot2::geom_point(size = 0.4, alpha = 0.6) +
    ggplot2::scale_colour_viridis_c(name = "Mitochondrial fraction (%)") +
    ggplot2::geom_hline(yintercept = as.numeric(q$min_genes),
                        colour = PAL$highlight, linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = as.numeric(q$max_genes),
                        colour = PAL$highlight, linetype = "dashed", linewidth = 0.4) +
    ggplot2::labs(x = unname(QC_LABELS["total_counts"]),
                  y = unname(QC_LABELS["n_genes_by_counts"]),
                  title = "QC thresholds (red = cut-offs)")
  save_fig(cfg, "02-01-02-unit1-qc-scatter-thresholds", p,
           width = W_SINGLE, height = mm(64))
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# 主流程（01_qc.py:219-410 同语义）
# ---------------------------------------------------------------------------
run_01_qc <- function(cfg) {
  ensure_dirs(cfg)
  set_seed(cfg)
  data_dir <- cfg$output$data_dir
  res_dir <- cfg$output$results_dir

  info <- read_json_or_none(file.path(data_dir, "dataset_info.json")) %||% list()
  organism <- info$organism %||% "Homo sapiens"

  raw_path <- file.path(data_dir, "raw.rds")
  if (!file.exists(raw_path)) {
    stop(sprintf("缺输入 %s —— 先跑 00_fetch.R", raw_path), call. = FALSE)
  }
  counts <- readRDS(raw_path)
  n0 <- nrow(counts); g0 <- ncol(counts)
  log_info(sprintf("读入 %d 细胞 x %d 基因", n0, g0))

  # ---- 1. 指标 ------------------------------------------------------------
  qc <- add_qc_metrics(counts, organism)
  metrics <- qc$metrics
  gene_sets <- qc$gene_sets
  log_info(sprintf("QC 指标已算（%s）",
                   if (length(qc$qc_vars)) paste(qc$qc_vars, collapse = ", ") else "无线粒体/核糖体基因命中"))
  if (length(gene_sets$legacy_prefix_extra_hb) || length(gene_sets$legacy_prefix_extra_ribo)) {
    log_info(sprintf("QC 基因集（%s）：旧前缀规则会多收 hb=[%s]、ribo=[%s] —— 本轮的指标不含它们",
                     gene_sets$source,
                     paste(gene_sets$legacy_prefix_extra_hb, collapse = ", "),
                     paste(gene_sets$legacy_prefix_extra_ribo, collapse = ", ")))
  }

  # ---- 2. 过滤前的图 ------------------------------------------------------
  # 先画再滤 —— 滤完再画就看不到"滤掉了什么"。
  plot_qc_figures(cfg, metrics)

  # ---- 3. 过滤（M23：两维度各自计数、名字与语义一致）----------------------
  # Python 版顺序：filter_cells(min_genes) -> filter_genes(min_cells)
  #   -> max_genes（按 n_genes 截）-> max_pct_mt。R 版同序。
  q <- cfg$qc
  # 细胞维度：min_genes（metrics 行 = 细胞，先行对齐免歧义）
  keep_cell <- metrics$n_genes_by_counts >= as.integer(q$min_genes)
  n_after_min_genes <- sum(keep_cell)
  g_after_filter_cells <- g0   # filter_cells 不改基因维度，但记下来
  counts <- counts[keep_cell, , drop = FALSE]
  metrics <- metrics[keep_cell, , drop = FALSE]

  # 基因维度：min_cells（counts 列 = 基因）
  gkeep <- Matrix::colSums(counts > 0) >= as.integer(q$min_cells)
  counts <- counts[, gkeep, drop = FALSE]
  g_after_min_cells <- ncol(counts)

  # 细胞维度续：max_genes / max_pct_mt（基因删了，total_counts 变 ——
  # 但 n_genes_by_counts / pct_counts_mt 语义按过滤前指标走，与 Python 版
  # 的 adata.obs 值一致：scanpy filter_genes 也不重算 obs 指标）
  if (!is.null(q$max_genes)) {
    keep_cell <- metrics$n_genes_by_counts < as.integer(q$max_genes)
    n_after_maxg <- sum(keep_cell)
    counts <- counts[keep_cell, , drop = FALSE]
    metrics <- metrics[keep_cell, , drop = FALSE]
  } else {
    n_after_maxg <- nrow(counts)
  }
  if (!is.null(q$max_pct_mt)) {
    keep_cell <- metrics$pct_counts_mt < as.numeric(q$max_pct_mt)
    n_after_mt <- sum(keep_cell)
    counts <- counts[keep_cell, , drop = FALSE]
    metrics <- metrics[keep_cell, , drop = FALSE]
  } else {
    n_after_mt <- nrow(counts)
  }

  n_final_pre_doublet <- nrow(counts)
  log_info(sprintf(paste0("过滤: 细胞 %d -> %d (min_genes) -> %d (max_genes) -> %d (max_pct_mt)；",
                          "基因 %d -> %d (min_cells，删掉 %d 个)"),
                   n0, n_after_min_genes, n_after_maxg, n_after_mt,
                   g0, g_after_min_cells, g0 - g_after_min_cells))

  # ---- 4. 双细胞 ----------------------------------------------------------
  db <- run_doublet(counts, metrics, cfg)
  n_after_db <- n_final_pre_doublet
  db_flag <- rep(FALSE, n_final_pre_doublet)
  if (identical(db$status, "ok")) {
    db_flag <- db$is_doublet
    if (db$n_predicted_doublets > 0) {
      counts <- counts[!db_flag, , drop = FALSE]
      metrics <- metrics[!db_flag, , drop = FALSE]
      n_after_db <- nrow(counts)
      log_info(sprintf("去除双细胞: %d -> %d", n_after_db + db$n_predicted_doublets, n_after_db))
    } else {
      log_info("scDblFinder 未标出双细胞")
    }
  } else {
    log_warn(sprintf("双细胞检测未执行: %s —— %s", db$status, db$reason %||% ""))
  }

  # ---- 5. ambient RNA -----------------------------------------------------
  amb <- assess_ambient_rna(metrics, cfg)
  if (!identical(amb$status, "ok")) {
    log_warn(sprintf("ambient RNA: %s —— %s", amb$status,
                     amb$reason %||% amb$note %||% ""))
  }

  if (nrow(counts) < 50) {
    stop(sprintf("过滤后只剩 %d 个细胞，无法继续分析", nrow(counts)), call. = FALSE)
  }

  # ---- 6. 落盘 ------------------------------------------------------------
  out <- file.path(data_dir, "qc_filtered.rds")
  saveRDS(list(counts = counts, metrics = metrics), out)
  log_info(sprintf("已写出 %s（%d 细胞 x %d 基因）", out, nrow(counts), ncol(counts)))

  # QC 汇总表：每个细胞一行，便于事后复查"到底滤掉了谁"
  qc_cells <- data.frame(
    cell = rownames(metrics), metrics, row.names = NULL,
    predicted_doublet = db_flag[seq_len(nrow(metrics))],
    doublet_score = if (identical(db$status, "ok")) db$doublet_score else NA_real_)
  utils::write.csv(qc_cells, file.path(res_dir, "qc_cells.csv"), row.names = FALSE)

  status <- list(
    dataset_id = cfg$dataset_id,
    organism = organism,
    n_cells_raw = n0,
    n_genes_raw = g0,
    filtering = list(
      min_genes = as.integer(q$min_genes),
      max_genes = q$max_genes,
      max_pct_mt = q$max_pct_mt,
      min_cells = as.integer(q$min_cells),
      n_after_min_genes = n_after_min_genes,
      n_after_max_genes = n_after_maxg,
      n_after_mt = n_after_mt,
      n_after_doublets = n_after_db,
      n_removed_total = n0 - nrow(counts),
      frac_removed = round((n0 - nrow(counts)) / max(n0, 1), 5),
      n_genes_before_filter = g0,
      n_genes_after_filter_cells = g_after_filter_cells,
      n_genes_after_min_cells = g_after_min_cells,
      n_genes_removed_by_min_cells = g0 - g_after_min_cells),
    qc_gene_sets = gene_sets,
    doublet_detection = db[!names(db) %in% c("is_doublet", "doublet_score")],
    ambient_rna = amb,
    n_cells_final = nrow(counts),
    n_genes_final = ncol(counts),
    status = "ok")
  write_json(file.path(res_dir, "qc_status.json"), status)
  status
}

# ---------------------------------------------------------------------------
# 入口（01_qc.py:413-427 同语义；E-56）
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L || identical(Sys.getenv("SCRNA_STEP_MAIN"), "01_qc")) {
  args <- parse_args()
  cfg <- load_config(args$config)
  t0 <- Sys.time()
  tryCatch({
    out <- run_01_qc(cfg)
    record_step(cfg, "qc", "ok", as.numeric(difftime(Sys.time(), t0, units = "secs")),
                result_status = result_status_of(out))
    out
  }, error = function(e) {
    record_step(cfg, "qc", "failed", as.numeric(difftime(Sys.time(), t0, units = "secs")),
                message = conditionMessage(e), result_status = "failed")
    stop(e)
  })
}
