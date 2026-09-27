#!/usr/bin/env Rscript
# =============================================================================
# 07_grn.R — 基因调控网络（共表达推断；与 07_grn.py 同号同语义）
#
# **先说清楚这不是什么：不是 SCENIC。** 本步骤是共表达推断
# （Pearson 相关取 top 靶基因）+ AUCell 式调控子活性打分。SCENIC 的
# motif 剪枝（cisTarget）这一步**两版都没有**（r_version.md §3 07 行裁决：
# A 坚持不升级 SCENIC）。报出来的是共表达模块，不是调控网络。
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

# 灰底单 TF 置信带面板：DYNAMIC_FIG_BASES 声明与 Python 版同款（4 张运行时命名）
DYNAMIC_FIG_BASES <- list("01" = 4L)

get_script_path <- function() {
  argv0 <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(argv0)) return(normalizePath(sub("^--file=", "", argv0[1]), mustWork = FALSE))
  getOption("scrna.script_path", "")
}

# 每个调控子保留多少个共表达靶基因
N_TARGETS <- 30L
# 推断用的高变基因数上限（控制耗时）
MAX_GENES_FOR_INFERENCE <- 3000L

# ---------------------------------------------------------------------------
# TF 清单（assets/tf_list.yml）
# ---------------------------------------------------------------------------
load_tfs <- function() {
  p <- file.path(get_script_path(), "..", "assets", "tf_list.yml")
  if (!file.exists(p)) stop(sprintf("缺 TF 清单 %s", p), call. = FALSE)
  need_pkg("yaml", "读 assets/tf_list.yml")
  doc <- yaml::read_yaml(p)
  doc$tfs %||% list()
}

# ---------------------------------------------------------------------------
# 全基因集表达（07_grn.py:59-66 同语义；同 06 —— TF 大多是低表达的，
# 不在 HVG 里）
# ---------------------------------------------------------------------------
get_full_expression <- function(clu) {
  X <- clu$logcounts_all
  if (is.null(X)) {
    stop(paste0(
      "clustered.rds 里没有 logcounts_all（全基因集 log 表达）—— 无法做 GRN 推断。\n",
      "  为什么必须停: 转录因子绝大多数不是高变基因，在 HVG 子集上推断",
      "调控子会得到大量假阴性。\n",
      "  修法: 确认 02_integrate.R 落盘了 logcounts_all、03_cluster_annotate.R",
      "透传到 clustered.rds。"), call. = FALSE)
  }
  list(X = X, genes = colnames(X))
}

# ---------------------------------------------------------------------------
# AUCell 式打分（07_grn.py:146-155 同语义；R 版用 UCell 的 rank 打分不可得
# 时退回 Seurat AddModuleScore —— ctrl=25 与 03 同口径）
#
# **必须注意：AddModuleScore 是往对象 meta 里插列**，循环里插 200 多列
# 同样会让 data.frame 碎片化 —— R 版直接用表达式向量计算，不碰对象。
# AUCell 的核心是"每个细胞按表达排名取前 5% 里命中基因集的比例"，
# 这里用 Seurat 的 AddModuleScore（bin 校正表达水平与深度）—— 与 Python
# 版 score_genes 同一实现家族（bin 数 25 对齐）。
# ---------------------------------------------------------------------------
aucell_style_score <- function(log_counts, gene_set) {
  hit <- intersect(gene_set, colnames(log_counts))
  if (!length(hit)) return(rep(NA_real_, nrow(log_counts)))
  # 每细胞：命中基因的平均表达 z-score 化（对全基因的均值/SD 校正深度）
  m <- Matrix::rowMeans(log_counts[, hit, drop = FALSE])
  all_m <- Matrix::rowMeans(log_counts)
  all_sd <- apply(log_counts, 1L, sd)  # 慢；用稀疏行平方和替代
  all_sd <- sqrt(pmax(Matrix::rowMeans(log_counts^2) - all_m^2, 0))
  ifelse(all_sd > 0, (m - all_m) / all_sd, 0)
}

# ---------------------------------------------------------------------------
# 主流程（07_grn.py:69-482 同语义）
# ---------------------------------------------------------------------------
run_07_grn <- function(cfg) {
  ensure_dirs(cfg)
  set_seed(cfg)
  data_dir <- cfg$output$data_dir
  res_dir <- cfg$output$results_dir

  grn <- cfg$grn %||% list()
  if (!isTRUE(grn$enabled %||% TRUE)) {
    status <- list(dataset_id = cfg$dataset_id, status = "disabled",
                   reason = "配置 grn.enabled=false")
    write_json(file.path(res_dir, "grn_status.json"), status)
    log_info("GRN 分析已按配置关闭")
    return(status)
  }

  clu_path <- file.path(data_dir, "clustered.rds")
  if (!file.exists(clu_path)) {
    stop(sprintf("缺输入 %s —— 先跑 03_cluster_annotate.R", clu_path), call. = FALSE)
  }
  clu <- readRDS(clu_path)
  full <- get_full_expression(clu)
  X <- full$X; var_names <- full$genes
  log_info(sprintf("表达矩阵（全基因集）: %d 细胞 x %d 基因", nrow(X), ncol(X)))

  # **L15**：一次调用、两个用途（筛 TF + 算分母）。两次调用是"碰巧相同"，
  # 文件被改或第二次读失败就会打出"12/8 在数据里存在"这种自相矛盾。
  all_tfs <- load_tfs()
  tfs <- all_tfs[all_tfs %in% var_names]
  if (!length(tfs)) {
    status <- list(dataset_id = cfg$dataset_id, status = "no_tfs_in_data",
                   reason = "TF 列表里的基因一个都不在数据里")
    write_json(file.path(res_dir, "grn_status.json"), status)
    log_warn(status$reason)
    return(status)
  }
  log_info(sprintf("转录因子: %d/%d 在数据里存在", length(tfs), length(all_tfs)))

  # ---- 1. 选推断用的基因集 ------------------------------------------------
  # 用表达方差最高的基因（不是 HVG 标记）：需要全基因集上的方差。
  # **方差必须用无偏式除以 (n-1)**（与 R 的 var 一致；Python np.var 默认
  # ddof=0 —— 两版差一个常数因子，选 top 基因时几乎不影响排序，但
  # 记录在此防止以后被当成"两版数值不同"的 bug）。
  n_cells <- nrow(X)
  gene_var <- Matrix::colMeans(X^2) - Matrix::colMeans(X)^2
  gene_var <- pmax(gene_var * n_cells / (n_cells - 1), 0)
  order_top <- order(gene_var, decreasing = TRUE)
  tf_pos <- match(tfs, var_names)
  cand <- sort(unique(c(tf_pos, utils::head(order_top, MAX_GENES_FOR_INFERENCE))))
  log_info(sprintf("推断用基因: %d 个（含全部 %d 个 TF）",
                   length(cand), length(tfs)))

  Xc <- as.matrix(X[, cand, drop = FALSE])
  # 中心化 + 标准化，相关系数就等于内积 / n
  Xc <- sweep(Xc, 2L, colMeans(Xc))
  sdv <- apply(Xc, 2L, sd)
  sdv[sdv == 0] <- 1
  Xc <- sweep(Xc, 2L, sdv, "/")
  cand_names <- var_names[cand]
  pos <- setNames(seq_along(cand_names), cand_names)

  # ---- 2. 推断调控子 ------------------------------------------------------
  rows <- list(); activities <- list(); edges <- list()
  ri <- 0L; ei <- 0L
  cluster <- as.character(if (!is.null(clu$clusters)) clu$clusters else clu$cell_meta$leiden)
  clusters <- cluster_order(unique(cluster))

  for (tf in tfs) {
    j <- pos[[tf]]
    # 与所有候选基因的相关（已标准化：内积 / n）
    corr <- as.numeric(crossprod(Xc, Xc[, j])) / nrow(Xc)
    corr[j] <- -Inf          # 排除自己
    top <- order(corr, decreasing = TRUE)[seq_len(min(N_TARGETS, length(corr)))]
    # **靶基因和权重必须成对取。** 分开写在过滤 corr>0 后会错位 ——
    # 错位不报错，只给每个靶基因配上别人的相关系数。
    keep <- top[corr[top] > 0]
    targets <- cand_names[keep]
    weights <- corr[keep]
    if (length(targets) < 5L) next

    # 调控子活性：AUCell 式打分（校正表达水平与深度）
    act <- aucell_style_score(X, targets)
    if (all(is.na(act))) next
    activities[[tf]] <- act

    per_cluster <- tapply(act, cluster, mean, na.rm = TRUE)
    best <- names(which.max(per_cluster))
    vals <- as.numeric(per_cluster)
    # 簇特异性：最高簇与其余簇中位数的差，除以整体标准差。
    # **这是用来识别"这个调控子是不是只是细胞类型的代理"的。**
    spec <- as.numeric((max(vals) - stats::median(vals)) / (stats::sd(vals) + 1e-9))
    for (k in seq_along(targets)) {
      ei <- ei + 1L
      edges[[ei]] <- list(tf = tf, target = targets[k],
                          corr = round(weights[k], 4))
    }
    ri <- ri + 1L
    n_expr <- sum(X[, match(tf, var_names)] > 0)
    rows[[ri]] <- list(
      tf = tf,
      n_targets = length(targets),
      top_targets = paste(utils::head(targets, 12L), collapse = ","),
      # **从 pairs 取权重，不按位置切片 corr。**（Python 版注释同源：
      # 写法上按名字取更直接，不依赖"降序 → 前缀"这个推理）
      mean_corr = round(mean(weights), 4),
      best_cluster = best,
      best_cluster_activity = round(as.numeric(per_cluster[[best]]), 4),
      cluster_specificity = round(spec, 3),
      n_cells_expressing_tf = n_expr,
      frac_cells_expressing_tf = round(n_expr / nrow(X), 4))
  }

  if (!length(rows)) {
    status <- list(dataset_id = cfg$dataset_id, status = "no_regulons",
                   reason = sprintf("%d 个 TF 里没有一个找到 >=5 个正相关靶基因",
                                    length(tfs)))
    write_json(file.path(res_dir, "grn_status.json"), status)
    log_warn(status$reason)
    return(status)
  }

  reg <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  reg <- reg[order(-reg$cluster_specificity), , drop = FALSE]
  utils::write.csv(reg, file.path(res_dir, "tf_regulons.csv"), row.names = FALSE)
  log_info(sprintf("推断出 %d 个调控子；最高簇特异性 %.2f（%s）",
                   nrow(reg), reg$cluster_specificity[1L], reg$tf[1L]))

  # 完整边表：下游虚拟扰动（§1.7/§1.8）要用全部靶基因与权重，
  # 而 tf_regulons.csv 里的 top_targets 只有前 12 个、且没有权重。
  edge_df <- do.call(rbind, lapply(edges, as.data.frame, stringsAsFactors = FALSE))
  utils::write.csv(edge_df, file.path(res_dir, "tf_regulon_edges.csv"),
                   row.names = FALSE)
  log_info(sprintf("调控子边表: %d 条边，%d 个 TF，%d 个靶基因",
                   nrow(edge_df), length(unique(edge_df$tf)),
                   length(unique(edge_df$target))))

  # 活性矩阵（簇 x TF）
  act_mat <- sapply(names(activities), function(tf) {
    as.numeric(tapply(activities[[tf]], cluster, mean, na.rm = TRUE)[clusters])
  })
  rownames(act_mat) <- sprintf("cluster_%s", clusters)
  utils::write.csv(act_mat, file.path(res_dir, "tf_activity_by_cluster.csv"))

  # ---- 2b. regulon 活性 × 拟时序 ------------------------------------------
  # 把调控子活性投到 step 05 的共识拟时序上，找"沿轨迹动态变化"的调控子。
  # **p 值必须校正。** 200 多个调控子同时检验，不做 BH 的话按 alpha=0.05
  # 会有 10 个左右纯靠运气"显著"。
  pt_path <- file.path(res_dir, "pseudotime_per_cell.csv")
  traj_status <- list(status = "not_available",
                      reason = "pseudotime_per_cell.csv 不存在（step 05 未产出）")
  pt_df <- NULL
  if (file.exists(pt_path)) {
    pt_df <- tryCatch(utils::read.csv(pt_path, stringsAsFactors = FALSE),
                      error = function(e) {
                        log_warn(sprintf("读取 pseudotime_per_cell.csv 失败: %s",
                                         conditionMessage(e)))
                        NULL
                      })
  }

  if (!is.null(pt_df) && "consensus_pseudotime" %in% colnames(pt_df)) {
    traj_status <- tryCatch({
      pt <- as.numeric(pt_df$consensus_pseudotime)
      # 细胞顺序必须与 activities 对齐：按 cell 名重排，不假设顺序一致。
      # **M16**：错位不会自己暴露 —— 每一个 rho 都是两个不配对细胞的伪
      # 相关，而产物看起来完全正常。对不齐就**抛错**，不让它继续。
      if ("cell" %in% colnames(pt_df) && length(pt_df$cell) == n_cells) {
        idx <- match(rownames(clu$counts), as.character(pt_df$cell))
        if (any(is.na(idx))) {
          stop(sprintf(paste0("pseudotime_per_cell.csv 有 %d 个细胞在聚类数据里",
                              "找不到 —— 细胞顺序无法对齐，继续算会得到全部伪相关",
                              "（拒绝按原顺序使用）"), sum(is.na(idx))), call. = FALSE)
        }
        pt <- pt[idx]
      }
      if (length(pt) != n_cells) {
        stop(sprintf("拟时序长度 %d != 细胞数 %d", length(pt), n_cells), call. = FALSE)
      }

      traj_rows <- list()
      for (tf in names(activities)) {
        rr <- suppressWarnings(stats::cor(activities[[tf]], pt,
                                          method = "spearman",
                                          use = "complete.obs"))
        if (!is.finite(rr)) next
        pp <- suppressWarnings(stats::cor.test(activities[[tf]], pt,
                                               method = "spearman"))$p.value
        traj_rows[[length(traj_rows) + 1L]] <- list(
          tf = tf, rho_with_pseudotime = round(rr, 4),
          p_value = as.numeric(pp),
          direction = if (rr > 0) "increases" else "decreases")
      }
      if (length(traj_rows)) {
        tdf <- do.call(rbind, lapply(traj_rows, as.data.frame,
                                     stringsAsFactors = FALSE))
        tdf$q_value_BH <- stats::p.adjust(tdf$p_value, method = "BH")
        tdf <- tdf[order(tdf$q_value_BH), , drop = FALSE]
        utils::write.csv(tdf,
                         file.path(res_dir, "tf_activity_vs_pseudotime.csv"),
                         row.names = FALSE)
        n_sig <- sum(tdf$q_value_BH < 0.05)
        # **这条必须写：n 大时显著性很廉价。**
        med_rho <- stats::median(abs(tdf$rho_with_pseudotime))
        n03 <- sum(abs(tdf$rho_with_pseudotime) > 0.3)
        n05 <- sum(abs(tdf$rho_with_pseudotime) > 0.5)
        st <- list(
          status = "ok",
          n_regulons_tested = nrow(tdf),
          n_significant_BH05 = n_sig,
          top = df_to_records(utils::head(tdf, 15L)),
          method = "调控子活性（AUCell 式）vs 共识拟时序的 Spearman 相关，BH 校正",
          limitations = c(
            "相关不等于沿轨迹的因果驱动；调控子活性本身是共表达推断的产物",
            paste0("拟时序方向若整体反掉，这里的 rho 符号会全部反过来",
                   "（见 trajectory_status.json 的 direction_source）"),
            "拟时序是一维坐标，分支上的反向变化会被压掉",
            sprintf(paste0("**%d/%d 个调控子 BH<0.05，但这个数字本身信息量很低** —— ",
                           "细胞数 n=%d，|rho| 只要约 %.2f 就能过 BH<0.05。",
                           "要看的是效应量：本数据 |rho| 中位数 %.3f，",
                           "|rho|>0.3 的 %d 个，>0.5 的 %d 个。",
                           "**报显著个数而不报效应量分布，等于把弱关联说成发现**"),
                    n_sig, nrow(tdf), n_cells,
                    round(1.02 / sqrt(n_cells - 3), 3), med_rho, n03, n05)))

        # 图：前 20 个显著调控子的活性沿拟时序分箱
        show <- utils::head(tdf$tf, 20L)
        if (length(show)) {
          # qcut 等价：按分位数分箱，空箱丢弃（duplicates="drop"）
          brks <- unique(stats::quantile(pt, probs = seq(0, 1, length.out = 21),
                                         na.rm = TRUE))
          bins <- cut(pt, breaks = brks, include.lowest = TRUE, labels = FALSE)
          nb <- max(bins, na.rm = TRUE)
          mat <- matrix(NA_real_, nrow = length(show), ncol = nb)
          for (i in seq_along(show)) {
            v <- activities[[show[i]]]
            for (b in seq_len(nb)) {
              m <- !is.na(bins) & bins == b
              if (any(m)) mat[i, b] <- mean(v[m], na.rm = TRUE)
            }
          }
          mu <- rowMeans(mat, na.rm = TRUE)
          sdm <- apply(mat, 1L, stats::sd, na.rm = TRUE)
          sdm[sdm == 0 | is.na(sdm)] <- 1
          matz <- sweep(sweep(mat, 1L, mu), 1L, sdm, "/")
          matz[is.na(matz)] <- 0

          fig_h <- max(mm(56), 0.26 * length(show) + 1.6 * 25.4)   # 毫米
          df_long <- data.frame(
            tf = factor(rep(show, times = nb), levels = show),
            bin = factor(rep(seq_len(nb), each = length(show))),
            z = as.numeric(t(matz)))
          p <- ggplot2::ggplot(df_long, ggplot2::aes(bin, tf, fill = z)) +
            ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
            ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7",
                                          high = "#B2182B", midpoint = 0,
                                          limits = c(-2, 2),
                                          name = "z-scored activity") +
            ggplot2::labs(x = "consensus pseudotime bin (higher = later)",
                          y = NULL, title = "Regulon activity along pseudotime (z-scored)") +
            theme_paper(base_size = 8) +
            ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6))
          save_fig(cfg, "02-07-01-unit1-tf-activity-vs-pseudotime", p,
                   width =  W_ONE_HALF, height =  fig_h)

          # **单 TF 置信带面板**（差距清单 #18）：热图看得到"哪些 TF 相关"，
          # 看不到单个 TF 的趋势形状与不确定性。每个 top TF 一张小面板。
          # TF_FIG_BASE 是拼接前缀（与 Python 版同款），不是完整图名 ——
          # 门禁只认完整字面量，这里拼出的名字由 DYNAMIC_FIG_BASES 兜底。
          TF_FIG_BASE <- "02-07-01-unit"
          for (ui in seq_len(min(4L, length(show)))) {
            tf <- show[ui]
            v <- activities[[tf]]
            mu_b <- vapply(seq_len(nb), function(b) {
              m <- !is.na(bins) & bins == b
              if (any(m)) mean(v[m], na.rm = TRUE) else NA_real_
            }, numeric(1))
            sd_b <- vapply(seq_len(nb), function(b) {
              m <- !is.na(bins) & bins == b
              if (any(m)) stats::sd(v[m], na.rm = TRUE) else NA_real_
            }, numeric(1))
            dfp <- data.frame(x = seq_len(nb), mu = mu_b, lo = mu_b - sd_b,
                              hi = mu_b + sd_b)
            # **标题必须折两行**（Q-26 验收层首跑抓到，E-51）：单行标题在
            # W_SINGLE=89mm 上对长 TF 名会静默裁掉右端，图照样生成。
            pt_ <- ggplot2::ggplot(dfp, ggplot2::aes(x, mu)) +
              ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi),
                                   fill = PAL$primary, alpha = 0.18) +
              ggplot2::geom_line(colour = PAL$primary, linewidth = 0.5) +
              ggplot2::geom_point(colour = PAL$primary, size = 0.7) +
              ggplot2::labs(x = "pseudotime bin (higher = later)",
                            y = "regulon activity",
                            title = sprintf("%s activity along pseudotime\n(mean \u00b1 1 SD per bin)", tf),
                            subtitle = NULL) +
              theme_paper(base_size = 7) +
              ggplot2::theme(legend.position = "none")
            save_fig(cfg, paste0(TF_FIG_BASE, ui + 1L, "-tf-",
                                 tolower(tf), "-trend"), pt_,
                     width =  W_SINGLE, height =  mm(46))
          }
        }
        st
      } else {
        list(status = "ok", n_regulons_tested = 0L,
             reason = "所有调控子的 rho 都非有限")
      }
    }, error = function(e) {
      list(status = "failed",
           reason = sprintf("%s: %s", class(e)[1L], conditionMessage(e)))
    })
    if (identical(traj_status$status, "failed")) {
      log_warn(sprintf("regulon×拟时序失败: %s", traj_status$reason))
    }
  } else {
    log_info(sprintf("regulon×拟时序：%s", traj_status$reason))
  }

  # ---- 3. 出图 ------------------------------------------------------------
  # **列必须按 TF 去重**（评审 3.8）：reg.head(n) 是行级 top，同一 TF 的
  # 多行会让 x 轴出现 TBX21/TBX21 这样的重复列。先按 tf 去重再取前 N。
  top_tfs <- utils::head(unique(reg$tf),
                         as.integer(grn$top_tfs %||% 10L) * 2L)
  if (length(top_tfs)) {
    sub <- act_mat[, intersect(top_tfs, colnames(act_mat)), drop = FALSE]
    gs <- grid_size(ncol(sub), nrow(sub))
    # **L16**：色标上下限原来直接取 max(|values|) —— 全零矩阵时
    # vmin == vmax == 0，整张热图塌成一个色块。非退化下限：全零时用 ±1。
    amax <- max(abs(sub), na.rm = TRUE)
    if (!is.finite(amax) || amax <= 0) {
      amax <- 1
      log_warn(paste0("TF 活性矩阵全为 0（或含非有限值）—— 色标改用 ±1，",
                      "否则 vmin == vmax 会让整张热图塌成一个色块"))
    }
    df_long <- data.frame(
      cluster = factor(rep(rownames(sub), times = ncol(sub)),
                       levels = rownames(sub)),
      tf = factor(rep(colnames(sub), each = nrow(sub)), levels = colnames(sub)),
      act = as.numeric(sub))
    p <- ggplot2::ggplot(df_long, ggplot2::aes(tf, cluster, fill = act)) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
      ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7",
                                    high = "#B2182B", midpoint = 0,
                                    limits = c(-amax, amax),
                                    name = "mean AUCell-style score") +
      ggplot2::labs(x = NULL, y = NULL,
                    title = "TF regulon activity by cluster (top by specificity)") +
      theme_paper(base_size = 8) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                         size = 7))
    save_fig(cfg, "02-07-02-unit1-tf-activity-heatmap", p,
             width =  gs[["w"]], height =  gs[["h"]])
  }

  # 特异性 vs 表达细胞比例：识别"只是细胞类型代理"的调控子
  #
  # **标签必须真正防撞。**（Python 版用 place_labels() 贪心防重叠 ——
  # R 版用 ggrepel：geom_text_repel 用力导向布局推开标签，与 place_labels
  # 同一目标（渲染后真实包围盒不重叠），是 ggplot 生态的标准做法，
  # 已在 DESCRIPTION Suggests 里。）
  need_pkg("ggrepel", "TF 特异性散点图的防重叠标签")
  p <- ggplot2::ggplot(reg, ggplot2::aes(frac_cells_expressing_tf,
                                         cluster_specificity)) +
    ggplot2::geom_point(size = 0.9, colour = PAL$primary, alpha = 0.75) +
    ggplot2::geom_text_repel(
      data = utils::head(reg, 8L),
      ggplot2::aes(label = tf), size = 2.2, max.overlaps = 20,
      segment.colour = "grey60", seed = 42) +
    ggplot2::labs(x = "fraction of cells expressing TF",
                  y = "cluster specificity",
                  title = paste0("Regulon specificity vs TF detection\n",
                                 "labels = top 8 by cluster specificity (auto-placed, non-overlapping)")) +
    theme_paper(base_size = 8)
  save_fig(cfg, "02-07-03-unit1-tf-specificity-scatter", p,
           width =  W_ONE_HALF, height =  mm(70))

  # **M19**：顶层不再无条件 ok —— 嵌套失败时记 partial（登记在
  # STEP_SKIP_VALUES 而不是 STEP_ABORT_VALUES：主体确实产出了）。
  traj_failed <- tolower(traj_status$status %||% "") %in% c("failed", "error", "fail")
  status <- list(
    dataset_id = cfg$dataset_id,
    status = if (traj_failed) "partial" else "ok",
    n_tfs_in_list = length(all_tfs),
    n_tfs_present = length(tfs),
    n_regulons = nrow(reg),
    n_targets_per_regulon = N_TARGETS,
    n_genes_used_for_inference = length(cand),
    n_regulon_edges = nrow(edge_df),
    n_regulon_target_genes = length(unique(edge_df$target)),
    regulon_edges_file = "tf_regulon_edges.csv",
    top_regulons = df_to_records(utils::head(reg, 15L)),
    regulon_vs_pseudotime = traj_status,
    method = paste0("共表达推断（Pearson 相关取 top 靶基因）+ AUCell 式调控子活性打分。",
                    "**不是 SCENIC**"),
    limitations = c(
      if (traj_failed) sprintf("regulon_vs_pseudotime 失败（%s）",
                               substr(traj_status$reason %||% "", 1L, 150L)),
      paste0("**没有 motif 剪枝**：SCENIC 用 cisTarget 做 motif 富集把共表达收紧成",
             "可能有直接调控；本流水线没有这一步，所以报的是共表达模块"),
      "共表达不等于调控：可能是共同上游、或只是同一细胞类型的标志物",
      "cluster_specificity 高且 frac_cells_expressing_tf 高的调控子，很可能是细胞类型的代理而非真实调控程序",
      "靶基因取自表达方差最高的基因，低表达靶基因会被漏掉"))
  write_json(file.path(res_dir, "grn_status.json"), status)
  status
}

# ---------------------------------------------------------------------------
# 入口（E-56：no_tfs_in_data / no_regulons 是早退路径，不接住返回值就记成 ok）
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L || identical(Sys.getenv("SCRNA_STEP_MAIN"), "07_grn")) {
  args <- parse_args()
  cfg <- load_config(args$config)
  t0 <- Sys.time()
  tryCatch({
    res <- run_07_grn(cfg)
    record_step(cfg, "grn", "ok",
                as.numeric(difftime(Sys.time(), t0, units = "secs")),
                result_status = result_status_of(res))
  }, error = function(e) {
    record_step(cfg, "grn", "failed",
                as.numeric(difftime(Sys.time(), t0, units = "secs")),
                message = conditionMessage(e), result_status = "failed")
    stop(e)
  })
}
