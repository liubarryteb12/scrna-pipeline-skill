#!/usr/bin/env Rscript
# =============================================================================
# 06_communication.R — 细胞间通讯（配体-受体；与 06_communication.py 同号同语义）
#
# **两条路都跑，都如实标注。**
#   1. **liana rank_aggregate**（Bioconductor/GitHub，装得上就跑）→ liana_results.csv
#   2. **自建共表达打分**（始终跑）→ cell_communication.csv
#      内置配体-受体对 + 置换检验 + BH 校正。
#
# 两个都产出并**量化一致性**（Spearman + top-N 重叠）—— 不一致本身就是
# 发现，不是谁错了。R 版的 liana 走 liana_wrap()（GitHub 装包，
# references/r_version.md §3 06 行裁决：等级 A/B）。
#
# 必须说清楚的局限（两版相同）：共表达不等于通讯；稳态丰度不等于蛋白
# 水平；自建打分是启发式。
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

# 脚本定位：与 01_qc.R 同款 —— 独立 Rscript 用 --file=；被 main_analysis.R
# source 时由入口先设置 options(scrna.script_path=...)。
get_script_path <- function() {
  argv0 <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(argv0)) return(normalizePath(sub("^--file=", "", argv0[1]), mustWork = FALSE))
  getOption("scrna.script_path", "")
}

N_PERMUTATIONS <- 200L
LIANA_N_PERMS <- 100L

# ---------------------------------------------------------------------------
# 内置配体-受体对（assets/ligand_receptor.yml）
# ---------------------------------------------------------------------------
load_lr_pairs <- function() {
  p <- file.path(dirname(get_script_path()), "..", "assets", "ligand_receptor.yml")
  if (!file.exists(p)) stop(sprintf("缺内置配体-受体库 %s", p), call. = FALSE)
  need_pkg("yaml", "读 assets/ligand_receptor.yml")
  doc <- yaml::read_yaml(p)
  doc$pairs %||% list()
}

# ---------------------------------------------------------------------------
# 全基因集表达（06_communication.py:57-79 同语义）
#
# **本步骤最容易踩的坑，实测踩过。** 配体/受体基因大多是低表达的细胞因子/
# 趋化因子，几乎不进 HVG —— 实测 HVG 子集 38 对只有 3 对可用，全基因集
# 27 对，差 9 倍。缺全基因集 = "没找到显著通讯"变成假阴性，而假阴性
# 看起来和"真的没有通讯"一模一样。缺失时**明确报错**，不悄悄退回 HVG。
# ---------------------------------------------------------------------------
get_full_expression <- function(clu) {
  X <- clu$logcounts_all
  if (is.null(X)) {
    stop(paste0(
      "clustered.rds 里没有 logcounts_all（全基因集 log 表达）—— 无法做通讯分析。\n",
      "  为什么必须停: 配体/受体基因绝大多数不是高变基因，在 HVG 子集上做",
      "通讯分析会得到大量假阴性（实测 3/38 vs 27/38）。\n",
      "  修法: 确认 02_integrate.R 落盘了 logcounts_all、03_cluster_annotate.R",
      "透传到 clustered.rds。"), call. = FALSE)
  }
  # **logcounts_all 布局 = 基因 x 细胞**（02_integrate.R:284 落盘 t(data 层)，
  # 03 透传日志 ncol=基因数同口径）。基因名在 rownames，不在 colnames ——
  # 之前写 colnames(X) 拿到的是细胞名，配体/受体基因全部 miss，
  # 主链会静默退化成 0 分（正是本函数要防的假阴性形态）。
    # **logcounts_all 实际布局 = 细胞 x 基因**：02_integrate.R:284 对 Seurat data 层
  #（基因x细胞）做 t()，得到的是 Python adata.raw 语义（细胞x基因）。
  # run17 实锤：直接 rownames(X) 拿到的是**细胞名**（"2574 基因 x 13714 细胞"），
  # 配体/受体/TF 全部 miss → no_pairs_in_data / no_tfs_in_data。
  # 这里统一转成 **基因 x 细胞** 再交出去（下游消费方 b4ee9a6 起全部按基因x细胞写）。
  X <- Matrix::t(X)
  list(X = X, genes = rownames(X))
}

# 给定基因集与细胞下标，返回平均表达（log 后）—— 语义对齐
# 06_communication.py:82-89；R 版主链用 score_combo（子矩阵索引），本函数
# 保留供外部核查（与 Python 版同名同语义的对照物）。
mean_expression <- function(genes, idx_cells, X) {
  hit <- genes[genes %in% rownames(X)]
  if (!length(hit)) return(0)
  mean(X[hit, idx_cells, drop = FALSE])
}

# ---------------------------------------------------------------------------
# liana（R 版主工具；06_communication.py:108-167 同语义）
#
# **任何异常都吞掉并写进 info** —— liana 跑不动不该让整步失败，
# 但必须让人看见它跑不动。
# ---------------------------------------------------------------------------
try_liana <- function(clu, group_labels, cfg) {
  # 注意：LIANA_N_PERMS 只是**登记意图**（Python 侧 n_perms=100 的对照），
  # liana_wrap 顶层没有 n_perms 参数，实际走 liana 默认值 —— 不谎报为已生效。
  info <- list(attempted = TRUE, status = NULL, reason = "",
               n_perms_intended = LIANA_N_PERMS,
               n_perms_actual = "liana 默认")
  has <- requireNamespace("liana", quietly = TRUE)
  if (!has) {
    info$status <- "package_missing"
    info$reason <- "liana 未安装（GitHub 装包，references/r_version.md §3 06 行）"
    log_warn(sprintf("liana 不可用：%s", info$reason))
    return(list(res = NULL, info = info))
  }
  info$version <- tryCatch(as.character(utils::packageVersion("liana")),
                           error = function(e) "unknown")
  out <- tryCatch({
    # **liana 硬性要求 SCE 同时有 counts + logcounts 两个 assay**（run29 实锤
    # `liana expects counts and logcounts to be present in the SCE object`）。
    # SCE assays 必须同尺寸（run20 实锤）→ 只能用 HVG 对（counts/logcounts
    # 同为 HVG 子集）；logcounts_all 全基因集配不出同尺寸 counts。
    # 代价：LR 基因不在 HVG 时 liana 会低估 —— 记进 info$gene_scope，
    # 不静默（Python 版 liana 走 raw 全基因集，这是两版的真实差距）。
    # 落盘方向：counts/logcounts 都是 细胞x基因（run17 实锤），t() 成 SCE
    # 约定的 基因x细胞。
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = Matrix::t(clu$counts),
                    logcounts = Matrix::t(clu$logcounts)))
    info$gene_scope <- "HVG 子集（liana 要求 counts+logcounts 同尺寸，
                        全基因集 logcounts_all 配不出同尺寸原始计数）"
    # colLabels / reducedDims 是 **SingleCellExperiment 包**的导出（run25 实锤
    # 'colLabels<-' is not an exported object from 'namespace:SummarizedExperiment'
    # —— SummarizedExperiment 只提供 colData 等， SCE 的簇标签/降维槽要挂 SCE 命名空间）。
    # 簇标签按 SCE **列名（细胞名）显式匹配重排** —— clusters 向量的名序与
    # t(logcounts_all) 的列序来自不同落盘对象，集合相同不代表顺序相同；
    # colLabels<- 按位置赋值，顺序错 = 标签错配且不报错（比崩更危险）。
    cells_sce <- colnames(sce)
    miss <- sum(!cells_sce %in% names(group_labels))
    if (miss > 0L) {
      stop(sprintf("liana SCE 有 %d 个细胞在簇标签向量里找不到 —— 无法对齐", miss),
           call. = FALSE)
    }
    labels_aligned <- factor(group_labels[cells_sce])
    names(labels_aligned) <- NULL
    SingleCellExperiment::colLabels(sce) <- labels_aligned
    # reducedDim 的行（细胞）同样按列名对齐，缺行名的 PCA 直接不给（liana 只在
    # 需要表达方式时才用 reducedDim，缺了退默认 logcounts，不硬塞错位矩阵）。
    if (!is.null(rownames(clu$pca)) &&
        setequal(rownames(clu$pca), cells_sce)) {
      SingleCellExperiment::reducedDims(sce) <- list(PCA = clu$pca[cells_sce, , drop = FALSE])
    }
    # liana_wrap 顶层没有 n_perms/seed 参数（那是 natmi 等各 method 的内部
    # 参数，liana_wrap 的 ... 不做透传登记，多余参数直接 unused argument ——
    # 与 Read10X var.names.repair 同类的 scanpy 语义误置）。
    res <- liana::liana_wrap(sce, idents_col = "label",
                             resource = "Consensus")
    agg <- liana::liana_aggregate(res)
    info$status <- "ok"
    info$n_rows <- nrow(agg)
    info$columns <- sort(colnames(agg))
    log_info(sprintf("liana_wrap 完成：%d 行，版本 %s", nrow(agg), info$version))
    list(res = agg, info = info)
  }, error = function(e) {
    info$status <- "failed"
    info$reason <- sprintf("%s: %s", class(e)[1L], conditionMessage(e))
    log_warn(sprintf("liana 跑失败，退回自建共表达打分：%s", info$reason))
    list(res = NULL, info = info)
  })
  out
}

# ---------------------------------------------------------------------------
# 自建打分 vs liana 的一致性（06_communication.py:170-236 同语义）
#
# 两种方法给出不同排序是常态：liana 的 rank aggregate 综合了 7 种方法的
# 秩，自建打分只是表达量乘积。报 Spearman 与 top-N 重叠，让"差多少"
# 变成数字而不是印象。
# ---------------------------------------------------------------------------
compare_with_liana <- function(liana_res, own) {
  out <- list(compared = FALSE)
  if (is.null(liana_res) || is.null(own) || nrow(own) == 0L) return(out)
  tryCatch({
    score_col <- intersect(c("aggregate_rank", "consensus_score",
                             "magnitude_rank", "expr_prod", "lr_means"),
                           colnames(liana_res))[1L]
    if (is.na(score_col)) {
      out$reason <- sprintf("liana 结果里没有可识别的分数列（有 %s）",
                            paste(utils::head(colnames(liana_res), 8), collapse = ", "))
      return(out)
    }
    lr <- liana_res
    # liana 返回列名经 data.frame 自动 make.names 变成 ligand.complex /
    # receptor.complex（run34 实锤：代码写 snake_case → paste(NULL,NULL)=
    # character(0) → "replacement has 0 rows, data has 212"）。两种都认。
    pick <- function(df, candidates) {
      for (cn in candidates) if (cn %in% colnames(df)) return(cn)
      NA_character_
    }
    lig_col <- pick(lr, c("ligand_complex", "ligand.complex"))
    rec_col <- pick(lr, c("receptor_complex", "receptor.complex"))
    src_col <- pick(lr, c("source", "target"))
    tgt_col <- pick(lr, c("target", "receiver"))
    if (is.na(lig_col) || is.na(rec_col) || is.na(src_col) || is.na(tgt_col)) {
      out$reason <- sprintf("liana 结果缺少必要列（有 %s）",
                            paste(utils::head(colnames(lr), 10), collapse = ", "))
      return(out)
    }
    lr$.pair <- paste(lr[[lig_col]], lr[[rec_col]], sep = "^")
    lr$.combo <- paste(lr$.pair, lr[[src_col]], lr[[tgt_col]], sep = "|")
    lr_score <- tapply(lr[[score_col]], lr$.combo, mean)

    ow <- own
    ow$.pair <- paste(ow$ligand, ow$receptor, sep = "^")
    ow$.combo <- paste(ow$.pair, ow$sender, ow$receiver, sep = "|")
    ow_score <- tapply(ow$score, ow$.combo, mean)

    # 自建打分越大越强；liana 的 *_rank 越小越强，要翻向
    ascending <- grepl("_rank$", score_col)
    if (ascending) lr_score <- -lr_score

    common <- intersect(names(ow_score), names(lr_score))
    out$compared <- TRUE
    out$score_column_used <- score_col
    out$rank_direction <- if (ascending) "越小越强（已翻向）" else "越大越强"
    out$n_common_combinations <- length(common)
    out$n_own_only <- length(setdiff(names(ow_score), names(lr_score)))
    out$n_liana_only <- length(setdiff(names(lr_score), names(ow_score)))
    if (length(common) >= 5L) {
      # run35 实锤：自建 score 全 NA 时 cor(use="complete.obs") 抛
      # "no complete element pairs" —— 不是对比逻辑错，是上游打分全 NA 的
      # 连带症状。防御：完整对 <2 时如实记 reason，不抛错淹没根因。
      n_complete <- sum(is.finite(ow_score[common]) & is.finite(lr_score[common]))
      if (n_complete < 2L) {
        out$comparison_reason <- sprintf(
          "完整对 %d（自建侧非有限 %d/%d）—— 上游打分全 NA 的连带，先修打分根因",
          n_complete, sum(!is.finite(ow_score[common])), length(common))
      } else {
        rho <- stats::cor(ow_score[common], lr_score[common],
                          method = "spearman",
                          use = "complete.obs")
        out$spearman_rho <- round(as.numeric(rho), 4)
      }
    }
    for (n in c(10L, 25L, 50L)) {
      if (length(common) >= n) {
        a <- utils::head(names(sort(ow_score[common], decreasing = TRUE)), n)
        b <- utils::head(names(sort(lr_score[common], decreasing = TRUE)), n)
        out[[sprintf("top%d_overlap", n)]] <- length(intersect(a, b))
        out[[sprintf("top%d_overlap_frac", n)]] <- round(length(intersect(a, b)) / n, 3)
      }
    }
    log_info(sprintf("自建 vs liana：共同组合 %d，Spearman rho=%s，top25 重叠 %s/25",
                     out$n_common_combinations,
                     out$spearman_rho %||% "NA",
                     out$top25_overlap %||% "NA"))
  }, error = function(e) {
    # **不用 `<<-`**（b431 教训：error handler 里的超赋值会在全局静默造
    # 变量 —— M8 家族全局污染形态）；out 是本函数闭包内的 list，普通
    # `<-` 在 handler 里同样改的是这个绑定（tryCatch 的 handler 是在
    # 调用帧求值的），直接返回 out 由调用方取值。
    out$reason <- sprintf("对比失败：%s: %s", class(e)[1L], conditionMessage(e))
    log_warn(out$reason)
  })
  out
}

# ---------------------------------------------------------------------------
# 主流程（06_communication.py:239-542 同语义）
# ---------------------------------------------------------------------------
run_06_communication <- function(cfg) {
  ensure_dirs(cfg)
  set_seed(cfg)
  data_dir <- cfg$output$data_dir
  res_dir <- cfg$output$results_dir

  comm <- cfg$communication %||% list()
  if (!isTRUE(comm$enabled %||% TRUE)) {
    status <- list(dataset_id = cfg$dataset_id, status = "disabled",
                   reason = "配置 communication.enabled=false")
    write_json(file.path(res_dir, "communication_status.json"), status)
    log_info("细胞通讯分析已按配置关闭")
    return(status)
  }

  clu_path <- file.path(data_dir, "clustered.rds")
  if (!file.exists(clu_path)) {
    stop(sprintf("缺输入 %s —— 先跑 03_cluster_annotate.R", clu_path), call. = FALSE)
  }
  clu <- readRDS(clu_path)
  cell_meta <- clu$cell_meta
  group_key <- if ("celltype" %in% colnames(cell_meta)) "celltype" else "leiden"
  clusters_vec <- if (group_key == "celltype") cell_meta$celltype else clu$clusters
  clusters_vec <- as.character(clusters_vec)
  groups <- sort(unique(clusters_vec))
  if (length(comm$clusters)) {
    groups <- intersect(groups, as.character(comm$clusters))
  }
  if (length(groups) < 2L) {
    status <- list(dataset_id = cfg$dataset_id, status = "not_applicable",
                   reason = sprintf("只有 %d 个簇/细胞类型，通讯分析需要至少 2 个",
                                    length(groups)))
    write_json(file.path(res_dir, "communication_status.json"), status)
    log_warn(status$reason)
    return(status)
  }

  pairs <- load_lr_pairs()
  full <- get_full_expression(clu)
  X <- full$X; var_names <- full$genes
  lookup <- var_names
  log_info(sprintf("表达矩阵（全基因集）: %d 基因 x %d 细胞", nrow(X), ncol(X)))

  # **L7（R-03 裁决）**：score 的量纲由上游标准化决定，target_sum 不落盘
  # 读者无法复现量级。从 integration_status.json 读**实际用的**值。
  integ <- read_json(file.path(res_dir, "integration_status.json")) %||% list()
  target_sum <- integ$target_sum
  norm_note <- if (!is.null(target_sum)) {
    sprintf("表达量来自 logcounts_all（log1p(counts / 每细胞总计数 * target_sum)，target_sum=%s）",
            as.character(target_sum))
  } else {
    "**未能从 integration_status.json 读到 target_sum** —— 本表的 score 量级无法复现"
  }

  present <- Filter(function(pr) {
    any(pr$ligand %in% lookup) && any(pr$receptor %in% lookup)
  }, pairs)
  present <- lapply(present, function(pr) {
    lig <- pr$ligand[pr$ligand %in% lookup]
    rec <- pr$receptor[pr$receptor %in% lookup]
    if (length(lig) && length(rec)) {
      list(name = pr$name, ligand = lig, receptor = rec,
           ligand_missing = setdiff(pr$ligand, lig),
           receptor_missing = setdiff(pr$receptor, rec))
    } else NULL
  })
  present <- Filter(Negate(is.null), present)
  log_info(sprintf("配体-受体对: %d/%d 的基因在数据里存在",
                   length(present), length(pairs)))

  if (!length(present)) {
    status <- list(dataset_id = cfg$dataset_id, status = "no_pairs_in_data",
                   reason = sprintf("数据库里 %d 对，没有一对的配体与受体基因同时出现在数据里",
                                    length(pairs)))
    write_json(file.path(res_dir, "communication_status.json"), status)
    log_warn(status$reason)
    return(status)
  }

  masks <- lapply(groups, function(g) which(clusters_vec == g))
  names(masks) <- groups
  n_cells <- ncol(X)
  seed <- as.integer(cfg$analysis$seed)

  # **只抽涉及的配体/受体基因行**（06_communication.py:310-317 同优化；
  # X = 基因x细胞，基因按行取）。masks 的细胞下标做**列**下标。
  needed <- sort(unique(unlist(lapply(present, function(pr) c(pr$ligand, pr$receptor)))))
  small <- X[needed, , drop = FALSE]
  log_info(sprintf("置换用子矩阵: %d 个配体/受体基因 x %d 细胞",
                   nrow(small), ncol(small)))

  score_combo <- function(lig_genes, rec_genes, s_idx, d_idx) {
    li <- intersect(lig_genes, rownames(small))
    ri <- intersect(rec_genes, rownames(small))
    if (!length(li) || !length(ri)) return(0)
    a <- mean(small[li, s_idx, drop = FALSE])
    b <- mean(small[ri, d_idx, drop = FALSE])
    val <- as.numeric(a * b)
    # Python 版 NaN<=0 判 False 静默续跑；R if(NaN) 直接 fatal（run19 实锤
    # missing value where TRUE/FALSE needed）。这里把非有限显式传出去，
    # 调用端如实记 NA + note，不让 NaN 静默变 0（谎报）也不让步骤崩。
    if (length(val) != 1L || !is.finite(val)) {
      # run35 诊断实锤：CD274 行 @x 无 NA、全细胞均值有限（0.0041），但组合
      # 仍 NA —— a/b 谁是 NA 必须当场拆开。列子集 mean 出 NA 的两个候选：
      # 列子集里混进显式 NA（@x 路径）或 mean 分派异常。首例打印全部输入。
      if (n_nonfinite == 0L) {
        la <- small[li, s_idx, drop = FALSE]
        lb <- small[ri, d_idx, drop = FALSE]
        log_warn(sprintf(
          paste0("首个非有限打分拆解: pair=%s a_NA=%s b_NA=%s ",
                 "a_anyNA_x=%s a_len_x=%d a_sum=%s a_len=%d ",
                 "b_anyNA_x=%s b_len_x=%d b_sum=%s b_len=%d"),
          pr_name_cur, is.na(a), is.na(b),
          anyNA(la@x), length(la@x), sprintf("%.6g", sum(la@x)), length(la),
          anyNA(lb@x), length(lb@x), sprintf("%.6g", sum(lb@x)), length(lb)))
      }
      return(NA_real_)
    }
    val
  }

  rows <- list(); ri_ <- 0L
  pr_name_cur <- ""   # 诊断用：当前 pair 名（score_combo 拆解日志引用）
  n_evaluated <- 0L    # 真正做过置换检验的组合数
  n_zero_dropped <- 0L # 打分恒为 0、未做置换的组合数
  n_nonfinite <- 0L    # 打分非有限（NA/NaN）的组合数 —— 如实入账，不静默归 0
  p_res <- round(1.0 / (N_PERMUTATIONS + 1L), 8)
  for (pr in present) {
    for (src in groups) {
      for (dst in groups) {
        if (identical(src, dst)) next
        pr_name_cur <- pr$name
        obs <- score_combo(pr$ligand, pr$receptor, masks[[src]], masks[[dst]])
        if (is.na(obs)) {
          # 非有限打分：第一次打诊断日志（定位 NaN 源头用），组合如实记 NA。
          # run34 实锤 1134/1134 全 NA 但 small_NA=0 —— is.na(稀疏矩阵) 查的是
          # **显式存储**的 NA；Matrix 的 mean(sparseMatrix) = mean(as(x,"sparseVector"))
          # = sum(x@x)/n，@x 里**任何一个显式 NA 毒死整个 mean**。所以诊断必须
          # 直接查 @x 槽（anyNA(x@x)），并给出 mean 的稀疏实现展开值。
          if (n_nonfinite == 0L) {
            sub <- small[intersect(pr$ligand, rownames(small)), , drop = FALSE]
            log_warn(sprintf(
              paste0("首个非有限打分诊断: pair=%s %s->%s ligand=[%s] receptor=[%s] ",
                     "n_src=%d n_dst=%d small_NA=%d anyNA_at_x=%s n_at_x=%d ",
                     "mean_direct=%s sum_at_x=%s len_at_x=%d"),
              pr$name, src, dst, paste(pr$ligand, collapse = ","),
              paste(pr$receptor, collapse = ","),
              length(masks[[src]]), length(masks[[dst]]),
              sum(is.na(sub)), anyNA(sub@x), length(sub@x),
              sprintf("%.6g", sum(sub@x) / length(sub)),
              sprintf("%.6g", sum(sub@x)), length(sub@x)))
          }
          n_nonfinite <- n_nonfinite + 1L
          ri_ <- ri_ + 1L
          rows[[ri_]] <- list(pair = pr$name, sender = src, receiver = dst,
                              ligand = paste(pr$ligand, collapse = ","),
                              receptor = paste(pr$receptor, collapse = ","),
                              score = NA, null_mean = NA,
                              p_value = NA, n_permutations = 0L, tested = FALSE,
                              note = "obs 非有限（NA/NaN）—— 未检验，不得当 0 分")
          next
        }
        if (obs <= 0) {
          # **零分组合仍要进 rows（审计 S4 / 台账 E-58）。**
          # 掉它们 = BH 分母凭空小 43%，p_adj_bh 系统性偏小。
          # 零分组合的 p=1.0 是它在这个检验家庭里的正确取值。
          n_zero_dropped <- n_zero_dropped + 1L
          ri_ <- ri_ + 1L
          rows[[ri_]] <- list(pair = pr$name, sender = src, receiver = dst,
                              ligand = paste(pr$ligand, collapse = ","),
                              receptor = paste(pr$receptor, collapse = ","),
                              score = 0, null_mean = NA,
                              p_value = 1.0, n_permutations = 0L, tested = FALSE)
          next
        }
        # 置换：打乱细胞标签，看这个分数有多容易随机出现
        n_s <- length(masks[[src]]); n_d <- length(masks[[dst]])
        null <- numeric(N_PERMUTATIONS)
        set.seed(seed)   # 每 combo 重置 —— 与 Python 版 rng per-permutation 对齐
        for (i in seq_len(N_PERMUTATIONS)) {
          perm <- sample.int(n_cells)
          null[i] <- score_combo(pr$ligand, pr$receptor,
                                 perm[seq_len(n_s)],
                                 perm[n_s + seq_len(n_d)])
        }
        p <- (sum(null >= obs) + 1L) / (N_PERMUTATIONS + 1L)
        n_evaluated <- n_evaluated + 1L
        ri_ <- ri_ + 1L
        rows[[ri_]] <- list(pair = pr$name, sender = src, receiver = dst,
                            ligand = paste(pr$ligand, collapse = ","),
                            receptor = paste(pr$receptor, collapse = ","),
                            score = round(obs, 5),
                            null_mean = round(mean(null), 5),
                            p_value = round(p, 5),
                            n_permutations = N_PERMUTATIONS,
                            p_value_resolution = p_res,
                            p_value_at_floor = abs(p - p_res) < 1e-12,
                            tested = TRUE)
      }
    }
  }

  if (!length(rows)) {
    status <- list(dataset_id = cfg$dataset_id, status = "no_signal",
                   reason = "所有配体-受体对在所有簇对上的打分都为 0")
    write_json(file.path(res_dir, "communication_status.json"), status)
    log_warn(status$reason)
    return(status)
  }

  res <- do.call(rbind, lapply(rows, function(r) {
    r$null_mean <- if (is.na(r$null_mean)) NA_real_ else r$null_mean
    as.data.frame(r, stringsAsFactors = FALSE)
  }))
  res <- res[order(-res$score, na.last = TRUE), , drop = FALSE]
  # BH 校正：**检验家庭 = 全部组合，含零分那些**（审计 S4）。
  # NA p（非有限 obs）不能进 p.adjust —— NA 会传染整列全 NA（p.adjust(na.rm 不存在)）。
  # 语义：NA p 的组合**没有参与检验**，从检验家庭剔除；p_adj_bh 记 NA + note 列已说明。
  has_p <- !is.na(res$p_value)
  res$p_adj_bh <- NA_real_
  res$p_adj_bh[has_p] <- stats::p.adjust(res$p_value[has_p], method = "BH")
  res$p_adj_bh <- round(res$p_adj_bh, 5)
  utils::write.csv(res, file.path(res_dir, "cell_communication.csv"),
                   row.names = FALSE, na = "")
  n_sig <- sum(res$p_adj_bh < 0.05, na.rm = TRUE)
  log_info(sprintf("通讯打分: %d 个组合（%d 做了置换，%d 打分恒为 0，%d 打分非有限未检验），BH 后 %d 个 p<0.05",
                   nrow(res), n_evaluated, n_zero_dropped, n_nonfinite, n_sig))

  # 热图：配体-受体对 x 接收细胞类型 的总分
  # **按 pair 选前 20，不是按三元组选**（head(20) 拿到 20 个组合，
  # pivot 后塌成 3 个 pair —— 实测图 8/8 空白）。行序按总分从高到低。
  pair_score <- tapply(res$score, res$pair, sum, na.rm = TRUE)   # NA score 组合不拖垮 pair 总分
  pair_score <- sort(pair_score, decreasing = TRUE)
  n_pairs <- min(20L, length(pair_score))
  top_pairs <- utils::head(names(pair_score), n_pairs)
  top <- res[res$pair %in% top_pairs, , drop = FALSE]
  mat <- tapply(top$score, list(pair = top$pair, receiver = top$receiver), sum, na.rm = TRUE)
  mat[is.na(mat)] <- 0
  mat <- mat[intersect(top_pairs, rownames(mat)), , drop = FALSE]

  # 高度按实际行数算（不按请求数），宽夹到 W_ONE_HALF（上限若用 W_DOUBLE
  # 中间值算出非标宽，实测 153.7mm）；标题写**实际画出来的**行数。
  fig_h_mm <- min(W_DOUBLE, max(66, round(4.1 * nrow(mat) + 48)))
  p <- ggplot2::ggplot(
    as.data.frame(as.table(mat)),
    ggplot2::aes(x = receiver, y = pair, fill = Freq)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::scale_fill_viridis_c(name = "score") +
    ggplot2::labs(x = "receiver", y = "ligand-receptor pair",
                  title = sprintf("Top %d ligand-receptor pairs by summed score",
                                  nrow(mat))) +
    theme_paper(base_size = 8) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_fig(cfg, "02-06-01-unit1-communication-heatmap", p,
           width =  W_ONE_HALF, height =  fig_h_mm)
  log_info(sprintf("通讯热图: %d 个 pair x %d 个接收类型（%d 组合 -> %d pair，取前 %d）",
                   nrow(mat), ncol(mat), nrow(res), length(pair_score), n_pairs))

  # ---- liana（文档 §2.7 指定的主工具）------------------------------------
  # **必须传 named 版本**（名字=细胞名）：try_liana 里 SCE 列名要按名字显式
  # 匹配重排。run28 实锤：L258 的 as.character() 已剥掉名字，裸值传进去 →
  # 「liana SCE 有 2574 个细胞在簇标签向量里找不到」。03 落盘的 clu$clusters
  # 本身是 named（名字=细胞名，03_cluster_annotate.R:424）—— 这里直接传原
  # 对象；clusters_vec（裸值）留给上面的自建打分（masks 只用值）。
  liana_out <- try_liana(clu, clu$clusters, cfg)
  liana_info <- liana_out$info
  liana_cmp <- list(compared = FALSE)
  if (!is.null(liana_out$res)) {
    utils::write.csv(as.data.frame(liana_out$res),
                     file.path(res_dir, "liana_results.csv"), row.names = FALSE)
    liana_cmp <- compare_with_liana(as.data.frame(liana_out$res), res)
    if (isTRUE(liana_cmp$compared)) {
      liana_cmp$interpretation <- paste0(
        "liana 的 rank aggregate 综合了多种方法的秩，自建打分只是表达量乘积",
        " —— **两者排序不同是预期内的，不是谁错了**。这里报出来是为了让",
        "『差多少』有数字。")
    }
  }

  status <- list(
    dataset_id = cfg$dataset_id,
    status = "ok",
    group_key = group_key,
    n_groups = length(groups),
    n_pairs_in_database = length(pairs),
    n_pairs_usable = length(present),
    n_genes_used = ncol(X),
    target_sum = target_sum,
    score_scale_note = norm_note,
    score_definition = paste0("`mean(配体基因在 sender 的表达) × ",
                              "mean(受体基因在 receiver 的表达)` —— **两个均值相乘，",
                              "不是几何平均、不是最大值**；量纲随 target_sum 整体缩放"),
    gene_set_note = paste0("用的是**全基因集**（logcounts_all），不是 HVG 子集。",
                           "配体/受体基因大多是低表达的细胞因子/趋化因子，几乎不进 HVG",
                           " —— 实测在 HVG 子集上 38 对里只有 3 对可用，全基因集 27 对，",
                           "差 9 倍"),
    n_combinations = nrow(res),
    n_combinations_tested = n_evaluated,
    n_combinations_zero_score = n_zero_dropped,
    n_significant_bh = n_sig,
    bh_family_note = sprintf(paste0("**BH 的检验家庭 = 全部组合（%d 个），含打分为 0 的 %d 个**",
                                    "（它们的 p 记为 1.0）。旧实现把零分组合 continue 掉，",
                                    "分母只剩 %d 个，p_adj_bh 系统性偏小（审计 S4）"),
                             nrow(res), n_zero_dropped, n_evaluated),
    n_permutations = N_PERMUTATIONS,
    p_value_resolution = p_res,
    p_value_note = sprintf(paste0("置换检验的 p 只能取 1/(N+1) 的整数倍 —— N=%d 时最小非零值",
                                  "是 %.6f，共 %d 个离散档位。p_value_at_floor=TRUE 表示该组合的",
                                  "观测分数一次都没被置换超越，它是**分辨率上限**而不是精确估计"),
                           N_PERMUTATIONS, p_res, N_PERMUTATIONS + 1L),
    top_pairs = df_to_records(utils::head(res, 15L)),
    liana = liana_info,
    liana_vs_builtin = liana_cmp,
    method = paste0(
      "**两条路都跑了**：(1) 自建数据库驱动的配体-受体共表达打分 + 簇标签置换检验 + BH 校正 → ",
      "cell_communication.csv；(2) ",
      if (identical(liana_info$status, "ok"))
        sprintf("liana_wrap rank_aggregate v%s → liana_results.csv", liana_info$version)
      else
        sprintf("liana **未能运行**（%s：%s）", liana_info$status,
                substr(as.character(liana_info$reason), 1L, 120L))),
    primary_tool_per_spec = "liana（文档 §2.7；R 版经 liana_wrap，Python 版经 li.mt.rank_aggregate）",
    limitations = c(
      "共表达不等于通讯：没有空间信息时，只能说两类细胞分别表达了配体和受体",
      "表达量是稳态丰度，不等于蛋白水平，也不等于分泌量",
      "自建打分是启发式，不是 liana 的 consensus rank aggregate",
      "置换检验打乱的是细胞标签，保留了每种细胞类型的细胞数",
      sprintf(paste0("**置换检验的 p 值有分辨率下限**：N=%d 时最小非零值是 1/%d ≈ %.6f，",
                     "全表只有 %d 个可能的 p 值。p_value_at_floor 标出的组合是**分辨率上限**，",
                     "不要当成精确估计"), N_PERMUTATIONS, N_PERMUTATIONS + 1L,
              p_res, N_PERMUTATIONS + 1L),
      sprintf("**内置库只有 %d 对**（免疫为主），远少于 CellChatDB 的数千对；覆盖不全时",
              length(pairs)),
      if (identical(liana_info$status, "ok"))
        "liana 用的是它自带的 consensus 资源（CellChatDB + CellPhoneDB 等），与内置库不是同一套配体-受体对"
      else
        "liana 未运行，本轮只有自建打分这一条路",
      "CellChat（文档 §2.7 的辅工具）在 R 版**可装**（Bioconductor），本轮按 §3 裁决未接 —— 两条路的自建/liana 主链已覆盖；如需 CellChat 由 R-17 评估")
  )
  write_json(file.path(res_dir, "communication_status.json"), status)
  status
}

# ---------------------------------------------------------------------------
# 入口（E-56：no_pairs_in_data / no_signal 是早退路径，不接住返回值就记成 ok）
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L || identical(Sys.getenv("SCRNA_STEP_MAIN"), "06_communication")) {
  args <- parse_args()
  cfg <- load_config(args$config)
  t0 <- Sys.time()
  tryCatch({
    out <- run_06_communication(cfg)
    record_step(cfg, "communication", "ok",
                as.numeric(difftime(Sys.time(), t0, units = "secs")),
                result_status = result_status_of(out))
  }, error = function(e) {
    record_step(cfg, "communication", "failed",
                as.numeric(difftime(Sys.time(), t0, units = "secs")),
                message = conditionMessage(e), result_status = "failed")
    stop(e)
  })
}
