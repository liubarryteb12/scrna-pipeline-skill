#!/usr/bin/env Rscript
# =============================================================================
# 08_virtual_perturbation.R — 虚拟扰动（与 08_virtual_perturbation.py 同号同语义）
#
# **先说清楚这不是什么：不是 PerturbNet。** 本步骤两个引擎：
#   - 一阶引擎：沿 07 的调控子边传 Δz（线性、单跳）—— 便宜、可解释、
#     但它给出的"影响"严格说是"网络 expresses 这个扰动"，不是因果预测。
#   - Tenifold 引擎：调 scripts/lib/tenifoldKnk.R（既有资产，等级 D
#     零差异——两版跑**同一个** R 脚本，数值可比）。它自己在 R 里建
#     基因调控网络再做流形距离比较。
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

get_script_path <- function() {
  argv0 <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(argv0)) return(normalizePath(sub("^--file=", "", argv0[1]), mustWork = FALSE))
  getOption("scrna.script_path", "")
}

DEFAULT_TOP_N <- 20L
MIN_TARGETS <- 5L
DEFAULT_OVEREXPRESS_SD <- 1.0
ENGINES <- c("tenifold", "first_order", "both")
DEFAULT_TENIFOLD <- list(max_genes = 600L, n_net = 10L, n_cells = 500L,
                         n_comp = 3L, td_k = 3L, ma_n_dim = 2L, n_cores = 1L,
                         timeout_sec = 1500L)
# 网络基因数是唯一要紧的旋钮：scTenifoldNet 运行时间表 1000 基因 x 5000 细胞
# ≈ 2h52min，全基因集绝对跑不动；CI 上限 45min。
TENIFOLD_R <- "tenifold_knk.R"

# ---------------------------------------------------------------------------
# 全基因集表达 + 原始计数（08_virtual_perturbation.py:167-403 同语义）
#
# **不能用 clustered.rds 里的 counts**：02_integrate.R 存进 clustered.rds
# 的 counts 已缩到 HVG 子集。全基因集原始计数只在 qc_filtered.rds
# （01_qc 落盘、未 normalize）。细胞顺序不一致时显式对齐 + WARN ——
# 错位不会报错，只会把每个细胞的基因表达配给别的细胞。
# ---------------------------------------------------------------------------
open_full_counts <- function(cfg, clu) {
  qc_path <- file.path(cfg$output$data_dir, "qc_filtered.rds")
  if (!file.exists(qc_path)) {
    stop(sprintf("缺 qc_filtered.rds —— 08 的全基因集原始计数只能来自 step 01 的落盘（clustered.rds 的 counts 已缩到 HVG 子集）"),
                 call. = FALSE)
  }
  qc <- readRDS(qc_path)
  qc_counts <- qc$counts
  qc_cells <- rownames(qc$metrics)
  clu_cells <- rownames(clu$counts)
  if (!identical(qc_cells, clu_cells)) {
    n_miss <- sum(!clu_cells %in% qc_cells)
    if (n_miss > 0L) {
      stop(sprintf("%d 个聚类数据里的细胞在 qc_filtered.rds 找不到 —— 无法对齐",
                   n_miss), call. = FALSE)
    }
    log_warn("qc_filtered 与聚类数据的细胞顺序不一致 —— 按细胞名显式重排对齐")
    idx <- match(clu_cells, qc_cells)
    qc_counts <- qc_counts[idx, , drop = FALSE]
  }
  qc_counts
}

# ---------------------------------------------------------------------------
# 目标基因（load_targets 同语义）：优先 Part 1 交接 CSV，否则内部调控子
# ---------------------------------------------------------------------------
load_targets <- function(cfg, data_dir, res_dir) {
  pert <- (cfg$perturbation %||% list())
  top_n <- as.integer(pert$top_n %||% DEFAULT_TOP_N)
  targets_csv <- pert$targets_csv
  if (!is.null(targets_csv)) {
    p <- targets_csv
    if (!grepl("^([A-Za-z]:[\\\\/]|/)", p)) p <- file.path(data_dir, p)
    if (!file.exists(p)) {
      return(list(df = NULL, target_source = sprintf("configured_csv_missing:%s", p),
                  cross_language = NULL))
    }
    df <- utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
    cn <- names(df)
    cnl <- tolower(cn)
    gcol <- cn[match(c("gene", "target", "symbol"), cnl)][1L]
    if (is.na(gcol)) {
      return(list(df = NULL, target_source = sprintf("configured_csv_no_gene_column:%s", p),
                  cross_language = NULL))
    }
    lcol <- cn[grepl("^log2?fc(_|-)?", cnl)][1L]
    kept <- c(gcol, if (!is.na(lcol)) lcol)
    lost <- setdiff(cn, kept)
    record_cross_language(cfg, tool = "part1_de_to_part8",
                          direction = "in", keys = kept,
                          lost = list(lost),
                          note = "Part 1 差异表达交接：只保留 gene/logfc，其余列丢弃")
    out <- data.frame(gene = as.character(df[[gcol]]),
                      logfc = if (!is.na(lcol)) as.numeric(df[[lcol]]) else NA_real_,
                      stringsAsFactors = FALSE)
    out <- out[!is.na(out$gene) & out$gene != "", , drop = FALSE]
    out <- out[!duplicated(out$gene), , drop = FALSE]
    list(df = out, target_source = sprintf("part1_csv:%s", p),
         cross_language = list(n_rows = nrow(out), kept = kept, lost = list(lost)))
  } else {
    reg_path <- file.path(res_dir, "tf_regulons.csv")
    if (!file.exists(reg_path)) return(list(df = NULL, target_source = "no_tf_regulons",
                                            cross_language = NULL))
    reg <- utils::read.csv(reg_path, stringsAsFactors = FALSE)
    # **先按 TF 去重再取前 N**：tf_regulons.csv 是行级排名，同一 TF 可能
    # 有多行（按 best_cluster 拆）——不去重的话 top_n 会虚报。
    reg <- reg[!duplicated(reg$tf), , drop = FALSE]  # 已按特异性降序，保首行
    out <- data.frame(gene = utils::head(reg$tf, top_n), logfc = NA_real_,
                      stringsAsFactors = FALSE)
    list(df = out, target_source = sprintf("internal_regulons:top_%d", top_n),
         cross_language = NULL)
  }
}

# ---------------------------------------------------------------------------
# 一阶引擎（run_first_order_engine 同语义）
# ---------------------------------------------------------------------------
run_first_order_engine <- function(targets, edges, Z, cluster, delta_sd) {
  cells <- rownames(Z)
  clusters <- sort(unique(cluster))
  rows <- list(); skipped <- list(); ri <- 0L
  # 边查表：源 TF -> (靶基因下标, 权重)
  edge_map <- split(edges, edges$tf)

  for (g in targets$gene) {
    gi <- match(g, colnames(Z))
    if (is.na(gi)) {
      skipped[[length(skipped) + 1L]] <- list(gene = g, reason = "不在表达矩阵里")
      next
    }
    eg <- edge_map[[g]]
    if (is.null(eg)) {
      skipped[[length(skipped) + 1L]] <- list(gene = g, reason = "调控子无出边")
      next
    }
    tgt_idx <- match(eg$target, colnames(Z))
    ok <- !is.na(tgt_idx)
    w <- as.numeric(eg$corr)[ok]
    tgt_idx <- tgt_idx[ok]
    if (length(tgt_idx) < MIN_TARGETS) {
      skipped[[length(skipped) + 1L]] <- list(
        gene = g, reason = sprintf("可用靶基因只有 %d 个（<%d）", length(tgt_idx), MIN_TARGETS))
      next
    }
    # **network_sensitivity 与表达无关**，单独报：ko_magnitude = |Δz|×sens
    # 把"基因表达高"和"网络连接强"混在一个数里 —— 排在前面的可能只是表达高。
    sens <- sqrt(sum(w^2))
    for (ct in clusters) {
      sel <- cluster == ct
      if (sum(sel) < 3L) next
      mean_z_g <- mean(Z[sel, gi])
      mean_z_t <- colMeans(Z[sel, tgt_idx, drop = FALSE])
      # 敲除：把该基因的表达归零 = Δz = -z̄；过表达：+delta_sd
      dz_g <- -mean_z_g
      dz_t <- w * dz_g
      rows[[length(rows) + 1L]] <- list(
        gene = g, celltype = ct, n_cells = sum(sel),
        mean_z_in_celltype = round(mean_z_g, 4),
        network_sensitivity = round(sens, 4),
        ko_magnitude = round(abs(dz_g) * sens, 4),
        oe_magnitude = round(delta_sd * sens, 4),
        mean_abs_delta_z = round(mean(abs(dz_t)), 4),
        max_abs_delta_z = round(max(abs(dz_t)), 4),
        n_targets = length(tgt_idx))
    }
  }
  list(rows = do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE)),
       skipped = skipped)
}

# ---------------------------------------------------------------------------
# Tenifold 引擎（run_tenifold_engine 同语义）：R 版直接 source 既有资产，
# 不走 subprocess —— 同一个 R 进程，同一份代码（等级 D 零差异）。
# ---------------------------------------------------------------------------
select_tenifold_genes <- function(var_names, cand_genes, max_genes, variance) {
  cand_genes <- unique(intersect(cand_genes, var_names))
  # ① 全部候选必须在内（包对不在网络里的 gKO 是 stop 不是跳过）
  # ② 其余名额给方差最高的基因 ③ **输出按 var_names 原顺序**（网络断言
  # 行名逐位相符；顺序错位不报错只算垃圾网络）
  in_cand <- var_names %in% cand_genes
  rest <- setdiff(seq_along(var_names)[!in_cand],
                  seq_along(var_names)[in_cand])
  fill <- rest[order(variance[rest], decreasing = TRUE)]
  n_fill <- max(0L, max_genes - sum(in_cand))
  sel <- sort(unique(c(which(in_cand), utils::head(fill, n_fill))))
  sel
}

# 6 位**有效数字**（不是小数位）：signif(v, 6)。round(v, 6) 曾把 1.48e-08
# 归 0、摧毁排序（K-01b 真 bug）—— 十进制有效位用 signif 才对。
sig6 <- function(v) signif(v, 6L)

run_tenifold_engine <- function(cfg, X_counts, var_names, targets, tk) {
  tmp_dir <- file.path(cfg$output$results_dir, "_tenifold_tmp")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  res <- list(ok = FALSE)
  tryCatch({
    need_pkg("Matrix", "Tenifold 引擎读稀疏计数")
    cand <- targets$gene
    gi <- match(cand, var_names)
    if (any(is.na(gi))) {
      return(list(ok = FALSE, reason = sprintf("%d 个候选基因不在全基因集里",
                                               sum(is.na(gi)))))
    }
    cvar <- Matrix::colMeans(X_counts^2) - Matrix::colMeans(X_counts)^2
    sel <- select_tenifold_genes(var_names, cand, tk$max_genes, cvar)
    sel_names <- var_names[sel]

    # **方向见证（gene_sums）**：矩阵写文件前后各算一次行和 —— 方阵转置
    # 维度查不出转置错误，行和能（细胞 sum ≫ 基因 sum）。
    sub <- as.matrix(X_counts[, sel, drop = FALSE])
    sub <- t(sub)                                   # genes x cells（关键方向）
    gene_sums_before <- rowSums(sub)

    mat_file <- file.path(tmp_dir, "expr_matrix.txt")
    # 表头必须是**基因名**（n_genes 个）—— read_counts 把第 1 行当行名挂到
    # genes x cells 矩阵上（与 Python 08_virtual_perturbation.py:589 同契约）。
    # 曾错写 colnames(sub)（细胞名 2574 个）→ dimnames 长度不等直接崩。
    writeLines(c(paste(rownames(sub), collapse = "\t"),
                 apply(sub, 1L, paste, collapse = "\t")), mat_file)
    # meta 字段名与 tenifold_knk.R assert_orientation 的契约一致：
    # genes（有序全量）/ gene_sums（每基因总计数，方向见证），不是别名。
    meta <- list(n_genes = nrow(sub), n_cells = ncol(sub),
                 genes = sel_names,
                 gene_sums = as.numeric(gene_sums_before))
    write_json(file.path(tmp_dir, "meta.json"), meta)

    # targets 文件是**纯文本每行一个基因**（readLines 读，L243）——
    # 不能写 CSV 表头 `gene`，否则第一个"基因"就是字面量 "gene"（Python
    # 08_virtual_perturbation.py:599 同契约 write_text("\n".join(cand))）。
    writeLines(cand, file.path(tmp_dir, "targets.txt"))

    rscript <- file.path(dirname(get_script_path()), "lib", TENIFOLD_R)
    if (!file.exists(rscript)) {
      return(list(ok = FALSE, reason = sprintf("缺 %s", rscript)))
    }
    out_csv <- file.path(tmp_dir, "dist.csv")
    rmeta_json <- file.path(tmp_dir, "run_meta.json")
    cmd <- sprintf(paste(
      '"%s" "%s" --input "%s" --meta "%s" --targets "%s" --out_dir "%s"',
      "--n_net %d --n_cells %d --n_comp %d --td_k %d --ma_n_dim %d --n_cores %d"),
      file.path(R.home("bin"), "Rscript"), rscript, mat_file,
      file.path(tmp_dir, "meta.json"), file.path(tmp_dir, "targets.txt"),
      tmp_dir, tk$n_net, tk$n_cells, tk$n_comp, tk$td_k, tk$ma_n_dim, tk$n_cores)
    log_info(sprintf("Tenifold 引擎: %d 基因 x %d 细胞, n_net=%d（CI 预算内）",
                     nrow(sub), ncol(sub), tk$n_net))
    st <- system(cmd, timeout = tk$timeout_sec)
    if (!identical(st, 0L)) {
      return(list(ok = FALSE,
                  reason = sprintf("tenifold_knk.R 退出码 %s（矩阵与 meta 保留在 %s 供排查）",
                                   st, tmp_dir)))
    }
    if (!file.exists(out_csv)) {
      return(list(ok = FALSE, reason = "tenifold_knk.R 未产出 dist.csv"))
    }
    dist_df <- utils::read.csv(out_csv, stringsAsFactors = FALSE, check.names = FALSE)
    rmeta <- tryCatch(read_json(file.path(tmp_dir, "run_meta.json")),
                      error = function(e) list())
    empty_ko <- as.character(unlist(rmeta$empty_knockout_genes %||% list()))
    # 成功才清理临时目录（失败时矩阵是「R 吃到了什么」的唯一证据）
    unlink(tmp_dir, recursive = TRUE)
    res$ok <- TRUE
    res$dist <- dist_df
    res$empty_knockout_genes <- empty_ko
    res$n_genes_scored <- rmeta$n_genes %||% nrow(sub)
    res$engine_version <- rmeta$engine_version %||% TENIFOLD_R
    res
  }, error = function(e) {
    list(ok = FALSE, reason = sprintf("%s: %s", class(e)[1L], conditionMessage(e)))
  })
  res
}

# 距离矩阵 -> 每基因一行（compress_tenifold_distances 同语义）
compress_tenifold_distances <- function(dist_df, res) {
  # dist.csv：第 1 列基因名，其余列同名基因（对称矩阵）
  m <- as.matrix(dist_df[, -1L, drop = FALSE])
  rownames(m) <- dist_df[[1L]]
  rows <- list()
  for (g in sort(unique(c(rownames(m), colnames(m))))) {
    v <- c(m[g, ], m[, g])
    # 两个向量拼起来有重复（对称位计两次）——先去重保名称，再排
    v <- v[!duplicated(names(v))]
    v <- v[names(v) != g]
    v <- v[is.finite(v)]
    if (!length(v)) {
      rows[[length(rows) + 1L]] <- list(
        gene = g, tenifold_n_genes_scored = 0L,
        mean_distance = NA_real_, max_distance = NA_real_,
        top_gene = NA_character_, top_distance = NA_real_,
        target_outdegree = 0L,
        empty_knockout = TRUE,   # 结构证据（rmeta 的 empty_knockout_genes），不从整行 NA 反推
        top_distance_real = NA_real_)
      next
    }
    o <- order(v, decreasing = TRUE)
    rows[[length(rows) + 1L]] <- list(
      gene = g, tenifold_n_genes_scored = length(v),
      mean_distance = sig6(mean(v)), max_distance = sig6(max(v)),
      top_gene = names(o)[1L], top_distance = sig6(v[o[1L]]),
      target_outdegree = length(v), empty_knockout = FALSE,
      top_distance_real = sig6(v[o[1L]]))
  }
  out <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  out <- out[order(-out$mean_distance), , drop = FALSE]  # 与一阶 ko_magnitude 降序同向
  out
}

# 两引擎对比：**只比排名不比数值**（量纲不同：z 分数 L2 范数 vs 流形欧氏距离）
compare_engines <- function(fo, td) {
  if (is.null(fo) || !nrow(fo) || is.null(td) || !nrow(td)) return(NULL)
  fo_best <- fo[order(-fo$ko_magnitude), , drop = FALSE]
  fo_best <- fo_best[!duplicated(fo_best$gene), , drop = FALSE]
  # 空敲除基因显式剔除记账（悄悄 dropna 会少点不报）
  excl <- td$gene[td$empty_knockout]
  td2 <- td[!td$empty_knockout, , drop = FALSE]
  genes <- intersect(fo_best$gene, td2$gene)
  if (length(genes) < 3L) {
    return(list(status = "insufficient_overlap",
                n_shared = length(genes),
                excluded_empty_knockout = as.list(excl),
                excluded_reason = "这些候选在网络里没有可用边（空敲除），tenifold 距离全 NA",
                note = "共享基因不足 3 个，两引擎排名相关不报（两点永能连线）"))
  }
  a <- fo_best$ko_magnitude[match(genes, fo_best$gene)]
  b <- td2$mean_distance[match(genes, td2$gene)]
  ct <- suppressWarnings(stats::cor.test(a, b, method = "spearman"))
  top5_fo <- utils::head(genes[order(a, decreasing = TRUE)], 5L)
  top5_td <- utils::head(genes[order(b, decreasing = TRUE)], 5L)
  list(status = "ok",
       n_shared = length(genes),
       excluded_empty_knockout = as.list(excl),
       spearman_rho = round(as.numeric(ct$estimate), 4),
       p_value = as.numeric(ct$p.value),
       top5_overlap = length(intersect(top5_fo, top5_td)),
       note = paste0("**只比排名不比数值** —— 两引擎量纲不同（z 分数 L2 范数 vs 流形欧氏距离）。",
                     "一致性高也可能只是共享同一偏差（都从 07 的调控子边出发）"))
}

# ---------------------------------------------------------------------------
# 主流程（run_08_virtual_perturbation 同语义）
# ---------------------------------------------------------------------------
run_08_virtual_perturbation <- function(cfg) {
  ensure_dirs(cfg)
  set_seed(cfg)
  data_dir <- cfg$output$data_dir
  res_dir <- cfg$output$results_dir
  # 规则 14：开跑前删自己上一轮的 status（旧状态残留会被验收层当成本轮产物）
  unlink(file.path(res_dir, "virtual_perturbation_status.json"))

  pert <- (cfg$perturbation %||% list())
  if (!isTRUE(pert$enabled %||% TRUE)) {
    status <- list(dataset_id = cfg$dataset_id, status = "disabled",
                   reason = "配置 perturbation.enabled=false")
    write_json(file.path(res_dir, "virtual_perturbation_status.json"), status)
    return(status)
  }

  edges_path <- file.path(res_dir, "tf_regulon_edges.csv")
  if (!file.exists(edges_path)) {
    status <- list(dataset_id = cfg$dataset_id, status = "not_configured",
                   reason = "缺 tf_regulon_edges.csv（step 07 未产出）—— 虚拟扰动需要调控子边做先验")
    write_json(file.path(res_dir, "virtual_perturbation_status.json"), status)
    return(status)
  }
  edges <- utils::read.csv(edges_path, stringsAsFactors = FALSE)

  engine <- pert$engine %||% "both"
  if (!engine %in% ENGINES) {
    # **拼错就 raise，不静默退回** —— 静默退回会让"想跑 tenifold"的配置
    # 悄悄变成只跑一阶，产物看不出来。
    stop(sprintf("perturbation.engine='%s' 不在 [%s]", engine,
                 paste(ENGINES, collapse = ", ")), call. = FALSE)
  }

  clu_path <- file.path(data_dir, "clustered.rds")
  if (!file.exists(clu_path)) {
    stop(sprintf("缺输入 %s —— 先跑 03_cluster_annotate.R", clu_path), call. = FALSE)
  }
  clu <- readRDS(clu_path)

  # 细胞类型列：三写法都认；退回 leiden 打 WARN（解读完全变了）
  meta <- clu$cell_meta
  ckey <- pert$celltype_key
  ckey_candidates <- c(ckey, "celltype", "cell_type", "celltype_label")
  ckey_used <- NULL
  for (k in ckey_candidates) {
    if (!is.null(k) && k %in% colnames(meta)) { ckey_used <- k; break }
  }
  ckey_degraded <- FALSE
  if (is.null(ckey_used)) {
    ckey_used <- "leiden"
    ckey_degraded <- TRUE
    log_warn(paste0("找不到细胞类型列（celltype_key 及 celltype/cell_type/",
                    "celltype_label 都不在 meta）——退回 leiden 簇标签。",
                    "**解读完全变了**：下面报告的是「扰动对簇的影响」，",
                    "不是「对细胞类型的影响」"))
  }
  cluster <- as.character(meta[[ckey_used]])

  # ---- 表达矩阵与 z-score -------------------------------------------------
  X_counts <- open_full_counts(cfg, clu)
  genes_all <- colnames(X_counts)
  cells_all <- rownames(X_counts)
  # 逐基因 z-score（Δz 与边权同尺度）
  mu <- Matrix::colMeans(X_counts)
  sdv <- sqrt(pmax(Matrix::colMeans(X_counts^2) - mu^2, 0))
  sdv[sdv == 0] <- 1
  Z <- as.matrix(sweep(sweep(X_counts, 2L, mu), 2L, sdv, "/"))
  rownames(Z) <- cells_all; colnames(Z) <- genes_all

  # ---- 目标基因 -----------------------------------------------------------
  tgt <- load_targets(cfg, data_dir, res_dir)
  targets <- tgt$df
  if (is.null(targets) || !nrow(targets)) {
    status <- list(dataset_id = cfg$dataset_id, status = "no_candidates",
                   target_source = tgt$target_source,
                   reason = "没有可扰动的目标基因")
    write_json(file.path(res_dir, "virtual_perturbation_status.json"), status)
    return(status)
  }
  log_info(sprintf("扰动目标: %d 个（来源 %s）", nrow(targets), tgt$target_source))

  delta_sd <- as.numeric(pert$overexpress_sd %||% DEFAULT_OVEREXPRESS_SD)

  # ---- 2a 一阶引擎 --------------------------------------------------------
  fo <- NULL; fo_skipped <- list()
  if (engine %in% c("first_order", "both")) {
    r1 <- run_first_order_engine(targets, edges, Z, cluster, delta_sd)
    fo <- r1$rows; fo_skipped <- r1$skipped
  }

  # ---- 2b Tenifold 引擎 ---------------------------------------------------
  td <- NULL; td_ok <- FALSE; td_reason <- NULL; td_res <- NULL
  if (engine %in% c("tenifold", "both")) {
    tools_probe <- list(tenifold_r = file.exists(
      file.path(dirname(get_script_path()), "lib", TENIFOLD_R)))
    if (!tools_probe$tenifold_r) {
      td_reason <- sprintf("缺 lib/%s", TENIFOLD_R)
    } else {
      tk <- utils::modifyList(DEFAULT_TENIFOLD,
                              as.list(pert$tenifold %||% list()))
      td_res <- run_tenifold_engine(cfg, X_counts, genes_all, targets, tk)
      if (isTRUE(td_res$ok)) {
        td_ok <- TRUE
        td <- compress_tenifold_distances(td_res$dist, td_res)
        utils::write.csv(td, file.path(res_dir, "virtual_perturbation_tenifold.csv"),
                         row.names = FALSE)
        # Python 版还落 _tenifold_top.csv（08_virtual_perturbation.py:655-656
        # head(10)）—— 图名集合/产物集合两版必须相等（r_version.md §5），
        # R 版补齐同款。
        utils::write.csv(utils::head(td, 10L),
                         file.path(res_dir, "virtual_perturbation_tenifold_top.csv"),
                         row.names = FALSE)
        utils::write.csv(td_res$dist,
                         file.path(res_dir, "tenifold_perturbation_distances.csv"),
                         row.names = FALSE)
      } else {
        td_reason <- td_res$reason
        log_warn(sprintf("Tenifold 引擎未产出: %s", td_reason))
      }
    }
  }

  if ((is.null(fo) || !nrow(fo)) && !td_ok) {
    status <- list(dataset_id = cfg$dataset_id, status = "no_candidates",
                   target_source = tgt$target_source,
                   skipped_first_order = fo_skipped,
                   reason = "两个引擎都没有产出任何可比较的结果")
    write_json(file.path(res_dir, "virtual_perturbation_status.json"), status)
    return(status)
  }

  # ---- 一阶落盘 -----------------------------------------------------------
  use_td <- FALSE
  if (!is.null(fo) && nrow(fo)) {
    fo <- fo[order(-fo$ko_magnitude), , drop = FALSE]
    utils::write.csv(fo, file.path(res_dir, "virtual_perturbation.csv"), row.names = FALSE)
    top <- fo[order(-fo$ko_magnitude), , drop = FALSE]
    top <- top[!duplicated(top$gene), , drop = FALSE]
    names(top)[names(top) == "celltype"] <- "most_affected_celltype"
    utils::write.csv(top, file.path(res_dir, "virtual_perturbation_top.csv"),
                     row.names = FALSE)
  }

  # ---- 出图：一图只画一引擎（量纲不同不能同 x 轴）------------------------
  # Python 版结构：if use_td / else 各自 build `fig`，**分支外统一一次
  # save_fig** —— 图名只有一个（同名 ≠ 两张图，运行时只走一个分支）。
  # R 版同构：分支里只构建 p，末尾统一 save。
  p <- NULL
  if (td_ok && !all(td$empty_knockout)) {
    use_td <- TRUE
    d <- td[!td$empty_knockout & is.finite(td$mean_distance), , drop = FALSE]
    d <- utils::head(d[order(-d$mean_distance), , drop = FALSE], 20L)
    if (nrow(d)) {
      # **画点不画柱**：距离跨 3~4 个数量级（1.1e-09 ~ 4.5e-06）。
      # 线性轴上弱者柱长为零；对数轴的柱长会被读者当成倍数。
      # 位置编码（对数轴上的点）+ 每点右标真实数值。
      d$label <- sprintf("%s (%s)", d$gene,
                         formatC(d$mean_distance, format = "g", digits = 2))
      d$y <- seq_len(nrow(d))
      rng <- range(d$mean_distance)
      p <- ggplot2::ggplot(d, ggplot2::aes(mean_distance, y)) +
        ggplot2::geom_point(size = 1.4, colour = PAL$primary) +
        ggplot2::geom_text(ggplot2::aes(label = label), hjust = -0.05,
                           size = 2.1, nudge_y = 0.25) +
        ggplot2::scale_x_log10(limits = c(rng[1L] * 0.45, rng[2L] * 4.5)) +
        ggplot2::scale_y_reverse(breaks = NULL) +
        ggplot2::labs(x = "mean manifold distance after knockout (log scale)",
                      y = NULL,
                      title = "Virtual knockout: predicted effect size\n(scTenifoldKnk tensor network, manifold distance; log axis)") +
        theme_paper(base_size = 8) +
        ggplot2::theme(legend.position = "none")
      fig_w <- min(W_ONE_HALF, max(W_SINGLE, 0.30 * nrow(d) + 3.0 * 25.4))
    }
  } else if (!is.null(fo) && nrow(fo)) {
    # 一阶退回图：barh + ‖Δz‖₂
    d <- utils::head(fo, 20L)
    d$label <- sprintf("%s (%s)", d$gene, d$celltype)
    d$y <- factor(seq_len(nrow(d)), levels = rev(seq_len(nrow(d))),
                  labels = rev(d$label))
    p <- ggplot2::ggplot(d, ggplot2::aes(ko_magnitude, y)) +
      ggplot2::geom_col(fill = PAL$primary, width = 0.7) +
      ggplot2::labs(x = "knockout magnitude |Δz| × √Σw²", y = NULL,
                    title = "Virtual knockout: predicted effect size\n(one-hop linear propagation on a co-expression GRN)") +
      theme_paper(base_size = 8) +
      ggplot2::theme(legend.position = "none")
    fig_w <- W_ONE_HALF
  }
  if (!is.null(p)) {
    # save_fig 的签名是 width=/height=（common.R L571）—— 不是
    # fig_width_mm/fig_height_mm（那是 Python 版 save_fig 的关键字）。
    # R 的 ... 会把错名参数静默吞掉，图会以默认尺寸写出 —— 无声漂移。
    save_fig(cfg, "02-08-01-unit1-virtual-perturbation-effect", p,
             width = fig_w, height = mm(70))
  }

  comp <- if (td_ok && !is.null(fo) && nrow(fo)) compare_engines(fo, td) else NULL

  # ---- status -------------------------------------------------------------
  method_bits <- c()
  if (!is.null(fo) && nrow(fo)) {
    method_bits <- c(method_bits, paste0(
      "一阶引擎：Δz 沿 07 调控子边单跳传播（线性、无网络重构）——",
      "它给出的是「网络表达了这一扰动」的响应，不是因果预测"))
  }
  if (td_ok) {
    method_bits <- c(method_bits, sprintf(paste0(
      "Tenifold 引擎（%s）：基因调控网络 + 流形距离比较，",
      "**不是 PerturbNet**"), td_res$engine_version))
  }
  if (!td_ok && !is.null(td_reason)) {
    method_bits <- c(method_bits, sprintf("Tenifold 引擎未产出：%s", td_reason))
  }
  n_empty <- if (td_ok) sum(td$empty_knockout) else 0L
  status <- list(
    dataset_id = cfg$dataset_id,
    status = "ok",
    engines_used = c(if (!is.null(fo) && nrow(fo)) "first_order",
                     if (td_ok) "tenifold"),
    engine_configured = engine,
    target_source = tgt$target_source,
    n_targets = nrow(targets),
    celltype_key = ckey_used,
    celltype_key_degraded_to_leiden = ckey_degraded,
    overexpress_sd = delta_sd,
    n_first_order_rows = if (!is.null(fo)) nrow(fo) else 0L,
    first_order_skipped = fo_skipped,
    tenifold_ran = td_ok,
    tenifold_figure_drawn = use_td,
    n_empty_knockouts = n_empty,
    method_bits = method_bits,
    engine_comparison = comp,
    top_by_effect = df_to_records(utils::head(
      if (!is.null(fo) && nrow(fo)) fo else
        data.frame(gene = td$gene, magnitude = td$mean_distance), 10L)),
    limitations = c(
      "一阶引擎是线性单跳传播：非因果、不做多跳",
      "ko_magnitude = |Δz| × √Σw² 把「表达高」和「网络连接强」混在一个数里 —— 排前可能只是表达高（network_sensitivity 已分开报）",
      if (td_ok) paste0(
        "Tenifold 网络没有细胞类型分辨率 —— 结果是全群体网络的扰动，",
        "**不按细胞类型重新加权**（造出一个「分细胞类型的 Tenifold」",
        "再借它的名字发出去，是方法学造假）"),
      if (td_ok) sprintf(paste0("网络基因是全基因集的子集（max_genes=%d，", DEFAULT_TENIFOLD$max_genes),
                         "运行时间约束）；输入是原始计数（qc=FALSE，step 01 已滤过）"),
      if (td_ok && n_empty > 0L) sprintf(paste0(
        "%d 个候选在网络里出度为 0（空敲除）——", n_empty),
        "不表示「无影响」，表示这个网络表达不了该扰动"),
      "两引擎一致性高也可能只是共享同一偏差（都从 07 的调控子边出发）"))
  write_json(file.path(res_dir, "virtual_perturbation_status.json"), status)
  log_info(sprintf("虚拟扰动完成: engines=%s",
                   paste(status$engines_used, collapse = "+")))
  status
}

# ---------------------------------------------------------------------------
# 入口（E-56：no_candidates / not_configured / disabled 都是早退路径）
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L || identical(Sys.getenv("SCRNA_STEP_MAIN"), "08_virtual_perturbation")) {
  args <- parse_args()
  cfg <- load_config(args$config)
  t0 <- Sys.time()
  tryCatch({
    res <- run_08_virtual_perturbation(cfg)
    record_step(cfg, "virtual_perturbation", "ok",
                as.numeric(difftime(Sys.time(), t0, units = "secs")),
                result_status = result_status_of(res))
  }, error = function(e) {
    record_step(cfg, "virtual_perturbation", "failed",
                as.numeric(difftime(Sys.time(), t0, units = "secs")),
                message = conditionMessage(e), result_status = "failed")
    stop(e)
  })
}
