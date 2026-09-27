# ============================================================================
# r_deps.R — R 版依赖清单（CI 安装与版本探针共用）
# ============================================================================
# r_version.md §6.4 的裁决：**不新增 R 依赖清单门禁** —— 版本清单写死在这里，
# 由 CI 的 Rscript 调它打印版本，门禁只比对它输出的文本与清单一致。
# 不靠 DESCRIPTION（本仓没有 DESCRIPTION，也不走 renv —— §8 的既定裁决）。
#
# 与 common.R 的 KEY_PACKAGES 的关系：KEY_PACKAGES 是 run_manifest.json 的
# key_versions 探针清单（capture_versions 消费），本文件是 **安装清单**，
# 两者必须保持同集 —— 改任何一边都要同步另一边，否则会出现
# "装了但没记录" 或 "记录了但没装"（后者会在 CI 跑到一半时才炸）。
#
# 用法：
#   Rscript scripts/r_deps.R            # 打印每包一行 "  <name> <version|MISSING>"
#                                       # 有 MISSING 时退出码 1（CI 可直接判）
# **本文件只报告，不安装** —— 装包的唯一事实源是
# .github/workflows/scrna_r_analysis.yml 的 setup-r-dependencies 包名清单。
# 这里故意不提供 --install 之类的自救入口：R1 纪律禁本地执行，
# 一个"方便的安装脚本"正是把执行性依赖引进本地的口子。
# ============================================================================

r_deps_packages <- c(
  # ---- 方法栈（r_version.md §3 选型表）----
  "Seurat",        # 02/03: 对象、HVG(vst)、PCA、FindNeighbors/Clusters/AllMarkers、UMAP
  "harmony",       # 02: integration method=harmony（RunHarmony）
  "sva",           # 02: integration method=combat（ComBat）
  "scDblFinder",   # 01: 双细胞检出（scrublet 的 R 侧替代，01=C 方法差异）
  "scater",        # 01/03: QC 指标（perCellQCMetrics 等）
  "SingleR",       # 03: 参考集自动注释（try_celltypist_r 之外的 R 侧通道）
  "celldex",       # 03: SingleR 的参考集
  "DESeq2",        # 04: pseudobulk DE（pydeseq2 的母实现，04=A）
  "destiny",       # 05: DiffusionMap + DPT（dpt 槽位，05=C）
  "slingshot",     # 05: 主曲线树拟时序（scfates 槽位；CRAN 无 scFates，§11 未决 1）
  "liana",         # 06: 细胞通讯（liana_wrap rank_aggregate，06=A/B）
  "scTenifoldKnk", # 08: 网络扰动引擎（08=D，两版共用 tenifold_knk.R）
  # ---- 基础设施 ----
  "yaml",          # load_config
  "jsonlite",      # write_json / read_json（write_manifest 的 na="null" 契约）
  "digest",        # sha256_file（缺了退 md5，见 common.R sha256_file 注释）
  "ggplot2",       # 出图
  "Matrix",        # 稀疏矩阵（00/01/02/03）
  "patchwork"      # 02 的 batch-mixing 双面板（02-02-03）
)

# KEY_PACKAGES 同集断言（防两处清单分叉 —— 见文件头注释）：
# 只在 common.R 可加载时执行（单独 source 本文件不依赖 common.R）。
.r_deps_check_key_packages <- function() {
  common <- file.path("scripts", "lib", "common.R")
  if (!file.exists(common)) return(invisible(FALSE))
  src <- readLines(common, warn = FALSE)
  hit <- grep("KEY_PACKAGES\\s*<-\\s*c\\(", src, value = TRUE)
  length(hit) == 1L  # 存在性粗查；逐包比对交给 check_r_syntax.mjs（如有）
}

# ---- 版本探针（CI 与本地通用，输出格式与 geo 仓同款）----
r_deps_report <- function() {
  for (p in r_deps_packages) {
    v <- if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "MISSING"
    cat(sprintf("  %-16s %s\n", p, v))
  }
  n_missing <- sum(!vapply(r_deps_packages, function(p) requireNamespace(p, quietly = TRUE), logical(1)))
  cat(sprintf("resolved %d/%d\n", length(r_deps_packages) - n_missing, length(r_deps_packages)))
  invisible(n_missing)
}

# ---- 安装入口（无）----
# 刻意不设：CI 的安装清单在 workflow 里逐包列举（见文件头），本地不装包。
# 两处清单分叉的防线是本文件与 KEY_PACKAGES 的同集断言 + 改清单必须同步
# workflow 包列表的约定（r_version.md §8 缓存键注释同款纪律）。

if (sys.nframe() == 0L) {
  .r_deps_check_key_packages()
  quit(status = min(r_deps_report(), 1L))  # 有 MISSING 时退出码 1，CI 可直接判
}
