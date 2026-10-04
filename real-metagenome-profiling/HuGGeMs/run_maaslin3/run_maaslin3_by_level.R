# =============================================================================
#  run_maaslin3_by_level.R
#  按协变量层 (L1–L5) + 地域层 (LG) 批量运行 MaAsLin3
# =============================================================================
#  输入：GMB/maaslin3_v2/maaslin_inputs/_design.tsv
#        以及同目录下 <PID>__<LEVEL>__data.tsv / __meta.tsv
#        （由 01_make_maaslin_inputs.py 生成，已完成 complete-case 筛选）
#  输出：GMB/maaslin3_v2/Maaslin3-result-v2/<PID>/<LEVEL>/
#
#  设计要点：
#    1. 样本的缺失值问题已在 Python 端处理完毕（每个层都做了 complete-case
#       筛选并生成对齐的 data/meta），所以这里不需要再处理 NA。
#    2. 每个 队列×层 独立输出目录，避免断点续传误跳过。
#    3. 对照组统一设为 Health。
#    4. 兼容 maaslin3 不同版本的参数名（formula / fixed_effects）。
# =============================================================================

suppressPackageStartupMessages(library(maaslin3))

# -----------------------------------------------------------------------------
# 【重要】Windows 中文路径兼容
# -----------------------------------------------------------------------------
# 本脚本以 UTF-8 保存（与原始 RUN.R 一致）。
# 若在 RStudio 中运行：locale 通常已是中文，无需任何处理，直接可用。
# 若通过命令行 Rscript 运行且报 "unable to translate ... to UTF-8"，
# 说明 R 落在纯 C locale（MBCS=FALSE），此处会自动尝试启用多字节 locale。
#
# 注意：**优先使用 UTF-8 系的中文 locale**，否则 locale 与脚本编码不一致
#       （如脚本 UTF-8 而 locale 是 GBK），反而会报 "file name conversion problem"。
if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()$MBCS)) {
  fallback <- NULL
  for (loc in c("Chinese_China.utf8", "zh_CN.UTF-8", "Chinese_China.65001",
                "Chinese", "chs", "Chinese_China.936")) {
    ok <- tryCatch({
      suppressWarnings(Sys.setlocale("LC_CTYPE", loc)); TRUE
    }, warning = function(w) FALSE, error = function(e) FALSE)
    if (!ok || !isTRUE(l10n_info()$MBCS)) next
    if (is.null(fallback)) fallback <- loc
    if (grepl("utf", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE)) {
      fallback <- loc
      break
    }
  }
  if (!is.null(fallback) && is.null(attr(fallback, "done"))) {
    suppressWarnings(Sys.setlocale("LC_CTYPE", fallback))
  }
  cat("locale LC_CTYPE =", Sys.getlocale("LC_CTYPE"),
      "| MBCS =", l10n_info()$MBCS, "\n")
  if (!isTRUE(l10n_info()$MBCS)) {
    warning("未能启用多字节 locale，中文路径可能读取失败。\n",
            "建议：1) 直接在 RStudio 中运行本脚本；",
            "2) 或把数据目录复制/链接到不含中文的路径后修改上面的路径配置。")
  }
}

# -----------------------------------------------------------------------------
# 参数区
# -----------------------------------------------------------------------------
INPUT_DIR  <- "D:/Desktop/刘双江-实验室/硕士发表文章/各种材料/4.真实宏基因组样本信息/真实宏基因组样本分析/4.疾病健康队列分析/GMB/maaslin3_v2/maaslin_inputs"
OUTPUT_DIR <- "D:/Desktop/刘双江-实验室/硕士发表文章/各种材料/4.真实宏基因组样本信息/真实宏基因组样本分析/4.疾病健康队列分析/GMB/maaslin3_v2/Maaslin3-result-v2"

MAX_SIG    <- 0.1     # 显著性阈值。原 RUN.R 用 0.1；若要与已发表 261+237
                      # 的口径对齐，下游汇总时按 qval_joint < 0.05 过滤即可
                      # （all_results.tsv 会保留全部检验结果，不受此值影响）
MIN_PREV   <- 0.05    # 流行率过滤 = 5%（与文件名 prev5pct 一致）
CORES      <- 4       # 并行核数
SAVE_MODEL <- TRUE    # 是否保存模型对象（占用空间较大，可设 FALSE）

REFERENCE_LEVEL <- "Health"   # 对照组名称

# 只跑某些层？留空 = 全跑。例：ONLY_LEVELS <- c("L1","L5")
ONLY_LEVELS <- c()

# 只跑某些队列？留空 = 全跑。例：ONLY_PROJECTS <- c("PRJEB1220")
ONLY_PROJECTS <- c()

# -----------------------------------------------------------------------------
# 初始化
# -----------------------------------------------------------------------------
design_path <- file.path(INPUT_DIR, "_design.tsv")
if (!file.exists(design_path)) {
  stop("找不到设计表: ", design_path,
       "\n请先运行 01_make_maaslin_inputs.py")
}

design <- read.delim(design_path, sep = "\t", header = TRUE,
                     check.names = FALSE, stringsAsFactors = FALSE)
design <- design[design$status == "OK", , drop = FALSE]

if (length(ONLY_LEVELS) > 0) {
  design <- design[design$level %in% ONLY_LEVELS, , drop = FALSE]
}
if (length(ONLY_PROJECTS) > 0) {
  design <- design[design$project_id %in% ONLY_PROJECTS, , drop = FALSE]
}

cat("=======================================================\n")
cat("待运行任务数:", nrow(design), "\n")
cat("MaAsLin3 版本:", as.character(packageVersion("maaslin3")), "\n")
cat("=======================================================\n")

# --- 兼容不同版本的参数名 -----------------------------------------------------
m3_args <- names(formals(maaslin3))
FORMULA_ARG <- if ("formula" %in% m3_args) "formula" else "fixed_effects"
names(FORMULA_ARG) <- FORMULA_ARG
cat("使用公式参数名:", FORMULA_ARG, "\n\n")

log_lines <- character(0)
logit <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  log_lines <<- c(log_lines, msg)
}

# =============================================================================
# 主循环
# =============================================================================
for (i in seq_len(nrow(design))) {
  row <- design[i, ]
  pid   <- row$project_id
  lvl   <- row$level
  form  <- row$formula
  n_smp <- row$n_samples

  pre       <- row$file_prefix
  data_path <- file.path(INPUT_DIR, paste0(pre, "__data.tsv"))
  meta_path <- file.path(INPUT_DIR, paste0(pre, "__meta.tsv"))
  out_dir   <- file.path(OUTPUT_DIR, pid, lvl)
  done_flag <- file.path(out_dir, "all_results.tsv")

  logit(sprintf("[%d/%d] %s | %s | %s | n=%s",
                i, nrow(design), pid, lvl, form, n_smp))

  # --- 断点续传：每个 队列×层 独立目录，不会互相干扰 ---
  if (file.exists(done_flag)) {
    logit("        >>> [跳过] 已完成")
    next
  }

  if (!file.exists(data_path) || !file.exists(meta_path)) {
    logit("        >>> [警告] 输入文件缺失，跳过: ", data_path)
    next
  }

  tryCatch({
    # ---- 读 metadata，ID 作为行名 ----
    meta_df <- read.delim(meta_path, sep = "\t", header = TRUE,
                          row.names = 1, check.names = FALSE,
                          stringsAsFactors = FALSE)

    # ---- 强制数值列 ----
    for (nc in c("age", "BMI")) {
      if (nc %in% colnames(meta_df)) {
        meta_df[[nc]] <- suppressWarnings(as.numeric(meta_df[[nc]]))
      }
    }
    # ---- group / sex / country 转因子，Health 设为参照 ----
    meta_df$group <- factor(meta_df$group)
    if (REFERENCE_LEVEL %in% levels(meta_df$group)) {
      meta_df$group <- relevel(meta_df$group, ref = REFERENCE_LEVEL)
    } else {
      logit("        >>> [警告] 无 ", REFERENCE_LEVEL, " 参照，使用默认")
    }
    for (fc in c("sex", "country")) {
      if (fc %in% colnames(meta_df)) {
        meta_df[[fc]] <- factor(meta_df[[fc]])
      }
    }

    # ---- 检查模型矩阵是否可逆（防止某个协变量在子集后变成常数） ----
    rhs <- sub("^~\\s*", "", form)
    mf  <- model.frame(as.formula(form), data = meta_df, na.action = na.pass)
    mm  <- model.matrix(as.formula(form), data = meta_df)
    if (qr(mm)$rank < ncol(mm)) {
      logit("        >>> [跳过] 设计矩阵秩不足（协变量可能共线/常数）")
      next
    }

    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    # ---- 组装参数 ----
    args_list <- list(
      input_data     = data_path,
      input_metadata = meta_df,
      output         = out_dir,
      normalization  = "TSS",
      transform      = "LOG",
      standardize    = TRUE,
      augment        = TRUE,
      min_prevalence = MIN_PREV,
      max_significance = MAX_SIG,
      cores          = CORES,
      save_models    = SAVE_MODEL
    )
    args_list[[FORMULA_ARG]] <- form
    # 参照水平（格式 "变量,水平"）
    if ("reference" %in% m3_args) {
      args_list$reference <- paste0("group,", REFERENCE_LEVEL)
    }

    do.call(maaslin3, args_list)
    logit("        >>> [成功]")

  }, error = function(e) {
    logit("        >>> [错误] ", conditionMessage(e))
  })
}

# =============================================================================
# 写日志
# =============================================================================
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
log_path <- file.path(OUTPUT_DIR, paste0("run_log_",
                                         format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
writeLines(log_lines, log_path)

cat("\n=======================================================\n")
n_done <- length(list.files(OUTPUT_DIR, pattern = "^all_results.tsv$",
                            recursive = TRUE))
cat("全部任务处理完毕\n")
cat("已产生 all_results.tsv 的数量:", n_done, "\n")
cat("日志:", log_path, "\n")
cat("=======================================================\n")
