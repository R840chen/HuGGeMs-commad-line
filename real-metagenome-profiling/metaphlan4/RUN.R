# ==============================================================================
# 批量 MaAsLin 3 分析脚本 (最终精简版)
# 功能：只分析 Group | 自动匹配 | 支持断点续传
# ==============================================================================

library(maaslin3)

# ------------------------------------------------------------------------------
# 1. 路径设置 (已更新为你指定的路径)
# ------------------------------------------------------------------------------

# Metadata 所在文件夹
meta_dir <- "C:\\Users\\cc\\Desktop\\文章\\各种材料\\真实宏基因组样本分析\\4.疾病健康队列分析\\MPA4\\metadata的信息\\1\\1"

# Taxonomy (相对丰度) 所在文件夹
tax_dir <- "C:\\Users\\cc\\Desktop\\文章\\各种材料\\真实宏基因组样本分析\\4.疾病健康队列分析\\MPA4\\相对丰度matrix\\1"

# 输出总文件夹
output_base_dir <- "C:\\Users\\cc\\Desktop\\文章\\各种材料\\真实宏基因组样本分析\\4.疾病健康队列分析\\MPA4\\Maaslin3-result"

# ------------------------------------------------------------------------------
# 【注意】文件名后缀设置
# ------------------------------------------------------------------------------
# 仍然沿用你之前的长后缀来提取 Project ID。
# 如果你的新文件夹里文件名变短了(比如变成了 PRJ123.tsv)，请把下面这行改为: tax_suffix <- ".tsv"
tax_suffix <- "_species_abundance_matrix.prev5pct.maaslin_ready.tsv"

# ------------------------------------------------------------------------------
# 2. 获取文件列表
# ------------------------------------------------------------------------------

tax_files <- list.files(tax_dir, pattern = "\\.tsv$", full.names = FALSE)
print(paste("检测到任务总数:", length(tax_files)))

dir.create(output_base_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 3. 批量循环
# ------------------------------------------------------------------------------

for (tax_file in tax_files) {
  
  # --- A. 基础信息提取 ---
  project_id <- gsub(tax_suffix, "", tax_file, fixed = TRUE)
  meta_file <- paste0(project_id, ".tsv")
  
  full_path_tax <- file.path(tax_dir, tax_file)
  full_path_meta <- file.path(meta_dir, meta_file)
  project_output_dir <- file.path(output_base_dir, paste0("Result_", project_id))
  
  # --- 断点续传检测 ---
  # 如果输出文件夹里已经有了 all_results.tsv，说明上次跑完了，直接跳过
  if (dir.exists(project_output_dir) && file.exists(file.path(project_output_dir, "all_results.tsv"))) {
    print(paste(">>> [跳过] 项目已存在且完成:", project_id))
    next 
  }
  
  # 检查 Metadata 是否存在
  if (!file.exists(full_path_meta)) {
    warning(paste("[警告] 找不到对应的 metadata，跳过:", meta_file))
    next
  }
  
  print(paste("--------------------------------------------------"))
  print(paste("正在处理:", project_id))
  
  # --- B. 智能处理 Metadata ---
  tryCatch({
    # 读取数据
    meta_df <- read.table(full_path_meta, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
    
    # 1. 规范化 group 列名
    cols <- colnames(meta_df)
    grp_idx <- grep("^group$", cols, ignore.case = TRUE)
    if(length(grp_idx) > 0) colnames(meta_df)[grp_idx[1]] <- "group"
    
    # 2. 设定公式 (已修改：只看 group)
    current_formula <- "~ group"
    print(paste("应用公式:", current_formula))
    
    # 3. 确定对照组 (Reference)
    # 依然需要确定谁是"健康人"，否则无法计算差异
    unique_groups <- unique(meta_df$group)
    ref_label <- NULL
    possible_controls <- c("Health", "Healthy", "Control", "H", "Normal", "control")
    
    for (ctrl in possible_controls) {
      if (ctrl %in% unique_groups) {
        ref_label <- paste0("group,", ctrl)
        break
      }
    }
    
    if (is.null(ref_label)) {
      warning(paste("未找到标准对照组名称，使用默认基准:", unique_groups[1]))
    } else {
      print(paste("设定对照组:", ref_label))
    }
    
    # --- C. 运行 MaAsLin 3 ---
    
    # 准备参数列表
    args_list <- list(
      input_data = full_path_tax,
      input_metadata = meta_df,
      output = project_output_dir,
      formula = current_formula, # 这里只包含 group
      normalization = "TSS",
      transform = "LOG",
      standardize = TRUE,
      augment = TRUE,
      min_prevalence = 0.05,
      max_significance = 0.1,
      cores = 4,
      save_models = TRUE
    )
    
    # 如果找到了对照组，加入参数；没找到则不加(使用默认)
    if (!is.null(ref_label)) {
      args_list$reference <- ref_label
    }
    
    # 调用函数
    do.call(maaslin3, args_list)
    
    print(paste("[成功] 项目完成:", project_id))
    
  }, error = function(e) {
    print(paste("[错误] 项目失败:", project_id))
    print(e$message)
  })
}

print("==================================================")
print("所有任务已更新完毕。")