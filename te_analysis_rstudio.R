#!/usr/bin/env Rscript
# ==============================================================================
# Stratified TE Analysis: CSV Output Only (Log2 Fold Enrichment)
# 功能：执行分层置换检验，计算 Log2FE 和 P-value，仅保存 CSV 结果，不绘图。
# ==============================================================================


suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(dplyr)
})


# ==============================================================================
# ⚙️ 1. 配置区域 (Configuration)
# ==============================================================================


# 🟢 文件路径
SUPERPOOL_FILE <- "superpool_25bp_example.bed"
REPEAT_FILE    <- "rmsk.sort.example.bed"


# 🟢 分析参数
OVERLAP_THRESHOLD <- 0.4   # 重叠比例阈值
N_PERMUTATIONS    <- 10000  # 置换次数 (建议 >= 1000)
RANDOM_SEED       <- 42    # 随机种子


# 🟢 目标类别 (只有这 4 类保留原名，其他归为 Other_Families)
TARGET_CLASSES <- c("DNA", "LINE", "SINE", "LTR")


# 🟢 分层列名 (必须与 Superpool 文件表头一致，支持大小写自动匹配)
STRATA_COLS <- c("GC_Bin", "Promoter_Bin", "Mappability_Bin", "Distance_Bin")


# ==============================================================================
# 📦 2. 环境检查与辅助函数
# ==============================================================================


install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    cat("📦 正在安装缺失的包:", paste(new_packages, collapse=", "), "\n")
    install.packages(new_packages, repos = "https://cloud.r-project.org")
    if(any(c("GenomicRanges", "IRanges") %in% new_packages)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(c("GenomicRanges", "IRanges"), update = FALSE, ask = FALSE)
    }
  }
  sapply(packages, library, character.only = TRUE, verbose = FALSE)
}


cat("🔧 检查环境依赖...\n")
required_packages <- c("data.table", "GenomicRanges", "IRanges", "dplyr")
install_if_missing(required_packages)
cat("✅ 环境准备就绪。\n\n")


# 计算 TE 组成比例
get_composition <- function(query_gr, repeat_gr, thresh) {
  hits <- findOverlaps(query_gr, repeat_gr, ignore.strand = TRUE)
  if (length(hits) == 0) return(list(class_prop = numeric(0), family_prop = numeric(0)))
  
  q_idx <- queryHits(hits)
  r_idx <- subjectHits(hits)
  q_w <- width(query_gr)[q_idx]
  ov_w <- width(pintersect(query_gr[q_idx], repeat_gr[r_idx]))
  ratios <- ov_w / q_w
  
  dt <- data.table(
    q_idx = q_idx, 
    ratio = ratios, 
    cls = mcols(repeat_gr)$RepClass_Clean[r_idx],
    fam = mcols(repeat_gr)$RepFamily_Clean[r_idx]
  )
  
  dt_filt <- dt[ratio >= thresh]
  if (nrow(dt_filt) == 0) return(list(class_prop = numeric(0), family_prop = numeric(0)))
  
  # 每个区域只保留重叠度最高的那个 TE
  dt_best <- dt_filt[, .SD[which.max(ratio)], by = q_idx]
  
  total <- length(query_gr)
  if (total == 0) return(list(class_prop = numeric(0), family_prop = numeric(0)))
  
  c_tab <- table(dt_best$cls)
  c_prop <- as.numeric(c_tab) / total
  names(c_prop) <- names(c_tab)
  
  f_tab <- table(dt_best$fam)
  f_prop <- as.numeric(f_tab) / total
  names(f_prop) <- names(f_tab)
  
  return(list(class_prop = c_prop, family_prop = f_prop))
}


# ==============================================================================
# 🚀 3. 主程序逻辑
# ==============================================================================


cat("🚀 开始执行 Stratified TE Analysis (CSV Only Mode)...\n")


# --- 检查文件 ---
if (!file.exists(SUPERPOOL_FILE)) stop(sprintf("❌ 错误：找不到背景文件 '%s'", SUPERPOOL_FILE))
if (!file.exists(REPEAT_FILE)) stop(sprintf("❌ 错误：找不到注释文件 '%s'", REPEAT_FILE))


# --- 查找目标文件 ---
input_files <- list.files(path = ".", pattern = "test\\.bed$", full.names = TRUE)
if (length(input_files) == 0) {
  all_beds <- list.files(path = ".", pattern = "\\.bed$", full.names = TRUE)
  input_files <- all_beds[!grepl("superpool|rmsk", basename(all_beds), ignore.case = TRUE)]
}
if (length(input_files) == 0) stop("❌ 错误：未找到任何目标文件 (*.robust.bed 或其他 .bed)。")


cat(sprintf("📂 找到 %d 个待处理文件。\n", length(input_files)))


# --- 加载背景数据 (一次性加载) ---
cat("\n📊 正在加载背景数据 (Superpool & Repeats)...\n")
tryCatch({
  # 1. 加载 Superpool
  sp_dt <- fread(SUPERPOOL_FILE, header = TRUE, sep = "\t", check.names = FALSE)
  
  # 列名标准化
  if ("chr" %in% names(sp_dt)) setnames(sp_dt, "chr", "Chr")
  if ("start" %in% names(sp_dt)) setnames(sp_dt, "start", "Start")
  if ("end" %in% names(sp_dt)) setnames(sp_dt, "end", "End")
  
  base_cols <- c("Chr", "Start", "End")
  if (!all(base_cols %in% names(sp_dt))) {
    stop(sprintf("Superpool 文件缺少基础列 (%s)。", paste(base_cols, collapse=", ")))
  }
  
  # 验证并修复分层列大小写
  missing_strata <- setdiff(STRATA_COLS, names(sp_dt))
  if (length(missing_strata) > 0) {
    lower_file <- tolower(names(sp_dt))
    lower_target <- tolower(STRATA_COLS)
    matches <- match(lower_target, lower_file)
    
    if (all(!is.na(matches))) {
      cat("   ⚠️ 检测到分层列大小写不一致，正在自动修正...\n")
      setnames(sp_dt, old = names(sp_dt)[matches], new = STRATA_COLS)
    } else {
      stop(sprintf("❌ Superpool 文件缺少分层列：%s\n当前可用列：%s", 
                   paste(missing_strata, collapse=", "), paste(names(sp_dt), collapse=", ")))
    }
  }
  
  # 创建分层 Key
  sp_dt[, strata_key := do.call(paste, c(.SD, sep = "_")), .SDcols = STRATA_COLS]
  sp_dt[, strata_key := as.factor(strata_key)]
  
  # 创建 GRanges
  gr_superpool <- makeGRangesFromDataFrame(sp_dt, 
                                           seqnames.field = "Chr", 
                                           start.field = "Start", 
                                           end.field = "End",
                                           keep.extra.columns = TRUE)
  n_sp <- length(gr_superpool)
  
  # 2. 加载 RepeatMasker
  rep_dt <- fread(REPEAT_FILE, header = FALSE, sep = "\t",
                  col.names = c("Chr", "Start", "End", "Strand", "Name", "RepClass", "RepFamily"))
  
  # 列名标准化
  if ("chr" %in% names(rep_dt)) setnames(rep_dt, c("chr","start","end","repclass","repfamily"), 
                                         c("Chr","Start","End","RepClass","RepFamily"))
  
  # --- 核心逻辑：严格的家族合并 ---
  cat("   正在应用家族过滤规则...\n")
  
  is_valid_class <- rep_dt$RepClass %in% TARGET_CLASSES
  is_uncertain <- grepl("\\?", rep_dt$RepFamily)
  
  rep_dt[, RepFamily_Clean := ifelse(
    !is_valid_class | is_uncertain, 
    "Other_Families", 
    RepFamily
  )]
  
  rep_dt[, RepClass_Clean := ifelse(
    RepClass %in% TARGET_CLASSES, 
    RepClass, 
    "TE_other"
  )]
  
  cat(sprintf("   ✅ 背景加载完成：%d 个 Superpool 区域，%d 个 Repeat 条目。\n", n_sp, nrow(rep_dt)))
  
  gr_repeats <- makeGRangesFromDataFrame(rep_dt, 
                                         seqnames.field = "Chr", 
                                         start.field = "Start", 
                                         end.field = "End",
                                         keep.extra.columns = TRUE)
  
}, error = function(e) {
  stop(sprintf("💥 加载背景数据失败：%s", e$message))
})


# --- 单个文件处理函数 ---
run_stats_only <- function(ctcf_file, gr_sp, gr_rep, n_sp, n_perms, thresh, seed) {
  set.seed(seed)
  
  base_name <- basename(ctcf_file)
  prefix <- tools::file_path_sans_ext(base_name)
  
  cat(sprintf("\n>> 正在处理：%s\n", base_name))
  
  # 1. 加载目标文件
  tryCatch({
    ctcf_dt <- fread(ctcf_file, header = FALSE, select = c(1, 2, 3))
    colnames(ctcf_dt) <- c("Chr", "Start", "End")
    
    if (is.character(ctcf_dt$Start[1])) {
      ctcf_dt <- fread(ctcf_file, header = TRUE, select = c(1, 2, 3))
      names(ctcf_dt)[1:3] <- c("Chr", "Start", "End")
    }
    
    if ("chr" %in% names(ctcf_dt)) setnames(ctcf_dt, c("chr","start","end"), c("Chr","Start","End"))
    
    gr_ctcf <- makeGRangesFromDataFrame(ctcf_dt, 
                                        seqnames.field = "Chr", 
                                        start.field = "Start", 
                                        end.field = "End",
                                        keep.extra.columns = FALSE)
    n_ctcf <- length(gr_ctcf)
    
    if (n_ctcf > n_sp) {
      cat("   ⚠️ 警告：目标区域数 > 背景池区域数。跳过置换检验。\n")
      return(FALSE) 
    }
  }, error = function(e) {
    cat(sprintf("   ❌ 读取文件失败：%s\n", e$message))
    return(FALSE)
  })
  
  # 2. 映射目标区域的分层信息
  hits_map <- findOverlaps(gr_ctcf, gr_sp, type = "any", ignore.strand = TRUE)
  
  if (length(hits_map) == 0) {
    cat("   ❌ 错误：目标区域与背景池没有任何重叠！请检查染色体命名。\n")
    return(FALSE)
  }
  
  q_idx <- as.integer(queryHits(hits_map))
  s_idx <- as.integer(subjectHits(hits_map))
  sp_keys <- mcols(gr_sp)$strata_key
  
  map_df <- data.frame(q = q_idx, s = s_idx, stringsAsFactors = FALSE)
  map_df <- map_df[!duplicated(map_df$q), ]
  
  target_keys <- rep(NA_character_, length(gr_ctcf))
  target_keys[map_df$q] <- as.character(sp_keys[map_df$s])
  
  na_count <- sum(is.na(target_keys))
  if (na_count > 0) {
    valid_idx <- which(!is.na(target_keys))
    gr_ctcf <- gr_ctcf[valid_idx]
    target_keys <- target_keys[valid_idx]
    
    if (length(gr_ctcf) == 0) {
      cat("   ❌ 错误：所有目标区域都被过滤掉了。\n")
      return(FALSE)
    }
  }
  
  mcols(gr_ctcf)$strata_key <- factor(target_keys)
  
  # 3. 计算观测值 (Observed)
  obs <- get_composition(gr_ctcf, gr_rep, thresh)
  obs_c <- obs$class_prop
  obs_f <- obs$family_prop
  
  # 4. 分层置换检验
  cat(sprintf("   正在运行 %d 次分层置换...\n", n_perms))
  
  sp_strata_groups <- split(1:n_sp, mcols(gr_sp)$strata_key)
  target_strata_counts <- table(mcols(gr_ctcf)$strata_key)
  valid_strata <- intersect(names(target_strata_counts), names(sp_strata_groups))
  
  if (length(valid_strata) == 0) {
    cat("   ❌ 错误：目标文件和背景池没有共同的分层。\n")
    return(FALSE)
  }
  
  null_c_list <- vector("list", n_perms)
  null_f_list <- vector("list", n_perms)
  
  for (i in 1:n_perms) {
    if (i %% 200 == 0) cat(sprintf("     进度：%d/%d\n", i, n_perms))
    
    sampled_indices <- c()
    for (st in valid_strata) {
      n_needed <- target_strata_counts[st]
      candidates <- sp_strata_groups[[st]]
      
      if (length(candidates) < n_needed) {
        sampled_indices <- c(sampled_indices, candidates)
      } else {
        sampled_indices <- c(sampled_indices, sample(candidates, n_needed, replace = FALSE))
      }
    }
    
    gr_samp <- gr_sp[sampled_indices]
    res <- get_composition(gr_samp, gr_rep, thresh)
    null_c_list[[i]] <- res$class_prop
    null_f_list[[i]] <- res$family_prop
  }
  
  # 5. 汇总统计量 (计算 Mean, Log2FE, Empirical P-value)
  agg_stats <- function(lst, obs_vec, n_perms) {
    all_names <- unique(c(names(obs_vec), unlist(lapply(lst, names))))
    if (length(all_names) == 0) return(data.frame())
    
    mat <- matrix(0, nrow = n_perms, ncol = length(all_names))
    colnames(mat) <- all_names
    
    for (i in 1:n_perms) {
      v <- lst[[i]]
      if (length(v) > 0) {
        common <- intersect(names(v), all_names)
        if (length(common) > 0) mat[i, common] <- v[common]
      }
    }
    
    means <- colMeans(mat)
    means[means == 0] <- 1e-9 
    
    # 计算 Empirical P-value
    p_vals <- sapply(names(obs_vec), function(nm) {
      if (!nm %in% colnames(mat)) return(1.0)
      obs_val <- obs_vec[nm]
      null_vals <- mat[, nm]
      (sum(null_vals >= obs_val) + 1) / (n_perms + 1)
    })
    
    # 计算 Fold Enrichment (Linear)
    fe_vals <- sapply(names(obs_vec), function(nm) {
      if (!nm %in% names(means)) return(NA)
      obs_vec[nm] / means[nm]
    })
    
    # 计算 Log2 Fold Enrichment
    log2_fe_vals <- log2(fe_vals)
    log2_fe_vals[is.infinite(log2_fe_vals) | is.na(log2_fe_vals)] <- NA 
    
    df <- data.frame(
      Name = names(obs_vec),
      Observed_Fraction = as.numeric(obs_vec),
      Expected_Mean = means[names(obs_vec)],
      Fold_Enrichment_Linear = as.numeric(fe_vals),
      Log2_Fold_Enrichment = log2_fe_vals,
      P_Value = as.numeric(p_vals),
      stringsAsFactors = FALSE
    )
    return(df)
  }
  
  res_c <- agg_stats(null_c_list, obs_c, n_perms)
  res_f <- agg_stats(null_f_list, obs_f, n_perms)
  
  if (nrow(res_c) > 0) res_c$Type <- "Class"
  if (nrow(res_f) > 0) res_f$Type <- "Family"
  
  # 6. 保存 CSV
  all_res <- rbind(res_c, res_f)
  
  # 确保列顺序一致
  col_order <- c("Type", "Name", "Observed_Fraction", "Expected_Mean", "Fold_Enrichment_Linear", "Log2_Fold_Enrichment", "P_Value")
  all_res <- all_res[, col_order]
  
  csv_out <- sprintf("%s_TE_Log2FE_Stats_0.4.csv", prefix)
  fwrite(all_res, csv_out) # 使用 fwrite 速度更快
  
  cat(sprintf("   💾 结果已保存：%s\n", csv_out))
  
  gc()
  return(TRUE)
}


# --- 批量执行 ---
success_count <- 0
for (f in input_files) {
  if (run_stats_only(f, gr_superpool, gr_repeats, n_sp, N_PERMUTATIONS, OVERLAP_THRESHOLD, RANDOM_SEED)) {
    success_count <- success_count + 1
  }
}


cat("\n==================================================\n")
cat(sprintf("🎉 全部完成！成功处理 %d/%d 个文件。\n", success_count, length(input_files)))
cat("输出文件说明:\n")
cat("   - *_TE_Log2FE_Stats_0.4.csv : 包含 Class 和 Family 的完整统计结果 (Log2FE, P-value 等)\n")