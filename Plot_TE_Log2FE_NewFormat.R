#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-


# ==============================================================================
# TE Enrichment: Grouped Class Barplot + Grouped Family Volcano Plot
# FIXED VERSION: Resolved data.table := assignment errors
# ==============================================================================


suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(stringr)
  library(ggrepel) 
})


# --- CONFIGURATION ---
INPUT_CSV <- "hg38.5.way.TE_Log2FE_Stats_0.4.csv" 
OUTPUT_PREFIX <- "hg38.5.way.TE_Analysis_Log2FE_Grouped_Fixed"


SIG_P_VAL <- 0.05
LOG_P_CUTOFF <- -log10(SIG_P_VAL) 
LOG2FC_CUTOFF <- 0 


# 🎨 颜色设置 (根据你的数据 Category 动态调整)
CAT_COLORS <- c(
  "5way"      = "#C05C5C", 
  "5way TB"   = "#E6A86C", 
  "5way HSTB" = "#7FB085" 
)


COL_ENRICHED  <- "Enriched"
COL_DEPLETED  <- "Depleted"
COL_NS        <- "Not_Sig"


# --- HELPER: Get Category Colors ---
get_category_colors <- function(cats, predefined_colors) {
  missing_cats <- setdiff(cats, names(predefined_colors))
  if (length(missing_cats) > 0) {
    cat("⚠️ 警告：以下 Category 未定义颜色，将使用默认灰色:", paste(missing_cats, collapse=", "), "\n")
    extra_colors <- setNames(rep("#999999", length(missing_cats)), missing_cats)
    return(c(predefined_colors, extra_colors))
  }
  return(predefined_colors)
}


# --- HELPER: Get Stars ---
get_stars <- function(p_val) {
  if (is.na(p_val)) return("")
  if (p_val < 0.0001) return("****")
  if (p_val < 0.001)  return("***")
  if (p_val < 0.01)   return("**")
  if (p_val < 0.05)   return("*")
  return("")
}


# --- MAIN LOGIC ---


cat("Loading data from:", INPUT_CSV, "\n")
if (!file.exists(INPUT_CSV)) stop("File not found.")


dt_raw <- fread(INPUT_CSV)
colnames(dt_raw) <- trimws(colnames(dt_raw))


# 检查新格式所需的列
req_cols <- c("Category", "Type", "Name", "Observed_Fraction", "Expected_Mean", "Log2_Fold_Enrichment", "P_Value")
missing_cols <- setdiff(req_cols, colnames(dt_raw))
if (length(missing_cols) > 0) {
  stop("Missing columns: ", paste(missing_cols, collapse=", "))
}


cat("Data loaded. Rows:", nrow(dt_raw), "\n")


# 动态获取颜色
all_cats <- unique(dt_raw$Category)
FINAL_CAT_COLORS <- get_category_colors(all_cats, CAT_COLORS)
CAT_ORDER <- names(FINAL_CAT_COLORS) 


# ==============================================================================
# 🔀 STEP 1: SPLIT BY TYPE
# ==============================================================================


dt_class_raw <- dt_raw[Type == "Class"]
dt_family_raw <- dt_raw[Type == "Family"]


cat("Split: Class=", nrow(dt_class_raw), ", Family=", nrow(dt_family_raw), "\n")


# ==============================================================================
# 🏭 PROCESS CLASS (Grouped Barplot) - FIXED SYNTAX
# ==============================================================================


process_class_bar <- function(dt_in) {
  if (is.null(dt_in) || nrow(dt_in) == 0) return(NULL)
  
  # 确保是 data.table
  dt <- as.data.table(dt_in)
  
  dt[, CleanName := as.character(Name)]
  
  # 1. 定义显著性状态 (使用标准 data.table 语法)
  dt[, Sig_Status := COL_NS]
  dt[!is.na(P_Value) & P_Value < SIG_P_VAL & Log2_Fold_Enrichment > 0, Sig_Status := COL_ENRICHED]
  dt[!is.na(P_Value) & P_Value < SIG_P_VAL & Log2_Fold_Enrichment <= 0, Sig_Status := COL_DEPLETED]
  
  # 2. 填充缺失的 Category/Name 组合
  all_names <- unique(dt[, Name])
  # 使用 CJ 创建模板
  template <- CJ(Category = CAT_ORDER, Name = all_names, Type = "Class", unique = TRUE)
  dt_complete <- merge(template, dt, by = c("Category", "Name", "Type"), all.x = TRUE)
  
  # 确保合并后也是 data.table
  setDT(dt_complete)
  
  # 3. 填充 NA (修复后的语法：每个操作独立一行)
  dt_complete[is.na(Sig_Status), Sig_Status := COL_NS]
  dt_complete[is.na(Observed_Fraction), Observed_Fraction := 0]
  dt_complete[is.na(Log2_Fold_Enrichment), Log2_Fold_Enrichment := 0]
  dt_complete[is.na(P_Value), P_Value := 1]
  
  # 4. 排序 Name
  scores <- dt_complete[, .(max_abs_fc = max(abs(Log2_Fold_Enrichment), na.rm=TRUE),
                            is_sig = any(Sig_Status != COL_NS)), by = Name]
  scores <- scores[order(-is_sig, -max_abs_fc, Name)]
  
  dt_complete[, Name := factor(Name, levels = scores[, Name])]
  dt_complete[, Category := factor(Category, levels = CAT_ORDER)]
  
  # 5. 添加星号
  dt_complete[, Star := sapply(P_Value, get_stars)]
  dt_complete[is.na(Star), Star := ""]
  
  return(dt_complete)
}


# ==============================================================================
# 🌋 PROCESS FAMILY (Volcano Data) - FIXED SYNTAX
# ==============================================================================


process_family_volcano <- function(dt_in) {
  if (is.null(dt_in) || nrow(dt_in) == 0) return(NULL)
  
  dt <- as.data.table(dt_in)
  dt[, CleanName := as.character(Name)]
  
  # 2. 计算 -log10(P)
  dt[, NegLog10P := -log10(P_Value + 1e-300)]
  dt[is.na(P_Value), NegLog10P := 0]
  
  # 3. 定义状态
  dt[, Sig_Status := COL_NS]
  dt[!is.na(P_Value) & P_Value < SIG_P_VAL & Log2_Fold_Enrichment > 0, Sig_Status := COL_ENRICHED]
  dt[!is.na(P_Value) & P_Value < SIG_P_VAL & Log2_Fold_Enrichment <= 0, Sig_Status := COL_DEPLETED]
  
  # 4. 🔥 核心过滤：只保留显著富集
  dt <- dt[Sig_Status == COL_ENRICHED]
  
  if (nrow(dt) == 0) {
    cat("Warning: No significantly enriched families found.\n")
    return(NULL)
  }
  
  # 5. 标注策略：每个 Category 标注最显著的前 3 个
  dt[, Label := ""]
  
  # 使用 data.table 的分组更新语法
  # 先找出需要标注的行索引
  top_hits <- dt[, head(.SD[order(P_Value)], 3), by = Category]
  
  # 将标签合并回主表
  if (nrow(top_hits) > 0) {
    # 创建映射表
    label_map <- top_hits[, .(Category, Name, Label = CleanName)]
    # 更新主表
    dt[label_map, on = .(Category, Name), Label := i.Label]
  }
  
  # 排序
  dt <- dt[order(Category, -NegLog10P, -Log2_Fold_Enrichment)]
  
  return(dt)
}


# Execute Processing
dt_class_final <- process_class_bar(dt_class_raw)
dt_family_final <- process_family_volcano(dt_family_raw)


# ==============================================================================
# 🎨 PLOTTING
# ==============================================================================


# --- Plot 1: Class Grouped Barplot ---
make_class_plot <- function(dt_data) {
  if (is.null(dt_data) || nrow(dt_data) == 0) return(NULL)
  
  # 创建填充颜色列
  dt_data[, FillColor := NA_character_]
  
  for (cat in CAT_ORDER) {
    dt_data[Category == cat & Sig_Status == COL_ENRICHED, FillColor := FINAL_CAT_COLORS[cat]]
    dt_data[Category == cat & Sig_Status == COL_DEPLETED, FillColor := "#555555"]
    dt_data[Category == cat & Sig_Status == COL_NS, FillColor := "#DDDDDD"]
  }
  
  max_obs <- max(dt_data[, Observed_Fraction], na.rm = TRUE)
  if (is.infinite(max_obs) || is.na(max_obs)) max_obs <- 0.2
  offset_val <- max_obs * 0.05
  
  dt_sig <- dt_data[Star != "", ]
  
  p <- ggplot(dt_data, aes(x = .data[["Name"]], y = .data[["Observed_Fraction"]], fill = .data[["FillColor"]])) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, na.rm = TRUE) +
    geom_text(data = dt_sig, aes(y = .data[["Observed_Fraction"]] + offset_val, label = .data[["Star"]]), 
              position = position_dodge(width = 0.8), vjust = 0, hjust = 0.5, 
              size = 3.5, fontface = "bold", color = "black", na.rm = TRUE) +
    scale_fill_identity(name = "Status", guide = "none") +
    labs(title = "TE Class Composition & Enrichment", 
         x = "TE Class", 
         y = "Observed Fraction") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", color = "#2C3E50"),
          axis.title = element_text(face = "bold", color = "#2C3E50"),
          legend.position = "right")
  
  return(p)
}


# --- Plot 2: Family Volcano Plot (Log2FC X-axis) ---
make_volcano_plot <- function(dt_data) {
  if (is.null(dt_data) || nrow(dt_data) == 0) return(NULL)
  
  p <- ggplot(dt_data, aes(x = .data[["Log2_Fold_Enrichment"]], y = .data[["NegLog10P"]], color = .data[["Category"]])) +
    geom_hline(yintercept = LOG_P_CUTOFF, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_vline(xintercept = LOG2FC_CUTOFF, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "solid", color = "gray80", linewidth = 0.5) +
    
    geom_point(size = 3.5, alpha = 0.9) +
    
    geom_text_repel(data = dt_data[Label != ""], 
                    aes(label = .data[["Label"]]), 
                    size = 4.5, 
                    box.padding = 0.6, 
                    point.padding = 0.6, 
                    segment.color = "grey50",
                    max.overlaps = 20,
                    fontface = "bold") +
    
    scale_color_manual(values = FINAL_CAT_COLORS, name = "Category") +
    labs(title = "Significantly Enriched TE Families",
         subtitle = sprintf("Filtered: P < %.2f & Log2FC > 0", SIG_P_VAL),
         x = expression(log[2]~"(Fold Enrichment)"),
         y = expression(-log[10]~italic(p)-value)) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", color = "#2C3E50"),
      plot.subtitle = element_text(hjust = 0.5, color = "gray50", size = 10),
      axis.title = element_text(face = "bold", color = "#2C3E50"),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    ) +
    coord_cartesian(xlim = c(-0.5, max(dt_data[, Log2_Fold_Enrichment], na.rm = TRUE) * 1.1)) 
  
  return(p)
}


cat("Generating plots...\n")
p_class <- make_class_plot(dt_class_final)
p_family <- make_volcano_plot(dt_family_final)


# Save
save_plot <- function(p, suffix) {
  if (is.null(p)) return()
  ggsave(paste0(OUTPUT_PREFIX, suffix, ".pdf"), p, width = 12, height = 8, dpi = 300)
  ggsave(paste0(OUTPUT_PREFIX, suffix, ".png"), p, width = 12, height = 8, dpi = 300)
  cat("Saved:", paste0(OUTPUT_PREFIX, suffix), "\n")
}


save_plot(p_class, "_Class_Bar")
save_plot(p_family, "_Family_Volcano")


# Combined Plot
if (!is.null(p_class) && !is.null(p_family)) {
  combined <- grid.arrange(p_class, p_family, ncol = 2, widths = c(5, 6))
  ggsave(paste0(OUTPUT_PREFIX, "_Combined.pdf"), combined, width = 18, height = 8, dpi = 300)
  cat("Saved Combined.\n")
}


cat("All done successfully!\n")