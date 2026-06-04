# ==============================================================================
# SCRIPT: TE Enrichment Robustness Plotter (No Legend)
# DESCRIPTION:
#   Reads a single merged CSV file and plots TE enrichment robustness.
#   Automatically installs missing packages.
#   NO LEGEND is displayed (colors/shapes are defined in config for consistency).
#
# KEY FEATURES:
#   1. Auto-Installation of packages (ggplot2, dplyr, stringr, scales).
#   2. Configurable Category & Subgroup mappings at the top.
#   3. Enhanced Significance Stars: **** (<0.0001), *** (<0.001), ** (<0.01), * (<0.05).
#   4. No Background Line (Superpool_Mean ignored).
#   5. NO LEGEND displayed on the plot.
# ==============================================================================




# ================= 1. PACKAGE MANAGEMENT =================
required_packages <- c("ggplot2", "dplyr", "stringr", "scales")


for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("[INFO] Installing missing package:", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
  library(pkg, character.only = TRUE)
}




# ================= 2. CONFIGURATION (EDIT HERE) =================


# --- File Paths ---
INPUT_FILE <- "combined_superpool.csv"
OUTPUT_IMG <- "Figure_SX_Robustness_NoLegend_v5.png"
OUTPUT_PDF <- "Figure_SX_Robustness_NoLegend_v5.pdf"


# --- Category & Color Mapping ---
# Define colors for each category found in your CSV.
colors_map <- c(
  "hs"        = "#D62728",   # Human-Specific (All)
  "hs TB"     = "#FF7F0E",   # Human-Specific (TB)
  "hs HSTB"   = "#FFBB78",   # Human-Specific (HSTB)
  "5way"      = "#1F77B4",   # Conserved (All)
  "5way TB"   = "#52A1D6",   # Conserved (TB)
  "5way HSTB" = "#9EC2E8"    # Conserved (HSTB)
)


# --- Subgroup Visual Styles ---
# Maps subgroup types (detected from category name) to line styles and shapes.
linestyles_map <- c("All" = "solid", "TB" = "dashed", "HSTB" = "dotdash")
shapes_map     <- c("All" = 16, "TB" = 15, "HSTB" = 17) # 16=Circle, 15=Square, 17=Triangle


# --- Significance Thresholds ---
SIG_STARS_4 <- 0.0001  # ****
SIG_STARS_3 <- 0.001   # ***
SIG_STARS_2 <- 0.01    # **
SIG_STARS_1 <- 0.05    # *




# ================= 3. DATA LOADING & CLEANING =================


cat("[INFO] Reading data from:", INPUT_FILE, "\n")


if (!file.exists(INPUT_FILE)) {
  stop("[ERROR] File not found. Please check the path: ", INPUT_FILE)
}


df <- read.csv(INPUT_FILE, stringsAsFactors = FALSE)


# Clean column names
colnames(df) <- trimws(colnames(df))


# Validate required columns
required_cols <- c("category", "Threshold", "Observed", "P_Value")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("[ERROR] Missing required columns: ", paste(missing_cols, collapse = ", "))
}


# Clean category strings
df$category <- trimws(as.character(df$category))




# ================= 4. DATA PREPROCESSING =================


df <- df %>%
  mutate(
    # Identify subgroup based on keywords
    subgroup = case_when(
      str_detect(category, "HSTB") ~ "HSTB",
      str_detect(category, "TB")   ~ "TB",
      TRUE                         ~ "All"
    ),
    # Factor for consistent ordering (uses order from colors_map)
    category_factor = factor(category, levels = names(colors_map)),
    # Convert P_Value to numeric
    p_val_num = as.numeric(P_Value)
  )


if (nrow(df) == 0) {
  stop("[ERROR] No valid data found after cleaning.")
}


# Prepare X-axis ticks
x_ticks_vals <- sort(unique(df$Threshold))
x_ticks_labels <- paste0(as.integer(x_ticks_vals * 100), "%")




# ================= 5. SIGNIFICANCE STAR CALCULATION =================


stars_df <- df %>%
  mutate(
    symbol = case_when(
      p_val_num < SIG_STARS_4 ~ "****",
      p_val_num < SIG_STARS_3 ~ "***",
      p_val_num < SIG_STARS_2 ~ "**",
      p_val_num < SIG_STARS_1 ~ "*",
      TRUE                    ~ ""
    ),
    y_offset = case_when(
      subgroup == "All" ~ 0.025,
      TRUE              ~ 0.02
    ),
    label_y = Observed + y_offset
  ) %>%
  filter(symbol != "")




# ================= 6. PLOTTING =================


p <- ggplot(df, aes(x = Threshold, y = Observed, color = category_factor, linetype = subgroup, shape = subgroup)) +
  
  # Draw lines and points
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5, stroke = 1.2) +
  
  # Add significance stars
  geom_text(
    data = stars_df,
    aes(x = Threshold, y = label_y, label = symbol, color = category_factor),
    size = 5, 
    fontface = "bold",
    vjust = 0
  ) +
  
  # Apply Manual Scales (Labels removed since there is no legend)
  scale_color_manual(
    values = colors_map
  ) +
  scale_linetype_manual(
    values = linestyles_map
  ) +
  scale_shape_manual(
    values = shapes_map
  ) +
  
  scale_x_continuous(
    breaks = x_ticks_vals,
    labels = x_ticks_labels
  ) +
  
  labs(
    title = "TE Enrichment Robustness Across Functional Subgroups",
    x = "Overlap Threshold (TE Coverage)",
    y = "Fraction of TE-associated CBSs"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 20)),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    # REMOVE LEGEND
    legend.position = "none", 
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )




# ================= 7. FINAL ADJUSTMENTS & SAVING =================


# Set Y-axis limits dynamically
global_max <- max(df$Observed, na.rm = TRUE)
y_limit <- min(global_max * 1.25, 1.0)
p <- p + ylim(0, y_limit)


cat(sprintf("[INFO] Setting Y-axis limit to %.2f\n", y_limit))
cat("[INFO] Legend removed as requested.\n")


# Save Outputs
ggsave(OUTPUT_IMG, plot = p, width = 12, height = 8, dpi = 300, limitsize = FALSE)
ggsave(OUTPUT_PDF, plot = p, width = 12, height = 8, device = cairo_pdf, limitsize = FALSE)


cat("\n[SUCCESS] Plot saved successfully:\n")
cat("  - Image:", OUTPUT_IMG, "\n")
cat("  - PDF:  ", OUTPUT_PDF, "\n")


# Display plot
print(p)