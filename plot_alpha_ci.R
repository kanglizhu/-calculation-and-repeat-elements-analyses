library(ggplot2)
library(dplyr)


# ================= Configuration Area =================
# Replace with your actual CSV file path
INPUT_FILE <- "alpha_data.csv" 
OUTPUT_FILE <- "alpha_bootstrap_updated.png"
# ======================================================


# 1. Data Loading
# Reads the dataset from a CSV file. 
# Ensure the file has columns: way_info, alpha, ci_lower, ci_upper, type, maf
cat("Reading data from:", INPUT_FILE, "\n")


if (!file.exists(INPUT_FILE)) {
  stop("Error: File not found. Please check the file path.")
}


data <- read.csv(INPUT_FILE, stringsAsFactors = FALSE)


# Verify required columns exist
required_cols <- c("way_info", "alpha", "ci_lower", "ci_upper", "type", "maf")
if (!all(required_cols %in% colnames(data))) {
  stop("Error: The CSV file is missing required columns. Please ensure it contains: ", 
       paste(required_cols, collapse = ", "))
}


# 2. Data Preprocessing
plot_data <- data %>%
  mutate(
    # Convert to factors to fix legend order
    type = factor(type, levels = c("intron", "synonymous")),
    maf = factor(maf, levels = c("10MAF", "25MAF")),
    # Ensure way_info follows the specific evolutionary context order
    way_info = factor(way_info, levels = c("hs", "hs_50KB", "5way", "5way_50KB"))
  )


# 3. Plotting
p <- ggplot(plot_data, aes(x = way_info, y = alpha, color = maf, shape = type)) +
  
  # Plot error bars
  # position_dodge ensures different groups at the same x-axis position are separated horizontally
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                width = 0.15,          # Width of the error bar caps
                size = 0.8,            # Line thickness
                position = position_dodge(width = 0.6)) +
  
  # Plot points
  geom_point(size = 4, 
             stroke = 1.2,             # Point border thickness
             position = position_dodge(width = 0.6)) +
  
  # Add y=0 reference line (Neutral evolution line)
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  # Manually set colors (Professional Blue/Red palette)
  scale_color_manual(values = c("10MAF" = "#2E86AB", "25MAF" = "#D62828"), 
                     name = "MAF Threshold") +
  
  # Manually set shapes (17=Triangle, 16=Circle)
  scale_shape_manual(values = c("intron" = 17, "synonymous" = 16), 
                     name = "Control Type",
                     labels = c("intron" = "Intron", "synonymous" = "Synonymous")) +
  
  # Labels and Title
  labs(
    title = "Proportion of Adaptation (α) by Evolutionary Context",
    subtitle = "Comparison of Intron vs. Synonymous Controls with 95% Bootstrap CIs",
    x = "Evolutionary Context (Way Info)",
    y = expression(paste(alpha, " Estimate"))
  ) +
  
  # Theme customization
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, face = "bold"), # X-axis labels horizontal
    panel.grid.major.x = element_blank(), # Remove vertical grid lines
    panel.grid.minor = element_blank(),
    legend.position = "top", 
    legend.box = "horizontal", # Arrange legend horizontally
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
    axis.title = element_text(face = "bold")
  )


# 4. Display Plot
print(p)


# 5. Save Plot
ggsave(OUTPUT_FILE, plot = p, width = 10, height = 6, dpi = 300)
cat("Plot successfully saved as:", OUTPUT_FILE, "\n")