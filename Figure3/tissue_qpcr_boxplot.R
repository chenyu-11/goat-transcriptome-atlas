# Load required packages
library(tidyverse)
library(ggplot2)
library(readxl)  # For reading Excel files

# Read data from the Excel file
data <- read_excel("qpcr-muscle4.xlsx")

# Inspect data structure and confirm column names
print(head(data))
print(str(data))
print(colnames(data))

# Define tissue order (as requested)
tissue_order <- c("Heart", "Longest dorsi muscle", "Quadriceps femoris muscle", "Liver", "Lung", "Spleen", "Skin")

# Convert Tissue column to a factor with the specified order
data$Tissue <- factor(data$Tissue, levels = tissue_order)

# Calculate the mean for each individual (Sample)
data_mean <- data %>%
  group_by(Tissue, Sample) %>%
  summarise(across(c(NMRK2, TNNI2, MYF6,ATF4), mean, na.rm = TRUE), .groups = "drop")

# Ensure Tissue in data_mean is also a factor with the same order
data_mean$Tissue <- factor(data_mean$Tissue, levels = tissue_order)

# View the processed data
print(data_mean)

# Get the tissue types actually present in the data
tissues_present <- unique(data_mean$Tissue)
print(paste("Tissues present in data:", paste(tissues_present, collapse = ", ")))

# Define color mapping (in the specified order)
color_mapping_all <- c(
  "Heart" = "#D35240",
  "Longest dorsi muscle" = "#daa520",
  "Quadriceps femoris muscle" = "#c8a600",
  "Liver" = "#00ff9f",
  "Lung" = "#51a4a6",
  "Spleen" = "#e96508",
  "Skin" = "#ff7110"
)

# Keep only colors for tissues present in the data, ordered by tissue_order
color_mapping <- color_mapping_all[names(color_mapping_all) %in% tissues_present]
# Reorder according to tissue_order
color_mapping <- color_mapping[tissue_order[tissue_order %in% names(color_mapping)]]

# Automatically assign colors for any new tissues not in the predefined colors
new_tissues <- tissues_present[!tissues_present %in% names(color_mapping_all)]
if(length(new_tissues) > 0) {
  library(RColorBrewer)
  new_colors <- brewer.pal(max(3, length(new_tissues)), "Set3")[1:length(new_tissues)]
  names(new_colors) <- new_tissues
  color_mapping <- c(color_mapping, new_colors)
  print(paste("Auto-assigned colors for:", paste(new_tissues, collapse = ", ")))
}

print("Color mapping used:")
print(color_mapping)

# Convert data to long format for plotting
data_long <- data_mean %>%
  pivot_longer(cols = c( MYF6,ATF4,NMRK2,TNNI2), 
               names_to = "Gene", 
               values_to = "Expression")

# Ensure Tissue in data_long is also a factor with the same order
data_long$Tissue <- factor(data_long$Tissue, levels = tissue_order)

# Draw boxplots (combined plot) - add full border
p <- ggplot(data_long, aes(x = Tissue, y = Expression, fill = Tissue)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = color_mapping) +
  theme_classic() +  # Use the classic theme, which has black axis lines
  theme(
    # Add a full border (top, bottom, left, right)
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),  # Add black border
    
    # Axis line settings (keep left and bottom axis lines)
    axis.line = element_line(color = "black", size = 0.8),
    axis.ticks = element_line(color = "black", size = 0.8),
    axis.ticks.length = unit(0.2, "cm"),
    
    # Axis text settings
    axis.text = element_text(color = "black", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    
    # Axis title settings
    axis.title = element_text(color = "black", size = 12, face = "bold"),
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    
    # Facet strip text settings
    strip.text = element_text(color = "black", size = 12, face = "bold"),
    strip.background = element_rect(fill = "gray90", color = "black", size = 0.5),
    strip.placement = "inside",  # Place facet titles inside the border
    
    # Legend settings
    legend.position = "none",
    
    # Background and grid
    panel.grid.major = element_line(color = "gray80", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    
    # Ensure complete facet borders
    panel.spacing = unit(0.5, "lines")  # Spacing between facets
  ) +
  labs(
    title = "Gene Expression by Tissue",
    x = "Tissue",
    y = "Relative mRNA expression"
  )

# Display the plot
print(p)

# Save the combined plot
ggsave("gene_expression_boxplot_combined.png", p, width = 10, height = 8, dpi = 300)
ggsave("gene_expression_boxplot_combined.pdf", p, width = 10, height = 8, dpi = 300)