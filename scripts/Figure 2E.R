# Load required libraries
library(ggplot2)
library(readxl)
library(dplyr)
library(scales)
library(ggrepel)  # for intelligent label placement

# Set working directory to the GO enrichment folder
setwd("GO")

# 1. Read data
data <- read_excel("go_results.xlsx")

# 2. Inspect data structure
str(data)
head(data)

# 3. Ensure data formats are correct
data$Shape <- as.factor(data$Shape)
data$Color <- as.character(data$Color)
data$Label <- as.character(data$Label)

# 4. Create mapping from label to color/shape
# Ensure each label has a unique color and shape
label_info <- data %>%
  distinct(Label, Color, Shape) %>%
  arrange(Label)

# Create color mapping (directly use colors from your data)
color_mapping <- setNames(label_info$Color, label_info$Label)

# Create shape mapping (directly use shapes from your data)
shape_mapping <- setNames(as.numeric(as.character(label_info$Shape)), label_info$Label)

# 5. Draw scatter plot with labels annotated directly on points
p <- ggplot(data, aes(x = logFDR, y = `Fold Enrichment`)) +
  
  # Add points, size determined by Genes, color and shape mapped directly from Label
  geom_point(aes(size = Genes, color = Label, shape = Label),
             alpha = 0.8, stroke = 1.2) +
  
  # Use colors specified in your data
  scale_color_manual(values = color_mapping, guide = "none") +
  
  # Use shapes specified in your data
  scale_shape_manual(values = shape_mapping, guide = "none") +
  
  # Adjust point size range
  scale_size_continuous(
    range = c(3, 10),  # adjust size range
    name = "Number of Genes",
    breaks = scales::pretty_breaks(n = 5),
    guide = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      label.position = "right"
    )
  ) +
  
  # Use ggrepel to add text labels (avoid overlap)
  geom_text_repel(
    aes(label = Label, color = Label),  # text color matches point color
    size = 3.5,                         # text size
    fontface = "bold",                  # bold text
    segment.color = "grey60",           # connector line color
    segment.size = 0.3,                 # connector line width
    segment.alpha = 0.6,                # connector line transparency
    box.padding = 0.5,                  # text box padding
    point.padding = 0.3,                # distance from point to text box
    min.segment.length = 0.2,           # minimum connector line length
    max.overlaps = Inf,                 # allow unlimited overlap handling
    max.time = 2,                       # maximum computation time (seconds)
    max.iter = 10000,                   # maximum iterations
    show.legend = FALSE                 # do not show in legend
  ) +
  
  # Add axis labels
  labs(
    x = "-log10(FDR)",
    y = "Fold Enrichment",
    title = "GO Enrichment Analysis",
    subtitle = "Point size represents number of genes"
  ) +
  
  # Beautify theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
    legend.position = "right",           # only show Genes size legend
    legend.title.align = 0.5,            # center legend title
    legend.box.margin = margin(0, 0, 0, 10),  # legend outer margin
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(20, 20, 20, 20)  # plot margin
  )

# Display the plot
print(p)

# Save the plot (high resolution)
ggsave("GO_enrichment_plot_with_labels.pdf",
       plot = p,
       width = 14,
       height = 10,
       dpi = 600,
       bg = "white")

# 6. Optional: create a validation plot to confirm color and shape mapping
cat("Color mapping validation:\n")
for(i in 1:nrow(label_info)) {
  cat(sprintf("%2d. %-40s: color=%s, shape=%s\n", 
              i, label_info$Label[i], 
              label_info$Color[i], label_info$Shape[i]))
}

# 7. If too many labels overlap severely, try adjusting ggrepel parameters
if (nrow(data) > 10) {
  message("Many data points (", nrow(data), "), optimizing label layout...")
  
  # Create optimized version
  p_optimized <- ggplot(data, aes(x = logFDR, y = `Fold Enrichment`)) +
    geom_point(aes(size = Genes, color = Label, shape = Label),
               alpha = 0.8, stroke = 1.2) +
    scale_color_manual(values = color_mapping, guide = "none") +
    scale_shape_manual(values = shape_mapping, guide = "none") +
    scale_size_continuous(range = c(3, 12), name = "Number of Genes") +
    
    # Optimized label parameters: smaller label size, larger spacing
    geom_text_repel(
      aes(label = Label, color = Label),
      size = 3.0,  # smaller label size
      fontface = "plain",  # not bold
      segment.color = "grey50",
      segment.size = 0.2,
      segment.alpha = 0.5,
      box.padding = 0.8,   # larger padding
      point.padding = 0.5,
      min.segment.length = 0.1,
      max.overlaps = 20,   # limit overlap handling
      force = 1.0,         # adjust repulsion force
      nudge_x = 0.15,
      nudge_y = 0.15,
      show.legend = FALSE
    ) +
    
    labs(x = "-log10(FDR)", y = "Fold Enrichment",
         title = "GO Enrichment Analysis (Optimized Labels)") +
    theme_minimal() +
    theme(legend.position = "right")
  
  ggsave("GO_enrichment_plot_optimized.png",
         plot = p_optimized,
         width = 14, height = 10, dpi = 600, bg = "white")
  
  print(p_optimized)
}