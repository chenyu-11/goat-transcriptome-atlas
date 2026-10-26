
rm(list = ls())
library(Rtsne)
library(ggplot2)

# Read input data (files should be in the working directory)
data <- read.csv("filtered_gene_expression.csv", row.names = 1)

# Read sample group information
groups <- read.table("group.txt", header = TRUE, sep = "\t")

# Log transformation
data <- log2(data + 1)

# Z-score standardization
data <- scale(data)
data <- t(scale(t(data), center = TRUE, scale = FALSE))

# Transpose
data <- t(data)

# Set random seed for reproducibility
set.seed(321)

# Perform t-SNE dimensionality reduction
tsne_out <- Rtsne(
  data,
  dims = 2,
  pca = TRUE,
  max_iter = 1000,
  theta = 0.4,
  perplexity = 20,
  verbose = FALSE
)

# Convert t-SNE results to a data frame
tsne_result <- as.data.frame(tsne_out$Y)
colnames(tsne_result) <- c("tSNE1", "tSNE2")

# Add sample group information
tsne_result$group <- groups$group[
  match(rownames(tsne_result), groups$sample)
]
tsne_result$group <- factor(
  tsne_result$group,
  levels = unique(groups$group)
)

# Plot t-SNE results, colored by group
ggplot(tsne_result, aes(x = tSNE1, y = tSNE2, color = group)) +
  geom_point(size = 1) +
  scale_color_manual(values = c(
    "RUM" = "#90D92D", "RTC" = "#D2E8AF", "OMS" = "#658238", "ABM" = "#8fbc8f",
    "RCT" = "#f9c9e1", "CCM" = "#dda0dd", "ILM" = "#FF64B9", "CLN" = "#C66D9E",
    "LDM" = "#daa520", "QFM" = "#c8a600", "PRF" = "#ccccff", "SCF" = "#a4a4d5",
    "IMF" = "#d2d1ea", "SKN" = "#ff7110", "SPL" = "#e96508", "THY" = "#ff4d40",
    "CBL" = "#465F89", "BRN" = "#4D77C1", "CCX" = "#91AFE2", "KDN" = "#ccff00",
    "BLDR" = "#dfff00", "LNG" = "#51a4a6", "TRC" = "#4798b3", "HRT" = "#D3524F",
    "PNS" = "#ab42c5", "ESO" = "#e63995", "SMG" = "#998F9B", "LVR" = "#00ff9f",
    "ART" = "#8b4513", "PIT" = "#533dd4", "THA" = "#4636a5",
    "PG" = "#5e4cc5", "ADR" = "#ffd700"
  )) +
  theme_minimal() +
  labs(title = "t-SNE visualization of sample expression data")

