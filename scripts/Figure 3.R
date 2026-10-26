# Load required packages
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(patchwork)
library(reshape2)
library(corrplot)
library(pheatmap)

# Set working directory
setwd("ASE-snpeff-annotation")

# Read data
data <- read_excel("ASE-snpeff_annotation_results.xlsx")

# View column names to confirm
colnames(data)

# Data preprocessing (correctly handle columns with % in names)
data <- data %>%
  mutate(
    # Wrap column names with special characters in backticks
    ASE_SNP_pct = as.numeric(gsub("%", "", `ASE_SNP_%`)),
    Gene_Annotation_Rate = as.numeric(gsub("%", "", Gene_Annotation_Rate)),
    ASE_Gene_ge3_Rate = as.numeric(gsub("%", "", ASE_Gene_ge3_Rate)),
    # If Avg_SNPs_per_Gene is already numeric, use directly
    Avg_SNPs_per_Gene = as.numeric(Avg_SNPs_per_Gene)
  )

# Rename tissue column (if first column has no name)
colnames(data)[1] <- "Tissue"

# Extract main statistical metrics
main_stats <- data %>%
  select(Tissue, Total_Tested_SNPs, ASE_SNPs, ASE_Genes, ASE_Genes_ge3, Genes_w_ge5, 
         Avg_SNPs_per_Gene, Gene_Annotation_Rate, ASE_Gene_ge3_Rate)

# Functional classification data (long format)
func_data <- reshape2::melt(
  data %>% select(Tissue, MISSENSE, NONSENSE, SILENT),
  id.vars = "Tissue",
  variable.name = "Function_Class",
  value.name = "Count"
)

# Region annotation data (long format)
region_data <- data %>%
  select(Tissue, DOWNSTREAM, EXON, INTERGENIC, INTRON,
         SPLICE_SITE_REGION, TRANSCRIPT, UPSTREAM, UTR_3_PRIME, UTR_5_PRIME) %>%
  pivot_longer(
    cols = -Tissue,            # convert all columns except Tissue to long format
    names_to  = "Region",
    values_to = "Count"
  )

# ==================== Figure 1: ASE SNPs vs Genes (bubble plot) ====================
p1 <- ggplot(main_stats, aes(x = ASE_Genes, y = ASE_SNPs, 
                             color = Gene_Annotation_Rate,
                             size = Avg_SNPs_per_Gene)) +
  geom_point(alpha = 0.8) +
  geom_text(aes(label = Tissue), vjust = -1.2, size = 2.5, show.legend = FALSE) +
  scale_color_gradient(low = "blue", high = "red", name = "Annotation Rate (%)") +
  scale_size_continuous(range = c(2, 8), name = "Avg SNPs/Gene") +
  labs(title = "ASE SNPs vs Genes by Tissue",
       x = "Number of ASE Genes", 
       y = "Number of ASE SNPs") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    # Add black axis lines
    axis.line = element_line(color = "black", size = 0.8),
    # Add tick marks (pointing inward)
    axis.ticks = element_line(color = "black", size = 0.6),
    axis.ticks.length = unit(0.15, "cm"),
    # Ensure tick labels are retained
    axis.text = element_text(color = "black", size = 10)
  )

p1

# ==================== Figure 2: Tissue ASE activity ranking (bar chart) ====================
p2 <- main_stats %>%
  mutate(Tissue = factor(Tissue, levels = Tissue[order(-ASE_SNPs)])) %>%
  ggplot(aes(x = Tissue, y = ASE_SNPs, fill = ASE_Genes)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(ASE_Genes, " genes")), 
            vjust = -0.3, size = 3) +
  labs(title = "ASE SNPs by Tissue",
       x = "", y = "Number of ASE SNPs",
       fill = "ASE Genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2

# ==================== Figure 3: Genes_ge3 proportion (horizontal bar chart) ====================
p3 <- main_stats %>%
  mutate(Ge3_Ratio = ASE_Genes_ge3 / ASE_Genes * 100) %>%
  mutate(Tissue = factor(Tissue, levels = Tissue[order(Ge3_Ratio)])) %>%
  ggplot(aes(x = Tissue, y = Ge3_Ratio, fill = Ge3_Ratio)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Genes with ≥3 ASE SNPs (%)", 
       x = "", y = "Percentage") +
  scale_fill_gradient(low = "#80b1d3", high = "#80b1d3", name = "Ratio (%)") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    # Add black axis lines
    axis.line = element_line(color = "black", size = 0.8),
    # Add black tick marks
    axis.ticks = element_line(color = "black", size = 0.6),
    axis.ticks.length = unit(0.15, "cm"),
    # Ensure tick labels are black
    axis.text = element_text(color = "black", size = 10),
    # Axis titles in black (if any)
    axis.title = element_text(color = "black", size = 12, face = "bold"),
    # Remove grid lines (optional, for a cleaner plot)
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

p3

# ==================== Figure 4: Gene annotation rate (horizontal bar chart) ====================
p4 <- main_stats %>%
  mutate(Tissue = factor(Tissue, levels = Tissue[order(Gene_Annotation_Rate)])) %>%
  ggplot(aes(x = Tissue, y = Gene_Annotation_Rate, fill = Gene_Annotation_Rate)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Gene Annotation Rate by Tissue", 
       x = "", y = "Annotation Rate (%)") +
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "Rate (%)") +
  theme_minimal() +
  theme(legend.position = "bottom")

p4

# ==================== Figure 5: Functional classification stacked bar chart (3 classes) ====================
p5 <- func_data %>%
  filter(Function_Class %in% c("MISSENSE", "NONSENSE", "SILENT")) %>%
  mutate(Function_Class = factor(Function_Class, 
                                 levels = c("NONSENSE", "MISSENSE", "SILENT"))) %>%  # reverse order: NONSENSE on top, SILENT at bottom
  ggplot(aes(x = Tissue, y = Count, fill = Function_Class)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Functional Classification of ASE SNPs",
       x = "", y = "Number of SNPs",
       fill = "Function Class") +
  scale_fill_manual(values = c("SILENT" = "#caa38e", 
                               "MISSENSE" = "#f3d4c4", 
                               "NONSENSE" = "#1b9e77")) +
  theme_minimal() +
  theme(
    # Rotate x-axis labels
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    # Add black axis lines
    axis.line = element_line(color = "black", size = 0.8),
    # Add black tick marks
    axis.ticks = element_line(color = "black", size = 0.6),
    axis.ticks.length = unit(0.15, "cm"),
    # Axis titles
    axis.title = element_text(color = "black", size = 12, face = "bold"),
    # Legend position (optional)
    legend.position = "right",
    # Remove grid lines (optional, for a cleaner plot)
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

p5

# ==================== Figure 6: Genomic region heatmap (main regions) ====================
p6 <- data %>%
  select(Tissue, UTR_3_PRIME, UTR_5_PRIME, DOWNSTREAM, UPSTREAM,
         INTERGENIC, INTRON, EXON, TRANSCRIPT) %>%
  pivot_longer(
    cols = -Tissue,            # convert all except Tissue to long format
    names_to  = "Region",
    values_to = "Count"
  ) %>%
  ggplot(aes(x = Tissue, y = Region, fill = Count)) +
  geom_tile() +
  scale_fill_gradient(low = "#fff5f0", high = "darkred", name = "Count") +
  labs(title = "Genomic Region Distribution",
       x = "", y = "Region") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p6


