# ============================================
# Step 1: Load required packages
# ============================================
library(readxl)
library(UpSetR)

# ============================================
# Step 2: Read the Excel file
# ============================================
file_path <- "12_groups_gene.xlsx"
raw_data <- read_excel(file_path, col_names = TRUE)

# View column names
names(raw_data)

# ============================================
# Step 3: Build gene lists directly (keep all original values, including date formats)
# ============================================

# Convert each column to a vector and remove NAs
gene_lists <- list(
  FST = na.omit(raw_data$FST),
  Tajima_pos = na.omit(raw_data$Tajima_pos),
  Tajima_neg = na.omit(raw_data$Tajima_neg),
  Muscle = na.omit(raw_data$Muscle),
  Adipose = na.omit(raw_data$Adipose),
  SKN = na.omit(raw_data$SKN),
  SPL = na.omit(raw_data$SPL),
  forestomach = na.omit(raw_data$forestomach),
  ABM = na.omit(raw_data$ABM),
  CLN = na.omit(raw_data$CLN),
  ILM = na.omit(raw_data$ILM),
  RCT = na.omit(raw_data$RCT)
)

# View the number of genes in each group
sapply(gene_lists, length)

# ============================================
# Step 4: Build a binary matrix
# ============================================

# Get all unique genes
all_genes <- unique(unlist(gene_lists))

# Build the matrix
binary_matrix <- matrix(0, nrow = length(all_genes), ncol = length(gene_lists))
colnames(binary_matrix) <- names(gene_lists)
rownames(binary_matrix) <- all_genes

# Fill the matrix
for (i in 1:length(gene_lists)) {
  binary_matrix[gene_lists[[i]], i] <- 1
}

# View the result
head(binary_matrix, 10)

# Check the number of genes in each group (validation)
colSums(binary_matrix)

# ============================================
# Step 5: Export the matrix (optional)
# ============================================
write.csv(binary_matrix, "binary_matrix_12groups.csv", row.names = TRUE)

# ============================================
# Step 6: Draw the UpSet plot
# ============================================
upset(fromList(gene_lists),
      nsets = 12,
      nintersects = 50, # Use a sufficiently large number to show all sets
      order.by = "freq",
      keep.order = TRUE,
      sets = c("FST", "Tajima_neg", "Tajima_pos", "Muscle", "Adipose", "SKN", "SPL", "forestomach", "ABM", "CLN", "ILM", "RCT"),
      main.bar.color = "#00304b",
      mainbar.y.label = "Intersection Size",
      sets.x.label = "Set Size",
      text.scale = 1.5,
      mb.ratio = c(0.6, 0.4))  # Adjust the ratio between the main plot and the set size plot