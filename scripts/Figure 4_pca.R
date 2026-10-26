
# Set working directory
setwd("PCA")

# Load functions
library(ggplot2)
library(readxl)

# Read files
pca <- read.csv("pca.csv", header = T, sep = ",")
name <- read.csv("pca.name.csv", header = T, sep = ",")

# Read eigenvalue file (contains percentages of each principal component)
# Assume the eigenvalue file has one percentage value per line, or contains percentage information
eigenvals <- scan("pca.eigenval")

# Method 1: if the file already contains percentage values
percentages <- eigenvals  # use the read values directly

# Sort by PCH column from small to large
pca_sorted <- pca[order(pca$PCH), ]
name_sorted <- name[order(name$PCH), ]

#========================================================PCA1-PCA2
# Create labels (directly use percentages read from the file)
# Note: ensure percentages has at least 3 values
if(length(percentages) >= 3) {
  xlab <- paste0("PC1 (", percentages[1], "%)")
  ylab <- paste0("PC2 (", percentages[2], "%)")
} else {
  # If fewer than 3 values, use default labels
  warning("Not enough principal component data in eigenvalue file, using default labels")
  xlab <- "PC1"
  ylab <- "PC2"
}

# Plot
plot(pca_sorted$PCA1, pca_sorted$PCA2, 
     pch = pca_sorted$PCH,      # use fillable symbols
     col = pca_sorted$COL,      # border color
     bg = pca_sorted$BG,        # fill color
     xlab = xlab, 
     ylab = ylab,
     cex = 1.5)          # appropriately increase point size

# Add legend
legend("bottomright", 
       legend = name_sorted$BREED, 
       pch = name_sorted$PCH, 
       col = name_sorted$COL, 
       pt.bg = name_sorted$BG,
       pt.cex = 1.0,      # adjust point size in legend (default 1.0)
       cex = 0.9,        # adjust legend text size
       bty = "n")         # remove legend border


#========================================================PCA2-PCA3
# Create labels (directly use percentages read from the file)
# Note: ensure percentages has at least 3 values
if(length(percentages) >= 3) {
  xlab <- paste0("PC2 (", percentages[2], "%)")
  ylab <- paste0("PC3 (", percentages[3], "%)")
} else {
  # If fewer than 3 values, use default labels
  warning("Not enough principal component data in eigenvalue file, using default labels")
  xlab <- "PC2"
  ylab <- "PC3"
}

# Plot
plot(pca_sorted$PCA2, pca_sorted$PCA3, 
     pch = pca_sorted$PCH,      # use fillable symbols
     col = pca_sorted$COL,      # border color
     bg = pca_sorted$BG,        # fill color
     xlab = xlab, 
     ylab = ylab,
     cex = 1.5)          # appropriately increase point size

# Add legend
legend("bottomleft", 
       legend = name_sorted$BREED,
       pch = name_sorted$PCH, 
       col = name_sorted$COL, 
       pt.bg = name_sorted$BG,
       pt.cex = 1.0,      # adjust point size in legend (default 1.0)
       cex = 0.9,        # adjust legend text size
       bty = "n")         # remove legend border


#========================================================PCA1-PCA3
# Create labels (directly use percentages read from the file)
# Note: ensure percentages has at least 3 values
if(length(percentages) >= 3) {
  xlab <- paste0("PC1 (", percentages[1], "%)")
  ylab <- paste0("PC3 (", percentages[3], "%)")
} else {
  # If fewer than 3 values, use default labels
  warning("Not enough principal component data in eigenvalue file, using default labels")
  xlab <- "PC1"
  ylab <- "PC3"
}

# Plot
plot(pca_sorted$PCA1, pca_sorted$PCA3, 
     pch = pca_sorted$PCH,      # use fillable symbols
     col = pca_sorted$COL,      # border color
     bg = pca_sorted$BG,        # fill color
     xlab = xlab, 
     ylab = ylab,
     cex = 1.5)          # appropriately increase point size

# Add legend
legend("bottomright", 
       legend = name_sorted$BREED,
       pch = name_sorted$PCH, 
       col = name_sorted$COL, 
       pt.bg = name_sorted$BG,
       pt.cex = 1.0,      # adjust point size in legend (default 1.0)
       cex = 0.9,        # adjust legend text size
       bty = "n")         # remove legend border