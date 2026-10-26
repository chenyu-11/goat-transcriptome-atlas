
rm(list=ls())

library(CMplot)

setwd("fst")

### Read main data file
data <- read.table("weir.plot.fst", header=T, sep="\t", stringsAsFactors = F)

# Rename columns
colnames(data) <- c("CHROM", "POS", "FST")
data$CHROM <- gsub("chr", "", data$CHROM) # remove "chr" characters
data$SNP <- paste0(data$CHROM, ":", data$POS)
data$CHROM <- as.character(data$CHROM)  # ensure chromosome is character
data$FST[data$FST < 0] <- 0 # replace negative values with 0
data <- data[, c("SNP", "CHROM", "POS", "FST")] # reorder columns

#### Determine threshold lines - top 5%
# Calculate the 95th percentile of the FST column
percentile_95 <- quantile(data[, 4], 0.95, na.rm=TRUE)
print(paste("95th percentile threshold:", percentile_95))

# Calculate the 99th percentile of the FST column
percentile_99 <- quantile(data[, 4], 0.99, na.rm=TRUE)
print(paste("99th percentile threshold:", percentile_99))

## Plotting
# Define the plotting function
plot_manhattan <- function(file_type) {
  tryCatch({
    CMplot(data[,1:4], plot.type="m", 
           type = "h",
           col=c("#80b1d3","#4a7c9c"),
           LOG10=FALSE, ylim=NULL, threshold=c(percentile_95),
           threshold.lty=c(2,2), threshold.lwd=c(2,2), 
           threshold.col=c("red","blue"), amplify=FALSE,
           signal.col=c("red","green"), 
           signal.cex=c(1,1), signal.pch=c(19,19),
           ylab="Fst",  # set ylab only once
           cex=0.5,     # set cex only once
           file = file_type,
           dpi = 300,
           file.output = TRUE,
           verbose = TRUE,
           width = 14,
           height = 6)
  }, error = function(e) {
    print(paste("Plotting error:", e$message))
    
    # Fallback: plot with the 99th percentile threshold instead
    print("Attempting to plot with the 99th percentile threshold...")
    CMplot(data[,1:4], plot.type="m", 
           #type = "h",
           col=c("#80b1d3","#4a7c9c"),
           LOG10=FALSE, ylim=NULL, threshold=c(percentile_99),
           threshold.lty=c(2,2), threshold.lwd=c(2,2), 
           threshold.col=c("red","blue"), amplify=FALSE,
           signal.col=c("red","green"), 
           signal.cex=c(1,1), signal.pch=c(19,19),
           ylab="Fst",
           cex=0.7,
           file = file_type,
           dpi = 300,
           file.output = TRUE,
           verbose = TRUE,
           width = 14,
           height = 6)
  })
}

#### Output PDF and JPG
# First output PDF
plot_manhattan("pdf")

# Second output JPG
plot_manhattan("jpg")