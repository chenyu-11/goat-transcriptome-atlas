# ================================================
# Gene enrichment analysis bubble plot (Fold enrichment version)
# X-axis: Fold enrichment   Color: -log10(PValue) / PValue
# ================================================

library(ggplot2)
library(dplyr)
library(scales)

# ==== Read input data ====
inputFile <- "go-bubble.csv"
rt <- read.table(inputFile, header = TRUE, sep = ",", check.names = FALSE)

# ==== Handle column names: rename "Fold enrichment" (with space) to FoldEnrichment ====
if ("Fold enrichment" %in% colnames(rt)) {
  rt <- rt %>% rename(FoldEnrichment = `Fold enrichment`)
}

# ==== If FoldEnrichment is not in the input, compute it from Count/ListTotal and PopHits/PopTotal ====
if (!"FoldEnrichment" %in% colnames(rt)) {
  # Requires List Total, Pop Hits, Pop Total columns
  need <- c("List Total", "Pop Hits", "Pop Total")
  if (all(need %in% colnames(rt))) {
    rt <- rt %>%
      mutate(FoldEnrichment = (Count / `List Total`) / (`Pop Hits` / `Pop Total`))
  } else {
    stop("❌ The input file has no FoldEnrichment column and is missing List Total / Pop Hits / Pop Total, so it cannot be computed.")
  }
}

# ==== Check whether required columns exist ====
required_cols <- c("Term", "FoldEnrichment", "Count", "PValue")
if (!all(required_cols %in% colnames(rt))) {
  stop("❌ The input file must contain columns: Term, FoldEnrichment (or computable columns), Count, PValue")
}

# ================================================
# Version 1: Color = -log10(PValue)
# ================================================
rt1 <- rt %>%
  mutate(log10PValue = -log10(PValue)) %>%
  arrange(FoldEnrichment) %>%
  mutate(Term = factor(Term, levels = Term))

p1 <- ggplot(rt1, aes(x = FoldEnrichment, y = Term)) +
  geom_point(aes(size = Count, color = log10PValue)) +
  scale_color_gradientn(
    colours = c("#4575B4", "#91BFDB", "#E0F3F8", "#FEE090", "#FC8D59", "#D73027"),
    name = expression(-log[10](PValue))
  ) +
  labs(size = "Count", x = "Fold enrichment", y = "Term") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10)
  )

print(p1)
outFile1 <- paste0("bubble_logP_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
ggsave(outFile1, plot = p1, width = 7, height = 5)
cat("✅ Image saved as: ", outFile1, "\n")

# ================================================
# Version 2: Color = PValue (red → yellow → blue; smaller is redder = more significant)
# ================================================
rt2 <- rt %>%
  arrange(FoldEnrichment) %>%
  mutate(Term = factor(Term, levels = Term))

p2 <- ggplot(rt2, aes(x = FoldEnrichment, y = Term)) +
  geom_point(aes(size = Count, color = PValue)) +
  scale_color_gradientn(
    colours = c("#D73027", "#FC8D59", "#FEE090", "#E0F3F8", "#91BFDB", "#4575B4"),
    values = rescale(c(0, 0.01, 0.05, 0.1, 0.2, 1)),
    name = "PValue"
  ) +
  labs(size = "Count", x = "Fold enrichment", y = "Term") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10)
  )

print(p2)
outFile2 <- paste0("bubble_P_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
ggsave(outFile2, plot = p2, width = 7, height = 5)
cat("✅ Image saved as: ", outFile2, "\n")