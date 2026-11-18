# ============================================================
# Introduction
# ------------------------------------------------------------
# This script draws ADMIXTURE ancestry composition bar plots
# for K = 2 to K = 9. For each K, it reads a genotype .fam file
# and the corresponding .Q file, reorders samples by a
# predefined population order, assigns a custom color to each
# ancestry component, and then plots one horizontal bar per
# individual with population boundaries and labels.
#
# Input files (renamed to hide project-specific prefixes):
#   genotype.fam  : sample and population information
#   K2.Q ... K9.Q : ADMIXTURE ancestry proportions for each K
#
# The population order is defined in the vector `pop` below.
# Colors are set manually for each K so that components are
# visually distinct and consistent across panels.
# ============================================================

rm(list = ls())
dev.off()

par(mfrow = c(8, 1))

####=====================================K=2
### k2
fam = read.table("genotype.fam", header = F, sep = "\t", stringsAsFactors = F)
k2 = read.table("K2.Q", header = F, sep = " ")
rownames(k2) = fam[, 2]

## reorder fam and k2 file
sub = NULL; fam1 = NULL
for (i in 1:length(pop)) {
  sub = fam[fam[, 1] %in% pop[i], ]
  fam1 = rbind(fam1, sub)
}
k2 = k2[fam1[, 2], ]

sub = NULL; name = NULL
popn = as.data.frame(table(fam1[, 1]))
for (i in 1:length(pop)) {
  sub = popn[popn[, 1] %in% pop[i], ]
  name = rbind(name, sub)
}

## get colors
#library(rcolors)
#col= sample(get_color(rcolors$t2m_29lev,4))

## Custom colors
col = c(
  "#b2df8a",
  "#d6604d"
)

# Set smaller margins
par(mar = c(1, 4, 1, 1))  # bottom, left, top, right margins

## Draw the plot
barplot(t(as.matrix(k2)),
        #        col = col[1:4],
        col = col,
        ylab = "K=2",
        ylim = c(0, 1),
        border = NA,
        space = 0,
        xaxt = "n")

## Add dividing lines
segments(cumsum(name[, 2]), 0, cumsum(name[, 2]), 1, cex = 0.9, col = "grey")

## Add tick marks at the end of each population
axis(1, at = cumsum(name[, 2]), labels = FALSE, cex = 7, line = 0, srt = 135, adj = 0.5, cex.axis = 0.6, cex.lab = 3, cex.main = 6.5)
text(x = cumsum(name[, 2]), y = par()$usr[3] - 0.05 * (par()$usr[4] - par()$usr[3]),
     labels = name[, 1], srt = 0, adj = 1, xpd = TRUE, cex = 0.8)

####=====================================K=3
### k3
fam = read.table("genotype.fam", header = F, sep = "\t", stringsAsFactors = F)
k3 = read.table("K3.Q", header = F, sep = " ")
rownames(k3) = fam[, 2]

## reorder fam and k3 file
sub = NULL; fam1 = NULL
for (i in 1:length(pop)) {
  sub = fam[fam[, 1] %in% pop[i], ]
  fam1 = rbind(fam1, sub)
}
k3 = k3[fam1[, 2], ]

sub = NULL; name = NULL
popn = as.data.frame(table(fam1[, 1]))
for (i in 1:length(pop)) {
  sub = popn[popn[, 1] %in% pop[i], ]
  name = rbind(name, sub)
}

## get colors
#library(rcolors)
#col= sample(get_color(rcolors$t2m_29lev,4))

## Custom colors
col = c(
  "#b2df8a",
  "#d6604d",
  "#d89700"
)

# Set smaller margins
par(mar = c(1, 4, 1, 1))  # bottom, left, top, right margins

## Draw the plot
barplot(t(as.matrix(k3)),
        #        col = col[1:4],
        col = col,
        ylab = "K=3",
        ylim = c(0, 1),
        border = NA,
        space = 0,
        xaxt = "n")

## Add dividing lines
segments(cumsum(name[, 2]), 0, cumsum(name[, 2]), 1, cex = 0.9, col = "grey")

## Add tick marks at the end of each population
axis(1, at = cumsum(name[, 2]), labels = FALSE, cex = 7, line = 0, srt = 135, adj = 0.5, cex.axis = 0.6, cex.lab = 3, cex.main = 6.5)
text(x = cumsum(name[, 2]), y = par()$usr[3] - 0.05 * (par()$usr[4] - par()$usr[3]),
     labels = name[, 1], srt = 0, adj = 1, xpd = TRUE, cex = 0.8)

####=====================================K=4
### k4
fam = read.table("genotype.fam", header = F, sep = "\t", stringsAsFactors = F)
k4 = read.table("K4.Q", header = F, sep = " ")
rownames(k4) = fam[, 2]

## reorder fam and k4 file
sub = NULL; fam1 = NULL
for (i in 1:length(pop)) {
  sub = fam[fam[, 1] %in% pop[i], ]
  fam1 = rbind(fam1, sub)
}
k4 = k4[fam1[, 2], ]

sub = NULL; name = NULL
popn = as.data.frame(table(fam1[, 1]))
for (i in 1:length(pop)) {
  sub = popn[popn[, 1] %in% pop[i], ]
  name = rbind(name, sub)
}

## get colors
#library(rcolors)
#col= sample(get_color(rcolors$t2m_29lev,4))

## Custom colors
col = c(
  "#b2df8a",
  "#d6604d",
  "#d89700",
  "#a6cee3"
)

# Set smaller margins
par(mar = c(1, 4, 1, 1))  # bottom, left, top, right margins

## Draw the plot
barplot(t(as.matrix(k4)),
        #        col = col[1:4],
        col = col,
        ylab = "K=4",
        ylim = c(0, 1),
        border = NA,
        space = 0,
        xaxt = "n")

## Add dividing lines
segments(cumsum(name[, 2]), 0, cumsum(name[, 2]), 1, cex = 0.9, col = "grey")

## Add tick marks at the end of each population
axis(1, at = cumsum(name[, 2]), labels = FALSE, cex = 7, line = 0, srt = 135, adj = 0.5, cex.axis = 0.6, cex.lab = 3, cex.main = 6.5)
text(x = cumsum(name[, 2]), y = par()$usr[3] - 0.05 * (par()$usr[4] - par()$usr[3]),
     labels = name[, 1], srt = 0, adj = 1, xpd = TRUE, cex = 0.6)

####=====================================K=5
### k5
fam = read.table("genotype.fam", header = F, sep = "\t", stringsAsFactors = F)
k5 = read.table("K5.Q", header = F, sep = " ")
rownames(k5) = fam[, 2]

## reorder fam and k5 file
sub = NULL; fam1 = NULL
for (i in 1:length(pop)) {
  sub = fam[fam[, 1] %in% pop[i], ]
  fam1 = rbind(fam1, sub)
}
k5 = k5[fam1[, 2], ]

sub = NULL; name = NULL
popn = as.data.frame(table(fam1[, 1]))
for (i in 1:length(pop)) {
  sub = popn[popn[, 1] %in% pop[i], ]
  name = rbind(name, sub)
}

## get colors
#library(rcolors)
#col= sample(get_color(rcolors$t2m_29lev,4))

## Custom colors
col = c(
  "#b2df8a",
  "#d6604d",
  "#d89700",
  "#a6cee3",
  "#6a3d9a"
)

# Set smaller margins
par(mar = c(1, 4, 1, 1))  # bottom, left, top, right margins

## Draw the plot
barplot(t(as.matrix(k5)),
        #        col = col[1:4],
        col = col,
        ylab = "K=5",
        ylim = c(0, 1),
        border = NA,
        space = 0,
        xaxt = "n")

## Add dividing lines
segments(cumsum(name[, 2]), 0, cumsum(name[, 2]), 1, cex = 0.9, col = "grey")

## Add tick marks at the end of each population
axis(1, at = cumsum(name[, 2]), labels = FALSE, cex = 7, line = 0, srt = 135, adj = 0.5, cex.axis = 0.6, cex.lab = 3, cex.main = 6.5)
text(x = cumsum(name[, 2]), y = par()$usr[3] - 0.05 * (par()$usr[4] - par()$usr[3]),
     labels = name[, 1], srt = 0, adj = 1, xpd = TRUE, cex = 0.6)

####=====================================K=6
### k6
fam = read.table("genotype.fam", header = F, sep = "\t", stringsAsFactors = F)
k6 = read.table("K6.Q", header = F, sep = " ")
rownames(k6) = fam[, 2]

## reorder fam and k6 file
sub = NULL; fam1 = NULL
for (i in 1:length(pop)) {
  sub = fam[fam[, 1] %in% pop[i], ]
  fam1 = rbind(fam1, sub)
}
k6 = k6[fam1[, 2], ]

sub = NULL; name = NULL
popn = as.data.frame(table(fam1[, 1]))
for (i in 1:length(pop)) {
  sub = popn[popn[, 1] %in% pop[i], ]
  name = rbind(name, sub)
}

## get colors
#library(rcolors)
#col= sample(get_color(rcolors$t2m_29lev,4))

## Custom colors
col = c(
  "#b2df8a",
  "#d6604d",
  "#d89700",
  "#a6cee3",
  "#6a3d9a",
  "#ff7f00"
)

# Set smaller margins
par(mar = c(1, 4, 1, 1))  # bottom, left, top, right margins

## Draw the plot
barplot(t(as.matrix(k6)),
        #        col = col[1:4],
        col = col,
        ylab = "K=6",
        ylim = c(0, 1),
        border = NA,
        space = 0,
        xaxt = "n")

## Add dividing lines
segments(cumsum(name[, 2]), 0, cumsum(name[, 2]), 1, cex = 0.9, col = "grey")

## Add tick marks at the end of each population
axis(1, at = cumsum(name[, 2]), labels = FALSE, cex = 7, line = 0, srt = 135, adj = 0.5, cex.axis = 0.6, cex.lab = 3, cex.main = 6.5)
text(x = cumsum(name[, 2]), y = par()$usr[3] - 0.05 * (par()$usr[4] - par()$usr[3]),
     labels = name[, 1], srt = 0, adj = 1, xpd = TRUE, cex = 0.6)

####=====================================K=7
### k7
fam = read.table("genotype.fam", header = F, sep = "\t", stringsAsFactors = F)
k7 = read.table("K7.Q", header = F, sep = " ")
rownames(k7) = fam[, 2]

## reorder fam and k7 file
sub = NULL; fam1 = NULL
for (i in 1:length(pop)) {
  sub = fam[fam[, 1] %in% pop[i], ]
  fam1 = rbind(fam1, sub)
}
k7 = k7[fam1[, 2], ]

sub = NULL; name = NULL
popn = as.data.frame(table(fam1[, 1]))
for (i in 1:length(pop)) {
  sub = popn[popn[, 1] %in% pop[i], ]
  name = rbind(name, sub)
}

## get colors
#library(rcolors)
#col= sample(get_color(rcolors$t2m_29lev,4))

## Custom colors
col = c(
  "#b2df8a",
  "#d6604d",
  "#d89700",
  "#a6cee3",
  "#6a3d9a",
  "#ff7f00",
  "#1f78b4"
)

# Set smaller margins
par(mar = c(1, 4, 1, 1))  # bottom, left, top, right margins

## Draw the plot
barplot(t(as.matrix(k7)),
        #col = col[1:4],
        col = col,
        ylab = "K=7",
        ylim = c(0, 1),
        border = NA,
        space = 0,
        xaxt = "n")

## Add dividing lines
segments(cumsum(name[, 2]), 0, cumsum(name[, 2]), 1, cex = 0.9, col = "grey")

## Add tick marks at the end of each population
axis(1, at = cumsum(name[, 2]), labels = FALSE, cex = 7, line = 0, srt = 135, adj = 0.5, cex.axis = 0.6, cex.lab = 3, cex.main = 6.5)
text(x = cumsum(name[, 2]), y = par()$usr[3] - 0.05 * (par()$usr[4] - par()$usr[3]),
     labels = name[, 1], srt = 0, adj = 1, xpd = TRUE, cex = 0.6)

####=====================================K=8
### k8
fam = read.table("genotype.fam", header = F, sep = "\t", stringsAsFactors = F)
k8 = read.table("K8.Q", header = F, sep = " ")
rownames(k8) = fam[, 2]

## reorder fam and k8 file
sub = NULL; fam1 = NULL
for (i in 1:length(pop)) {
  sub = fam[fam[, 1] %in% pop[i], ]
  fam1 = rbind(fam1, sub)
}
k8 = k8[fam1[, 2], ]

sub = NULL; name = NULL
popn = as.data.frame(table(fam1[, 1]))
for (i in 1:length(pop)) {
  sub = popn[popn[, 1] %in% pop[i], ]
  name = rbind(name, sub)
}

## get colors
#library(rcolors)
#col= sample(get_color(rcolors$t2m_29lev,4))

## Custom colors
col = c(
  "#b2df8a",
  "#d6604d",
  "#d89700",
  "#a6cee3",
  "#6a3d9a",
  "#ff7f00",
  "#1f78b4",
  "#fb9a99"
)

# Set smaller margins
par(mar = c(1, 4, 1, 1))  # bottom, left, top, right margins

## Draw the plot
barplot(t(as.matrix(k8)),
        #        col = col[1:4],
        col = col,
        ylab = "K=8",
        ylim = c(0, 1),
        border = NA,
        space = 0,
        xaxt = "n")

## Add dividing lines
segments(cumsum(name[, 2]), 0, cumsum(name[, 2]), 1, cex = 0.9, col = "grey")

## Add tick marks at the end of each population
axis(1, at = cumsum(name[, 2]), labels = FALSE, cex = 7, line = 0, srt = 135, adj = 0.5, cex.axis = 0.6, cex.lab = 3, cex.main = 6.5)
text(x = cumsum(name[, 2]), y = par()$usr[3] - 0.05 * (par()$usr[4] - par()$usr[3]),
     labels = name[, 1], srt = 0, adj = 1, xpd = TRUE, cex = 0.6)

####=====================================K=9
### k9
fam = read.table("genotype.fam", header = F, sep = "\t", stringsAsFactors = F)
k9 = read.table("K9.Q", header = F, sep = " ")
rownames(k9) = fam[, 2]

## reorder fam and k9 file
sub = NULL; fam1 = NULL
for (i in 1:length(pop)) {
  sub = fam[fam[, 1] %in% pop[i], ]
  fam1 = rbind(fam1, sub)
}
k9 = k9[fam1[, 2], ]

sub = NULL; name = NULL
popn = as.data.frame(table(fam1[, 1]))
for (i in 1:length(pop)) {
  sub = popn[popn[, 1] %in% pop[i], ]
  name = rbind(name, sub)
}

## get colors
#library(rcolors)
#col= sample(get_color(rcolors$t2m_29lev,4))

## Custom colors
col = c(
  "#b2df8a",
  "#d6604d",
  "#d89700",
  "#a6cee3",
  "#6a3d9a",
  "#ff7f00",
  "#1f78b4",
  "#fb9a99",
  "#b3a7b8"
)

# Set smaller margins
par(mar = c(1, 4, 1, 1))  # bottom, left, top, right margins

## Draw the plot
barplot(t(as.matrix(k9)),
        #        col = col[1:4],
        col = col,
        ylab = "K=9",
        ylim = c(0, 1),
        border = NA,
        space = 0,
        xaxt = "n")

## Add dividing lines
segments(cumsum(name[, 2]), 0, cumsum(name[, 2]), 1, cex = 0.9, col = "grey")

## Add tick marks at the end of each population
axis(1, at = cumsum(name[, 2]), labels = FALSE, cex = 7, line = 0, srt = 135, adj = 0.5, cex.axis = 0.6, cex.lab = 3, cex.main = 6.5)
text(x = cumsum(name[, 2]), y = par()$usr[3] - 0.05 * (par()$usr[4] - par()$usr[3]),
     labels = name[, 1], srt = 0, adj = 1, xpd = TRUE, cex = 0.6)
