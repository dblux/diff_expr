library(VennDiagram)
library(cowplot)
library(gridGraphics)
library(grid)
source("functions.R")
setwd("~/projects/phd/diff_expr/")

# Split subset of class A samples into two
split_adhoc <- function(df) {
  id <- colnames(df)
  col_fltr <- grepl("(1|2)$", id) | grepl("X[1-3]_.3", id)
  return(list(df[, col_fltr], df[, !col_fltr]))
}

# FALSE POSITIVE RATE -----------------------------------------------------
# MAQC data set
# Quantile normalised: Log2
# AFX_1_A1: Site 1 Sample A Replicate 1
# 6 different test sites
# A: UHRR, B: HBRR, C & D: Mixture, N: Normal, T: Tumour
# 5 technical replicates for each sample type
raw_classA <- read.table('data/MAQC-I/processed/classA_entrez.tsv',
                   header = T, row.names = 1)
exp_classA <- 2^raw_classA
a_ls <- split_adhoc(exp_classA)
a1 <- a_ls[[1]]
a2 <- a_ls[[2]]
a3 <- exp_classA[,1:15]
a4 <- exp_classA[,16:30]
a5 <- exp_classA[,1:5]
a6 <- exp_classA[,6:10]

raw_classB <- read.table('data/MAQC-I/processed/classB_entrez.tsv',
                         header = T, row.names = 1)
exp_classB <- 2^raw_classB
b_ls <- split_adhoc(exp_classB)
b1 <- b_ls[[1]]
b2 <- b_ls[[2]]
b3 <- exp_classB[,1:15]
b4 <- exp_classB[,16:30]
b5 <- exp_classB[,1:5]
b6 <- exp_classB[,6:10]

# T-test
fpr_ttest <- function(df1, df2) {
  pvalue_ttest <- row_ttest(df1, df2)
  fp_ttest <- sum(pvalue_ttest <= 0.05)
  FPR_ttest <- fp_ttest/length(pvalue_ttest)
  print(paste0(fp_ttest, "/", length(pvalue_ttest)))
  return(list(pvalue_ttest, FPR_ttest))
}

# Ad hoc split
fpr_ttest_A <- fpr_ttest(a1, a2)
print(fpr_ttest_A[[2]])
fpr_ttest_B <- fpr_ttest(b1, b2)
print(fpr_ttest_B[[2]])

# Split by batches
fpr_ttest_A_batch <- fpr_ttest(a3, a4)
print(fpr_ttest_A_batch[[2]])
fpr_ttest_B_batch <- fpr_ttest(b3, b4)
print(fpr_ttest_B_batch[[2]])

# Log fold change
fpr_logfc <- function(df1, df2) {
  logfc <- log_fc(df1, df2)
  fp_logfc <- sum(abs(logfc) >= 1)
  FPR_logfc <- fp_logfc/length(logfc)
  return(list(logfc, FPR_logfc))
}

# Ad hoc split
fpr_logfc_A <- fpr_logfc(a1, a2)
print(fpr_logfc_A[[2]])
fpr_logfc_B <- fpr_logfc(b1, b2)
print(fpr_logfc_B[[2]])

# Split by batches
fpr_logfc_A_batch <- fpr_logfc(a3, a4)
print(fpr_logfc_A_batch[[2]])
fpr_logfc_B_batch <- fpr_logfc(b3, b4)
print(fpr_logfc_B_batch[[2]])

# Volcano plot: FPR
par(mfrow = c(2,2))
plot(fpr_logfc_A[[1]], -log10(fpr_ttest_A[[1]]),
     col = 1 + (fpr_ttest_A[[1]] <= 0.05))

plot(fpr_logfc_A_batch[[1]], -log10(fpr_ttest_A_batch[[1]]),
     col = 1 + (fpr_ttest_A_batch[[1]] <= 0.05))

plot(fpr_logfc_B[[1]], -log10(fpr_ttest_B[[1]]),
     col = 1 + (fpr_ttest_B[[1]] <= 0.05))

plot(fpr_logfc_B_batch[[1]], -log10(fpr_ttest_B_batch[[1]]),
     col = 1 + (fpr_ttest_B_batch[[1]] <= 0.05))

volcano_plot <- recordPlot()
save_eps(volcano_plot, "dump/volcano.eps", fig_width = 8)

# Jaccard coefficient -----------------------------------------------------
jacc_coeff <- function(vec1, vec2) {
  # Generate overlap list
  overlap_list <- calculate.overlap(list(vec1,vec2))
  # Calculate venndiagram areas
  venn_area <- sapply(overlap_list, length)
  grid.newpage()
  venn_plot <- draw.pairwise.venn(venn_area[1], venn_area[2], venn_area[3],
                           category = c("D1", "D2"),
                           cex = 3, fontfamily = "sans",
                           cat.cex = 3, cat.fontfamily = "sans",
                           margin = 0.1)
  union <- (venn_area[1] + venn_area[2] - venn_area[3])
  print(unname(venn_area[3]/union))
  return(venn_plot)
}

# Ad hoc split
pvalue_ttest1 <- row_ttest(a1, b1)
DEG1 <- rownames(a1)[pvalue_ttest1 <= 0.05]
pvalue_ttest2 <- row_ttest(a2, b2)
DEG2 <- rownames(a1)[pvalue_ttest2 <= 0.05]
jacc_fig1 <- jacc_coeff(DEG1, DEG2)

# Split by batch
pvalue_ttest5 <- row_ttest(a3, b3)
DEG5 <- rownames(a1)[pvalue_ttest5 <= 0.05]
pvalue_ttest6 <- row_ttest(a4, b4)
DEG6 <- rownames(a1)[pvalue_ttest6 <= 0.05]
jacc_fig2 <- jacc_coeff(DEG5, DEG6)

pvalue_ttest7 <- row_ttest(a3, b4)
DEG7 <- rownames(a1)[pvalue_ttest7 <= 0.05]
pvalue_ttest8 <- row_ttest(a4, b3)
DEG8 <- rownames(a1)[pvalue_ttest8 <= 0.05]
jacc_fig3 <- jacc_coeff(DEG7, DEG8)

# Single batch
pvalue_ttest9 <- row_ttest(a5, b5)
DEG9 <- rownames(a1)[pvalue_ttest9 <= 0.05]
pvalue_ttest10 <- row_ttest(a6, b6)
DEG10 <- rownames(a1)[pvalue_ttest10 <= 0.05]
jacc_fig7 <- jacc_coeff(DEG9, DEG10)

pvalue_ttest11 <- row_ttest(a5, b6)
DEG11 <- rownames(a1)[pvalue_ttest11 <= 0.05]
pvalue_ttest12 <- row_ttest(a6, b5)
DEG12 <- rownames(a1)[pvalue_ttest12 <= 0.05]
jacc_fig8 <- jacc_coeff(DEG10, DEG11)

# LOGFC
# Ad hoc split
logfc_1 <- log_fc(a1, b1)
FCDEG1 <- rownames(a1)[abs(logfc_1) >= 1]
logfc_2 <- log_fc(a2, b2)
FCDEG2 <- rownames(a2)[abs(logfc_2) >= 1]
jacc_fig4 <- jacc_coeff(FCDEG1, FCDEG2)

# Split by batch
logfc5 <- log_fc(a3, b3)
FCDEG5 <- rownames(a1)[logfc5 <= 0.05]
logfc6 <- log_fc(a4, b4)
FCDEG6 <- rownames(a1)[logfc6 <= 0.05]
jacc_fig5 <- jacc_coeff(FCDEG5, FCDEG6)

logfc7 <- log_fc(a3, b4)
FCDEG7 <- rownames(a1)[logfc7 <= 0.05]
logfc8 <- log_fc(a4, b3)
FCDEG8 <- rownames(a1)[logfc8 <= 0.05]
jacc_fig6 <- jacc_coeff(FCDEG7, FCDEG8)

# Single batch
logfc9 <- log_fc(a5, b5)
FCDEG9 <- rownames(a1)[abs(logfc9) >= 1]
logfc10 <- log_fc(a6, b6)
FCDEG10 <- rownames(a1)[abs(logfc10) >= 1]
jacc_fig7 <- jacc_coeff(FCDEG9, FCDEG10)

logfc11 <- log_fc(a5, b6)
FCDEG11 <- rownames(a1)[abs(logfc11) >= 1]
logfc12 <- log_fc(a6, b5)
FCDEG12 <- rownames(a1)[abs(logfc12) >= 1]
jacc_fig8 <- jacc_coeff(FCDEG11, FCDEG12)

jacc_list <- list(jacc_fig1, jacc_fig2, jacc_fig3, jacc_fig4, jacc_fig5, jacc_fig6)
i <- 1
for(venn_plot in jacc_list) {
  grid.newpage()
  grid.draw(venn_plot)
  venn <- recordPlot()
  fpath <- sprintf("dump/venn%s.png", i)
  save_fig(venn, fpath, 600, 600)
  i <- i + 1
}

# Volcano plot: Many batch
par(mfrow = c(2,1))
plot(logfc5, -log10(pvalue_ttest5),
     col = 1 + (pvalue_ttest5 <= 0.05),
     main = "A3 vs. B3")
abline(v = c(-1,1))

plot(logfc7, -log10(pvalue_ttest7),
     col = 1 + (pvalue_ttest7 <= 0.05),
     main = "A3 vs. B4")
abline(v = c(-1,1))

DE_volcano <- recordPlot()
save_fig(DE_volcano, "dump/DE-volcano.eps",
         fig_height = 10)

# Volcano plot: Single batch
par(mfrow = c(2,1))
plot(logfc9, -log10(pvalue_ttest9),
     col = 1 + (pvalue_ttest9 <= 0.05),
     main = "A5 vs. B5")
abline(v = c(-1,1))

plot(logfc11, -log10(pvalue_ttest11),
     col = 1 + (pvalue_ttest11 <= 0.05),
     main = "A5 vs. B6")
abline(v = c(-1,1))

DE_single_volcano <- recordPlot()
save_fig(DE_single_volcano, "dump/DE_singlebatch-volcano.eps",
         fig_height = 10)
