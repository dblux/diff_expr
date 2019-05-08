source("functions.R")
setwd("~/projects/phd/diff_expr/")

# Split subset of class A samples into two
split_adhoc <- function(df) {
  id <- colnames(df)
  col_fltr <- grepl("(1|2)$", id) | grepl("X[1-3]_.3", id)
  return(list(df[, col_fltr], df[, -col_fltr]))
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
a_ls <- adhoc_split(exp_classA)
a1 <- a_ls[[1]]
a2 <- a_ls[[2]]

raw_classB <- read.table('data/MAQC-I/processed/classB_entrez.tsv',
                         header = T, row.names = 1)
exp_classB <- 2^raw_classB
b_ls <- adhoc_split(exp_classB)
b1 <- b_ls[[1]]
b2 <- b_ls[[2]]

# T-test
pvalue_ttest <- row_ttest(a1, a2)
fp_ttest <- sum(pvalue_ttest <= 0.05)
FPR_ttest <- fp_ttest/length(pvalue_ttest)

# pvalue_ttest <- row_ttest(b1, b2)
# fp_ttest <- sum(pvalue_ttest <= 0.05)

# Log fold change
logfc <- log_fc(a1, a2)
fp_logfc <- sum(abs(logfc) >= 1)
FPR_logfc <- fp_logfc/length(logfc)

# logfc <- log_fc(b1, b2)
# fp_logfc <- sum(abs(logfc) >= 1)

# Volcano plot
plot(logfc, -log10(pvalue_ttest))

# Jaccard coefficient -----------------------------------------------------
VENN_FPATH <- "dump/venn.tiff"
jacc_coeff <- function(vec1, vec2) {
  # Generate overlap list
  overlap_list <- calculate.overlap(list(vec1,vec2))
  # Calculate venndiagram areas
  venn_area <- sapply(overlap_list, length)
  grid.newpage()
  venn_plot <- draw.pairwise.venn(venn_area[1], venn_area[2], venn_area[3],
                           category = c("S1", "S2"),
                           cex = 1, fontfamily = "sans",
                           cat.cex = 1, cat.fontfamily = "sans",
                           margin = 0.01)
  print(venn_plot)
  union <- (venn_area[1]+venn_area[2]-venn_area[3])
  return(unname(venn_area[3]/union))
}

pvalue_ttest1 <- row_ttest(a1, b1)
DEG1 <- rownames(a1)[pvalue_ttest1 <= 0.05]
pvalue_ttest2 <- row_ttest(a2, b2)
DEG2 <- rownames(a1)[pvalue_ttest2 <= 0.05]
jacc_ttest <- jacc_coeff(DEG1, DEG2)
