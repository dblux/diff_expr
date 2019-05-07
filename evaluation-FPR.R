source("functions.R")
setwd("~/projects/phd/diff_expr/")

# FALSE POSITIVE RATE -----------------------------------------------------
# MAQC data set
# Quantile normalised: Log2
# AFX_1_A1: Site 1 Sample A Replicate 1
# 6 different test sites
# A: UHRR, B: HBRR, C & D: Mixture, N: Normal, T: Tumour
# 5 technical replicates for each sample type
maqc <- read.table('data/MAQC-I/processed/rma_original.tsv',
                   header = T, row.names = 1)

# Split subset of class A samples into two
a1_col <- c(which(endsWith(colnames(a_maqc),"A1")),
            which(endsWith(colnames(a_maqc),"A2")),
            which(grepl("[1-3]_A3", colnames(a_maqc))))
a2_col <- setdiff(1:30, a1_col)
