# Load packages
library(dplyr)
source("functions.R")
# Set working directory
setwd("~/projects/phd/diff_expr")

# DATA --------------------------------------------------------------------
# Ovarian cancer data set 1
# Not log2
ds1 <- read.table('data/ovarian_cancer/GSE18521/processed/processed_original.tsv',
                  header = T, row.names = 1)
control_ds1 <- ds1[54:63]
tumour_ds1 <- ds1[1:53]

# Ovarian cancer data set 2
# Pre-processed using RMA - Log2 and quantile normalised
ds2 <- read.table('data/ovarian_cancer/GSE26712/processed/processed_original.tsv',
                  header = T, row.names = 1)
control_ds2 <- ds2[1:10]
tumour_ds2 <- ds2[11:195]

# MAQC data set
maqc <- read.table('data/MAQC-I/processed/filtered_original.tsv', header = T, row.names = 1)
colnames(maqc)
a_maqc <- maqc[,c(1:5,21:25,41:45,61:65,81:85,101:105)]
b_maqc <- maqc[,c(6:10,26:30,46:50,66:70,86:90,106:110)]

a1_col <- c(which(endsWith(colnames(a_maqc),"A1")),
            which(endsWith(colnames(a_maqc),"A2")),
            which(grepl("[1-3]_A3", colnames(a_maqc))))
a2_col <- setdiff(1:30, a1_col)

# Breast cancer data set 1 (n rows x m columns)
# n_control <- 9
# n_patient <- 33
# 
# control <- read.table('data/yeoh_2002/processed/processed_normal.tsv', header = T, row.names = 1)
# patient <- read.table('data/yeoh_2002/processed/processed_TEL-AML1.tsv', header = T, row.names = 1)

# control <- read.table('data/yeoh_2002/processed/filtered_normal.tsv', header = T, row.names = 1,
#                       colClasses = c("character", rep(c("numeric","NULL"), n_control)))
# patient <- read.table('data/yeoh_2002/processed/filtered_TEL-AML1.tsv', header = T, row.names = 1,
#                       colClasses = c("character", rep(c("numeric","NULL"), n_patient)))

# Rename column names
# colnames(control) <- paste0('C', 1:n_control)
# colnames(patient) <- paste0('P', 1:n_patient)
# MAIN --------------------------------------------------------------------

gfs_control_ds1 <- GFS(control_ds1)
mean_gfs_control_ds1 <- apply(gfs_control_ds1, 1, mean)

highexpr_probes <- names(mean_gfs_control_ds1[mean_gfs_control_ds1 > 0.5])
highexpr_genes <- unname(sapply(highexpr_probes, function(x) probeset_annot[x,]))
gene_list <- highexpr_genes[!sapply(highexpr_genes, function(x) grepl("///", x))]
dup_gene <- highexpr_genes[sapply(highexpr_genes, function(x) grepl("///", x))]
dup_gene1 <- unlist(strsplit(dup_gene, " /// "))
final_list <- c(gene_list, dup_gene1)

# Evaluate CLT-Net
mean_diff <- calc_diff(tumour_ds2[1:20,], control_ds2[1:20,], 0.5)
par(mfrow=c(1,2))
plot(density(data.matrix(control_ds2[8,])))
plot(density(mean_diff[8,]))

# Data set is log2.
# Try on ALL - AML
# Focus on subnetwork

n <- 13
diff <- vector()
for (i in tumour_ds2[n,]) {
  for (j in control_ds2[n,]) {
    diff <- c(diff, i - j)
    
  }
}
plot(density(diff))

y <- sample(tumour_ds2[1,],10, replace=F)

# Performs row-wise t-test on original expression data
ds1_tt_pvalue <- row_ttest(control_ds1, tumour_ds1)
ds2_tt_pvalue <- row_ttest(control_ds2, tumour_ds2)

# Returns list of diff expr gene symbol
diff_probes_tt_ds1 <- row.names(control_ds1)[ds1_tt_pvalue < 0.05 & !is.na(ds1_tt_pvalue)]
diff_probes_tt_ds2 <- row.names(control_ds2)[ds2_tt_pvalue < 0.05 & !is.na(ds2_tt_pvalue)]
# diff_gene_symbol <- sapply(diff_probes, function(x) as.character(probeset_annot[x,][1]))

# Evaluate reproducibility of row-wise t-test
rowtt_intersect <- length(intersect(diff_probes_tt_ds1, diff_probes_tt_ds2))
rowtt_union <- length(union(diff_probes_tt_ds1, diff_probes_tt_ds2))
rowtt_jacc <- rowtt_intersect/rowtt_union

# Performs CLT-Test
x_bar_ds1 <- calc_diff(control_ds1, tumour_ds1)
pvalue_ds1 <- apply(x_bar_ds1, 1, ttest_onesample, 0)
diff_probes_clt_ds1 <- names(pvalue_ds1)[!is.na(pvalue_ds1) & pvalue_ds1 < 0.05]

x_bar_ds2 <- calc_diff(control_ds2, tumour_ds2)
pvalue_ds2 <- apply(x_bar_ds2, 1, ttest_onesample, 0)
diff_probes_clt_ds2 <- names(pvalue_ds2)[!is.na(pvalue_ds2) & pvalue_ds2 < 0.05]

# Evaluate reproducibility of CLT-Test
clt_intersect <- length(intersect(diff_probes_clt_ds1, diff_probes_clt_ds2))
clt_union <- length(union(diff_probes_clt_ds1, diff_probes_clt_ds2))
clt_jacc <- clt_intersect/clt_union

# Print results
length(diff_probes_tt_ds1)
length(diff_probes_tt_ds2)
length(diff_probes_clt_ds1)
length(diff_probes_clt_ds2)
clt_jacc
rowtt_jacc

# GFS transform
gfs_control_ds1 <- GFS(control_ds1)
gfs_tumour_ds1 <- GFS(tumour_ds1)
gfs_control_ds2 <- GFS(control_ds2)
gfs_tumour_ds2 <- GFS(tumour_ds2)

# GFS: Row-wise t-test
gfs_ds1_tt_pvalue <- row_ttest(gfs_control_ds1, gfs_tumour_ds1)
gfs_ds2_tt_pvalue <- row_ttest(gfs_control_ds2, gfs_tumour_ds2)

# Identify diff expr genes
diff_probes_gfs_tt_ds1 <- row.names(gfs_control_ds1)[gfs_ds1_tt_pvalue < 0.05 & !is.na(gfs_ds1_tt_pvalue)]
diff_probes_gfs_tt_ds2 <- row.names(gfs_control_ds2)[gfs_ds2_tt_pvalue < 0.05 & !is.na(gfs_ds2_tt_pvalue)]
# gfs_tt_diff_gene_symbol <- sapply(gfs_tt_diff_probesets, function(x) as.character(probeset_annot[x,][1]))

# Evaluate reproducibility of row-wise t-test
gfs_rowtt_intersect <- length(intersect(diff_probes_gfs_tt_ds1, diff_probes_gfs_tt_ds2))
gfs_rowtt_union <- length(union(diff_probes_gfs_tt_ds1, diff_probes_gfs_tt_ds2))
gfs_rowtt_jacc <- gfs_rowtt_intersect/gfs_rowtt_union

# CLT-Net
gfs_x_bar_ds1 <- calc_diff(gfs_control_ds1, gfs_tumour_ds1)
gfs_pvalue_ds1 <- apply(gfs_x_bar_ds1, 1, ttest_onesample, 0)
diff_probes_gfs_clt_ds1 <- names(gfs_pvalue_ds1)[!is.na(gfs_pvalue_ds1) & gfs_pvalue_ds1 < 0.05]

gfs_x_bar_ds2 <- calc_diff(gfs_control_ds2, gfs_tumour_ds2)
gfs_pvalue_ds2 <- apply(gfs_x_bar_ds2, 1, ttest_onesample, 0)
diff_probes_gfs_clt_ds2 <- names(gfs_pvalue_ds2)[!is.na(gfs_pvalue_ds2) & gfs_pvalue_ds2 < 0.05]

# Evaluate reproducibility of CLT-Test
gfs_clt_intersect <- length(intersect(diff_probes_gfs_clt_ds1, diff_probes_gfs_clt_ds2))
gfs_clt_union <- length(union(diff_probes_gfs_clt_ds1, diff_probes_gfs_clt_ds2))
gfs_clt_jacc <- gfs_clt_intersect/gfs_clt_union

# Print results
length(diff_probes_gfs_tt_ds1)
length(diff_probes_gfs_tt_ds2)
length(diff_probes_gfs_clt_ds1)
length(diff_probes_gfs_clt_ds2)
gfs_rowtt_jacc
gfs_clt_jacc
# Datasets are of different sizes!!!

df1 <- data.frame(a=1:10, b=21:30)
df2 <- cbind(df1,df1)

x <- all_diff(control_ds2, tumour_ds2)

y <- vector()
for (i in control_ds2[1,]) {
  for (j in tumour_ds2[1,]) {
    a <- j - i
    y <- c(y,a)
  }
}

# Bootstrap idea
z <- all_diff(cds2[1:100,], tds2[1:100,])
plot(density(z[5,]))
bootstrap <- replicate(1000, mean(sample(z[2,], 30, replace=T)))
plot(density(bootstrap))

ad <- all_diff(control_ds1, tumour_ds1)

# Distribution of probe intensities for a single gene
cds2 <- 2^control_ds2
tds2 <- 2^tumour_ds2

# Visualisation of x_bar
view_xbar <- function(a) {
  check <- list(gfs_patient, gfs_control, x_bar)
  print(lapply(check, '[', a,))
  hist(x_bar[a,], breaks=10)
}

gene_diff <- apply(x_bar, 1, mean)
gene_sd <- apply(x_bar, 1, sd)

gene_diff[1]
gene_sd[1]

# GFS Demonstration
X <- array(floor(rgamma(500, 9, 0.5)*10), c(100,5))
X
Y <- GFS(X, 0.05, 0.50, 4)
Y
hist(Y[,1], breaks = 30)

sum(Y[,4] == 0)

# VISUALISATION OF DATA

# TO-DO

# Dataset split into two / two datasets of the same disease (jacc)
# Create golden standard using intersection of DE genes from the two datasets
# Check in generation of subnetworks
# Evaluate normalisation techniques (e.g. discretised GFS) comparing the s.d. of GFS across samples against the s.d. of raw values, p-value, and log fold change
# No need to remove ambiguous probesets for one-one genes!!!
# How to evaluate?

# Simulated data
# Reproducibility
# Batch effects

# regularized t-test SAM: microarray. For one sample to one sample

# break up pathway into smaller. control or patient
# check whether One-net is equivalent to existing t-test

# ESSNet rotation test

# HYPGEOM -----------------------------------------------------------------
# Import ovarian cancer data set 1
# Not log2
ovarian_data1 <- read.table("data/ovarian_cancer/GSE18521/processed/GSE18521_entrez1.tsv",
                              header = T, row.names = 1)
control_data1 <- ovarian_data1[1:10]
tumour_data1 <- ovarian_data1[11:62]

# Ovarian cancer data set 2
# Pre-processed using RMA - Log2 and quantile normalised
ovarian_data2 <- read.table('data/ovarian_cancer/GSE26712/processed/GSE26712_entrez.tsv',
                            header = T, row.names = 1)
# Reverse the log2
ovarian_data2a <- 2^ovarian_data2
control_data2 <- ovarian_data2a[1:10]
tumour_data2 <- ovarian_data2a[11:195]

# Load subnetwork gene list
subnetwork_genelist <- read.table("data/subnetwork/hsa-nea/ovarian_cancer/subnet_genelist-breast-metastasis.tsv",
                                  header = T)
list_genelist <- split(subnetwork_genelist$gene, subnetwork_genelist$subnetwork)
MIN_SIZE <- 5
# Number of subnetworks with size smaller than min size
num_subnetwork <- sum(sapply(list_genelist, length) < MIN_SIZE)

# Significant genes for data set 1
pvalue_1 <- row_ttest(control_data1, tumour_data1)
significant_genes1 <- rownames(control_data1)[pvalue_1 <= 0.01]
# Significant genes for data set 2
pvalue_2 <- row_ttest(control_data2, tumour_data2)
significant_genes2 <- rownames(control_data2)[pvalue_2 <= 0.01]

num_genes_microarray <- nrow(control_data1)
# Additional arguments: num_genes_microarray
# Hypergeometric test
hypgeom <- function(gene_list, significant_genes) {
  num_significant_genes <- length(significant_genes)
  num_genes_subnetwork <- length(gene_list)
  num_sig_subnetwork <- sum(gene_list %in% significant_genes)
  pvalue <- 1 - phyper(num_sig_subnetwork - 1,
                       num_significant_genes,
                       num_genes_microarray - num_significant_genes,
                       num_genes_subnetwork)
  return(pvalue)
}
# Significant subnetworks in dataset 1
hypgeom_pvalue1 <- lapply(list_genelist, hypgeom, significant_genes1)
sig_subnetworks1 <- names(hypgeom_pvalue1)[hypgeom_pvalue1 <= 0.05]
# Significant subnetworks in dataset 2
hypgeom_pvalue2 <- lapply(list_genelist, hypgeom, significant_genes2)
sig_subnetworks2 <- names(hypgeom_pvalue2)[hypgeom_pvalue1 <= 0.05]


# Import breast metastasis data
metastasis_data <- read.table("data/breast_metastasis/GSE2034/processed/data_labelled.tsv",
                             header = T, row.names = 1)
# Import subnetworks
subnetworks <- read.table("data/subnetwork/hsa-nea/subnetworks-breast_metastasis.tsv",
                          header = T)
# Create gene list for each subnetwork
list_subnetworks <- split(subnetworks[2:3], subnetworks$subnetwork_name)
subnetwork_genelist <- lapply(list_subnetworks, function(x) sort(unique(unname(unlist(x)))))

# Save subnetwork gene list
subnetwork_matrix <- lapply(subnetwork_genelist, as.data.frame)
subnetwork_df <- do.call(rbind, subnetwork_matrix)
subnetwork_df["subnetwork"] <- sub("\\..*$", "", rownames(subnetwork_df))
colnames(subnetwork_df)[1] <- "gene"
subnetwork_df1 <- subnetwork_df[,2:1]
rownames(subnetwork_df1) <- NULL
write.table(subnetwork_df1,
            "data/subnetwork/hsa-nea/subnet_genelist-breast-metastasis.tsv",
            col.names = T, row.names = F, quote = F)
