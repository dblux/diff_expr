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
subnetwork_genelist <- read.table("data/subnetwork/nea-hsa/ovarian_cancer/subnet_genelist-breast-metastasis.tsv",
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