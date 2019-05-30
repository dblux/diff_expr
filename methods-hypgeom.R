# Initialise --------------------------------------------------------------
library(VennDiagram)
# Import data -------------------------------------------------------------
# Import ovarian cancer data set 1
data1_rpath <- "data/ovarian_cancer/GSE18521/processed/mas5_qnorm.tsv"
ovarian_data1 <- read.table(data1_rpath, header = T, row.names = 1)
ovarian1_classA <- ovarian_data1[1:10]
ovarian1_classB <- ovarian_data1[11:63]

# Ovarian cancer data set 2
data2_rpath <- 'data/ovarian_cancer/GSE26712/processed/mas5_qnorm.tsv'
ovarian_data2 <- read.table(data2_rpath, header = T, row.names = 1)
ovarian2_classA <- ovarian_data2[1:10]
ovarian2_classB <- ovarian_data2[11:195]

# Load subnetwork gene list
geneset_rpath <- "data/subnetwork/nea-hsa/ovarian_cancer/geneset-ovarian.tsv"
geneset_df <- read.table(geneset_rpath, header = T)
print(head(geneset_df))
geneset_list <- split(geneset_df$gene, geneset_df$subnetwork_id)
print(head(geneset_list))
MIN_SIZE <- 5
# Check for number of subnetworks with size smaller than min size
sum(sapply(list_genelist, length) < MIN_SIZE)

# # Filter out subnetworks that fall below min size
# fltr_genelist <- list_genelist[(sapply(list_genelist, length) < MIN_SIZE)]
# fltr_genelist

# Hypergeometric test -----------------------------------------------------
# Significant genes for data set 1
pvalue_1 <- row_ttest(ovarian1_classA, ovarian1_classB)
sum(pvalue_1 <= 0.01, na.rm = T)
significant_genes1 <- rownames(ovarian1_classA)[pvalue_1 <= 0.01]

# Significant genes for data set 2
pvalue_2 <- row_ttest(ovarian2_classA, ovarian2_classB)
sum(pvalue_2 <= 0.01, na.rm = T)
significant_genes2 <- rownames(ovarian2_classA)[pvalue_2 <= 0.01]

# Evaluate reproducibility of genes
hypgeom_genes_plot <- plot_venn2(significant_genes1, significant_genes2)

# Additional arguments: num_genes_microarray
# Hypergeometric test for one subnetwork
hypgeom <- function(gene_list, significant_genes, num_genes_microarray) {
  num_significant_genes <- length(significant_genes)
  subnetwork_size <- length(gene_list)
  num_sig_subnetwork <- sum(gene_list %in% significant_genes)
  pvalue <- 1 - phyper(num_sig_subnetwork - 1,
                       num_significant_genes,
                       num_genes_microarray - num_significant_genes,
                       subnetwork_size)
  return(pvalue)
}

# Number of genes in the microarray
num_genes_microarray <- nrow(ovarian1_classA)

# Significant subnetworks in dataset 1
hypgeom_pvalue1 <- lapply(geneset_list, hypgeom, significant_genes1, num_genes_microarray)
sig_subnetworks3 <- names(hypgeom_pvalue1)[hypgeom_pvalue1 <= 0.05]

# Significant subnetworks in dataset 2
hypgeom_pvalue2 <- lapply(geneset_list, hypgeom, significant_genes2, num_genes_microarray)
sig_subnetworks4 <- names(hypgeom_pvalue2)[hypgeom_pvalue2 <= 0.05]

# Plot venn diagram
# Arguments: 2 vectors
plot_venn2 <- function(vec1, vec2) {
  # Generate overlap list
  overlap_list <- calculate.overlap(list(vec1,vec2))
  # Calculate venndiagram areas
  venn_area <- sapply(overlap_list, length)
  grid.newpage()
  venn_plot <- draw.pairwise.venn(venn_area[1], venn_area[2], venn_area[3],
                                  category = c("D1", "D2"),
                                  cex = 2, fontfamily = "sans",
                                  cat.cex = 2, cat.fontfamily = "sans",
                                  margin = 0.1)
  union <- (venn_area[1] + venn_area[2] - venn_area[3])
  print(unname(venn_area[3]/union))
  return(venn_plot)
}
hypgeom_venn <- plot_venn2(sig_subnetworks3, sig_subnetworks4)

# Plot venn
grid.newpage()
grid.draw(hypgeom_genes_plot)
venn <- recordPlot()
fpath <- "dump/hypgeom-genes-venn.png"
save_fig(venn, fpath, 600, 600)
