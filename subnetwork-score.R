setwd("~/projects/phd/cs6216/")

# Naive row-wise two-sample t-test for every gene
# Does a t-test between every row of matrices a and b
# Returns a vector of p-values (length: nrow(a))
row_ttest <- function (a,b) {
  tt_pvalue <- numeric(nrow(a))
  names(tt_pvalue) <- rownames(a)
  for (i in 1:nrow(a)) {
    try(tt_pvalue[i] <- t.test(a[i,], b[i,])$p.value, silent = T)
  }
  return (tt_pvalue)
}

# Load data (179: No relapse, 107: Relapse)
data <- read.table('data/raw/data_entrez1.tsv', sep="\t",
                   header = T, row.names = 1, strip.white = T)
# Identify labels
y <- read.table('data/raw/metastasis_annot.tsv', sep="\t",
                   header = T, row.names = 1)
# Relabel columns as labels
colnames(data) <- sapply(colnames(data), function(x) y[x,])
# Reorder dataframe
reordered_index <- c(which(colnames(data) == 0), which(colnames(data) == 1))
reordered_data <- data[,reordered_index]

colnames(reordered_data) <- c(paste("0", 1:179, sep = "_"),
                              paste("1", 1:107, sep = "_"))

# Save re-ordered data
write.table(reordered_data, "data/raw/data_labelled.tsv",
            quote = F, sep = "\t", row.names = T, col.names = T)

# Load data (179: No relapse, 107: Relapse)
reordered_data <- read.table('data/raw/data_labelled.tsv', sep="\t",
                             header = T, row.names = 1, strip.white = T)

no_relapse <- reordered_data[,1:179]
relapse <- reordered_data[,180:286]

# Selects genes based on t-test values
p_value <- row_ttest(no_relapse, relapse)
# Select
NUM_GENES <- 20
selected_genes <- names(head(sort(p_value), NUM_GENES))
fltr_data <- reordered_data[rownames(reordered_data) %in% selected_genes,]
# heatmap(data.matrix(fltr_data))

# Select genes that are represented in the pathways
total_pathway_genes <- read.table("data/pathway/all_pathway_nodes.tsv",
                                  sep="\t", header = F, strip.white = T)
rep_genes <- unname(unlist(total_pathway_genes))
fltr_data <- reordered_data[rownames(reordered_data) %in% rep_genes,]

# Save filtered data
write.table(t(fltr_data), "data/gene/genes_1983.tsv",
            quote = F, sep = "\t", row.names = F, col.names = F)

# Save original data to numpy format
write.table(t(reordered_data), "data/gene/genes_all.tsv",
            quote = F, sep = "\t", row.names = F, col.names = F)

# PATHWAY SCORES ----------------------------------------------------------
GENELIST_DIR <- "data/pathway/gene_list/"
pathway_files <- list.files(GENELIST_DIR)

load_genelist <- function(fpath) {
  data <- read.table(paste0(GENELIST_DIR, fpath),
                     header = F, strip.white = T)
  return(unname(unlist(data)))
}

genelist_list <- sapply(pathway_files, load_genelist)
pathway_score <- function(gene_list) {
  fltr_data <- reordered_data[rownames(reordered_data) %in% gene_list,]
  # Return mean pathway value
  return(apply(fltr_data, 2, mean))
}
pathway_values <- sapply(genelist_list, pathway_score)
rownames(pathway_values)

# Save pathway data to numpy format
write.table(pathway_values, "data/pathway/X_pathway.tsv",
            quote = F, sep = "\t", row.names = F, col.names = F)

# EDA ---------------------------------------------------------------------
data <- read.table("data/raw/data_entrez1.tsv",
                   sep="\t", header=T, row.names=1, strip.white = T)

# Select genes with highest CV
row_mean <- apply(data, 1, mean)
row_sd <- apply(data, 1, sd)
row_cv <- row_sd/row_mean
ls_genes <- names(sort(row_cv)[1:50])
fltr_data <- data[rownames(data) %in% ls_genes,]

# Identify labels
y <- read.table('data/raw/metastasis_annot.tsv', sep="\t",
                header = T, row.names = 1)
# Relabel columns as labels
colnames(fltr_data) <- sapply(colnames(fltr_data), function(x) y[x,])

heatmap(data.matrix(fltr_data))
