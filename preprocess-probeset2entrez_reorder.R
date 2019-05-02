setwd("~/projects/phd/diff_expr/")
# DATA --------------------------------------------------------------------
DATA_PATH <- "data/ovarian_cancer/GSE26712/processed/filtered_original.tsv"
data_probeset <- read.table(DATA_PATH, header = T, row.names = 1)
# CHECK: What microarray platform is the data from?
# Imports probeset annotations
PROBESET_PATH <- "../info/microarray/HG-U133A/annot_entrez-GPL96.tsv"
probeset_annot <- read.table(PROBESET_PATH,
                             sep="\t", header=T, row.names=1,
                             stringsAsFactors=F, strip.white = T)
# Filters out ambiguous and AFFY probesets
fltr_annot <- probeset_annot[grepl("[0-9]_at", rownames(probeset_annot))
                             & !startsWith(rownames(probeset_annot), "A"), , drop=F]

METADATA_PATH <- "data/ovarian_cancer/GSE26712/README/annot.tsv"
metadata <- read.table(METADATA_PATH, header = T, row.names = 2)

# MAIN --------------------------------------------------------------------
# Returns entrez ID for all probe sets
entrez <- unname(sapply(rownames(data_probeset), function(x) probeset_annot[x,]))

# # First entrez ID selected for ambiguous probe sets
# correction <- sub(" ///.*$", "", entrez[grepl("///", entrez)])
# # Entrez ID to be substituted
# entrez[grepl("///", entrez)] <- correction

# Indices of ambiguous probe sets and probe sets with no corresponding entrez ID to be deleted
list_del <- which(grepl("///", entrez) | entrez == "")
# Identifies genes that have multiple probesets mapping to it
freq_gene <- table(entrez)
dup_genes <- names(freq_gene[freq_gene > 1])
for (i in dup_genes) {
  # Rows of dataframe with the same entrez ID
  same_rows <- data_probeset[entrez == i,]
  # Assign indices as rownames
  rownames(same_rows) <- which(entrez == i)
  # Rows that do not have the maximum sum are deleted
  row_del <- as.integer(rownames(same_rows[-which.max(apply(same_rows,1,sum)),]))
  # Concat with existing list of indices to be deleted
  list_del <- c(list_del, row_del)
}
# Rows are deleted
data_genes <- data_probeset[-list_del,]
fltr_entrez <- entrez[-list_del]
# Assigning entrez ID to df
rownames(data_genes) <- fltr_entrez

# Reorder with normal samples in front patients at back
col_names <- colnames(data_genes)
col_names1 <- sapply(col_names, function(x) metadata[x,])
# TODO: Decide order of samples
# reordered_data <- data_genes[c(66:75, 13:65)]
head(data_genes)
head(reordered_data)
# Save annotated file
OUTPUT_PATH <- "data/ovarian_cancer/GSE26712/processed/GSE26712_entrez.tsv"
write.table(data_genes, OUTPUT_PATH,
            sep = "\t", quote = F)
