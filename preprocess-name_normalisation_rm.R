# Initialisation ----------------------------------------------------------
library(dplyr)
library(ggplot2)
library(reshape2)
library(cowplot)
source("functions.R")
setwd("~/projects/phd/diff_expr/")

# Data --------------------------------------------------------------------
DATA_RPATH <- "data/ovarian_cancer/GSE26712/processed/mas5_original.tsv"
raw_data <- read.table(DATA_RPATH, header = T, row.names = 1)
colnames(raw_data)
# Sample name -------------------------------------------------------------
ANNOT_ID_RPATH <- "data/ovarian_cancer/GSE26712/README/annot.tsv"
annot_sample_id <- read.table(ANNOT_ID_RPATH,
                              sep = "\t", header = T,
                              row.names = 2, stringsAsFactors = F)
head(annot_sample_id, 5)
sample_id <- sapply(colnames(raw_data), function(x) annot_sample_id[x, 1])

# Mapped row names of data
unname(sample_id)

# Ordered row names of data
unname(sample_id)[c(66:75, 13:65)]
# Select and order raw data
selected_data <- raw_data[, c(66:75, 13:65)]
colnames(selected_data)

# Assigned short IDs
short_id <- c(paste0("A", 1:10), paste0("B", 1:185))
short_id
annot_accession <- data.frame(geo_accession = colnames(raw_data), short_id = short_id)
annot_accession
# Write annotation to file
ANNOT_ID_WPATH <- "data/ovarian_cancer/GSE26712/processed/annot-shortid.tsv"
write.table(annot_accession, ANNOT_ID_WPATH,
            sep = "\t", row.names = F, col.names = T, quote = F)

# Change colnames of data
colnames(raw_data) <- short_id

# Reorder columns ---------------------------------------------------------
# MAQC - Ref
sample_id <- colnames(raw_data)
order_index <- c(which(grepl("A[0-9]$", sample_id)), which(grepl("B[0-9]$", sample_id)))
ordered_data <- raw_data[, order_index]
colnames(ordered_data)

# Batch information -------------------------------------------------------
BATCH_INFO <- rep(1:2, c(10,185))
# BATCH_INFO <- rep(rep(1:6, each = 5), 2)
print(BATCH_INFO)

# Visualisation ----------------------------------------------------------
# Evaluation: Normalisation
# Assumes that dataframe has been log-transformed
plot_evaluation <- function(df, batch_info) {
  # Melt dataframe
  melt_df <- melt(df, variable.name = "ID")
  # Plot boxplot
  boxplot <- ggplot(melt_df, aes(x = ID, y = value, col = ID)) + 
              geom_boxplot(show.legend = F) +
              theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # Plot density curve
  pdf <- ggplot(melt_df, aes(x = value, col = ID)) + 
          geom_density(show.legend = T) +
          xlim(0, 16)
  # Total probe intensities for each chip
  mean_tibble <- melt_df %>% group_by(ID) %>%
                  summarise(mean = mean(value))
  # Scatter plot
  scatter <- ggplot(mean_tibble, aes(x = ID, y = mean, col = factor(batch_info))) +
              geom_point(show.legend = F) + 
              theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # Principal component analysis
  col_logical <- apply(t(df), 2, sum) != 0 & apply(t(df), 2, var) != 0
  pca_df <- t(df)[, col_logical]
  pca_obj <- prcomp(pca_df, center = T, scale. = T)
  top_pc <- as.data.frame(pca_obj$x[,1:4])
  pc1_pc2 <- ggplot(top_pc, aes(x = PC1, y = PC2, col = factor(batch_info))) +
              geom_point(size = 3, show.legend = F)
  pc3_pc4 <- ggplot(top_pc, aes(x = PC3, y = PC4, col = factor(batch_info))) +
              geom_point(size = 3, show.legend = F)
  # Plot all graphs
  pca <- plot_grid(pc1_pc2, pc3_pc4)
  multiplot <- plot_grid(scatter, boxplot, pdf, pca,
                        nrow = 4)
  return(multiplot)
}

before_norm <- plot_evaluation(log2_transform(raw_data), BATCH_INFO)
save_plot("dump/GSE26712-before.eps", before_norm,
          base_height = 12, base_width = 8)

# Normalisation -----------------------------------------------------------
colnames(raw_data)
classA <- raw_data[, 1:10]
classB <- raw_data[, 11:195]

norm_classA <- norm_quantile(classA)
norm_classB <- norm_quantile(classB)
norm_data <- cbind(norm_classA, norm_classB)
ncol(norm_data)

# DATA_WPATH <- "data/platinum_spike/GSE21344/processed/qnorm_mas5_010.tsv"
# write.table(norm_data, DATA_WPATH,
#             sep = "\t", quote = F, col.names = T, row.names = T)

# Remove probesets --------------------------------------------------------
# Filters out ambiguous and AFFY probesets from dataframe
# Rowname of affymetrix probesets
filter_probesets <- function(df) {
  logical_vec <- grepl("[0-9]_at", rownames(df)) & !startsWith(rownames(df), "AFFX")
  print(paste0("No. of ambiguous and AFFY probesets removed: ",
               nrow(df) - sum(logical_vec)))
  return(df[logical_vec, , drop=F])
}

# tail(rownames(platinum_spike))
fltr_data <- filter_probesets(norm_data)

# Map probesets to IDs ----------------------------------------------------
# Removes ambiguous probesets and probesets with no ID
# Selects maximum if two probesets match to same gene
# CHECK: What microarray platform is the data from?
ANNOT_PROBESET_RPATH <- "../info/microarray/HG-U133A/annot_entrez-GPL96.tsv"
annot_data <- affy2id(fltr_data, ANNOT_PROBESET_RPATH)

# # Filter for total set of spiked in cRNA
# LABEL_RPATH <- "data/platinum_spike/GSE21344/processed/probeset_foldchange_full.tsv"
# annot_label <- read.table(LABEL_RPATH, sep = "\t", header = F,
#                           row.names = 1, stringsAsFactors = F)
# # Total set of spiked in cRNA
# spikein_probeset <- rownames(annot_label)
# fltr_platinum <- norm_data[rownames(norm_data) %in% spikein_probeset,]

plot_fltr <- plot_evaluation(log2_transform(annot_data), BATCH_INFO)
save_plot("dump/GSE26712-after.eps", plot_fltr,
          base_height = 12, base_width = 8)

PROCESSED_DATA_WPATH <- sub("processed/.*$", "processed/mas5_qnorm.tsv", DATA_RPATH)
print(PROCESSED_DATA_WPATH)
write.table(annot_data, PROCESSED_DATA_WPATH,
            sep = "\t", quote = F, col.names = T, row.names = T)

# TODO: Normalise before filtering for probesets or vice versa???
