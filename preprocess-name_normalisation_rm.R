library(dplyr)
library(ggplot2)
library(reshape2)
library(cowplot)
source("functions.R")
setwd("~/projects/phd/diff_expr/")

# Data --------------------------------------------------------------------
DATA_RPATH <- "data/platinum_spike/GSE21344/processed/mas5_original010.tsv"
raw_data <- read.table(DATA_RPATH, header = T, row.names = 1)
# Sample name -------------------------------------------------------------
ANNOT_ID_RPATH <- "data/platinum_spike/GSE21344/processed/annot-sample_name.tsv"
annot_sample_id <- read.table(ANNOT_ID_RPATH,
                              sep = "\t", header = T,
                              row.names = 1, stringsAsFactors = F)
head(annot_sample_id, 15)
sample_id <- sapply(colnames(raw_data), function(x) annot_sample_id[x, 1])
# Mapped row names of data
unname(sample_id)
# Ordered row names of data
unname(sample_id)[c(66:75, 13:65)]
# Select and order raw data
selected_data <- raw_data[, c(66:75, 13:65)]
colnames(selected_data)
# Assigned short IDs
short_id <- c(paste0("A", 1:10), paste0("B", 1:53))
annot_accession <- data.frame(colnames(selected_data), short_id)
# Write annotation to file
ANNOT_ID_WPATH <- "data/ovarian_cancer/GSE18521/processed/annot-accession.tsv"
write.table(annot_accession, ANNOT_ID_WPATH,
            sep = "\t", row.names = F, col.names = F, quote = F)

# Change colnames of data
# colnames(selected_data) <- short_id
colnames(raw_data) <- sample_id

# Batch/class information
BATCH_INFO <- rep(1:6, each = 3)
print(BATCH_INFO)
# QC: Class A ------------------------------------------------------------
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
  col_logical <- apply(t(df), 2, sum) != 0
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
save_plot("dump/platinum_mas-before.eps", before_norm,
          base_height = 12, base_width = 8)

# Normalisation -----------------------------------------------------------
colnames(selected_data)
norm_classA <- norm_quantile(fltr_platinum[, 1:9])
norm_classB <- norm_quantile(fltr_platinum[, 10:18])
norm_data <- cbind(norm_classA, norm_classB)
after_norm <- plot_evaluation(log2_transform(norm_data), BATCH_INFO)
save_plot("dump/platinum_mas-afternorm.eps", after_norm,
          base_height = 12, base_width = 8)

DATA_WPATH <- "data/platinum_spike/GSE21344/processed/qnorm_mas5_010.tsv"
write.table(norm_data, DATA_WPATH,
            sep = "\t", quote = F, col.names = T, row.names = T)

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
fltr_data <- filter_probesets(selected_data)

# Map probesets to IDs ----------------------------------------------------
# # Removes ambiguous probesets and probesets with no ID
# # Selects maximum if two probesets match to same gene
# # CHECK: What microarray platform is the data from?
# ANNOT_PROBESET_RPATH <- "../info/microarray/HG-U133_Plus_2/annot_entrez-GPL570.tsv"
# annot_data <- affy2id(fltr_data, ANNOT_PROBESET_RPATH)

# Filter for total set of spiked in cRNA
LABEL_RPATH <- "data/platinum_spike/GSE21344/processed/probeset_foldchange_full.tsv"
annot_label <- read.table(LABEL_RPATH, sep = "\t", header = F,
                          row.names = 1, stringsAsFactors = F)
# Total set of spiked in cRNA
spikein_probeset <- rownames(annot_label)
fltr_platinum <- raw_data[rownames(raw_data) %in% spikein_probeset,]
# rm_platinum <- raw_data[!(rownames(raw_data) %in% spikein_probeset),]

plot_fltr <- plot_evaluation(log2_transform(fltr_platinum), BATCH_INFO)
save_plot("dump/platinum_mas5-fltr.eps", plot_fltr,
          base_height = 12, base_width = 8)

write.table(annot_data, "data/ovarian_cancer/GSE18521/processed/qnorm_data.tsv",
            sep = "\t", quote = F, col.names = T, row.names = T)

# TODO: Normalise before filtering for probesets or vice versa???
