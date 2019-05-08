library(dplyr)
library(ggplot2)
library(reshape2)
library(cowplot)
library(gridGraphics)
source("functions.R")
setwd("~/projects/phd/diff_expr/")

# Data --------------------------------------------------------------------
maqc <- read.table('data/MAQC-I/processed/rma_original.tsv',
                   header = T, row.names = 1)
# Substring colnames of MAQC
colnames(maqc) <- substring(colnames(maqc), 5)

# SUBSETS of MAQC
# All sample type: A (6 batches)
sample_a <- maqc[,grepl("A[0-9]", colnames(maqc))]
nrow(sample_a)
# All sample type: B (6 batches)
sample_b <- maqc[,grepl("B[0-9]", colnames(maqc))]

# QC: Class A ------------------------------------------------------------
# Evaluation: Normalisation
# Assumes that dataframe has been log-transformed
plot_evaluation <- function(df) {
  # Melt dataframe
  melt_df <- melt(df, variable.name = "ID")
  # Plot boxplot
  boxplot <- ggplot(melt_df, aes(x = ID, y = value, col = ID)) + 
              geom_boxplot(show.legend = F) +
              theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # Plot density curve
  pdf <- ggplot(melt_df, aes(x = value, col = ID)) + 
          geom_density(show.legend = T)
  # Total probe intensities for each chip
  mean_tibble <- melt_df %>% group_by(ID) %>%
                  summarise(mean = mean(value))
  # Scatter plot
  scatter <- ggplot(mean_tibble, aes(x = ID, y = mean, col = factor(rep(1:6, each = 5)))) +
              geom_point(show.legend = F) + 
              theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # Principal component analysis
  pca_obj <- prcomp(t(df), center = T, scale. = T)
  top_pc <- as.data.frame(pca_obj$x[,1:4])
  pc1_pc2 <- ggplot(top_pc, aes(x = PC1, y = PC2, col = factor(rep(1:6, each = 5)))) +
              geom_point(size = 3, show.legend = F)
  pc3_pc4 <- ggplot(top_pc, aes(x = PC3, y = PC4, col = factor(rep(1:6, each = 5)))) +
              geom_point(size = 3, show.legend = F)
  # Plot all graphs
  pca <- plot_grid(pc1_pc2, pc3_pc4)
  multiplot <- plot_grid(scatter, boxplot, pdf, pca,
                        nrow = 4)
  return(multiplot)  
}

before_norm <- plot_evaluation(sample_a)
save_plot("dump/a-before1.eps", before_norm,
          base_height = 12, base_width = 8)

# PCA plot
# Clustering
# RLE plot
# MA plot

# Normalisation -----------------------------------------------------------
qnorm_sample_a <- norm_quantile(sample_a)
after_norm <- plot_evaluation(qnorm_sample_a)
save_plot("dump/a-after.eps", after_norm,
          base_height = 12, base_width = 8)

gfs_sample <- GFS(sample_a)
gfs_norm <- plot_evaluation(gfs_sample)
save_plot("dump/norm_gfs.eps", gfs_norm,
          base_height = 12, base_width = 8)

# Remove probesets --------------------------------------------------------
# Filters out ambiguous and AFFY probesets
filter_probesets <- function(df) {
  logical_vec <- grepl("[0-9]_at", rownames(df)) & !startsWith(rownames(df), "A")
  print(paste0("No. of ambiguous and AFFY probesets removed: ",
               nrow(df) - sum(logical_vec)))
  return(df[logical_vec, , drop=F])  
}
fltr_qnorm <- filter_probesets(qnorm_sample_a)

# Map probesets to IDs ----------------------------------------------------
# Removes ambiguous probesets and probesets with no ID
# Selects maximum if two probesets match to same gene
# CHECK: What microarray platform is the data from?
ANNOT_FPATH <- "../info/microarray/HG-U133_Plus_2/annot_entrez-GPL570.tsv"
final_data <- affy2entrez(fltr_qnorm, ANNOT_FPATH)

# Checking probesets that are removed
rm_probeset <- final_data[[2]]
rm_plot <- plot_evaluation(rm_probeset)
save_plot("dump/a-remove.eps", rm_plot,
          base_height = 12, base_width = 8)

final_plot <- plot_evaluation(final_data)
save_plot("dump/a-final_plot.eps", final_plot,
          base_height = 12, base_width = 8)

# Save annotated file
OUTPUT_PATH <- "data/MAQC-I/processed/classB_entrez.tsv"
write.table(final_data, OUTPUT_PATH,
            sep = "\t", quote = F)
