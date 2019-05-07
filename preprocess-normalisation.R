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

# All sample type: B (6 batches)
sample_b <- maqc[,grepl("B[0-9]", colnames(maqc))]

# QC: Class A ------------------------------------------------------------
# Evaluation: Normalisation
# Assumes that dataframe has been log-transformed
plot_evaluation <- function(df) {
  # Melt dataframe
  melt_df <- melt(df, variable.name = "ID")
  # Plot boxplot
  boxplot <- ggplot(melt_df, aes(x = ID, y = value, color = ID)) + 
              geom_boxplot(show.legend = F) +
              theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # Plot density curve
  pdf <- ggplot(melt_df, aes(x = value, color = ID)) + 
          geom_density(show.legend = T)
  # Total probe intensities for each chip
  mean_tibble <- melt_df %>% group_by(ID) %>%
                  summarise(mean = mean(value))
  # Scatter plot
  scatter <- ggplot(mean_tibble, aes(x = ID, y = mean, col = ID)) +
    geom_point(show.legend = F) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # Plot all graphs
  multiplot <- plot_grid(scatter, boxplot, pdf,
                        nrow = 3)
  return(multiplot)  
}

before_norm <- plot_evaluation(sample_b)
save_plot("dump/b-before.eps", before_norm,
          base_height = 10, base_width = 8)

# PCA plot
# Clustering
# RLE plot
# MA plot

# Normalisation -----------------------------------------------------------
qnorm_sample_b <- norm_quantile(sample_b)
after_norm <- plot_evaluation(qnorm_sample_b)
save_plot("dump/b-after.eps", after_norm,
          base_height = 10, base_width = 8)

gfs_sample <- GFS(sample_a)
gfs_norm <- plot_evaluation(gfs_sample)
save_plot("dump/norm_gfs.eps", gfs_norm,
          base_height = 10, base_width = 8)

# Remove probesets --------------------------------------------------------
# Filters out ambiguous and AFFY probesets
filter_probesets <- function(df) {
  return(df[grepl("[0-9]_at", rownames(df)) & !startsWith(rownames(df), "A"), , drop=F])  
}
fltr_qnorm <- filter_probesets(qnorm_sample_b)

# Map probesets to IDs ----------------------------------------------------
# Removes ambiguous probesets and probesets with no ID
# Selects maximum if two probesets match to same gene
# CHECK: What microarray platform is the data from?
ANNOT_FPATH <- "../info/microarray/HG-U133_Plus_2/annot_entrez-GPL570.tsv"
final_data <- affy2entrez(fltr_qnorm, ANNOT_FPATH)

final_plot <- plot_evaluation(final_data)
save_plot("dump/b-final_plot.eps", final_plot,
          base_height = 10, base_width = 8)

# Save annotated file
OUTPUT_PATH <- "data/MAQC-I/processed/classB_entrez.tsv"
write.table(final_data, OUTPUT_PATH,
            sep = "\t", quote = F)
