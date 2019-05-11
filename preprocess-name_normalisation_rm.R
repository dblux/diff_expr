library(dplyr)
library(ggplot2)
library(reshape2)
library(cowplot)
source("functions.R")
setwd("~/projects/phd/diff_expr/")

# Data --------------------------------------------------------------------
platinum_spike <- read.table("data/platinum_spike/GSE21344/processed/rma_original.tsv",
                              header = T, row.names = 1)
# Sample name -------------------------------------------------------------
annot_sample_id <- read.table("data/platinum_spike/GSE21344/processed/annot-sample_name.tsv",
                              sep = "\t", header = T, row.names = 1, stringsAsFactors = F)
sample_id <- sapply(colnames(platinum_spike), function(x) annot_sample_id[x, 1])
colnames(platinum_spike) <- sample_id

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
  scatter <- ggplot(mean_tibble, aes(x = ID, y = mean, col = factor(rep(1:6, each = 3)))) +
              geom_point(show.legend = F) + 
              theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # Principal component analysis
  pca_obj <- prcomp(t(df), center = T, scale. = T)
  top_pc <- as.data.frame(pca_obj$x[,1:4])
  pc1_pc2 <- ggplot(top_pc, aes(x = PC1, y = PC2, col = factor(rep(1:6, each = 3)))) +
              geom_point(size = 3, show.legend = F)
  pc3_pc4 <- ggplot(top_pc, aes(x = PC3, y = PC4, col = factor(rep(1:6, each = 3)))) +
              geom_point(size = 3, show.legend = F)
  # Plot all graphs
  pca <- plot_grid(pc1_pc2, pc3_pc4)
  multiplot <- plot_grid(scatter, boxplot, pdf, pca,
                        nrow = 4)
  return(multiplot)  
}

before_norm <- plot_evaluation(platinum_spike)
save_plot("dump/spike-before.eps", before_norm,
          base_height = 12, base_width = 8)

# Normalisation -----------------------------------------------------------
qnorm_sample_a <- norm_quantile(sample_a)
after_norm <- plot_evaluation(qnorm_sample_a)
save_plot("dump/a-after.eps", after_norm,
          base_height = 12, base_width = 8)

# Remove probesets --------------------------------------------------------
# Filters out ambiguous and AFFY probesets from dataframe
filter_probesets <- function(df) {
  logical_vec <- grepl("[0-9]_at", rownames(df)) & !startsWith(rownames(df), "AFFX")
  print(paste0("No. of ambiguous and AFFY probesets removed: ",
               nrow(df) - sum(logical_vec)))
  return(df[logical_vec, , drop=F])  
}

tail(rownames(platinum_spike))
fltr_platinum <- filter_probesets(platinum_spike)

# Map probesets to IDs ----------------------------------------------------
# Removes ambiguous probesets and probesets with no ID
# Selects maximum if two probesets match to same gene
# CHECK: What microarray platform is the data from?
ANNOT_FPATH <- "../info/microarray/Drosophila_G2/annot_flybase.tsv"
flybase_platinum <- affy2id(fltr_platinum, ANNOT_FPATH)

# Filter for total set of spiked in cRNA
LABEL_FPATH <- "data/platinum_spike/GSE21344/processed/flybase_foldchange.tsv"
annot_label <- read.table(LABEL_FPATH, sep = "\t", header = T,
                          row.names = 1, stringsAsFactors = F)
# Total set of spiked in cRNA
spikein_flybase <- rownames(annot_label)
flybase_platinum[rownames(flybase_platinum) %in% spikein_flybase,]

sink("dump/flybase_gn.txt")
cat(rownames(flybase_platinum), sep = "\n")
sink()

spikein_flybase

final_plot <- plot_evaluation(final_data)
save_plot("dump/a-final_plot.eps", final_plot,
          base_height = 12, base_width = 8)

# Save annotated file
OUTPUT_PATH <- "data/MAQC-I/processed/classB_entrez.tsv"
write.table(final_data, OUTPUT_PATH,
            sep = "\t", quote = F)
