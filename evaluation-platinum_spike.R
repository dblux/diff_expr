# Initialise --------------------------------------------------------------
library(ggplot2)
library(cowplot)
library(reshape2)
library(samr)
# library(rgl)
setwd("~/projects/phd/diff_expr/")
source("functions.R")

# Import normalised platinum spike dataset
DATA_RPATH <- "data/platinum_spike/GSE21344/processed/qnorm_mas5.tsv"
qnorm_data <- read.table(DATA_RPATH, header = T,
                         row.names = 1, stringsAsFactors = F)
classA <- qnorm_data[, 1:9]
classB <- qnorm_data[, 10:18]

DATA_RPATH1 <- "data/platinum_spike/GSE21344/processed/mas5_no_norm.tsv"
qnorm_data1 <- read.table(DATA_RPATH1, header = T,
                         row.names = 1, stringsAsFactors = F)
classA1 <- qnorm_data1[, 1:9]
classB1 <- qnorm_data1[, 10:18]

LABEL_RPATH <- "data/platinum_spike/GSE21344/processed/ordered_foldchange_full.tsv"
# Import platinum spike labels
probeset_foldchange <- read.table(LABEL_RPATH,
                                  header = F, row.names = 1, stringsAsFactors = F)

# ordered_foldchange <- probeset_foldchange[order(rownames(probeset_foldchange)), , drop =][]][] F]
# write.table(ordered_foldchange, "data/platinum_spike/GSE21344/processed/ordered_foldchange_full.tsv",
#             sep = "\t", col.names = F, quote = F)

# Convert foldchange labels to vector
vec_foldchange <- probeset_foldchange[, 1]

names(vec_foldchange) <- rownames(probeset_foldchange)
# Creates binary labels with TRUE being differentially expressed
binary_label <- vec_foldchange != 1 & vec_foldchange != 0
sum(vec_foldchange == 1)

# Check that dataset is aligned with labels
print(paste("No. of mismatches between data and labels:",
            sum(!(rownames(qnorm_data) == names(vec_foldchange)))))

# Evaluation --------------------------------------------------------------
# Plots ROC and calculates AUC in a primitive fashion (i.e. ROC is step function)
# Does not resolve ties in the score
# Assumption: Lower score the better, ROC is step function
plot_roc <- function(score_list, label_vec,
                     is_bigger_better_vec = rep(F, length(score_list)),
                     name_vec = NULL) {
  # Function to plot a single ROC curve and calculate AUC
  ROC_AUC <- function(score_vec, is_bigger_better, color) {
    # Sort label vector according to score vector in ascending order
    sort_label <- label_vec[order(score_vec, decreasing = is_bigger_better)]
    # Dataframe of TPR and FPR
    df <- data.frame(TPR = cumsum(sort_label)/sum(sort_label),
                     FPR = cumsum(!sort_label)/sum(!sort_label))
    # Insert 0, 0 at the first row
    roc_df <- rbind(c(0,0), df)
    # Calculates change in FPR at each step
    dFPR <- c(0, diff(roc_df$FPR))
    # Sum of individual rectangle steps
    AUC <- sum(roc_df$TPR * dFPR)
    lines(roc_df$FPR, roc_df$TPR,
          col = color, lwd = 2)
    return(AUC)
  }
  # Initialise plot
  color_index <- (1:length(score_list)) + 1
  plot(NULL,
       xlab = "1 - Specificity", ylab = "Sensitivity",
       xlim = c(0,1), ylim = c(0,1),
       xaxs = "i", yaxs = "i")
  abline(0, 1, lty = 5)
  auc_vec <- mapply(ROC_AUC, score_list, is_bigger_better_vec,
                         color_index, SIMPLIFY = T)
  # If name_vec is not NULL display legend
  if (!is.null(name_vec)) {
    format_text <- function(name, auc) {
      return(sprintf("%s (%.3f)", name, auc))
    }
    legend_text <- mapply(format_text, name_vec, auc_vec)
    legend("bottomright", inset = 0.03, lty = 1, lwd = 2,
           legend = legend_text, col = color_index)
  }
  return(auc_vec)
}

# T-test
pvalue_ttest <- row_ttest(classA, classB)
predict_ttest <- pvalue_ttest <= 0.05
# Converts NA labels to F labels
predict_ttest[is.na(predict_ttest)] <- F
sum(predict_ttest)

pvalue_ttest <- row_ttest(classA, classB)
predict_ttest <- pvalue_ttest <= 0.05
# Converts NA labels to F labels
predict_ttest[is.na(predict_ttest)] <- F
sum(predict_ttest)

# Expression data of genes that were not predicted correctly
false_negative <- qnorm_data[!(predict_ttest == binary_label) & binary_label,]
false_positive <- qnorm_data[!(predict_ttest == binary_label) & !binary_label,]

# LogFC
# Depends on whether mean or median is used
foldchange_median <- calc_logfc(classA, classB, func = median)
abs_fc <- abs(foldchange_median)
# sum(is.infinite(abs_fc))
predict_logfc <- abs_fc >= log2(1.2)
# Converts NA to F labels
sum(predict_logfc)

plot(foldchange_median, -log10(pvalue_ttest))

# Evaluation
score_list <- list(pvalue_ttest, abs_fc)
roc_auc <- plot_roc(score_list, binary_label, 
                    c(F,T), c("T-test", "Log fold-change"))
plot_roc <- recordPlot()
save_fig(plot_roc, "dump/platinum-roc.eps", fig_height = 8)

eval_ttest <- evaluation_report(predict_ttest, binary_label)
eval_logfc <- evaluation_report(predict_logfc, binary_label)

abline(h = 0.87, v = 0.62)


# SAMR --------------------------------------------------------------------
# Data has to be log2 transformed before input
log_data <- log2_transform(qnorm_data)
# Create samr data file
samr_data <- list(x = data.matrix(qnorm_data),
                  y = rep(1:2, c(9,9)),
                  geneid = row.names(qnorm_data),
                  genenames = row.names(qnorm_data),
                  logged2 = T)
# Create samr data object
samr_obj <- samr(samr_data, resp.type = "Two class unpaired", nperms = 100)
# Delta represents the threshold distance in the qq-plot
delta_table <- samr.compute.delta.table(samr_obj)
# From the delta table, decide threshold threshold based on FDR
delta_table

delta <- 0.4
samr.plot(samr_obj, delta)
deg_table <- samr.compute.siggenes.table(samr_obj, delta, samr_data, delta_table)
str(deg_table)
up_genes <- deg_table$genes.up[,2]
down_genes <- deg_table$genes.lo[,2]
# DE genes identified by SAMR
deg_samr <- union(up_genes, down_genes)
# Prediction vector by SAMR
predict_samr <- rownames(qnorm_data) %in% deg_samr
evaluation_report(predict_samr, binary_label)

# Robust statistics -------------------------------------------------------
foldchange_mean <- calc_logfc(classA, classB, func = mean)
foldchange_median <- calc_logfc(classA, classB, func = median)

# Plot scatterplot
par(mfrow = c(2,1))
plot(log2(vec_foldchange), foldchange_mean)
abline(a = 0, b = 1)
plot(log2(vec_foldchange), foldchange_median)
abline(a = 0, b = 1)

# Plot violin plot
numeric_logfc <- as.numeric(levels(factor(logfc_label)))
factor_logfc <- factor(numeric_logfc)
ref_point <- data.frame(factor_logfc, numeric_logfc)

# Evaluate estimation of fold-change -------------------------------------
logfc_label <- log2(vec_foldchange)
median_logfc_factor <- data.frame(logfc_label, logfc_median)
violin_logfc <- ggplot(median_logfc_factor, aes(x = factor(logfc_label), y = logfc_median)) + 
                  geom_violin(aes(fill = factor(logfc_label))) +
                  geom_point(data = ref_point, aes(x = factor_logfc, y = numeric_logfc)) +
                  labs(y = "Measured log fold change (median)") +
                  theme(legend.position="none",
                        axis.title.x=element_blank(),
                        axis.text.x = element_text(angle = 90, hjust = 1))
violin_logfc
ggsave("dump/violin_logfc.eps", violin_logfc, width = 10, height = 10)

# MA Plot -----------------------------------------------------------------
calc_A <- function(df1, df2) {
  mean_vec1 <- apply(df1, 1, mean)
  mean_vec2 <- apply(df2, 1, mean)
  return(0.5*(log2(mean_vec1) + log2(mean_vec2)))
}

M <- log_fc(classA, classB)
A <- calc_A(classA, classB)

no_spikein <- vec_foldchange == 0
par(mfrow = c(1,1))
plot(A, M, col = 1 + no_spikein)

plot_df <- data.frame(A, M)
plot_nospike <- ggplot(plot_df[no_spikein,], aes(x = A, y = M)) +
                  geom_point(pch = 1) + stat_density2d() + 
                  geom_vline(xintercept = 3.5)

plot_spike <- ggplot(plot_df[!no_spikein,], aes(x = A, y = M)) +
  geom_point(pch = 1)

ggsave("dump/spike-scatter.eps", plot_spike)

# Background noise --------------------------------------------------------
plot_density <- function(df) {
  # Melt dataframe
  melt_df <- melt(df, variable.name = "ID")
  # Plot density curve
  density_plot <- ggplot(melt_df, aes(x = value, col = ID)) + 
    geom_density(show.legend = F)
  return(density_plot)
}

no_spikein <- vec_foldchange == 0

plot_spike_nospike <- function(df) {
  # Log transform data
  log_df <- log2_transform(df)
  nospike_df <- log_df[no_spikein,]
  spike_df <- log_df[!no_spikein,]
  density_nospike <- plot_density(nospike_df) + xlim(0, 7)
  density_spike <- plot_density(spike_df) + xlim(0, 16)
  density_compare <- plot_grid(density_nospike, density_spike)
  return(density_compare)
}

mas5_010 <- plot_spike_nospike(qnorm_data1)
mas5_010
save_plot("dump/mas5_010.eps", mas5_010, base_width = 10)

# Normalisation: Across genes ---------------------------------------------
### EXPERIMENT
no_spikein <- vec_foldchange == 0
log_df <- log2_transform(qnorm_data)
nospike_df <- log_df[no_spikein,]
spike_df <- log_df[!no_spikein,]

# No spike probesets showing signal
length(nospike_df[,1])
sum(nospike_df[,1] != 0)
sum(nospike_df[,1] != 0)/length(nospike_df[,1])

# Spike-in probesets showing no signal
length(spike_df[,1])
sum(spike_df[,1] == 0)
sum(spike_df[,1] == 0)/length(spike_df[,1])

# Histogram of probeset intensities for one sample
par(mfrow = c(1,2))
hist(nospike_df[,1], breaks = 30)
hist(spike_df[,1], breaks = 30)
plot_hist <- recordPlot()
save_fig(plot_hist, "dump/hist-mas5.eps",
         width = 10)

split_data <- function(df, vec_fc) {
  no_spikein <- vec_fc == 0
  nospike_df <- df[no_spikein,]
  spike_df <- df[!no_spikein,]
  spike_A <- spike_df[, 1:9]
  spike_B <- spike_df[, 10:18]
  nospike_A <- nospike_df[, 1:9]
  nospike_B <- nospike_df[, 10:18]
  return(list(spike_A, spike_B, nospike_A, nospike_B))
}

freq_nonzero <- function(df) {
  freq_nonzero <- table(rowSums(df != 0))
  graphics::barplot(freq_nonzero, space = 0)
  text(0.5:9.5, 1000, freq_nonzero)
  return(freq_nonzero)
}

par(mfrow=c(2,1))
data_ls <- split_data(qnorm_data, vec_foldchange)
list_004 <- data_ls[c(1,3)]
lapply(list_004, freq_nonzero)

data1_ls <- split_data(qnorm_data1, vec_foldchange)
list_010 <- data1_ls[c(1,3)]
lapply(list_010, freq_nonzero)

plot_010 <- recordPlot()
save_fig(plot_010, "dump/010-num_nonzero.eps",
         width = 6, height = 12)

### LOOK AT ONE PIECE OF INFO
spike_1 <- spike_classA[rowSums(spike_classA != 0) == 1,]
nospike_1 <- nospike_classA[rowSums(nospike_classA != 0) == 1,]
nospike_1_max <- apply(2^nospike_1, 1, max)
spike_1_max <- apply(2^spike_1, 1, max)

par(mfrow = c(1,2))
hist(spike_1_max, breaks = 15)
hist(nospike_1_max, breaks = 15)
plot_one_nonzero <- recordPlot()
save_fig(plot_one_nonzero, "dump/004-one_nonzero.eps")

length(nospike_1_max)
NOISE_THRESHOLD <- 20
sum(nospike_1_max < NOISE_THRESHOLD)/length(nospike_1_max)
sum(spike_1_max < NOISE_THRESHOLD)/length(spike_1_max)

# Inspect 1
spike_1 <- data1_ls[[1]][rowSums(data1_ls[[1]]!= 0) == 1,]
nospike_1 <- data1_ls[[3]][rowSums(data1_ls[[3]] != 0) == 1,]
nospike_1_max <- apply(nospike_1, 1, max)
spike_1_max <- apply(spike_1, 1, max)

hist(spike_1_max, breaks = 15)
hist(nospike_1_max, breaks = 15)

length(nospike_1_max)
NOISE_THRESHOLD <- 20
sum(nospike_1_max < NOISE_THRESHOLD)/length(nospike_1_max)
sum(spike_1_max < NOISE_THRESHOLD)/length(spike_1_max)

# Inspect 2
spike_2 <- data1_ls[[1]][rowSums(data1_ls[[1]]!= 0) == 2,]
nospike_2 <- data1_ls[[3]][rowSums(data1_ls[[3]] != 0) == 2,]
spike_4 <- data1_ls[[1]][rowSums(data1_ls[[1]]!= 0) == 4,]
nospike_4 <- data1_ls[[3]][rowSums(data1_ls[[3]] != 0) == 4,]
spike_7 <- data1_ls[[1]][rowSums(data1_ls[[1]]!= 0) == 7,]
nospike_7 <- data1_ls[[3]][rowSums(data1_ls[[3]] != 0) == 7,]
spike_8 <- data1_ls[[1]][rowSums(data1_ls[[1]]!= 0) == 8,]
nospike_8 <- data1_ls[[3]][rowSums(data1_ls[[3]] != 0) == 8,]
spike_9 <- data1_ls[[1]][rowSums(data1_ls[[1]]!= 0) == 9,]
nospike_9 <- data1_ls[[3]][rowSums(data1_ls[[3]] != 0) == 9,]
par(mfrow=c(2,2))
boxplot(t(spike_7), ylim = c(0,120), main = "Spike-in: 7")
boxplot(t(nospike_7), ylim = c(0,120), main = "No spike-in: 7")
boxplot(t(spike_8), ylim = c(0,120), main = "Spike-in: 8")
boxplot(t(nospike_8), ylim = c(0,120), main = "No spike-in: 8")

boxplot(t(spike_2), ylim = c(0,120), main = "Spike-in: 2")
boxplot(t(nospike_2), ylim = c(0,120), main = "No spike-in: 2")
boxplot(t(spike_4), ylim = c(0,120), main = "Spike-in: 4")
boxplot(t(nospike_4), ylim = c(0,120), main = "No spike-in: 4")

NOISE_THRESHOLD <- 40
sum(apply(spike_8, 1, mean) < NOISE_THRESHOLD)
sum(apply(nospike_8, 1, mean) < NOISE_THRESHOLD)
sum(apply(nospike_9, 1, mean) < NOISE_THRESHOLD)

# Mixture models
# Non-spikein probesets
index <- which.max(density(qnorm_nospike[,1])$y)
max_classA <- density(qnorm_nospike[,1])$x[index]
index <- which.max(density(qnorm_nospike[,18])$y)
max_classB <- density(qnorm_nospike[,18])$x[index]

# MOM estimator of shape and rate parameter
gamma_shape <- (mean(qnorm_nospike[,1])^2)/var(qnorm_nospike[,1])
gamma_rate <- mean(qnorm_nospike[,1])/var(qnorm_nospike[,1])

(657/5370)*((10975+659)/18707)/(5370/18707)
