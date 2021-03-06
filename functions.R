#!/usr/bin/env Rscript
library(dplyr)

# Used in GFS function. Bins score with range [0,1] into intervals
# E.g. 4 Intervals: Binned into 0.2, 0.4, 0.6, 0.8
bin <- function(score, num_intervals) {
  for (i in 1:num_intervals) {
    if (score <= i/num_intervals) {
      return (i/(num_intervals+1))
    }
  }
}

# Problem: When there are too many zeros and fewer values are assigned to be 0
# .. No zeroes will be assigned
# Gene Fuzzy Scoring function transforms gene expression values
# Wilson Goh's paper
# Dense rank is used
norm_gfs <- function(A, upper=0.05, lower=0.15, num_intervals=0) {
  cat(sprintf("Top %.2f of expressed genes are assigned GFS scores of 1\n", upper))
  cat(sprintf("Genes below the top %.2f of expressed genes are assigned GFS scores of 0\n", lower))
  # Rank function ranks largest value as 1 [-A is used]
  # Handle NaN?
  ranked_A <- apply(-A, 2, dense_rank)
  rownames(ranked_A) <- rownames(A)
  # Returns [1,] = upper, [2,] = lower
  qtile <- apply(ranked_A, 2, quantile, probs=c(upper, lower), names=F)
  
  if (num_intervals <= 0) {
    for (c in 1:ncol(ranked_A)) {
      # Calculate qtile range
      q_range <- qtile[2,c] - qtile[1,c]
      for (r in 1:nrow(ranked_A)) {
        if (ranked_A[r,c] <= qtile[1,c]) {
          # Assign 1s
          ranked_A[r,c] <- 1
        } else if (ranked_A[r,c] > qtile[2,c]){
          # Assign 0s
          ranked_A[r,c] <- 0
        } else {
          # Assign score
          score <- (qtile[2,c] - ranked_A[r,c]) / q_range
          ranked_A[r,c] <- score
        }
      }
    }
  } else {
    # Discrete intervals
    for (c in 1:ncol(ranked_A)) {
      # Calculate qtile range
      q_range <- qtile[2,c] - qtile[1,c]
      for (r in 1:nrow(ranked_A)) {
        if (ranked_A[r,c] <= qtile[1,c]) {
          # Assign 1s
          ranked_A[r,c] <- 1
        } else if (ranked_A[r,c] > qtile[2,c]){
          # Assign 0s
          ranked_A[r,c] <- 0
        } else {
          # Assign score
          score <- (qtile[2,c] - ranked_A[r,c]) / q_range
          # Round off score
          ranked_A[r,c] <- bin(score, num_intervals)
        }
      }
    }
  }
  return (as.data.frame(ranked_A))
}

# Quantile normalisation
# Takes in df where columns are samples and rows are genes
norm_quantile <- function(df) {
  sort_arr <- apply(df, 2, sort)
  # Creates reference distribution
  ref_distr <- apply(sort_arr, 1, mean)
  rank_arr <- apply(df, 2, rank, ties.method = "random")
  qnorm_arr <- apply(rank_arr, c(1,2), function(x) ref_distr[x])
  rownames(qnorm_arr) <- rownames(df)
  qnorm_df <- as.data.frame(qnorm_arr)
  return(qnorm_df)
}

# Function that calculates the matrix of mean differences from the patient and control matrix
# Arguments: a <- larger matrix, b <- smaller matrix
# CLT: As n -> infinity, approximates more accurately a normal distribution.
# If underlying distribution is normal, then no need for n to be large
# Returns: mean_diff matrix with dimension of smaller matrix
old_calc_diff <- function(a, b) {
  a_colsum <- apply(a, 1, sum)
  a_size <- ncol(a)
  mean_diff <- array(numeric(), dim(b))
  rownames(mean_diff) <- rownames(b)
  
  for (r in 1:nrow(b)) {
    for (c in 1:ncol(b)) {
      mean_diff[r,c] <- (a_colsum[r]/a_size) - b[r,c]
    }
  }
  return (mean_diff)
}

# Function that calculates the matrix of mean differences from the patient and control matrix
# Arguments: a <- larger matrix, b <- smaller matrix
# CLT: As n -> infinity, approximates more accurately a normal distribution.
# If underlying distribution is normal, then no need for n to be large
# Returns: mean_diff matrix with dimension of smaller matrix
# pct= 0.75 smallest=4. pct=0.8 smallest=5
calc_diff <- function(a, b, sample_pct) {
  a_size <- ncol(a)
  sample_size <- floor(sample_pct * a_size)
  if (sample_size == a_size) {print("Sample size too small..")}
  
  mean_diff <- array(numeric(), dim(b))
  rownames(mean_diff) <- rownames(b)
  
  for (r in 1:nrow(b)) {
    for (c in 1:ncol(b)) {
      subsample <- unlist(sample(a[r,], sample_size, replace = F))
      mean_diff[r,c] <- mean(subsample) - b[r,c]
    }
  }
  return (mean_diff)
}

all_diff <- function(a, b) {
  all_diff <- array(numeric(), c(nrow(a), ncol(a)*ncol(b)))
  rownames(all_diff) <- rownames(a)
  for (r in 1:nrow(a)) {
    print(r)
    for (col_a in 1:ncol(a)) {
      for (col_b in 1:ncol(b)) {
        c <- col_b + (col_a - 1) * ncol(b)
        all_diff[r,c] <- b[r,col_b] - a[r,col_a]
      }
    }
  }
  return (all_diff)
}

# Calculate one-sample t-test statistic
# Using sample sd as estimator of population sd
# Use t-distribution to calculate p value
ttest_onesample <- function(vector, mu) {
  n <- length(vector)
  t_stat <- (mean(vector) - mu) / (sd(vector)/sqrt(n))
  if (is.infinite(t_stat)|is.na(t_stat)) {
    return (NaN)
  }
  if (t_stat < 0) {
    p_value <- pt(t_stat, n-1) * 2
  } else if (t_stat > 0) {
    p_value <- (1-pt(t_stat, n-1)) * 2
  } else {
    # t_stat == 0
    p_value <- 1
  }
  return (p_value)
}

# Naive row-wise two-sample t-test for every probe
# Does a t-test between every row of matrices a and b
# Returns a vector of p-values or tstats (length: nrow(a))
row_ttest <- function(df1, df2, flag = "pvalue") {
  ttest_vec <- numeric(nrow(df1))
  names(ttest_vec) <- rownames(df1)
  if (flag == "pvalue") {
    for (i in 1:nrow(df1)) {
      try(ttest_vec[i] <- t.test(df1[i,], df2[i,])$p.value, silent = T)
    }
    ttest_vec[is.na(ttest_vec)] <- 1
  } else if (flag == "tstat") {
    for (i in 1:nrow(df1)) {
      try(ttest_vec[i] <- t.test(df1[i,], df2[i,])$statistic, silent = T)
    }
    ttest_vec[is.na(ttest_vec)] <- 0
  } else {
    stop("Flag not in options pvalue or tstat..")
  }
  return(ttest_vec)
}

# Arguments: 2 dataframes that are not log-transformed
# Log-fold change (class1/class2)
calc_logfc <- function(df1, df2, prior_value = 0.1, func = mean) {
  vec1 <- apply(df1, 1, func)
  vec2 <- apply(df2, 1, func)
  # log2(0) = -Inf; log2(Inf) = -Inf; log2(0/0) = NaN
  # Reassigns 0s with prior_expr
  vec1[vec1 == 0] <- prior_value
  vec2[vec2 == 0] <- prior_value
  fc <- vec1/vec2
  # Reassigns NaN (0/0) with 1
  fc[is.nan(fc)] <- 1
  return(log2(fc))
}

# Argument: Recorded plot
# Save figure as file format indicated
save_fig <- function(recorded_plot, fpath, width = 8, height = 5) {
  if (endsWith(fpath, ".eps")) {
    setEPS()
    postscript(fpath, width = width, height = height)
    replayPlot(recorded_plot)
    dev.off()
  } else if (endsWith(fpath, ".pdf")) {
    pdf(fpath, width = width, height = height)
    replayPlot(recorded_plot)
    dev.off()
  } else if (endsWith(fpath, ".png")) {
    png(fpath, width = width, height = height)
    replayPlot(recorded_plot)
    dev.off()
  } else if (endsWith(fpath, ".jpg")) {
    jpeg(fpath, width = width, height = height)
    replayPlot(recorded_plot)
    dev.off()
  } else {
    stop("File extension not supported...")
  }
}

# Arguments: Dataframe, probeset annotation filepath
# Maps affy probesets to ID
# Removes ambiguous probesets and probesets with no ID
# Selects maximum if two probesets match to same gene
affy2id <- function(df, annot_fpath) {
  probeset_annot <- read.table(annot_fpath,
                               sep="\t", header=T, row.names=1,
                               stringsAsFactors=F, strip.white = T)
  # Filters out ambiguous and AFFY probesets from annot
  fltr_annot <- probeset_annot[grepl("[0-9]_at", rownames(probeset_annot))
                               & !startsWith(rownames(probeset_annot), "A"), , drop=F]
  # Returns entrez ID for all probe sets
  id <- unname(sapply(rownames(df), function(x) probeset_annot[x,]))

  # Indices of ambiguous probe sets and probe sets with no corresponding entrez ID to be deleted
  list_del <- which(grepl("///", id) | id == "")
  print(paste0("No. of probesets mapping to multiple IDs removed: ", sum(grepl("///", id))))
  print(paste0("No. of probesets with no ID removed: ", sum(id == "")))
  # Identifies genes that have multiple probesets mapping to it
  freq_gene <- table(id)
  dup_genes <- names(freq_gene[freq_gene > 1])
  for (i in dup_genes) {
    # Rows of dataframe with the same entrez ID
    same_rows <- df[id == i,]
    # Assign indices as rownames
    rownames(same_rows) <- which(id == i)
    # Rows that do not have the maximum sum are deleted
    row_del <- as.integer(rownames(same_rows[-which.max(apply(same_rows,1,sum)),]))
    # Concat with existing list of indices to be deleted
    list_del <- c(list_del, row_del)
  }
  # Rows are deleted
  df_genes <- df[-list_del,]
  fltr_id <- id[-list_del]
  # Assigning entrez ID to df
  rownames(df_genes) <- fltr_id
  # # CONCEPT CHECK: Deleted rows
  # df_genes_del <- df[list_del,]
  # entrez_del <- entrez[list_del]
  print(paste0("Total no. of probesets removed (incl. probesets mapping to same gene): ",
               length(list_del)))
  return(df_genes)
}

# Evaluates DE analysis method
# Arguments: 2 logical vectors - predicted and truth labels
# Returns: Accuracy, Sensitivity, Precision, Specificity, NPR
evaluation_report <- function(predict_vec, label_vec) {
  accuracy <- sum(predict_vec == label_vec)/length(predict_vec)
  TP <- sum(predict_vec & label_vec)
  TN <- sum(!predict_vec & !label_vec)
  sensitivity_1 <- TP/sum(label_vec)
  precision_1 <- TP/sum(predict_vec)
  specificity_1 <- TN/sum(!label_vec)
  negative_predictive_rate_1 <- TN/sum(!predict_vec)
  Sensitivity <- c(sensitivity_1, specificity_1)
  Precision <- c(precision_1, negative_predictive_rate_1)
  df <- cbind(Sensitivity, Precision)
  rownames(df) <- c("DE", "Not DE")
  print(df)
  print(paste("Accuracy:", accuracy))
  metrics <- c(accuracy, sensitivity_1, precision_1,
               specificity_1, negative_predictive_rate_1)
  names(metrics) <- c("Accuracy", "Sensitivity", "Precision",
                      "Specificity", "NPR")
  return(metrics)
}

# Log2 transforms data and handles -Inf values
log2_transform <- function(df) {
  log2_df <- log2(df)
  logical_df <- is.infinite(data.matrix(log2_df))
  log2_df[logical_df] <- 0
  return(log2_df)
}

# Trimmed mean scaling
# Non-log values
norm_scaling <- function(df, target_mean = 500, trim_percentage = 0.02) {
  trimmed_mean <- apply(df, 2, mean, trim = trim_percentage)
  scaling_factor <- target_mean/trimmed_mean
  return(mapply(function(a,b) a*b, df, scaling_factor))
}

# Substring without n last characters
substring_head <- function(string, n) {
  return(substring(string, 1, nchar(string) - n))
}

# Plot venn diagram
jacc_coeff <- function(vec1, vec2) {
  # Generate overlap list
  overlap_list <- calculate.overlap(list(vec1,vec2))
  # Calculate venndiagram areas
  venn_area <- sapply(overlap_list, length)
  grid.newpage()
  venn_plot <- draw.pairwise.venn(venn_area[1], venn_area[2], venn_area[3],
                                  category = c("D1", "D2"),
                                  cex = 3, fontfamily = "sans",
                                  cat.cex = 3, cat.fontfamily = "sans",
                                  margin = 0.1)
  union <- (venn_area[1] + venn_area[2] - venn_area[3])
  print(unname(venn_area[3]/union))
  return(venn_plot)
}

# Side effects: Writes KEGG pathway dataframes (Entrez)
kegg_df <- function(kegg_fpath) {
  pathway_id <- substring(kegg_fpath, 29, 36)
  wpath <- sprintf("../info/KEGG/hsa_df/%s.tsv", pathway_id)
  # Parses file into dataframe
  df_hsa <- parseKGML2DataFrame(kegg_fpath)
  # Convert KEGG ID to Entrez ID
  df_hsa$from <- substring(df_hsa$from, 5)
  df_hsa$to <- substring(df_hsa$to, 5)
  write.table(df_hsa, wpath,
              quote = F, sep = "\t", row.names = F)
}

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

# Returns: ggplot2 of PCA plot
plot_pca <- function(df, batch_info) {
  # Principal component analysis
  # Removes columns with all zeroes
  col_logical <- apply(t(df), 2, var) != 0
  pca_df <- t(df)[, col_logical]
  pca_obj <- prcomp(pca_df, center = T, scale. = T)
  top_pc <- as.data.frame(pca_obj$x[,1:2])
  pc1_pc2 <- ggplot(top_pc, aes(x = PC1, y = PC2, col = factor(batch_info))) +
    geom_point(size = 2, show.legend = F) +
    geom_label(label = rownames(pca_obj$x),
               nudge_x = 1, nudge_y = 2, size = 4,
               show.legend = F)
  return(pc1_pc2)
}

# Min-max scaling function
# Returns: Scaled vector with range [0,1]
norm_minmax <- function(vec) {(vec-min(vec))/(max(vec)-min(vec))}
