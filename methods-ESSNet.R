# Initialisation ----------------------------------------------------------
library(parallel)
library(VennDiagram)
source("functions.R")
NUM_CORES <- detectCores()

# Constructs a matrix of all combinations
# Rows: Unique genes in the consolidated gene list
# Columns: All the possible combinations per gene (m*n)
construct_combi_arr <- function(geneset_genes, df1, df2, func = `-`) {
  # IMPT: If not indexing of rows will be wrong
  geneset_genes <- as.character(geneset_genes)
  # Driver function for each gene ID
  # For each gene_id, subset df row according to gene name (CHAR)
  # Outer product of two vectors to get all combinations
  # Arguments: Gene ID, function (e.g. subtract, divide)
  construct_combi_vec <- function(gene_id, df1, df2, func) {
    # Convert outer product matrix to vector
    as.vector(outer(unlist(df1[gene_id,]),
                    unlist(df2[gene_id,]),
                    func))
  }
  combi_arr <- t(simplify2array(mclapply(geneset_genes,
                                         construct_combi_vec,
                                         df1, df2, func,
                                         mc.cores = NUM_CORES
  )))
  rownames(combi_arr) <- geneset_genes
  return(combi_arr)
}

# Calculates t-statistic for ONE gene set
subnetwork_tstat <- function(gene_list, combi_arr) {
  # Calculate one-sample t-test statistic
  # Using sample sd as estimator of population sd
  calc_tstat <- function(vector, mu) {
    n <- length(vector)
    t_stat <- (mean(vector) - mu) / (sd(vector)/sqrt(n))
    return(ifelse (is.infinite(t_stat)|is.na(t_stat), NaN, t_stat))
  }
  # Construct subnetwork vector for each gene list
  subnetwork_vec <- as.vector(combi_arr[as.character(gene_list),])
  return(calc_tstat(subnetwork_vec, 0))
}

# Takes in original dfs
# Returns permutation null distribution for all subnetworks
# Rows: No. of subnetworks, Cols: No. of permutations
generate_permutation_null <- function(df1, df2, geneset_list, geneset_genes, max_perm) {
  # Combine both df
  df_stack <- cbind(df1, df2)
  n1 <- ncol(df1); n2 <- ncol(df2)
  num_perm <- choose(n1+n2, n1)
  cat(sprintf("No. of combinations (%d choose %d): %d\n",
              n1+n2, n1, num_perm))
  if (num_perm > max_perm) {
    # Monte Carlo permutation tests
    montecarlo_permutation <- function(df_stack) {
      permutation_indx <- sample(n1+n2, n1)
      df_p1 <- df_stack[, permutation_indx]
      df_p2 <- df_stack[, -permutation_indx]
      # Takes in original dfs
      perm_combi_arr <- construct_combi_arr(geneset_genes, df_p1, df_p2)
      perm_tstat_vec <- unlist(mclapply(geneset_list,
                                        subnetwork_tstat,
                                        perm_combi_arr,
                                        mc.cores = NUM_CORES))
      # Modifies global variable
      prog_indx <<- prog_indx + 1
      setTxtProgressBar(prog_bar, prog_indx)
      return(perm_tstat_vec)
    }
    cat(sprintf("Generating null distribution (%d Monte Carlo permutations):\n", max_perm))
    prog_indx <- 0
    # Starts progress bar
    prog_bar <- txtProgressBar(min = 0, max = max_perm,
                               width = 50, style = 3)
    montecarlo_null_arr <- replicate(max_perm, montecarlo_permutation(df_stack))
    close(prog_bar)
    return(montecarlo_null_arr)
    
  } else {
    # Standard permutation tests
    # Generates array with all possible combinations as columns
    combination_indx_arr <- combn(n1+n2, n1)
    # Returns t-statistic vector for SINGLE permutation
    standard_permutation <- function(permutation_indx, df_stack) {
      df_p1 <- df_stack[, permutation_indx]
      df_p2 <- df_stack[, -permutation_indx]
      perm_combi_arr <- construct_combi_arr(geneset_genes, df_p1, df_p2)
      perm_tstat_vec <- unlist(mclapply(geneset_list,
                                        subnetwork_tstat,
                                        perm_combi_arr,
                                        mc.cores = NUM_CORES))
      # Modifies global variable 
      prog_indx <<- prog_indx + 1
      setTxtProgressBar(prog_bar, prog_indx)
      return(perm_tstat_vec)
    }
    cat(sprintf("Generating null distribution (%d standard permutations):\n", num_perm))
    prog_indx <- 0
    # Starts progress bar
    prog_bar <- txtProgressBar(min = 0, max = num_perm,
                               width = 50, style = 3)
    permutation_null_arr <- apply(combination_indx_arr, 2, standard_permutation, df_stack)
    close(prog_bar)
    return(permutation_null_arr)
  }
}

# Arguments: Original t-statistic vector and null distr arr
# Calculate two-sided p-value from null distribution
empirical_pvalue <- function(tstat_vec, null_arr) {
  num_perm <- ncol(null_arr)
  # Driver function to calculate pvalue for SINGLE subnetwork
  calc_pvalue <- function(subnetwork_null_vec, tstat_scalar) {
    pvalue <- sum(abs(subnetwork_null_vec) >= abs(tstat_scalar))/num_perm
    return(pvalue)
  }
  # Tranposes such that columns are subnetworks for mapply
  null_df <- data.frame(t(null_arr))
  pvalue_vec <- mapply(calc_pvalue, null_df, tstat_vec)
  return(pvalue_vec)
}

# Arguments: T-statistic vector and combined sample size
# Returns theoretical p-value of two-sided paired t-test
theoretical_pvalue <- function(tstat_vec, n) {
  calc_pvalue <- function(t_stat, n) {
    # Handles infinity and NA values
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
    return(p_value)
  }
  return(sapply(tstat_vec, calc_pvalue, n))
}

# Arguments: List of gene sets, expression data (df1 and df2)
# Returns: List of pvalues for each gene set
essnet <- function(df1, df2, geneset_rpath, min_size = 5, max_perm = 1000) {
  geneset_df <- read.table(geneset_rpath, header = T)
  geneset_list <- split(geneset_df[,2], geneset_df[,1])
  print(head(geneset_list))
  # Check for number of subnetworks with size smaller than min size
  cat(sprintf("No. of subnetworks with size smaller than min. size: %d\n", 
              sum(sapply(geneset_list, length) < min_size)))
  # Remove subnetworks smaller than min_size if present
  if (sum(sapply(geneset_list, length) < min_size) > 0) {
    # Filter out subnetworks that fall below min size
    geneset_list <- geneset_list[(sapply(geneset_list, length) < min_size)]
  }
  geneset_genes <- sort(unique(geneset_df[,2]))
  # Name geneset_genes for functions
  names(geneset_genes) <- geneset_genes
  # IMPT: If not indexing of rows will be wrong
  geneset_genes <- as.character(geneset_genes)
  
  # Constructions original matrix of all combinations
  original_combi_arr <- construct_combi_arr(geneset_genes, df1, df2)
  # Calculates t-statistic for every gene set
  original_tstat_vec <- unlist(mclapply(geneset_list,
                                        subnetwork_tstat,
                                        original_combi_arr,
                                        mc.cores = NUM_CORES))
  permutation_null_df <- generate_permutation_null(df1, df2, geneset_list,
                                                   geneset_genes, max_perm)
  pvalue_vec <- empirical_pvalue(original_tstat_vec, permutation_null_df)
  return(pvalue_vec)
}

# Import data -------------------------------------------------------------
# Import ovarian cancer data set 1
data1_rpath <- "data/ovarian_cancer/GSE18521/processed/mas5_qnorm.tsv"
ovarian_data1 <- read.table(data1_rpath, header = T, row.names = 1)
ovarian1_classA <- ovarian_data1[1:5]
ovarian1_classB <- ovarian_data1[11:15]

# Import ovarian cancer data set 1
data1_rpath <- "data/ovarian_cancer/GSE18521/processed/mas5_qnorm.tsv"
ovarian_data1 <- read.table(data1_rpath, header = T, row.names = 1)
ovarian1_classA <- ovarian_data1[5:10]
ovarian1_classB <- ovarian_data1[43:49]

# Ovarian cancer data set 2
data2_rpath <- 'data/ovarian_cancer/GSE26712/processed/mas5_qnorm.tsv'
ovarian_data2 <- read.table(data2_rpath, header = T, row.names = 1)
ovarian2_classA <- ovarian_data2[1:5]
ovarian2_classB <- ovarian_data2[11:15]

geneset_rpath <- "data/subnetwork/nea-hsa/ovarian_cancer/geneset-ovarian.tsv"
pvalue1 <- essnet(ovarian1_classA, ovarian1_classB, geneset_rpath)
pvalue2 <- essnet(ovarian2_classA, ovarian2_classB, geneset_rpath)

sig_subnetworks1 <- names(pvalue1)[pvalue1 <= 0.05]
sig_subnetworks2 <- names(pvalue2)[pvalue2 <= 0.05]

length(sig_subnetworks1)

gold_std <- intersect(sig_subnetworks1, sig_subnetworks2)
# Save mapping from pathway names to ID
sink("dump/sig_intersection.txt", append = F)
cat(gold_std, fill = 1)
sink()

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
                                  cex = 3, fontfamily = "sans",
                                  cat.cex = 3, cat.fontfamily = "sans",
                                  margin = 0.1)
  union <- (venn_area[1] + venn_area[2] - venn_area[3])
  print(unname(venn_area[3]/union))
  return(venn_plot)
}

essnet_venn <- plot_venn2(sig_subnetworks1, sig_subnetworks2)
# Plot venn
grid.newpage()
grid.draw(essnet_venn)
venn <- recordPlot()
fpath <- "dump/essnet-venn.png"
save_fig(venn, fpath, 600, 600)

# Investigate -------------------------------------------------------------
# Import ovarian cancer data set 1
data1_rpath <- "data/ovarian_cancer/GSE18521/processed/mas5_qnorm.tsv"
ovarian_data1 <- read.table(data1_rpath, header = T, row.names = 1)
ovarian1_classA <- ovarian_data1[5:10]
ovarian1_classB <- ovarian_data1[43:49]

geneset_rpath <- "data/subnetwork/nea-hsa/ovarian_cancer/geneset-ovarian.tsv"
geneset_df <- read.table(geneset_rpath, header = T)
geneset_list <- split(geneset_df[,2], geneset_df[,1])
print(head(geneset_list))
# Check for number of subnetworks with size smaller than min size
cat(sprintf("No. of subnetworks with size smaller than min. size: %d\n", 
            sum(sapply(geneset_list, length) < 5)))
# Remove subnetworks smaller than min_size if present
if (sum(sapply(geneset_list, length) < 5) > 0) {
  # Filter out subnetworks that fall below min size
  geneset_list <- geneset_list[(sapply(geneset_list, length) < 5)]
}
geneset_genes <- sort(unique(geneset_df[,2]))
# Name geneset_genes for functions
names(geneset_genes) <- geneset_genes
# IMPT: If not indexing of rows will be wrong
geneset_genes <- as.character(geneset_genes)

# Constructions original matrix of all combinations
original_combi_arr <- construct_combi_arr(geneset_genes, ovarian1_classA, ovarian1_classB)
# Calculates t-statistic for every gene set
original_tstat_vec <- unlist(mclapply(geneset_list,
                                      subnetwork_tstat,
                                      original_combi_arr,
                                      mc.cores = NUM_CORES))
permutation_null_df <- generate_permutation_null(ovarian1_classA, ovarian1_classB, geneset_list,
                                                 geneset_genes, 1000)
pvalue_vec <- empirical_pvalue(original_tstat_vec, permutation_null_df)

theoretical_pvalue_vec <- theoretical_pvalue(original_tstat_vec, 20)

sum(theoretical_pvalue_vec < 0.05)

length(pvalue_vec)
sum(pvalue_vec <= 0.05)
names(pvalue_vec)[pvalue_vec < 0.05]

# No. of DE subnetworks (from same batch) = 2470

# Visualise null distributions of t-statistic
n <- 120
par(mfrow = c(2,2), mai = c(0.4, 0.4, 0.4, 0.2))
for (i in 1:n) {
  title <- sprintf("%s (%.4f)", rownames(permutation_null_df)[i], pvalue_vec[i])
  hist(permutation_null_df[i,], breaks = 20, main = title)
  # Plot t-statistic
  abline(v = c(-1,1) * original_tstat_vec[i], col = "red")
  if (i %in% c(4*1:n%/%4, n)) {
    fig <- recordPlot()
    fig_wpath <- sprintf("dump/null_distr%d.eps", i)
    save_fig(fig, fig_wpath, width = 10, height = 8)
  }
}

# Investigate DE subnetworks
gold_std <- read.table("dump/sig_intersection.txt", header = F)
gold_std <- unname(unlist(gold_std))
center_nodes <- substring(gold_std, 10)
dup_centernodes <- table(center_nodes)[table(center_nodes) > 1]
tail(sort(dup_centernodes, decreasing = T))

# length(unlist(gold_std)) = 1953
# sum(dup_centernodes) = 1493
# No. of genes represented by NEA-KEGG pathways = 2451

subnetworks_nea <- read.table("data/subnetwork/nea-hsa/ovarian_cancer/subnetworks-ovarian.tsv",
                              header = T, sep = "\t")
list_subnetworks <- split(subnetworks_nea[,2:3], subnetworks_nea[,1])
CENTER_STR <- "8897"
subset_subnetworks <- list_subnetworks[endsWith(names(list_subnetworks), CENTER_STR)]
subset_matrix <- lapply(subset_subnetworks, data.matrix)

# # Copy default Rgraphviz attributes
# graph_attr <- getDefaultAttrs()
# # Alter Rgraphviz attributes
# graph_attr$graph$bgcolor <- "white"
# graph_attr$node$fontsize <- 14

n <- length(subset_matrix)
par(mfrow = c(1,2), mai = c(0.4, 0.2, 0.4, 0.2))
for (i in 1:n) {
  # title <- sprintf("%s (%.4f)", rownames(null_AnB_arr)[i], geneset_AnB_pvalue[i])
  try(plot(ftM2graphNEL(subset_matrix[[i]]), attrs = graph_attr))
  # Plot t-statistic
  if (i %in% c(2*1:n%/%2, n)) {
    fig <- recordPlot()
    fig_wpath <- sprintf("dump/subnetwork_%s_%d.eps", CENTER_STR, i)
    save_fig(fig, fig_wpath, width = 10, height = 8)
  }
}

edges_subnetworks <- lapply(list_subnetworks, nrow)
gold_std_size <- unlist(edges_subnetworks[names(edges_subnetworks) %in% gold_std])
substd_size <- unlist(edges_subnetworks[!(names(edges_subnetworks) %in% gold_std)])
hist(gold_std_size, breaks = 30)
hist(substd_size, breaks = 30)

surv_data <- read.table("data/ovarian_cancer/GSE18521/README/survival.txt", sep = "\t",
                        stringsAsFactors = F)
raw_df <- t(surv_data)[14:66,]
rownames(raw_df) <- raw_df[,1]
survival_vec <- unlist(raw_df[,2])
survival_vec <- substring(survival_vec, 12)
survival_vec[endsWith(survival_vec, "(A)")] <- substring_head(survival_vec[endsWith(survival_vec, "(A)")], 4)
survival_vec






