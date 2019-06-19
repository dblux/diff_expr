# Initialisation ----------------------------------------------------------
library(igraph)
library(parallel)
library(KEGGgraph)
source("functions.R")
NUM_CORES <- detectCores()
setwd("~/projects/phd/diff_expr/")
# Rgraphviz is built on graph package: graphNEL format

# Import pathwayAPI as dataframe and removes all duplicates and loops
# Dataframe has three columns: pathway_id, from, to
# Returns: NAMED list of pathwayAPI pathway dataframes
import_df_pathway <- function(df_rpath) {
  pathway_df <- read.table(df_rpath, header=T,
                           sep="\t", stringsAsFactors = F)
  print(head(pathway_df))
  list_df <- split(pathway_df[,2:3], pathway_df[,1])
  
  # Removes duplicate edges and loops in df
  remove_dupl_loops <- function(df) {
    # Removes loops in graph
    processed_df <- df[!(df[,1] == df[,2]),]
    # Remove duplicated edges
    concat_str <- paste0(processed_df[,1], processed_df[,2])
    freq <- table(concat_str)
    dup_str <- names(freq[freq > 1])
    if (length(dup_str) == 1) {
      row_no <- which(concat_str == dup_str)[-1]
      # Returns df with row removed
      processed_df <- processed_df[-row_no,]
    } else if (length(dup_str) > 1) {
      all_rows <- vector()
      for (i in dup_str) {
        row_no <- which(concat_str == i)[-1]
        all_rows <- c(all_rows, row_no)
      }
      # Returns df with row removed
      processed_df <- processed_df[-all_rows,]
    }
    return(processed_df)
  }
  list_fltr_df <- lapply(list_df, remove_dupl_loops)
  nrow_removed <- sum(sapply(list_df, nrow)) - sum(sapply(list_fltr_df, nrow))
  cat(sprintf("Total no. of duplicated edges or loops removed: %s\n", nrow_removed))
  return(list_fltr_df)
}

# Arguments: KEGG df directory path
# Returns: NAMED list of KEGG pathway dataframes
# Duplicate edges and loops removed
import_kegg <- function(kegg_rdirpath) {
  kegg_df_rpath <- list.files(substring_head(kegg_rdirpath, 1), full.names = T)
  names(kegg_df_rpath) <- sapply(kegg_df_rpath, substring, 21, 28)
  # Returns: Processed df of SINGLE KEGG pathway
  process_kegg <- function(kegg_fpath) {
    
    df_hsa <- read.table(kegg_fpath, header = T, sep = "\t")
    # If dataframe is not NULL
    if (nrow(df_hsa) != 0) {
      # Removes loops in graph
      processed_df <- df_hsa[!(df_hsa$from == df_hsa$to), 1:2]
      # Remove duplicated edges
      concat_str <- with(processed_df, paste0(from, to))
      freq <- table(concat_str)
      dup_str <- names(freq[freq > 1])
      if (length(dup_str) == 1) {
        row_no <- which(concat_str == dup_str)[-1]
        # Returns df with row removed
        processed_df <- processed_df[-row_no,]
      } else if (length(dup_str) > 1) {
        all_rows <- vector()
        for (i in dup_str) {
          row_no <- which(concat_str == i)[-1]
          all_rows <- c(all_rows, row_no)
        }
        # Returns df with row removed
        processed_df <- processed_df[-all_rows,]
      }
      return(processed_df)
    }
  }
  list_raw <- lapply(kegg_df_rpath, process_kegg)
  # Removes NULL values in list
  list_kegg_df <- list_raw[!sapply(list_raw, is.null)]
  return(list_kegg_df)  
}

# Arguments: List of pathway dfs, list of highexpr genes, subnetwork name
# Returns list of gene sets
# Side effects: Saves edgelist and geneset tsv
generate_geneset <- function(list_pathway_df, highexpr_genes, dataset_name,
                             pathwaydb_name, class_id, min_size) {
  # Arguments: dataframe of one pathway, list of highexpr genes and min size
  # Returns a list of subnetwork edgelist df. Subnetworks are above min size
  subnetwork_pfsnet <- function(df, gene_list, min_size) {
    list_subnetwork <- list()
    # Only retain edges where both nodes are in list of genes
    fltr_df <- df[df$from %in% gene_list & df$to %in% gene_list,]
    # Skip df with too few rows to form subnetwork of minimum size
    # Does not ensure that subnetworks are above minimum size
    if (nrow(fltr_df) >= (min_size- 1)) {
      # Convert numeric matrix to string matrix in order to avoid vertex id
      char_array <- apply(data.matrix(fltr_df), c(1, 2), as.character)
      sub_graphs <- graph_from_edgelist(char_array)
      # Only returns graphs that have at least minimum size
      connected_list <- decompose.graph(sub_graphs, min.vertices = min_size)
      for (graph in connected_list) {
        subnetwork <- get.edgelist(graph)
        list_subnetwork <- append(list_subnetwork, list(subnetwork))
      }
    }
    # Only returns list if it is not empty
    if (length(list_subnetwork) > 0) {
      return(list_subnetwork)
    }
  }
  # Generates list of subnetwork edgelist df
  subnetworks_list_list <- lapply(list_pathway_df, subnetwork_pfsnet, highexpr_genes, min_size)
  subnetworks_list_list <- subnetworks_list_list[!sapply(subnetworks_list_list, is.null)]
  
  # Names the pathways and unlists and stacks subnetworks into a df
  vstack_df <- function(subnetwork_list, pathway_id) {
    # Expect max of 99 subnetworks for each pathway
    num_id <- sprintf("%.2d", 1:length(subnetwork_list))
    # length(subnetwork_names): No. of subnetworks in each pathway
    subnetwork_names <- paste(pathway_id, num_id, sep = "_")
    num_row <- lapply(subnetwork_list, nrow)
    subnetwork_id <- rep(subnetwork_names, num_row)
    subnetwork_df <- do.call(rbind, subnetwork_list)
    return(cbind(subnetwork_id, subnetwork_df))
  }
  named_subnetworks <- mapply(vstack_df, subnetworks_list_list, names(subnetworks_list_list))
  # Stacked edgelist df of all the subnetworks
  edgelist_vstack <- do.call(rbind, named_subnetworks)
  colnames(edgelist_vstack)[2:3] <- c("from", "to")
  # Add dataset_name and class_id to pathway id
  edgelist_vstack[,1] <- sprintf("%s-%s_%s", edgelist_vstack[,1],
                                 dataset_name, class_id)
  
  # Save subnetwork edgelist df
  edgelist_wpath <- sprintf("data/subnetwork/pfsnet/edgelist-%s_%s_%s.tsv",
                            pathwaydb_name, dataset_name, class_id)
  write.table(edgelist_vstack, edgelist_wpath,
              sep = "\t", col.names = T, row.names = F, quote = F)
  
  # Create gene list for each subnetwork split on COLUMN
  subnetworks_list <- split(edgelist_vstack[,2:3], edgelist_vstack[,1])
  # List of genes in each subnetwork
  geneset_list <- lapply(subnetworks_list, function(x) sort(unique(x)))
  
  # Convert list of genes to df in order to save
  subnetwork_matrix <- lapply(geneset_list, as.data.frame)
  genelist_df <- do.call(rbind, subnetwork_matrix)
  genelist_df$subnetwork_id <- sub("\\..*$", "", rownames(genelist_df))
  colnames(genelist_df)[1] <- "gene_id"
  rownames(genelist_df) <- NULL
  geneset_df <- genelist_df[,2:1]
  geneset_wpath <- sprintf("data/subnetwork/pfsnet/geneset-%s_%s_%s.tsv",
                           pathwaydb_name, dataset_name, class_id)
  write.table(geneset_df, geneset_wpath,
              sep = "\t", col.names = T, row.names = F, quote = F)
  
  return(geneset_list)
}

# Arguments: List of genes for each subnetwork, class_id of subnetwork
# Returns: T-statistic scalar
calc_subnetwork_tstat <- function(gene_list, class_id, gene_weight_df, df_AnB) {
  gene_list <- as.character(gene_list)
  subnetwork_data <- df_AnB[gene_list,]
  subnetwork_weight <- gene_weight_df[gene_list,]
  indx1 <- class_id
  # Assign opposite index based on class_id
  indx2 <- ifelse(indx1 == "A", "B", "A")
  # Calculate paired difference of subnetwork scores for each patient
  indv_diff <- apply(subnetwork_data, 2,
                     function(col) sum(col * subnetwork_weight[indx1]) - sum(col * subnetwork_weight[indx2]))
  
  # Calculate tstat based on vector of subnetwork scores
  calc_tstat <- function(vector, mu) {
    n <- length(vector)
    t_stat <- (mean(vector) - mu) / (sd(vector)/sqrt(n))
    return(ifelse (is.infinite(t_stat)|is.na(t_stat), NaN, t_stat))
  }
  return(calc_tstat(indv_diff, 0))
}

# Takes in original dfs and geneset lists
# Environment: gfs_data
# Returns permutation null distribution for all subnetworks as arr
# Rows: No. of subnetworks, Cols: No. of permutations
generate_permutation_null <- function(df_A, df_B, geneset_list_A, geneset_list_B, max_perm) {
  # Combine both df
  df_hstack <- cbind(df_A, df_B)
  print(head(df_hstack))
  n1 <- ncol(df_A); n2 <- ncol(df_B)
  num_perm <- choose(n1+n2, n1)
  cat(sprintf("No. of possible combinations (%d choose %d): %d\n",
              n1+n2, n1, num_perm))
  if (num_perm > max_perm) {
    # Monte Carlo permutation tests
    montecarlo_permutation <- function(df_hstack) {
      permutation_indx <- sample(n1+n2, n1)
      # Calculate gene weight matrix of permuted class A
      gene_weight_A <- apply(df_hstack[, permutation_indx], 1, mean)
      # Calculate gene weight matrix of permuted class B
      gene_weight_B <- apply(df_hstack[, -permutation_indx], 1, mean)
      gene_weight_df <- data.frame(A = gene_weight_A, B = gene_weight_B)
      
      subnetwork_A_tstat <- unlist(mclapply(geneset_list_A, calc_subnetwork_tstat,
                                            "A", gene_weight_df, df_hstack,
                                            mc.cores = NUM_CORES))
      subnetwork_B_tstat <- unlist(mclapply(geneset_list_B, calc_subnetwork_tstat,
                                            "B", gene_weight_df, df_hstack,
                                            mc.cores = NUM_CORES))
      
      # Modifies global variable
      prog_indx <<- prog_indx + 1
      setTxtProgressBar(prog_bar, prog_indx)
      
      # Concatenate both null tstat of genesets from class A and B
      return(c(subnetwork_A_tstat, subnetwork_B_tstat))
    }
    
    cat(sprintf("Generating null distribution (%d Monte Carlo permutations):\n", max_perm))
    prog_indx <- 0
    # Initialise progress bar
    prog_bar <- txtProgressBar(min = 0, max = max_perm,
                               width = 50, style = 3)
    
    # Generate null distribution array
    montecarlo_null_arr <- replicate(max_perm, montecarlo_permutation(df_hstack))
    close(prog_bar)
    cat(sprintf("Successful permutations: %d\n", prog_indx))
    return(montecarlo_null_arr)
  } else {
    # Standard permutation tests
    # Generates array with all possible combinations as columns
    combination_indx_arr <- combn(n1+n2, n1)
    # Returns t-statistic vector for SINGLE permutation
    standard_permutation <- function(permutation_indx, df_hstack) {
      # Calculate gene weight matrix of permuted class A
      gene_weight_A <- apply(df_hstack[, permutation_indx], 1, mean)
      # Calculate gene weight matrix of permuted class B
      gene_weight_B <- apply(df_hstack[, -permutation_indx], 1, mean)
      gene_weight_df <- data.frame(A = gene_weight_A, B = gene_weight_B)
      
      subnetwork_A_tstat <- unlist(mclapply(geneset_list_A, calc_subnetwork_tstat,
                                            "A", gene_weight_df, df_hstack,
                                            mc.cores = NUM_CORES))
      subnetwork_B_tstat <- unlist(mclapply(geneset_list_B, calc_subnetwork_tstat,
                                            "B", gene_weight_df, df_hstack,
                                            mc.cores = NUM_CORES))

      # Modifies global variable 
      prog_indx <<- prog_indx + 1
      setTxtProgressBar(prog_bar, prog_indx)
      
      # Concatenate both null tstat of genesets from class A and B
      return(c(subnetwork_A_tstat, subnetwork_B_tstat))
    }
    
    cat(sprintf("Generating null distribution (%d standard permutations):\n", num_perm))
    prog_indx <- 0
    # Initialise progress bar
    prog_bar <- txtProgressBar(min = 0, max = num_perm,
                               width = 50, style = 3)
    
    permutation_null_arr <- apply(combination_indx_arr, 2, standard_permutation, df_hstack)
    close(prog_bar)
    cat(sprintf("Successful permutations: %d\n", prog_indx))
    return(permutation_null_arr)
  }
}

# Arguments: Original t-statistic vector and null distr arr
# Calculate p-value from null distribution (One-tailed)
empirical_pvalue <- function(tstat_vec, null_arr) {
  num_perm <- ncol(null_arr)
  # Driver function to calculate pvalue for SINGLE subnetwork
  calc_pvalue <- function(subnetwork_null_vec, tstat_scalar) {
    # One-tailed t-test does not need to abs(null_vec)
    pvalue <- sum(subnetwork_null_vec >= tstat_scalar)/num_perm
    return(pvalue)
  }
  # Tranposes such that columns are subnetworks for mapply
  null_df <- data.frame(t(null_arr))
  # print(data.frame(null = colnames(null_df),
  #                  tstat = names(tstat_vec)))
  pvalue_vec <- mapply(calc_pvalue, null_df, tstat_vec)
  return(pvalue_vec)
}

# Arguments: T-statistic vector and combined sample size
# Returns theoretical p-value of one-sided paired t-test
theoretical_pvalue <- function(tstat_vec, n) {
  return(1 - sapply(tstat_vec, pt, df = n-1))
}

# Be careful of wpath for edgelist and genelist of subnetwork
pfsnet <- function(df_A, df_B, pathway_rpath, dataset_name, pathwaydb_name, min_size = 5, max_perm = 1000) {
  gfs_df_A <- norm_gfs(df_A)
  gfs_df_B <- norm_gfs(df_B)
  gfs_df_AnB <- cbind(gfs_df_A, gfs_df_B)
  # Create list of highly expressed genes
  gene_weight_A <- apply(gfs_df_A, 1, mean)
  highexpr_genes_A <- names(gene_weight_A[gene_weight_A > 0.5])
  gene_weight_B <- apply(gfs_df_B, 1, mean)
  highexpr_genes_B <- names(gene_weight_B[gene_weight_B > 0.5])
  # Gene weight matrix for class A and B
  original_gene_weight_df <- data.frame(A = gene_weight_A, B = gene_weight_B)
  
  # Create subnetworks for both classA and classB
  list_pathway <- import_df_pathway(pathway_rpath)
  
  geneset_A <- generate_geneset(list_pathway, highexpr_genes_A, dataset_name, pathwaydb_name, "A", min_size)
  geneset_B <- generate_geneset(list_pathway, highexpr_genes_B, dataset_name, pathwaydb_name, "B", min_size)
  
  geneset_A_tstat <- unlist(mclapply(geneset_A, calc_subnetwork_tstat,
                                     "A", original_gene_weight_df, gfs_df_AnB,
                                     mc.cores = NUM_CORES))
  geneset_B_tstat <- unlist(mclapply(geneset_B, calc_subnetwork_tstat,
                                     "B", original_gene_weight_df, gfs_df_AnB,
                                     mc.cores = NUM_CORES))
  
  null_AnB_arr <- generate_permutation_null(gfs_df_A, gfs_df_B, geneset_A, geneset_B, max_perm)
  
  # Concatenate both tstat of genesets from class A and B
  geneset_AnB_tstat <- c(geneset_A_tstat, geneset_B_tstat)
  geneset_AnB_pvalue <- empirical_pvalue(geneset_AnB_tstat, null_AnB_arr)
  return(geneset_AnB_pvalue)
}

# MAIN --------------------------------------------------------------------
# Imports data and returns two lists of top genes and gene weight matrix
ovarian_data1 <- read.table("data/ovarian_cancer/GSE18521/processed/mas5_qnorm.tsv",
                            header = T, row.names = 1)
ovarian1_A <- ovarian_data1[,5:10]
ovarian1_B <- ovarian_data1[,53:58]

PWAPI_FPATH <- "../info/pathwayAPI/pwapi_id_human-filtered_entrez.tsv"
ovarian1_pwapi_pvalue <- pfsnet(ovarian1_A, ovarian1_B, PWAPI_FPATH, "GSE18521", "PWAPI")

# No. of identified gene sets for class A and B
length(ovarian1_pwapi_pvalue)
# No. of significant subnetworks
ovarian1_pwapi_pvalue[ovarian1_pwapi_pvalue < 0.05]

# Investigate -------------------------------------------------------------
# Imports data and returns two lists of top genes and gene weight matrix
ovarian_data1 <- read.table("data/ovarian_cancer/GSE18521/processed/mas5_qnorm.tsv",
                            header = T, row.names = 1)
ovarian1_A <- ovarian_data1[,5:10]
ovarian1_B <- ovarian_data1[,53:58]

# Batch information
scan_dates_df <- read.table("data/ovarian_cancer/GSE26712/README/scan_dates.tsv", sep = "\t")
rownames(scan_dates_df) <- substring(rownames(scan_dates_df), 1, 9)
id_annot <- read.table("data/ovarian_cancer/GSE26712/processed/annot-shortid.tsv",
                       header = T, row.names = 1)
simple_id <- as.character(sapply(rownames(scan_dates_df), function(x) id_annot[x,]))
scan_dates <- cbind(simple_id, scan_dates_df)
scan_dates

# Ovarian cancer data set 2
data2_rpath <- 'data/ovarian_cancer/GSE26712/processed/mas5_qnorm.tsv'
ovarian_data2 <- read.table(data2_rpath, header = T, row.names = 1)
ovarian2_classA <- ovarian_data2[,1:6]
ovarian2_classB <- ovarian_data2[,c(23,25,26,28,29,31)]
colnames(ovarian2_classB)

gfs_df_A <- norm_gfs(ovarian2_classA )
gfs_df_B <- norm_gfs(ovarian2_classB)
gfs_df_AnB <- cbind(gfs_df_A, gfs_df_B)
# Create list of highly expressed genes
gene_weight_A <- apply(gfs_df_A, 1, mean)
highexpr_genes_A <- names(gene_weight_A[gene_weight_A > 0.5])
gene_weight_B <- apply(gfs_df_B, 1, mean)
highexpr_genes_B <- names(gene_weight_B[gene_weight_B > 0.5])
# Gene weight matrix for class A and B
original_gene_weight_df <- data.frame(A = gene_weight_A, B = gene_weight_B)

# Generate subnetworks
# KEGGDF_RDIRPATH <- "../info/KEGG/hsa_df/"
# list_kegg <- import_kegg(KEGGDF_RDIRPATH)
# geneset_A_kegg <- generate_geneset(list_kegg, highexpr_genes_A, "GSE26712", "KEGG", "A", 5)
# geneset_B_kegg <- generate_geneset(list_kegg, highexpr_genes_B, "GSE26712", "KEGG", "B", 5)

GSA_RPATH <-"data/subnetwork/pfsnet/geneset-KEGG_GSE26712_A.tsv"
geneset_A_df <- read.table(GSA_RPATH, header = T)
geneset_A <- split(geneset_A_df[,2], geneset_A_df[,1])

GSB_RPATH <-"data/subnetwork/pfsnet/geneset-KEGG_GSE26712_B.tsv"
geneset_B_df <- read.table(GSB_RPATH, header = T)
geneset_B <- split(geneset_B_df[,2], geneset_B_df[,1])

null_AnB_arr <- generate_permutation_null(gfs_df_A, gfs_df_B, geneset_A, geneset_B, 1000)

geneset_A_tstat <- unlist(mclapply(geneset_A, calc_subnetwork_tstat,
                                   "A", original_gene_weight_df, gfs_df_AnB,
                                   mc.cores = NUM_CORES))

geneset_B_tstat <- unlist(mclapply(geneset_B, calc_subnetwork_tstat,
                                   "B", original_gene_weight_df, gfs_df_AnB,
                                   mc.cores = NUM_CORES))

geneset_AnB_tstat <- c(geneset_A_tstat, geneset_B_tstat)
geneset_AnB_pvalue <- empirical_pvalue(geneset_AnB_tstat, null_AnB_arr)

names(geneset_AnB_pvalue)[geneset_AnB_pvalue <= 0.05]

pvalue1 <- theoretical_pvalue(geneset_AnB_tstat, 20)
names(pvalue1)[pvalue1 <= 0.05]

# Null distribution -------------------------------------------------------
# Visualise null distributions of t-statistic
k <- length(geneset_AnB_pvalue)
par(mfrow = c(2,2), mai = c(0.4, 0.4, 0.4, 0.2))
for (i in 1:k) {
  title <- sprintf("%s (%.4f)", rownames(null_AnB_arr)[i], geneset_AnB_pvalue[i])
  hist(null_AnB_arr[i,], breaks = 20, main = title)
  # Plot t-statistic
  abline(v = geneset_AnB_tstat[i], col = "red")
  if (i %in% c(4*1:k%/%4, k)) {
    fig <- recordPlot()
    fig_wpath <- sprintf("dump/null_distr%d.eps", i)
    save_fig(fig, fig_wpath, width = 10, height = 8)
  }
}

unname(geneset_AnB_pvalue)

# Visualise subnetworks ---------------------------------------------------
# Plot subnetworks generated from pfsnet
library(Rgraphviz)
df_A <- read.table("data/subnetwork/pfsnet/edgelist-KEGG_GSE18521_A.tsv", header = T)
list_subnetworks_A <- split(df_A[,2:3], df_A[,1])
list_arr_A <- lapply(list_subnetworks_A, data.matrix)
list_arr_A[1]

# Copy default Rgraphviz attributes
graph_attr <- getDefaultAttrs()
# Alter Rgraphviz attributes
graph_attr$graph$bgcolor <- "white"
graph_attr$node$fontsize <- 14

n <- length(list_arr_A)
par(mfrow = c(3,4), mai = c(0.4, 0.2, 0.4, 0.2))
for (i in 1:n) {
  # title <- sprintf("%s (%.4f)", rownames(null_AnB_arr)[i], geneset_AnB_pvalue[i])
  plot(ftM2graphNEL(list_arr_A[[i]]), attrs = graph_attr)
  # Plot t-statistic
  if (i %in% c(12*1:n%/%12, n)) {
    fig <- recordPlot()
    fig_wpath <- sprintf("dump/subnetworkA_%d.eps", i)
    save_fig(fig, fig_wpath, width = 10, height = 8)
  }
}

df_B <- read.table("data/subnetwork/pfsnet/edgelist-KEGG_GSE18521_B.tsv", header = T)
list_subnetworks_B <- split(df_B[,2:3], df_B[,1])
list_arr_B <- lapply(list_subnetworks_B, data.matrix)

n <- length(list_arr_B)
par(mfrow = c(3,4), mai = c(0.4, 0.2, 0.4, 0.2))
for (i in 1:n) {
  # title <- sprintf("%s (%.4f)", rownames(null_AnB_arr)[i], geneset_AnB_pvalue[i])
  plot(ftM2graphNEL(list_arr_B[[i]]), attrs = graph_attr)
  # Plot t-statistic
  if (i %in% c(12*1:n%/%12, n)) {
    fig <- recordPlot()
    fig_wpath <- sprintf("dump/subnetworkB_%d.eps", i)
    save_fig(fig, fig_wpath, width = 10, height = 8)
  }
}

