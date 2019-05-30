# Initialisation ----------------------------------------------------------
library(parallel)
library(VennDiagram)
source("functions.R")
NUM_CORES <- detectCores()
# Arguments: List of gene sets, expression data (df1 and df2)
# Returns: List of pvalues for each gene set
# Import data -------------------------------------------------------------
# Import ovarian cancer data set 1
data1_rpath <- "data/ovarian_cancer/GSE18521/processed/mas5_qnorm.tsv"
ovarian_data1 <- read.table(data1_rpath, header = T, row.names = 1)
ovarian1_classA <- ovarian_data1[1:5]
ovarian1_classB <- ovarian_data1[11:15]

# Ovarian cancer data set 2
data2_rpath <- 'data/ovarian_cancer/GSE26712/processed/mas5_qnorm.tsv'
ovarian_data2 <- read.table(data2_rpath, header = T, row.names = 1)
ovarian2_classA <- ovarian_data2[1:5]
ovarian2_classB <- ovarian_data2[11:15]

# Load subnetwork gene list
geneset_rpath <- "data/subnetwork/nea-hsa/ovarian_cancer/geneset-ovarian.tsv"
geneset_df <- read.table(geneset_rpath, header = T)
print(head(geneset_df))
geneset_list <- split(geneset_df$gene, geneset_df$subnetwork_id)
print(head(geneset_list))

MIN_SIZE <- 5
# Check for number of subnetworks with size smaller than min size
print(paste("No. of subnetworks with size smaller than min. size:", 
            sum(sapply(geneset_list, length) < MIN_SIZE)))

# # Filter out subnetworks that fall below min size
# fltr_genelist <- list_genelist[(sapply(list_genelist, length) < MIN_SIZE)]
# fltr_genelist

# ESSNet ------------------------------------------------------------------
essnet <- function(df1, df2, geneset_df) {
  print(head(geneset_df))
  geneset_list <- split(geneset_df$gene, geneset_df$subnetwork_id)
  print(head(geneset_list))
  print(length(geneset_list))
  
  geneset_genes <- sort(unique(geneset_df[,2]))
  # Name geneset_genes for functions
  names(geneset_genes) <- geneset_genes
  print(tail(geneset_genes))
  # IMPT: If not indexing of rows will be wrong
  geneset_genes <- as.character(geneset_genes)
  
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
    t_stat <- calc_tstat(subnetwork_vec, 0)
    return(t_stat)
  }
  
  # Constructions original matrix of all combinations
  original_combi_arr <- construct_combi_arr(geneset_genes, df1, df2)
  # Calculates t-statistic for every gene set
  original_tstat_vec <- unlist(mclapply(geneset_list,
                                        subnetwork_tstat,
                                        original_combi_arr,
                                        mc.cores = NUM_CORES))
  
  # Takes in original dfs
  # Returns permutation null distribution for all subnetworks
  # Rows: No. of subnetworks, Cols: No. of permutations
  generate_permutation_null <- function(df1, df2, geneset_list, geneset_genes, n = 1000) {
    # Combine both df
    df_stack <- cbind(df1, df2)
    n1 <- ncol(df1); n2 <- ncol(df2)
    num_perm <- choose(n1+n2, n1)
    cat(sprintf("No. of combinations (%d choose %d): %d\n",
                n1+n2, n1, num_perm))
    if (num_perm > n) {
      # Monte Carlo permutation tests
      montecarlo_permutation <- function(df_stack) {
        permutation_indx <- sample(n1+n2, n1)
        df_p1 <- df_stack[, permutation_indx]
        df_p2 <- df_stack[, -permutation_indx]
        # Takes in original dfs
        print(system.time({
          perm_combi_arr <- construct_combi_arr(geneset_genes, df_p1, df_p2)
        }))
        print(system.time({
          perm_tstat_vec <- unlist(mclapply(geneset_list,
                                            subnetwork_tstat,
                                            perm_combi_arr,
                                            mc.cores = NUM_CORES))
        }))
        # Modifies global variable
        prog_indx <<- prog_indx + 1
        setTxtProgressBar(prog_bar, prog_indx)
        return(perm_tstat_vec)
      }
      cat("Generating null distribution (Monte Carlo permutations):\n")
      prog_indx <- 0
      # Starts progress bar
      prog_bar <- txtProgressBar(min = 0, max = n,
                                 width = 50, style = 3)
      montecarlo_null_arr <- replicate(n, montecarlo_permutation(df_stack))
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
        print(system.time({
          perm_combi_arr <- construct_combi_arr(geneset_genes, df_p1, df_p2)
        }))
        print("Calculating t-stat!")
        print(system.time({
          perm_tstat_vec <- unlist(mclapply(geneset_list,
                                            subnetwork_tstat,
                                            perm_combi_arr,
                                            mc.cores = NUM_CORES))
        }))
        # Modifies global variable 
        prog_indx <<- prog_indx + 1
        setTxtProgressBar(prog_bar, prog_indx)
        return(perm_tstat_vec)
      }
      cat("Generating null distribution (standard permutations):\n")
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
  # Calculate p-value from null distribution
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
  
  permutation_null_df <- generate_permutation_null(df1, df2,
                                                   geneset_list = geneset_list,
                                                   geneset_genes = geneset_genes)
  
  emp_pvalue_vec <- empirical_pvalue(original_tstat_vec, permutation_null_df)
  return(emp_pvalue_vec)
}

pvalue1 <- essnet(ovarian1_classA, ovarian1_classB, geneset_df)
pvalue2 <- essnet(ovarian2_classA, ovarian2_classB, geneset_df)

sig_subnetworks1 <- names(pvalue1)[pvalue1 <= 0.05]
sig_subnetworks2 <- names(pvalue2)[pvalue2 <= 0.05]

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

# Takes a matrix with rows n and columns p
# Rotate matrix from each class separately
# Code reproduced from gsearot
rotation <- function(Y, method = 1) { 
  n <- dim(Y)[1]
  n <- 1000
  #Generate matrix with random standard normally distributed numbers
  #See: Rotation tests by Øyvind Langsrud
  W <- matrix(rnorm(n^2), nrow=n, ncol=n)
  
  if(method == 1) {
    # method 1 for matrices with population mean=0
    QRdata <- qr(W)
    Qm <- qr.Q(QRdata)
    
  } else if(method == 2) {
    # method 2 for matrices with population mean!=0
    E<-scale(W,scale=F)
    QRdata <- qr(E)
    W <- qr.Q(QRdata)[,1:(n-1)]
    E <- matrix(rnorm((n-1)^2),(n-1))
    QRdata <- qr(E)
    R.star <- qr.Q(QRdata)
    Qm <- n^(-1)*matrix(1,n,n)+(W%*%R.star)%*%t(W)
    
  } else {
    stop("method must be either 1 or 2 \n")    
  }
  
  #Multiply rotation matrix with original data matrix
  return(Qm%*%Y)
}
