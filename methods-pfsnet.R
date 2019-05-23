library(igraph)
source("functions.R")
# detach(package:Rgraphviz)
# Rgraphviz is built on graph package: graphNEL format
setwd("~/projects/phd/diff_expr/")

# FUNCTIONS ---------------------------------------------------------------
# Arguments: dataframe of one pathway, list of highly expr genes and min size
# Returns a list of subnetwork dataframes
generate_subnetwork <- function(df, gene_list, MIN_SIZE) {
  list_subnetwork <- list()
  # Only retain edges where both nodes are in list of genes
  fltr_df <- df[df$from %in% gene_list
                & df$to %in% gene_list,]
  # Removes loops in graph
  processed_df <- fltr_df[!(fltr_df$from == fltr_df$to),]
  # Remove duplicated edges
  concat_str <- with(processed_df, paste0(from, to))
  freq <- table(concat_str)
  dup_str <- names(freq[freq > 1])
  if (length(dup_str) == 1) {
    row_no <- which(concat_str == dup_str)[-1]
    # Returns df with row removed
    processed_df <- processed_df[-row_no,]
  } else if (length(dup_str) > 1) {
    for (i in dup_str) {
      row_no <- which(concat_str == i)[-1]
      # Returns df with row removed
      processed_df <- processed_df[-row_no,]
    }
  }
  # TODO: Check logic of tree
  if (nrow(processed_df) >= (MIN_SIZE - 1)) {
    # Convert numeric matrix to string matrix in order to avoid vertex id
    char_array <- apply(data.matrix(processed_df), c(1, 2), as.character)
    sub_graphs <- graph_from_edgelist(char_array)
    connected_list <- decompose.graph(sub_graphs, min.vertices = MIN_SIZE)
    for (graph in connected_list) {
      subnetwork <- get.edgelist(graph)
      list_subnetwork <- append(list_subnetwork, list(subnetwork))
    }
  }
  return(list_subnetwork)
}

# Arguments: List of subnetwork, subnetwork id, file path
# Saves list of subnetwork as one dataframe
save_subnetwork <- function(list_subnetwork, subnetwork_name, fpath) {
  num_id <- sprintf("%.2d", 1:length(list_subnetwork))
  subnetwork_id <- paste0(subnetwork_name, num_id)
  # Assign subnetwork ids
  names(list_subnetwork) <- subnetwork_id
  # Create column of subnetwork names
  num_row <- lapply(list_subnetwork, nrow)
  subnetwork_id <- rep(names(num_row), num_row)
  subnetwork_df <- do.call(rbind, list_subnetwork)
  all_subnetworks <- cbind(subnetwork_id, subnetwork_df)
  # Assign column names
  colnames(all_subnetworks)[2:3] <- c("from", "to")
  write.table(all_subnetworks, fpath,
              sep = "\t", row.names = F, quote = F)
}

# MAIN --------------------------------------------------------------------
# Import pathwayAPI
PWAPI_FPATH <- "../info/pathwayAPI/pwapi_id_human-filtered_entrez.tsv"
pwapi_df <- read.table(PWAPI_FPATH, header=T,
                       sep="\t", stringsAsFactors = F)
list_pwAPI <- split(pwapi_df[,2:3], pwapi_df$pathway)

# Import ovarian cancer data set 1
# Not log2
ovarian_data1 <- read.table("data/ovarian_cancer/GSE18521/processed/GSE18521_entrez1.tsv",
                            header = T, row.names = 1)
gfs_ovarian_data1 <- GFS(ovarian_data1)
classA_data1 <- gfs_ovarian_data1[,1:10]
classB_data1 <- gfs_ovarian_data1[,11:62]

# Create list of highly expressed genes
mean_beta_classA <- apply(classA_data1, 1, mean)
highexpr_genes_classA <- names(mean_beta_classA[mean_beta_classA > 0.5])
mean_beta_classB <- apply(classB_data1, 1, mean)
highexpr_genes_classB <- names(mean_beta_classB[mean_beta_classB > 0.5])

# Create subnetworks for both classA and classB
all_subnetwork_classA <- list()
all_subnetwork_classB <- list()
for (df_pwapi in list_pwAPI) {
  subnetwork_single_pathway_classA <- generate_subnetwork(df_pwapi, highexpr_genes_classA, 5)
  all_subnetwork_classA <- append(all_subnetwork_classA, subnetwork_single_pathway_classA)
  subnetwork_single_pathway_classB <- generate_subnetwork(df_pwapi, highexpr_genes_classB, 5)
  all_subnetwork_classB <- append(all_subnetwork_classB, subnetwork_single_pathway_classB)
}

# Save subnetworks for GSE18521 classA
fpath_classA <- "data/subnetwork/pfsnet/subnetworks-pfsnet_pwapi_GSE18521_A.tsv"
subnetwork_name_classA <- "GSE18521_A"
save_subnetwork(all_subnetwork_classA, subnetwork_name_classA, fpath_classA)
# Save subnetworks for GSE18521 classBs
fpath_classB <- "data/subnetwork/pfsnet/subnetworks-pfsnet_pwapi_GSE18521_B.tsv"
subnetwork_name_classB <- "GSE18521_B"
save_subnetwork(all_subnetwork_classB, subnetwork_name_classB, fpath_classB)

# # Ovarian cancer data set 2
# # Pre-processed using RMA - Log2 and quantile normalised
# ovarian_data2 <- read.table('data/ovarian_cancer/GSE26712/processed/GSE26712_entrez.tsv',
#                             header = T, row.names = 1)
# # Reverse the log2
# ovarian_data2a <- 2^ovarian_data2
# gfs_ovarian_data2 <- GFS(ovarian_data2)
# classA_data2 <- ovarian_data2a[,1:10]
# tumour_data2 <- ovarian_data2a[,11:195]
