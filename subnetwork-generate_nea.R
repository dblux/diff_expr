library(KEGGgraph)
library(Rgraphviz)
library(igraph)
source("functions.R")
# detach(package:Rgraphviz)
# Rgraphviz is built on graph package: graphNEL format
setwd("~/projects/phd/diff_expr/")

# FUNCTIONS ---------------------------------------------------------------
# Generate NEA regulator-effector style
# Select for subnetworks >= min number of edges!
# Some subnetworks have less than size of 5! (Directed edges)
# Allow for both directions reg -> effector and effector -> reg
# Receives a df that has two columns from and to
subnetwork_nea <- function(pathway_df, list_microarray, min_size = 5) {
  # Selects edges where both nodes are present in the microarray
  fltr_pathway <- pathway_df[pathway_df$from %in% list_microarray & pathway_df$to %in% list_microarray,]
  # Assumes there are no loops in the graph
  subnetwork_nodes <- unique(unlist(fltr_pathway))
  names(subnetwork_nodes) <- subnetwork_nodes
  # Generates list of subnetwork df
  subnetwork_list <- lapply(subnetwork_nodes,
                            function(x) fltr_pathway[(fltr_pathway$from == x | fltr_pathway$to == x),])
  # Calculate subnetwork size according to number of nodes
  subnetwork_size <- sapply(subnetwork_list, function(x) length(unique(unlist(x))))
  # Subset only subnetworks with size above minimum
  return(subnetwork_list[subnetwork_size >= min_size])
}

# NEA - KEGG --------------------------------------------------------------
# Select the intersect of probesets between two datasets OR
# Select dataset with the smallest number of probesets
# Group of probesets has to be a subset of the other dataset
data <- read.table('data/ovarian_cancer/GSE26712/processed/mas5_qnorm.tsv',
                  header = T, sep = "\t", row.names = 1)
GENES_MICROARRAY <- rownames(data)

# List of kegg pathway files
kegg_rdirpath <- "../info/KEGG/human_pathways/"
all_keggxml <- list.files(substring_head(kegg_rdirpath, 1), full.names = T)

# Generates subnetworks NEA style for one KEGG pathway
# Returns a dataframe of vertically stacked subnetworks for ONE KEGG pathway
subnetwork_nea_kegg <- function(kegg_fpath, genes_microarray, min_size) {
  path_name <- substring(kegg_fpath, 29, 36)
  # Parses file into dataframe
  df_hsa <- parseKGML2DataFrame(kegg_fpath)
  # If dataframe is not NULL
  if (nrow(df_hsa) != 0) {
    # Convert KEGG ID to Entrez ID
    df_hsa$from <- substring(df_hsa$from, 5)
    df_hsa$to <- substring(df_hsa$to, 5)
    # unwanted_edges <- c("phosphorylation", "indirect effect", "dephosphorylation",
    #                     "ubiquitination", "methylation", "state change", "missing interaction")
    # Removes unwanted edges and loops in graph
    processed_df <- df_hsa[!(df_hsa$from == df_hsa$to), 1:2]
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
    
    # Create subnetworks
    subnetwork_ls <- subnetwork_nea(processed_df, genes_microarray, min_size)
    if (length(subnetwork_ls) != 0) {
      # Generate single dataframe of all subnetworks for each pathway
      all_subnetwork <- do.call(rbind, subnetwork_ls)
      # Get subnetwork names from names of list
      subnetwork_name <- paste(path_name, names(subnetwork_ls), sep = "_")
      # Get size of each subnetwork
      subnetwork_size <- sapply(subnetwork_ls, nrow)
      subnetwork_id <- rep(subnetwork_name, subnetwork_size)
      df_subnetwork <- cbind(subnetwork_id, all_subnetwork)
      rownames(df_subnetwork) <- NULL
      print(kegg_fpath)
      return(df_subnetwork)
    }
  }
}

# NEA subnetwork generation on KEGG
kegg_nea_raw <- lapply(all_keggxml, subnetwork_nea_kegg, GENES_MICROARRAY, 5)
# Omit NULL values from list (from pathways without subnetworks)
kegg_nea_list <- kegg_nea_raw[!sapply(kegg_nea_raw, is.null)]

# Vertically stack list of dataframes to consolidate
df_subnetworks <- do.call(rbind, kegg_nea_list)

# Print total number of subnetworks
subnetwork_size <- table(df_subnetworks[,1])
print(paste0("Total no. of subnetworks before: ", length(subnetwork_size)))

# Removes duplicate subnetworks in dataframe
# Returns a dataframe
remove_duplicate_subnetworks <- function(df_subnetworks) {
  subnetwork_size <- table(df_subnetworks[,1])
  # Check whether there are any duplicate center nodes
  center_nodes <- substring(names(subnetwork_size), 10)
  center_freq <- table(center_nodes)
  duplicated_nodes <- center_freq[center_freq > 1]
  
  # Check for duplicate subnetworks in NEA
  duplicate_subnetworks <- list()
  # LOOP through each duplicated node
  for (node in names(duplicated_nodes)) {
    same_centernode <- subnetwork_size[endsWith(names(subnetwork_size), paste0("_", node))]
    length_freq <- table(same_centernode)
    fltr_length_freq <- length_freq[length_freq > 1]
    # LOOP through each length
    for (length in as.numeric(names(fltr_length_freq))) {
      list_identical <- list()
      # LOOP through pairwise comparison between every subnetwork with same length (~n^2/2)
      same_length <- names(same_centernode[same_centernode == length])
      for (i in 1:(length(same_length)-1)) {
        for (j in 1:i) {
          df1 <- df_subnetworks[df_subnetworks[1] == same_length[i+1], 2:3]
          df2 <- df_subnetworks[df_subnetworks[1] == same_length[j], 2:3]
          # Sort by both columns
          df1 <- df1[order(df1$from, df1$to),]
          df2 <- df2[order(df2$from, df2$to),]
          # Remove row names
          row.names(df1) <- NULL -> row.names(df2)
          # If dataframes are indentical
          if (identical(df1, df2)) {
            identical_pair <- c(same_length[i+1], same_length[j])
            # LOOP through list for this particular pairwise comparison
            if (length(list_identical) == 0) {
              # No need to check if pair is present in previous sets
              list_identical <- c(list_identical, list(identical_pair))
            } else {
              for (k in 1:length(list_identical)) {
                # If any of the pair is present in previous sets
                if (sum(identical_pair %in% list_identical[[k]]) >= 1) {
                  # Update previous set
                  list_identical[[k]] <- unique(c(list_identical[[k]], identical_pair))
                  break
                } else if (k == length(list_identical)) {
                  # Not present in any previous set
                  # Executes at last cycle of loop
                  list_identical <- c(list_identical, list(identical_pair))
                }
              }
            }
          }
        }
      }
      # Concat list to overall list
      if (length(list_identical) > 0) {
        duplicate_subnetworks <- c(duplicate_subnetworks, list_identical)
      }
    }
  }
  # Returns list of duplicate subnetworks: duplicate_subnetworks
  # Subnetworks to be removed (retain first subnetwork of each vector)
  rm_subnetwork <- unlist(lapply(duplicate_subnetworks, "[", -1))
  # Remove rows in df corresponding to subnetworks in rm_subnetwork
  fltr_df_subnetworks <- df_subnetworks[!(df_subnetworks$subnetwork_id %in% rm_subnetwork),]
  return(fltr_df_subnetworks)
}

# Remove duplicate subnetworks
fltr_df <- remove_duplicate_subnetworks(df_subnetworks)
num_fltr_subnetwork <- length(unique(fltr_df$subnetwork_id))
print(paste0("Total no. of subnetworks after: ", num_fltr_subnetwork))

# Save filtered subnetworks
write.table(fltr_df, "data/subnetwork/nea-hsa/ovarian_cancer/subnetworks-ovarian.tsv",
            quote = F, sep = "\t", row.names = F, col.names = T)

# NEA - PWAPI -------------------------------------------------------------
ds1 <- read.table('data/ovarian_cancer/GSE18521/processed/GSE18521_entrez1.tsv',
                  header = T, sep = "\t", row.names = 1)
GENES_MICROARRAY <- rownames(ds1)

PWAPI_FPATH <- "../info/pathwayAPI/pwapi_id_human-filtered_entrez.tsv"
pwapi_df <- read.table(PWAPI_FPATH, header=T,
                       sep="\t", stringsAsFactors = F)
list_pwAPI <- split(pwapi_df[,2:3], pwapi_df$pathway)
id_pwAPI <- names(list_pwAPI)

# Initialise list of dataframes
list_df <- list()
for (j in 1:length(list_pwAPI)) {
  df_pwapi <- list_pwAPI[[j]]
  pathway_id <- id_pwAPI[j]
  # If dataframe is not NULL
  if (nrow(df_pwapi) != 0) {
    # Removes loops in graph
    processed_df <- df_pwapi[!(df_pwapi$from == df_pwapi$to),]
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
    # Create subnetworks
    MIN_SIZE <- 5
    subnetwork_ls <- subnetwork_nea(processed_df, GENES_MICROARRAY, MIN_SIZE)
    # Remove subnetworks with less than 5 nodes
    subnetwork_size <- sapply(subnetwork_ls, function(x) length(unique(unlist(x))))
    subnetwork_ls1 <- subnetwork_ls[subnetwork_size >= MIN_SIZE]
    if (length(subnetwork_ls1) != 0) {
      # Generate single dataframe of all subnetworks for each pathway
      all_subnetwork <- do.call(rbind, subnetwork_ls1)
      # # Alternative foolproof method to get subnetwork name
      # subnetwork_name <- paste0("hsa04152_", names(hsa_subnetwork_ls1))
      # subnetwork_size <- sapply(hsa_subnetwork_ls1, nrow)
      # df_name <- rep(subnetwork_name, subnetwork_size)
      # Extract subnetwork name
      rname <- rownames(all_subnetwork)
      center_gene <- sub("\\..*$", "", rname)
      subnetwork_name <- paste(pathway_id, center_gene, sep = "_")
      df_subnetwork <- cbind(subnetwork_name, all_subnetwork)
      print(j)
      # Add df to list
      list_df <- append(list_df, list(df_subnetwork))
    }
  }
}
# Create df of all subnetworks
df_all_subnetwork <- do.call(rbind, list_df)
OUTPUT_FPATH <- "data/subnetwork/nea-pwapi/ovarian_cancer/raw_nea_pwapi.tsv"
write.table(df_all_subnetwork, OUTPUT_FPATH, quote = F, sep = "\t", row.names = F)

# Check for duplicate subnetworks in NEA
RAW_SUBNETWORK_PATH <- "data/subnetwork/nea-pwapi/ovarian_cancer/raw_nea_pwapi.tsv"
df_subnetworks <- read.table(RAW_SUBNETWORK_PATH, sep = "\t",
                              header = T, strip.white = T)
subnetwork_size <- table(df_subnetworks$subnetwork_name)
# Print total number of subnetworks
num_subnetworks <- length(subnetwork_size)
print(paste0("Total no. of subnetworks: ", num_subnetworks))

# Check whether there are any duplicate center nodes
# Ignore "pwapi001"
center_nodes <- substring(names(subnetwork_size), 10)
center_freq <- table(center_nodes)
duplicated_nodes <- center_freq[center_freq > 1]
# LOOP through each duplicated node
for (node in names(duplicated_nodes)) {
  # Paste underscore to avoid error with endsWith short node name
  same_centernode <- subnetwork_size[endsWith(names(subnetwork_size), paste0("_", node))]
  length_freq <- table(same_centernode)
  fltr_length_freq <- length_freq[length_freq > 1]
  # LOOP through each length
  for (length in as.numeric(names(fltr_length_freq))) {
    list_identical <- list()
    # LOOP through pairwise comparison between every subnetwork with same length (~n^2/2)
    same_length <- names(same_centernode[same_centernode == length])
    for (i in 1:(length(same_length)-1)) {
      for (j in 1:i) {
        df1 <- df_subnetworks[df_subnetworks[1] == same_length[i+1], 2:3]
        df2 <- df_subnetworks[df_subnetworks[1] == same_length[j], 2:3]
        # Sort by both columns
        df1 <- df1[order(df1$from, df1$to),]
        df2 <- df2[order(df2$from, df2$to),]
        # Remove row names
        row.names(df1) <- NULL -> row.names(df2)
        # If dataframes are indentical
        if (identical(df1, df2)) {
          identical_pair <- c(same_length[i+1], same_length[j])
          # LOOP through list for this particular pairwise comparison
          if (length(list_identical) == 0) {
            # No need to check if pair is present in previous sets
            list_identical <- c(list_identical, list(identical_pair))
          } else {
            for (k in 1:length(list_identical)) {
              # If any of the pair is present in previous sets
              if (sum(identical_pair %in% list_identical[[k]]) >= 1) {
                # Update previous set
                list_identical[[k]] <- unique(c(list_identical[[k]], identical_pair))
                break
              } else if (k == length(list_identical)) {
                # Not present in any previous set
                # Executes at last cycle of loop
                list_identical <- c(list_identical, list(identical_pair))
              }
            }
          }
        }
      }
    }
    # Concat list to overall list
    if (length(list_identical) > 0) {
      list_list <- c(list_list, list_identical)
    }
  }
}

# Saves list of duplicated subnetworks in rds format
print(list_list)
LIST_FPATH <- "data/subnetwork/nea-pwapi/ovarian_cancer/duplicated_subnetworks.rds"
saveRDS(list_list, LIST_FPATH)
# list_list <- readRDS(LIST_FPATH)

# List of subnetworks to be removed
rm_subnetwork <- unlist(lapply(list_list, "[", -1))
# Filtered df
fltr_df_subnetworks <- df_subnetworks[!(df_subnetworks$subnetwork_name %in% rm_subnetwork),]
num_fltr_subnetwork <- length(unique(fltr_df_subnetworks$subnetwork_name))
print(num_fltr_subnetwork)
# Save filtered subnetworks
SUBNETWORK_FPATH <- "data/subnetwork/nea-pwapi/ovarian_cancer/subnetworks-nea_pwapi_ovarian.tsv"
write.table(fltr_df_subnetworks, SUBNETWORK_FPATH,
            quote = F, sep = "\t", row.names = F, col.names = T)
