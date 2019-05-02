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
subnetwork_nea <- function(pathway_df, list_microarray, min_size) {
  # Selects edges where both nodes are present in the microarray
  fltr_pathway <- pathway_df[pathway_df$from %in% list_microarray & pathway_df$to %in% list_microarray,]
  # Generates frequency table of nodes in flattened df
  # Assumes there are no loops in the graph
  node_freq <- table(unname(unlist(fltr_pathway)))
  subnetwork_edges <- node_freq[node_freq >= min_size]
  subnetwork_nodes <- names(subnetwork_edges)
  names(subnetwork_nodes) <- names(subnetwork_edges)
  # Generates list of subnetwork df
  subnetwork_df <- lapply(subnetwork_nodes,
                          function(x) fltr_pathway[(fltr_pathway$from == x | fltr_pathway$to == x),])
  return(subnetwork_df)
}

# DATA --------------------------------------------------------------------
# Load breast metastasis data set
ds1 <- read.table('data/ovarian_cancer/GSE18521/processed/processed_entrez.tsv',
                  header = T, sep = "\t", row.names = 1)
genes_microarray <- rownames(ds1)

# TODO: Preserve largest connected portion of graph

# Load KEGG pathway
# kegg <- vector(mode = "list", length = 100)
hsa <- parseKGML2Graph("../info/KEGG/human_pathways/hsa04152.xml")
df_hsa <- parseKGML2DataFrame("../info/KEGG/human_pathways/hsa01521.xml")

# Single subnetwork generation
df_hsa$from <- substring(df_hsa$from, 5)
df_hsa$to <- substring(df_hsa$to, 5)
unwanted_edges <- c("phosphorylation", "indirect effect", "dephosphorylation",
                    "ubiquitination", "methylation", "state change",
                    "missing interaction", "binding/association")
# Remove unwanted edges and loops in graph
processed_df <- df_hsa[!(df_hsa$subtype %in% unwanted_edges) & !(df_hsa$from == df_hsa$to), 1:2]

# KEGG hsa04664 subnetworks
hsa_subnetwork_ls <- subnetwork_nea(processed_df, genes_microarray, 5)
pwapi_subnetwork_ls <- subnetwork_nea(all_pwAPI[[2]], genes_microarray, 5)
pwapi_subnetwork_ls

# Visualise subnetworks
# Plot multiple
# png("subnetworks.png", width = 1200, height = 1000)
par(mfrow = c(6,5))
for (i in 1:length(hsa_subnetwork_ls)) {
  # plot(ftM2graphNEL(data.matrix(i)))
}
# dev.off()

# Plot single subnetwork
graph <- ftM2graphNEL(data.matrix(hsa_subnetwork_ls[[1]]))
plot(graph)

networks <- read.table("data/subnetwork/hsa-nea/breast_metastasis/subnetworks-breast_metastasis.tsv",
                        header = T, sep = "\t")

df <- networks[1:6,]
length(unique(unlist(df)))
x <- as.list(1:10)

# NEA subnetwork generation on KEGG
all_keggxml <- list.files("../info/KEGG/human_pathways/")
for (fname in all_keggxml) {
  fpath <- paste0("../info/KEGG/human_pathways/", fname)
  path_name <- substring(fname, 1, 8)
  # Parses file into dataframe
  df_hsa <- parseKGML2DataFrame(fpath)
  # If dataframe is not NULL
  if (nrow(df_hsa) != 0) {
    # Convert KEGG ID to Entrez ID
    df_hsa$from <- substring(df_hsa$from, 5)
    df_hsa$to <- substring(df_hsa$to, 5)
    unwanted_edges <- c("phosphorylation", "indirect effect", "dephosphorylation",
                        "ubiquitination", "methylation", "state change", "missing interaction")
    # Removes unwanted edges and loops in graph
    processed_df <- df_hsa[!(df_hsa$subtype %in% unwanted_edges) & !(df_hsa$from == df_hsa$to), 1:2]
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
    # Create subnetworks
    MIN_SIZE <- 5
    subnetwork_ls <- subnetwork_nea(processed_df, genes_microarray, MIN_SIZE)
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
      subnetwork_name <- paste(path_name, sub("\\..*$", "", rname), sep = "_")
      df_subnetwork <- cbind(subnetwork_name, all_subnetwork)
      output_fpath <- paste0("data/subnetwork/hsa-nea/", path_name, ".tsv")
      write.table(df_subnetwork, output_fpath, quote = F, sep = "\t", row.names = F)
      print(output_fpath)
    }
  }
}

# Check for duplicate subnetworks in NEA
# identical_df <- data.frame()
list_list <- list()
all_hsa_nea <- list.files("data/subnetwork/hsa-nea/raw")
all_fpath <- paste0("data/subnetwork/hsa-nea/raw/", all_hsa_nea)
all_subnetworks <- lapply(all_fpath, read.table, sep = "\t", header = T, strip.white = T)
# Subset column 1
all_col1 <- lapply(all_subnetworks, "[", 1)
all_freq <- lapply(all_col1, table)
all_size <- do.call(c, all_freq)
df_subnetworks <- do.call(rbind, all_subnetworks)
# Print total number of subnetworks
num_subnetworks <- length(all_size)
print(paste0("Total no. of subnetworks: ", num_subnetworks))

# Check whether there are any duplicate center nodes
center_nodes <- substring(names(all_size), 10)
center_freq <- table(center_nodes)
duplicated_nodes <- center_freq[center_freq > 1]
# LOOP through each duplicated node
for (node in names(duplicated_nodes)) {
  same_centernode <- all_size[endsWith(names(all_size), paste0("_", node))]
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

list_list

# Saves list into RDS
saveRDS(list_list, "data/subnetwork/duplicated_list.rds")
list_list <- readRDS("data/subnetwork/duplicated_list.rds")

# Remove first element of each vector
rm_subnetwork <- unlist(lapply(list_list, "[", -1))
fltr_df_subnetworks <- df_subnetworks[!(df_subnetworks$subnetwork_name %in% rm_subnetwork),]
num_fltr_subnetwork <- length(unique(fltr_df_subnetworks$subnetwork_name))

# Save filtered subnetworks
write.table(fltr_df_subnetworks, "data/subnetwork/hsa-nea/all_subnetworks.tsv",
            quote = F, sep = "\t", row.names = F, col.names = T)

# TODO
# Print list of lists to file
write(c("hello","world"), file = "output1.txt", sep = "\t")

# MAIN --------------------------------------------------------------------
# PFSNet style generation
# Create list of highly expressed genes
gfs_control_ds1 <- GFS(control_ds1, 0.10, 0.90)
mean_gfs_control_ds1 <- apply(gfs_control_ds1, 1, mean)
highexpr_probes <- names(mean_gfs_control_ds1[mean_gfs_control_ds1 > 0.5])
highexpr_genes <- unname(sapply(highexpr_probes, function(x) probeset_annot[x,]))
gene_list <- highexpr_genes[!sapply(highexpr_genes, function(x) grepl("///", x))]
dup_gene <- highexpr_genes[sapply(highexpr_genes, function(x) grepl("///", x))]
dup_gene1 <- unlist(strsplit(dup_gene, " /// "))
final_list <- c(gene_list, dup_gene1)

# Filter pathway with this list
fltr_df <- hsa_df[hsa_df$from %in% final_list & hsa_df$to %in% final_list,]
fltr_graph <- ftM2graphNEL(data.matrix(fltr_df))
fltr_list <- connComp(fltr_graph)
subnetworks_pfs <- fltr_list[lapply(fltr_list, length) >= 5]
connComp()
fltr_df
plot(fltr_graph)
subnetworks_pfs

nodes_high <- paste0("hsa:", final_list[final_list %in% substring(hsa@nodes, 5)])

# Alternative
subGraph(nodes_high, hsa)

# subset_nodes <- sample(nodes(hsa1_g), 25)

# sub_g1 <- igraph.from.graphNEL(sub_g)
# dcomp_sub <- decompose(sub_g1)
# 
# for (i in 1:length(dcomp_sub)){
#   print(names(V(dcomp_sub[[i]])))
# }
# 
# plot(dcomp_sub[[3]])

# PLOT --------------------------------------------------------------------

# Copy default Rgraphviz attributes
graph_attr <- getDefaultAttrs()
# Alter Rgraphviz attributes
graph_attr$graph$bgcolor <- "white"
graph_attr$node$fontsize <- 50
# Make node attributes
node_attr <- makeNodeAttrs(hsa04664, shape = "ellipse", fillcolor = "#e0e0e0")
# Plot with attributes
plot(hsa, attrs=graph_attr, nodeAttrs = node_attr)

# Plotting using igraph. Continuous vector of edges
pw250 <- as.character(as.vector(t(data.matrix(all_pwAPI[[250]]))))
pw251 <- as.character(as.vector(t(data.matrix(all_pwAPI[[251]]))))
g2 <- graph(pw250, directed = F)
plot.igraph(g2)

# Plotting using RGraphviz
g250 <- ftM2graphNEL(data.matrix(all_pwAPI[[250]]))
plot(g250)
g251 <- ftM2graphNEL(data.matrix(all_pwAPI[[251]]))
plot(g251)

plot(hsa)

# Plot star graph using igraph
star31 <- make_star(6, mode = "undirected")
plot.igraph(star31, edge.width = 3, edge.color = "black",
            vertex.label=NA, vertex.size=20)

# TODO Choose the max! if multiple probes map to same gene. Upp bound.
# TODO Remove probes that map to multiple genes! Check the number

# Load pathwayAPI
raw_pwAPI <- read.table("../info/pathwayAPI/pathwayAPI_filtered_human_pathways_entrez_id.tsv",
                        sep = "\t", header=T, strip.white = T, stringsAsFactors = F)
all_pwAPI <- split(raw_pwAPI[2:3], raw_pwAPI$pathway)
head(all_pwAPI)

# CS6216 ------------------------------------------------------------------

data <- read.table('../diff_expr/data/breast_metastasis/GSE2034/processed/data_entrez.tsv',
                    header = T, sep = "\t", row.names = 1)
genes_microarray <- rownames(data)

# Convert KEGG to edgeList
kegg_fnames <- list.files("../info/KEGG/human_pathways/")
for (fname in kegg_fnames) {
  fpath <- paste0("../info/KEGG/human_pathways/", fname)
  path_name <- substring(fname, 1, 8)
  df_hsa <- parseKGML2DataFrame(fpath)
  # If df is not NULL
  if (nrow(df_hsa) != 0) {
    # Convert KEGG ID to Entrez ID
    df_hsa$from <- substring(df_hsa$from, 5)
    df_hsa$to <- substring(df_hsa$to, 5)
    unwanted_edges <- c("phosphorylation", "indirect effect", "dephosphorylation",
                        "ubiquitination", "methylation", "state change", "missing interaction")
    # Removes unwanted edges and loops in graph
    processed_df <- df_hsa[!(df_hsa$subtype %in% unwanted_edges) & !(df_hsa$from == df_hsa$to), 1:2]
    # Only preserves nodes present in gene microarray
    df <- processed_df[processed_df$from %in% genes_microarray & processed_df$to %in% genes_microarray]
    # Removes duplicated edges
    graph <- tryCatch(graph <- ftM2graphNEL(data.matrix(df)),
                      error = function(e) {
                        # Delete duplicate edges before loading
                        concat_str <- with(df, paste0(from, to))
                        freq <- table(concat_str)
                        dup_str <- names(freq[freq > 1])
                        if (length(dup_str) == 1) {
                          row_no <- which(concat_str == dup_str)[-1]
                          # Returns a graph of the matrix with row removed
                          return (ftM2graphNEL(data.matrix(df[-row_no,])))
                        } else {
                          for (i in dup_str) {
                            # Row indices to remove
                            row_no <- which(concat_str == i)[-1]
                            # Returns a graph of the matrix with row removed
                            df <- df[-row_no,]
                          }
                          return (ftM2graphNEL(data.matrix(df)))
                        }
                      })
    # Returns a list of connected components
    ls_components <- connComp(graph)
    component_nodes <- unlist(ls_components[which.max(sapply(ls_components, length))])
    largest_subgraph <- subGraph(component_nodes, graph)
    export_graph <- igraph.from.graphNEL(largest_subgraph, name = T, weight = F, unlist.attrs = F)
    output_fpath <- paste0("info/kegg_human-edgelist/", path_name, ".tsv")
    # Save graph as edge list
    write_graph(export_graph, output_fpath, format = "ncol")
    print(output_fpath)
  }
}

df <-read.table("info/kegg_human-edgelist/edited/hsa05332.tsv",
                sep = "\t")

graph <- tryCatch(graph <- ftM2graphNEL(data.matrix(df)),
                  error = function(e) {
                    # Delete duplicate edges before loading
                    # TODO: Change column names! V1, vV2
                    concat_str <- with(df, paste0(V1,V2))
                    freq <- table(concat_str)
                    dup_str <- names(freq[freq > 1])
                    if (length(dup_str) == 1) {
                      row_no <- which(concat_str == dup_str)[-1]
                      # Returns a graph of the matrix with row removed
                      return (ftM2graphNEL(data.matrix(df[-row_no,])))
                    } else {
                      for (i in dup_str) {
                        row_no <- which(concat_str == i)[-1]
                        # Returns a graph of the matrix with row removed
                        df <- df[-row_no,]
                      }
                      return (ftM2graphNEL(data.matrix(df)))
                    }
                  })

# log <- file("log.txt", "at")  # open an output file connection
# try(ftM2graphNEL(processed_df), outFile = log)
# close(log)

