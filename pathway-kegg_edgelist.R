setwd("~/projects/phd/diff_expr/")
library(KEGGgraph)
library(Rgraphviz)
library(igraph)
library(RColorBrewer)
source("functions.R")

data <- read.table('data/raw/data_entrez1.tsv',
                   header = T, sep = "\t", row.names = 1)
genes_microarray <- rownames(data)

# Convert KEGG to edgeList
kegg_fnames <- list.files("info/kegg_human-pathways/")
for (fname in kegg_fnames) {
  fpath <- paste0("info/kegg_human-pathways/", fname)
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
    df <- processed_df[processed_df$from %in% genes_microarray & processed_df$to %in% genes_microarray,]
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
                            concat_str <- concat_str[-row_no]
                          }
                          return (ftM2graphNEL(data.matrix(df)))
                        }
                      })
    # If graph is empty
    if (length(graph@nodes) == 0) {
      next
    }
    # Returns a list of connected components
    ls_components <- connComp(graph)
    component_nodes <- unlist(ls_components[which.max(sapply(ls_components, length))])
    largest_subgraph <- subGraph(component_nodes, graph)
    export_graph <- igraph.from.graphNEL(largest_subgraph, name = T, weight = F, unlist.attrs = F)
    output_fpath <- paste0("/home/dblux/projects/phd/cs6216/info/kegg_human-edgelist/",
                           path_name, ".tsv")
    # Save graph as edge list
    write_graph(export_graph, output_fpath, format = "ncol")
    print(output_fpath)
  }
}

# options(error=recover)
# options(error=NULL)
