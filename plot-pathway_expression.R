setwd("~/projects/phd/diff_expr/")
library(KEGGgraph)
library(Rgraphviz)
library(igraph)
library(RColorBrewer)
source("functions.R")

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

# Arguments: Named intensity values for each node
# Returns list for nAttrs in Rgraphviz
generate_node_attr <- function(intensity_vec, graph_nodes) {
  # Subset intensities of nodes present in graph
  node_intensities <- intensity_vec[graph_nodes]
  # Nodes not present in microarray (with NA values) assigned red font
  node_fontcolour <- node_intensities
  node_fontcolour[!is.na(node_fontcolour)] <- "black"
  node_fontcolour[is.na(node_fontcolour)] <- "red"
  
  # Nodes not present in microarray will be white in colour
  node_intensities[is.na(node_intensities)] <- 0
  # Minmax symmetrical limits for scaling
  top_limit <- max(max(intensity_vec), abs(min(intensity_vec)))
  # Scale intensity values according to top_limit (minmax)
  scaled_vec <- (node_intensities + top_limit)/(top_limit + top_limit)
  colour_gradient <- colorRamp(c("red", "white", "green"))
  node_colour <- rgb(colour_gradient(scaled_vec),
                     maxColorValue = 255)
  names(node_colour) <- graph_nodes
  names(node_fontcolour) <- graph_nodes
  node_attr <- list(fillcolor = node_colour, fontcolor = node_fontcolour)
  return(node_attr)
}

# Plots single KEGG pathway with overlaid diffexpr data
plot_kegg_expr <- function(kegg_df, kegg_id, diffexpr_vec) {
  kegg_graph <- ftM2graphNEL(data.matrix(kegg_df))
  node_attr <- generate_node_attr(diffexpr_vec, nodes(kegg_graph))
  plot(kegg_graph, attrs = graph_attr,
       nodeAttrs = node_attr, main = kegg_id)
  kegg_plot <- recordPlot()
  file_wpath <- sprintf("dump/%s.pdf", kegg_id)
  save_fig(kegg_plot, file_wpath, 20, 10)
}

ovarian_data <- read.table("data/ovarian_cancer/GSE18521/processed/mas5_qnorm.tsv",
                           header = T, row.names = 1)
ovarian_A <- ovarian_data[,5:10]
ovarian_B <- ovarian_data[,53:58]
tstat_vec <- row_ttest(ovarian_A, ovarian_B, "tstat")

# Copy default Rgraphviz attributes
graph_attr <- getDefaultAttrs()
# Alter Rgraphviz attributes
graph_attr$graph$bgcolor <- "white"
graph_attr$node$fontsize <- 45

list_kegg <- import_kegg("../info/KEGG/hsa_df/")
pathway_size <- sapply(list_kegg, nrow)
sort(unname(pathway_size))
subset_kegg <- list_kegg[pathway_size <= 1000]

# Have to pass in tstat_vec as list of one element
# Otherwise tstat_vec will be accessed one by one
# Saves plot of pathways
mapply(plot_kegg_expr, subset_kegg, names(subset_kegg), list(tstat_vec))

