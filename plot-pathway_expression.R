setwd("~/projects/phd/diff_expr/")
library(VennDiagram)
library(ggplot2)
library(reshape2)
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

num_genes <- length(unique(unlist(lapply(list_kegg, unlist))))
paste0("No. of genes represented by KEGG pathways: ", num_genes)

# Have to pass in tstat_vec as list of one element
# Otherwise tstat_vec will be accessed one by one
# Saves plot of pathways
mapply(plot_kegg_expr, subset_kegg, names(subset_kegg), list(tstat_vec))

# DATA EXPLORATION --------------------------------------------------------
ovarian_data1 <- read.table("data/ovarian_cancer/GSE18521/processed/mas5_qnorm.tsv",
                           header = T, row.names = 1)
ovarian_A1 <- ovarian_data1[,5:10]
ovarian_B1 <- ovarian_data1[,53:58]

ovarian_data2 <- read.table("data/ovarian_cancer/GSE26712/processed/mas5_qnorm.tsv",
                           header = T, row.names = 1)
ovarian_A2 <- ovarian_data2[,1:6]
ovarian_B2 <- ovarian_data2[,43:48]

# Argument: Non-log2 transformed df with multiple patient columns
plot_pdf <- function(df) {
  # Melt dataframe
  melt_df <- melt(log2_transform(df), variable.name = "ID")
  # fltr_df <- melt_df[melt_df$value > 1,]
  # Plot density curve
  pdf <- ggplot(melt_df, aes(x = value, col = ID)) + 
    geom_density(show.legend = T) +
    xlim(0, 16)
  return(pdf)
}

norm_cdf <- function(df) {
  for (c in 1:ncol(df)) {
    notzero <- df[,c] != 0
    df[,c][notzero] <- rank(df[,c][notzero])
    df[,c] <- df[,c]/sum(notzero)
  }
  return(df)
}

generate_subnetworks <- function(df_A, df_B, filepath_or_dirpath, flag = "dataframe", min_size = 5) {
  plot_coloured_graph <- function(igraph) {
    graphNEL <- igraph.to.graphNEL(igraph)
    graph_nodes <- nodes(graphNEL)
    node_colour <- rep("#FFFFFF", length(graph_nodes))
    names(node_colour) <- graph_nodes
    node_colour[graph_nodes %in% hionly_A] <- "#99CCCC"
    node_colour[graph_nodes %in% hionly_B] <- "#FFCCCC"
    node_colour[graph_nodes %in% both_hi_DE] <- "#66CC00"
    node_colour[graph_nodes %in% both_hi_notDE] <- "#CCFFCC"
    font_colour <- ifelse(graph_nodes %in% array_genes, "black", "red")
    names(font_colour) <- graph_nodes
    node_attr <- list(fillcolor = node_colour, fontcolor = font_colour)
    # Plot subnetwork
    plot(graphNEL, attrs = graph_attr, nodeAttrs = node_attr)
  }
  
  subset_data <- cbind(df_A, df_B)
  norm_subset <- norm_cdf(subset_data)
  # Highly expressed: All non-zeros
  log_df <- norm_subset != 0
  present_A <- norm_subset[,1:6][rowSums(log_df[,1:6]) == 6,]
  present_B <- norm_subset[,7:12][rowSums(log_df[,7:12]) == 6,]
  highexpr_A <- rownames(present_A)[rowMeans(present_A) > 0.5]
  highexpr_B <- rownames(present_B)[rowMeans(present_B) > 0.5]
  # jacc_coeff(highexpr_A, highexpr_B)
  
  hionly_A <- setdiff(highexpr_A, highexpr_B)
  hionly_B <- setdiff(highexpr_B, highexpr_A)
  both_hi <- intersect(highexpr_A, highexpr_B)
  array_genes <- rownames(ovarian_data)
  both_lo <- intersect(setdiff(array_genes, highexpr_A),
                       setdiff(array_genes, highexpr_B))
  
  log_fc <- calc_logfc(subset_data[both_hi,1:6], subset_data[both_hi,7:12])
  both_hi_DE <- names(log_fc)[abs(log_fc) > 1]
  both_hi_notDE <- setdiff(both_hi, both_hi_DE)
  
  # # Import pathwayAPI
  # pwapi_df <- read.table("../info/pathwayAPI/pwapi_id_human-filtered_entrez.tsv",
  #                        sep = "\t", header = T)
  # list_pwapi <- split(pwapi_df[,2:3], pwapi_df[,1])
  # # 5,6,7
  # N <- 7
  # pwapi_graph <- ftM2graphNEL(data.matrix(list_pwapi[[N]]))
  
  # Copy default Rgraphviz attributes
  graph_attr <- getDefaultAttrs()
  # Alter Rgraphviz attributes
  graph_attr$graph$bgcolor <- "white"
  graph_attr$node$fontsize <- 45
  
  # Import KEGG
  list_kegg <- import_kegg("../info/KEGG/hsa_df/")
  # pathway_size <- sapply(list_kegg, nrow)
  # sort(unname(pathway_size))
  # subset_kegg <- list_kegg[pathway_size <= 1000]

  # Generate subnetworks
  array_genes <- rownames(ovarian_data)
  # Open file connection
  if (flag == "dataframe") allsubnetwork_fileconn <- file(filepath_or_dirpath, "w")
  # showConnections(all = T)
  # closeAllConnections()
  for (i in 1:length(list_kegg)) {
    pathway_id <- names(list_kegg)[i]
    # igraph: if edgelist is numeric type will be treated as sequential ID
    kegg_matrix <- data.matrix(list_kegg[[i]])
    mode(kegg_matrix) <- "character"
    kegg_igraph <- graph_from_edgelist(kegg_matrix)
    lolo_vertices <- V(kegg_igraph)$name[V(kegg_igraph)$name %in% both_lo]
    # Remove all LL vertices
    subgraph_1 <- delete_vertices(kegg_igraph, lolo_vertices)
    list_conncomp_1 <- decompose(subgraph_1, min.vertices = min_size)
    j <- 0
    for (conncomp_1 in list_conncomp_1) {
      # Missing vertices that have either in-degree or out-degree = 0 are recursively removed
      repeat {
        in_degree <- degree(conncomp_1, mode = "in")
        out_degree <- degree(conncomp_1, mode = "out")
        filter_logvec <- !names(in_degree) %in% array_genes & (in_degree == 0 | out_degree == 0)
        peripheral_missing <- names(in_degree)[filter_logvec]
        # print(peripheral_missing)
        if (length(peripheral_missing) == 0) {
          break
        }
        print("Deleting vertices..")
        conncomp_1 <- delete_vertices(conncomp_1, peripheral_missing)
      }
      if (is_connected(conncomp_1)) {
        # If resulting subgraph is less than minimum size: NEXT iteration
        if (gorder(conncomp_1) < min_size) next
        j <- j + 1
        if (flag == "dataframe") {
          # Save subnetwork
          subnetwork_name <- sprintf("%s_%s", pathway_id, j)
          subnetwork_matrix <- cbind(subnetwork_name, as_edgelist(conncomp_1))
          write.table(subnetwork_matrix, allsubnetwork_fileconn,
                      quote = F, append = T, sep = "\t", row.names = F, col.names = F)
        } else if (flag == "plot") {
          print(gorder(conncomp_1))
          if (gorder(conncomp_1) < 200) {
            # Plot subnetwork
            plot_coloured_graph(conncomp_1)
            kegg_plot <- recordPlot()
            file_wpath <- sprintf("%s%s_%s.pdf", filepath_or_dirpath, pathway_id, j)
            save_fig(kegg_plot, file_wpath, 20, 10)
          }
        }
      } else {
        # All decomposed subgraphs are larger than min size
        list_conncomp_2 <- decompose(conncomp_1, min.vertices = min_size)
        for (conncomp_2 in list_conncomp_2) {
          j <- j + 1
          if (flag == "dataframe") {
            # Save subnetwork
            subnetwork_name <- sprintf("%s_%s", pathway_id, j)
            subnetwork_matrix <- cbind(subnetwork_name, as_edgelist(conncomp_2))
            write.table(subnetwork_matrix, allsubnetwork_fileconn,
                        quote = F, append = T, sep = "\t", row.names = F, col.names = F)
          } else if (flag == "plot") {
            print(gorder(conncomp_2))
            if (gorder(conncomp_2) < 200) {
              # Plot subnetwork
              plot_coloured_graph(conncomp_2)
              kegg_plot <- recordPlot()
              file_wpath <- sprintf("%s%s_%s.pdf", filepath_or_dirpath, pathway_id, j)
              save_fig(kegg_plot, file_wpath, 20, 10)
            }
          }
        }
      }
    }
  }
  # Close file connection
  if (flag == "dataframe") close(allsubnetwork_fileconn)
}

generate_subnetworks(ovarian_A2, ovarian_B2,
                     "dump/subnetwork/kegg_GSE26712/",
                     flag = "plot")

# Evaluation --------------------------------------------------------------

# Plotting ----------------------------------------------------------------
hist(log2_transform(ovarian_data[,5]), breaks = 70)
hist(rowMeans(present_A), breaks = 50)
hist(apply(present_subset, 1, median), breaks = 100)
hist(rowMeans(present_subset), breaks = 100)
rowMeans(present_subset)
rownames(log_df)[rowSums(log_df) == 6]

# Violin plot of mean normalised values against no. of zeros
df_plot <- as.data.frame(cbind(rowSums(log_df), rowMeans(norm_subset[,1:6])))
colnames(df_plot) <- c("num_zero", "avg")
point <- data.frame(a = factor(0:6), b = 0:6/6)

ggplot(df_plot, aes(factor(num_zero), avg)) +
  geom_violin(aes(fill = factor(num_zero))) +
  geom_point(data = point, aes(x = a, y = b))

# No. of non-expressed genes
num_zero <- apply(log2_transform(subset_data), 2, function(col) sum(col < 1))

density_plot <- plot_pdf(subset_data) +
  ggtitle("Density plot (expression values < 1 removed): GSE18521")

ggsave("dump/density_removelt1-GSE18521.png", density_plot)

curve(pnorm, xlim = c(-5,5))
x <- abs(c(rnorm(10000, mean = 5), rnorm(5000, mean = 10, sd = 2)))
df <- data.frame(a = x, b = 0)

hist(x, breaks = 100)
head(df)
y <- norm_cdf(df)
hist(y[,1], breaks = 100)

# Simulated data
rand_int <- sample(0:100, 1000, replace = T, prob = c(0.5, rep(0.005, 100)))
sim_data <- matrix(rand_int, nrow = 100, ncol = 10)

hist(log2_transform(ovarian_data[,1]), breaks = 100)
tstat_vec <- row_ttest(ovarian_A, ovarian_B, "tstat")

pvalue <- row_ttest(subset_data[1:6], subset_data[7:12])
tstat <- row_ttest(subset_data[1:6], subset_data[7:12], "tstat")
all_logfc <- calc_logfc(subset_data[1:6], subset_data[7:12])
hist(pvalue, breaks = 50)
plot(density(all_logfc))
plot(density(tstat), xlim = c(-5,5))

