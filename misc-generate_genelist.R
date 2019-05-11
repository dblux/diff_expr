setwd("~/projects/phd/cs6216/")

EDGELIST_DIR <- "info/kegg_human-edgelist/large/"
selected_pathways <- list.files(EDGELIST_DIR)

num_nodes <- function(fpath) {
  data <- read.table(paste0(EDGELIST_DIR, fpath), sep=" ",
                     header = F, strip.white = T)
  return(unique(unlist(data)))
}

# Returns vector of unique nodes for each pathway
pathway_nodes <- sapply(selected_pathways, num_nodes)

GENELIST_DIR <- "data/pathway/gene_list/"
fnames <- names(pathway_nodes)
for (i in 1:length(pathway_nodes)) {
  # Save gene list for each pathway
  write(pathway_nodes[[i]],
        file = paste0(GENELIST_DIR, fnames[i]),
        ncolumns = 1)
}

# Total unique nodes in all pathways
pathway_total_nodes <- unique(unname(unlist(pathway_nodes)))
write(pathway_total_nodes,
      file = "data/pathway/all_pathway_nodes.tsv",
      ncolumns = 1)

# Total unique nodes in all pathways: 1983
length(pathway_total_nodes)

# pathway_num_nodes <- sapply(selected_pathways, num_nodes)
# save(pathway_num_nodes, file = "data/pathway_num_nodes.rda")
