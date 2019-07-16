library(igraph)

df <- data.frame(from = 1:5, to = 2:6, type = LETTERS[1:5])
g1 <- graph_from_data_frame(df)
g2 <- add.edges(g1, c(1,2), attr = list(type="F"))

set.seed(1)
# Assign colour to vertices
V(g2)$color <- rep(c("green", "red"), 3)
# Plots curved edges to prevent overlap of multiple edges (layout)
plot(g2,
     layout = layout_nicely(g2),
     vertex.size = 20,
     label.cex = 1,
     edge.label = unlist(edge_attr(g2)),
     edge.arrow.size = 0.5,
     edge.width = 2)

# Save igraph as df
as_data_frame(g2)

V(g1)
edge_attr(g1)
vertex_attr(g1)

# Decompose graphs

# Parse KGML and compare with KEGG png

# Parse processed KEGG df
processed_kegg <- read.table("../info/KEGG/hsa_df/hsa04012.tsv",
                             sep = "\t", header = T)
g_processed <- graph_from_data_frame(processed_kegg)
plot(g_processed,
     layout = layout_as_tree(g_processed, mode = "all"),
     vertex.size = 5,
     vertex.label.size = 1,
     edge.color = "black",
     edge.arrow.size = 0.1,
     edge.width = 1)

# Inspect t-statistics of genes
# Inspect no. of genes represented by graphs
qnorm(0.95)
curve(dnorm, xlim = c(-3,3))
abline(v = -1.64)
pnorm(-1.64)

# Finalise idea in generating subnetworks
# Preliminary exploration on a few example pathways
# Differentiate between metabolic and signalling pathways


# Preliminary study -------------------------------------------------------
# Imports data and returns two lists of top genes and gene weight matrix
ovarian_data1 <- read.table("data/ovarian_cancer/GSE18521/processed/mas5_qnorm.tsv",
                            header = T, row.names = 1)
ovarian1_A <- ovarian_data1[,5:10]
ovarian1_B <- ovarian_data1[,53:58]

kegg_allpathway <- read.table("../info/KEGG/kegg-human_allpathway_genes.tsv",
                              sep = "\t", header = T)
kegg_genes <- unique(kegg_allpathway$Entrez.Gene.ID)
length(kegg_genes)
microarray_genes <- rownames(ovarian_data1)
length(intersect(microarray_genes, kegg_genes))

genes_freq <- table(kegg_allpathway$Entrez.Gene.ID)
sum(genes_freq == 1)


