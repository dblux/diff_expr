library(Rgraphviz)
library(igraph)
setwd("~/projects/phd/cs6216/")

# Check for pairs of directed edges
df <- read.table("info/kegg_human-edgelist/large/hsa04917.tsv",
                 header = F, sep = " ")
dir_one <- with(df, paste0(V1, V2))
dir_two <- with(df, paste0(V2, V1))
freq <- table(c(dir_one, dir_two))
print(sum(freq > 1))
print(sum(df$V1 == df$V2))
print(sum(table(dir_one) > 1))

G <- ftM2graphNEL(data.matrix(df))
plot(G)
