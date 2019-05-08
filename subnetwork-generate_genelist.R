# Import breast metastasis data
metastasis_data <- read.table("data/breast_metastasis/GSE2034/processed/data_labelled.tsv",
                              header = T, row.names = 1)
# Import subnetworks
subnetworks <- read.table("data/subnetwork/hsa-nea/subnetworks-breast_metastasis.tsv",
                          header = T)
# Create gene list for each subnetwork
list_subnetworks <- split(subnetworks[2:3], subnetworks$subnetwork_name)
subnetwork_genelist <- lapply(list_subnetworks, function(x) sort(unique(unname(unlist(x)))))

# Save subnetwork gene list
subnetwork_matrix <- lapply(subnetwork_genelist, as.data.frame)
subnetwork_df <- do.call(rbind, subnetwork_matrix)
subnetwork_df["subnetwork"] <- sub("\\..*$", "", rownames(subnetwork_df))
colnames(subnetwork_df)[1] <- "gene"
subnetwork_df1 <- subnetwork_df[,2:1]
rownames(subnetwork_df1) <- NULL
write.table(subnetwork_df1,
            "data/subnetwork/hsa-nea/subnet_genelist-breast-metastasis.tsv",
            col.names = T, row.names = F, quote = F)