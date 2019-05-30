# Command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Converts dataframe of subnetworks to dataframe of gene sets
subnetwork2geneset <- function(subnetworks_rpath, geneset_wpath) {
  # Import subnetworks
  subnetworks <- read.table(subnetworks_rpath, header = T)
  print(head(subnetworks))
  # Create gene list for each subnetwork split on COLUMN
  subnetworks_list <- split(subnetworks[2:3], subnetworks$subnetwork_id)
  geneset_list <- lapply(subnetworks_list, function(x) sort(unique(unname(unlist(x)))))
  # Save subnetwork gene list
  subnetwork_matrix <- lapply(geneset_list, as.data.frame)
  subnetwork_df <- do.call(rbind, subnetwork_matrix)
  subnetwork_df$subnetwork_id <- sub("\\..*$", "", rownames(subnetwork_df))
  # Change name of first column
  colnames(subnetwork_df)[1] <- "gene"
  # Remove rownames of df
  rownames(subnetwork_df) <- NULL
  # Save reordered df
  write.table(subnetwork_df[,2:1], geneset_wpath,
              col.names = T, row.names = F, quote = F)
}

# Run main function
subnetwork2geneset(args[1], args[2])
