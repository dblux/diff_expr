# Create flybase ID to fold change labels ---------------------------------
raw_clone <- read.table("data/platinum_spike/GSE21344/README/GSE21344_allclones_fcvalue.txt",
                         header = F, skip = 3, stringsAsFactors = F)
# Clone ID to flybase ID annotation
clone_flybase <- read.table("data/platinum_spike/GSE21344/README/annot-clone_flybaseclone_flybasegene.txt",
                            header = T, stringsAsFactors = F, strip.white = T)

colnames(raw_clone) <- c("flybase_id", "fold_change")

# 13 duplicates of clones exist in the labels
# Fold change of 1.7 was erroneously added
find_multiple <- function(vec) {
  return(names(table(vec))[table(vec) > 1])
}
duplicated_clones <- find_multiple(raw_clone$flybase_id)
# Output to file duplicated clones
sink("dump/platinum_spike_duplicated_clones.txt")
cat(duplicated_clones)
sink()

# Remove duplicated clones with fold change of 1.7
get_row_duplicate <- function(string) {
  return(which(raw_clone$flybase_id == string
               & raw_clone$fold_change == 1.7))
}

rownum_duplicate <- sapply(duplicated_clones, get_row_duplicate)
# Remove duplicates with fold change 1.7
raw_clone1 <- raw_clone[-rownum_duplicate,]

# Some flybase clone ID does not map to any flybase gene!
rm_clone <- clone_flybase$submitted_id[clone_flybase$converted_id == "-"]
print(paste0("No. of clones removed with no corresponding flybase gene: ",
             length(rm_clone)))
# Remove clones that have no corresponding flybase gene
raw_clone2 <- raw_clone1[!raw_clone1$flybase_id %in% rm_clone,]

# Some flybase clone ID has maps to multiple flybase genes!
# Add new rows if multiple flybase genes present
multiple_clone <- find_multiple(clone_flybase$submitted_id)
single_raw <- raw_clone2[!raw_clone2$flybase_id %in% multiple_clone,]

# TODO

raw_clone1[,1] <- flybase_id

write.table(raw_clone1, "data/platinum_spike/GSE21344/processed/flybase_gene-foldchange.tsv",
            sep = "\t", quote = F, row.names = F)

# Create GEO accession to sample title annot ------------------------------
annot_sample <- read.table("data/platinum_spike/GSE21344/README/annot-sample_title_geo.tsv",
                           row.names = 1, stringsAsFactors = F)
annot_sample1 <- t(annot_sample)
sample_name <- paste0(rep(c("A", "B"), each = 9),
                      rep(rep(1:3, each = 3), 2),
                      rep(c("X", "Y", "Z"), 6))
annot_arr <-cbind(annot_sample1, sample_name)
annot_arr1 <- annot_arr[, c(2,3,1)]
colnames(annot_arr1) <- c("sample_geo", "sample_name", "sample_title")

write.table(annot_arr1, "data/platinum_spike/GSE21344/processed/annot-sample_name.tsv",
            sep = "\t", quote = F, row.names = F)
