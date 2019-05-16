# Create flybase ID to fold change labels ---------------------------------
raw_clone <- read.table("data/platinum_spike/GSE21344/README/GSE21344_allclones_fcvalue.txt",
                         header = F, skip = 3, stringsAsFactors = F)
# Clone ID to flybase ID annotation
clone_flybase <- read.table("data/platinum_spike/GSE21344/README/annot-clone_flybaseclone_flybasegene.txt",
                            header = T, stringsAsFactors = F, strip.white = T)

colnames(raw_clone) <- c("submitted_id", "fold_change")

# 13 duplicates of clones exist in the labels
# Fold change of 1.7 was erroneously added
find_multiple <- function(vec) {
  return(names(table(vec))[table(vec) > 1])
}

duplicated_clones <- find_multiple(raw_clone$submitted_id)
print(raw_clone[raw_clone$submitted_id %in% duplicated_clones,])

# Output to file duplicated clones
# sink("dump/platinum_spike_duplicated_clones.txt")
# cat(duplicated_clones)
# sink()

# Remove duplicated clones with fold change of 1.7
get_row_duplicate <- function(string) {
  return(which(raw_clone$submitted_id == string
               & raw_clone$fold_change == 1.7))
}

rownum_duplicate <- sapply(duplicated_clones, get_row_duplicate)
# Remove duplicates with fold change 1.7
raw_clone1 <- raw_clone[-rownum_duplicate,]

# Some flybase clone ID has maps to multiple flybase genes!
# Perform a right outer join to preserve all foldchange values
raw_gene <- merge(clone_flybase[,c(1,3)], raw_clone1,
                  by = "submitted_id", all.y = T)
# Some flybase clone ID does not map to any flybase gene!
flybase_gene_fc <- raw_gene[!(raw_gene$converted_id == "-"), c(2,3)]
print(paste0("No. of clones with no corresponding gene ID: ",
             sum(raw_gene$converted_id == "-")))

find_multiple(flybase_gene_fc$converted_id)

write.table(flybase_gene_fc, "data/platinum_spike/GSE21344/processed/flybase_gene-foldchange.tsv",
            sep = "\t", quote = F, row.names = F)

# Above steps are redundant!!!
# Additional file provided with probeset labels
raw_foldchange <- read.table("data/platinum_spike/GSE21344/README/Zhu_2010-Platinum_Spike_add5.txt",
                            header = F, row.names = 1, stringsAsFactors = F, strip.white = T)
# Should not filter out ambiguous and AFFY probesets from dataframe!
# Rowname of affymetrix probesets
probeset_mask <- list("MC", "MF")
probeset_foldchange <- raw_foldchange[!(raw_foldchange[,1] %in% probeset_mask), , drop = F]

write.table(probeset_foldchange, "data/platinum_spike/GSE21344/processed/probeset_foldchange_full.tsv",
            sep = "\t", quote = F, col.names = F, row.names = T)
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
