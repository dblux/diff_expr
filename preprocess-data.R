# Pre-process data

# Download KEGG database
# hsa_pathways <- read.table("../info/KEGG/kegg-human_pathways-name.txt", header = T, sep = "\t")
# 
# for (i in hsa_pathways$Pathway) {
#   pw <- substr(i, 6, 13)
#   retrieveKGML(pw, "hsa", paste0("~/projects/phd/info/KEGG/human_pathways/", pw, ".xml"))
# }

# OVARIAN DATASET2 --------------------------------------------------------
ovarian_dataset2 <- read.table('data/ovarian_cancer/GSE26712/processed/filtered_original.tsv', header = T, row.names = 1)

ds2_annot <- read.table("data/ovarian_cancer/GSE26712/README/annot.tsv",
                        header=T, row.names=2, stringsAsFactors=F)

sample_names <- sapply(colnames(ovarian_dataset2), function(x) as.character(ds2_annot[x,]))
colnames(ovarian_dataset2) <- sample_names
ovarian_dataset2
write.table(df1, "data/ovarian_cancer/GSE18521/processed/processed_original1.tsv", sep="\t", col.names = NA)

# ALL DATASET -------------------------------------------------------------
patient <- read.table('data/yeoh_2002/processed/filtered_TEL-AML1.tsv', header = T, row.names = 1)

# Replace rows with 'A' and 'M' with 0s
for (i in 1:(ncol(control)/2)) {
  i <- i*2
  control[(control[,i] == 'A' | control[,i] == 'M'), i-1] <- 0
}

for (i in 1:(ncol(patient)/2)) {
  i <- i*2
  patient[(patient[,i] == 'A' | patient[,i] == 'M'), i-1] <- 0
}

# Removes gene call columns
control <- control[,seq(1,18,2)]
patient <- patient[seq(1,66,2)]

# Save control and patients dataframe
write.table(control, "data/yeoh_2002/processed/processed_control.tsv", sep="\t")
write.table(patient, "data/yeoh_2002/processed/processed_TEL-AML1.tsv", sep="\t")
