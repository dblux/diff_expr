setwd("~/projects/phd/diff_expr/")
ovarian_data1 <- read.table("data/ovarian_cancer/GSE18521/processed/GSE18521_entrez.tsv",
                            header = T, row.names = 1)
# Remove GSM461376
OUTPUT_PATH <- "data/ovarian_cancer/GSE18521/processed/GSE18521_entrez1.tsv"
write.table(ovarian_data1[,-39], OUTPUT_PATH,
            sep = "\t", quote = F)
