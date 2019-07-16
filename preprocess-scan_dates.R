# Match scan dates to own annotation
id_annot <- read.table("data/ovarian_cancer/GSE26712/processed/annot-shortid.tsv",
                       header = T, sep = "\t", row.names = 1)
scan_dates <- read.table("data/ovarian_cancer/GSE26712/README/scan_dates.tsv",
                         skip = 1, header = F, sep = "\t", row.names = 1)

geo_id <- unname(sapply(rownames(scan_dates), substring, 1, 9))
rownames(scan_dates) <- geo_id
scan_dates[,1] <- sapply(scan_dates[,1], substring, 1, 8)
rownames(scan_dates) <- id_annot[rownames(scan_dates),]
write.table(scan_dates, "data/ovarian_cancer/GSE26712/processed/scan_dates.tsv",
            quote = F, sep = "\t", row.names = T, col.names = F)
