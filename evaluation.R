
# FALSE POSITIVE RATE -----------------------------------------------------

# MAQC data set
maqc <- read.table('data/MAQC-I/processed/filtered_original.tsv', header = T, row.names = 1)
colnames(maqc)
a_maqc <- maqc[,c(1:5,21:25,41:45,61:65,81:85,101:105)]
b_maqc <- maqc[,c(6:10,26:30,46:50,66:70,86:90,106:110)]

a1_col <- c(which(endsWith(colnames(a_maqc),"A1")),
            which(endsWith(colnames(a_maqc),"A2")),
            which(grepl("[1-3]_A3", colnames(a_maqc))))
a2_col <- setdiff(1:30, a1_col)