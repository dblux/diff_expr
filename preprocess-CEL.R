library(affy)
setwd("~/projects/phd/diff_expr/")

# Microarray selection ----------------------------------------------------
CEL_DIRPATH <- "data/platinum_spike/GSE21344/raw/"
fpaths <- list.files(substring(CEL_DIRPATH, 1, nchar(CEL_DIRPATH) - 1),
                     full.names = T)
print(fpaths)
# Only select sample types A and sample types B
# subset_fpaths <- fpaths[!(grepl("C[0-9]", fpaths) | grepl("D[0-9]", fpaths))]
raw_data <- ReadAffy(filenames = fpaths)

# Quality control ---------------------------------------------------------
# Visualise each microarray
# Remove chips with too many artefacts
image(raw_data)
# Plot distribution of probe intensities for each microarray
par(mar = c(10,4,4,4),
    cex.axis = 0.8)
boxplot(raw_data, las = 2)
# hist(raw_data)

# Processing --------------------------------------------------------------
# Use mas5 algorithm to process CEL
# Includes all components of mas5
mas5_data <- mas5(raw_data)

# Use expresso function to mix and match algorithm
# Example: RMA without quantile normalisation
rma_data <- expresso(raw_data,
                    bgcorrect.method = "rma",
                    normalize = F,
                    pmcorrect.method = "pmonly",
                    summary.method = "medianpolish")

# Extract expression data -------------------------------------------------
expr_data <- exprs(mas5_data)
colnames(expr_data) <- substring(colnames(expr_data), 1, 9)
colnames(expr_data)

# mas5 detection calls
raw_mas5 <- mas5calls(raw_data)
mas5_call <- exprs(raw_mas5)
mas5_pvalue <- assayData(raw_mas5)[["se.exprs"]]
colnames(mas5_call) <- substring(colnames(mas5_call), 1, 9)
colnames(mas5_call)

# M call: 0.04 < p-value <= 0.06
# If call == "P" cell <- 1
mas5_call_numeric <- (mas5_call == "M") * 1
# Select p-value <= 0.05
mas5_filter <- (mas5_pvalue <= 0.1) * 1

# mas5 data preserving only "P" calls
# mas5_confident <- expr_data * mas5_call_numeric
mas5_confident <- expr_data * mas5_filter

OUTPUT_FPATH <- paste0(substring(CEL_DIRPATH, 1, nchar(CEL_DIRPATH) - 4),
                "processed/mas5_original010.tsv")
print(OUTPUT_FPATH)
# Save expression data
write.table(mas5_confident, OUTPUT_FPATH,
            sep = "\t", quote = F)
