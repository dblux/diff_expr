library(affy)
setwd("~/projects/phd/diff_expr/")

# Microarray selection ----------------------------------------------------
CEL_DIRPATH <- "data/MAQC-I/raw"
fpaths <- list.files(CEL_DIRPATH, full.names = T)
# Only select sample types A and sample types B
subset_fpaths <- fpaths[!(grepl("C[0-9]", fpaths) | grepl("D[0-9]", fpaths))]

raw_data <- ReadAffy(filenames = subset_fpaths)
attributes(raw_data)

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
# Extract expression data to dataframe
expr_data <- exprs(rma_data)
colnames(expr_data) <- substring(colnames(expr_data), 1, 8)
RMA_FPATH <- "data/MAQC-I/processed/rma_original.tsv"
# Save expression data
write.table(expr_data, RMA_FPATH,
            sep = "\t", quote = F)
