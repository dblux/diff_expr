# Initialisation ----------------------------------------------------------
library(affy)
setwd("~/projects/phd/diff_expr/")

# Microarray selection ----------------------------------------------------
CEL_DIRPATH <- "data/yeoh_2002/raw/"
fpaths <- list.files(substring_head(CEL_DIRPATH, 1), full.names = T)
print(fpaths)
# # Only select sample types A and sample types B
# subset_fpaths <- fpaths[!(grepl("C[0-9]", fpaths) | grepl("D[0-9]", fpaths))]
# subset_fpaths
raw_data <- ReadAffy(filenames = fpaths)

# Quality control ---------------------------------------------------------
# Visualise each microarray
# Remove chips with too many artefacts
DIR_WPATH <- sub("raw", "img", CEL_DIRPATH)
print(DIR_WPATH)

# # Creates new directory
# dir.create(DIR_WPATH)

# Saves each microarray image in folder
for (i in 1:length(fpaths)) {
  wpath <- sprintf("%s%03d.jpg", DIR_WPATH, i)
  jpeg(wpath)
  image(raw_data[, i])
  dev.off()
}

# # Plot boxplots of probe intensities for each microarray
# par(mar = c(10,4,4,4),
#     cex.axis = 0.8)
# boxplot(raw_data, las = 2)

# # Plot density curves of microarrays
# hist(raw_data)

# Processing CEL ----------------------------------------------------------
# MAS5 without trimmed mean scaling
mas5_data_obj <- mas5(raw_data, normalize = F)

# MAS5 detection call object
mas5_calls_obj <- mas5calls(raw_data)

# # RMA without quantile normalisation
# # Use expresso function to mix and match algorithm
# rma_data <- expresso(raw_data,
#                     bgcorrect.method = "rma",
#                     normalize = F,
#                     pmcorrect.method = "pmonly",
#                     summary.method = "medianpolish")

# Extract expression data -------------------------------------------------
expr_data <- exprs(mas5_data_obj)
# Modifies column names
colnames(expr_data)
new_names <- unname(sapply(colnames(expr_data),
                           function(x) substring(x, 1, nchar(x) - 4) ))
new_names
colnames(expr_data) <- new_names

# Assigns detection calls based on default threshold
# M call: 0.04 < p-value <= 0.06
mas5_call <- exprs(mas5_calls_obj)

# If call == "P" cell <- 1
mas5_call_numeric <- (mas5_call == "P") * 1
mas5_call[1:5, 1:6]
mas5_call_numeric[1:5, 1:6]

# # Extracts pvalues from MAS5 call object
# mas5_pvalue <- assayData(raw_mas5)[["se.exprs"]]
# # Select p-value <= 0.05
# mas5_filter <- (mas5_pvalue <= 0.1) * 1

# mas5 data preserving only "P" calls
mas5_confident <- expr_data * mas5_call_numeric

OUTPUT_FPATH <- paste0(substring(CEL_DIRPATH, 1, nchar(CEL_DIRPATH) - 8),
                "processed/mas5_TEL.tsv")
print(OUTPUT_FPATH)
# Save expression data
write.table(mas5_confident, OUTPUT_FPATH,
            sep = "\t", quote = F)
