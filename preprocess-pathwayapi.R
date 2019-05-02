# PATHWAYAPI --------------------------------------------------------------
# Load pathwayAPI
raw_pwAPI <- read.table("../info/pathwayAPI/pathwayAPI_human_pathways_entrez_id.txt", header=T, sep="\t")

splt_raw_pwAPI <- split(raw_pwAPI, raw_pwAPI$pathway)
ls_raw_pwAPI <- lapply(splt_raw_pwAPI, droplevels)

# Remove pathways with genes that are not in Entrez ID format
index <- vector()
for (i in 1:length(ls_raw_pwAPI)) {
  tryCatch(
    {
      as.numeric(levels(ls_raw_pwAPI[[i]]$gene1))
      as.numeric(levels(ls_raw_pwAPI[[i]]$gene2))
    }, warning = function(w) {
      index <<- c(index,i)
    }
  )
}

fltr_pwAPI <- ls_raw_pwAPI[-index]
new_pwAPI <- data.frame(Pathway=factor(), Gene1=factor(), Gene2=factor())
for (i in fltr_pwAPI) {
  new_pwAPI <- rbind(new_pwAPI, i)
}

write.table(new_pwAPI, "../info/pathwayAPI/pathwayAPI_filtered_human_pathways_entrez_id.tsv",
            sep = "\t", row.names=F, quote=F)

# Convert pathway name to ID
raw_pwAPI <- read.table("../info/pathwayAPI/pathwayAPI_filtered_human_pathways_entrez_id.tsv",
                        header=T, sep="\t", stringsAsFactors = F)
list_pwAPI <- split(raw_pwAPI[,2:3], raw_pwAPI$pathway)
id_pwAPI <- paste0("pwapi", sprintf("%.3d", 1:length(list_pwAPI)))
# Save mapping from pathway names to ID
sink("dump/pwapi_id.txt", append = F)
cat(paste(names(list_pwAPI), id_pwAPI, sep = "\t"),
    fill = 2)
sink()
# Change pathway names to ID
names(list_pwAPI) <- id_pwAPI
df_id <- do.call(rbind, list_pwAPI)
df_id$pathway_id <- substring(rownames(df_id), 1, 8)
df_pwapi <- df_id[,c(3,1,2)]
PWAPI_FPATH <- "../info/pathwayAPI/pwapi_id_human-filtered_entrez.tsv"
write.table(df_pwapi, PWAPI_FPATH, sep = "\t",
            row.names = F, quote = F)
