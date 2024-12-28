library("data.table")

#Declare the path of the global CoRE-BED elements directory
source_dir <- "/home/bettimj/gamazon_rotation/mod_core-bed/primed_promoter_whole_blood"
out_dir <- "/home/bettimj/gamazon_rotation/mod_core-bed/primed_promoter_whole_blood"

tissues <- c("hsc_and_b_cell")

#Loop through the gene-annotated bivalent promoter BED for each tissue, outputting only the unique genes as new files
for (tissue in tissues) {
	file <- fread(paste(source_dir, paste0("genes_", tissue, "_class4.bed"), sep = "/"), header = FALSE, quote = FALSE, sep = "\t")
	df <- as.data.frame(file)
	df <- unique(df[,11])
	write.table(df, file = paste(out_dir, paste0("unique_genes_", tissue, "_class4.bed"), sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

