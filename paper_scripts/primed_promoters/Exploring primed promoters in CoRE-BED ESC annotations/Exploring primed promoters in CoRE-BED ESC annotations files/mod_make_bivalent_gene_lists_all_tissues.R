library("data.table")

#Declare the path of the global CoRE-BED elements directory
source_dir <- "/home/bettimj/gamazon_rotation/mod_core-bed/bivalent_promoters/cat_all"
out_dir <- "/home/bettimj/gamazon_rotation/mod_core-bed/bivalent_promoters/cat_all"

tissues <- c("adipose", "blood", "bone", "brain", "cancer", "digestive", "endocrine", "endothelial", "epithelial", "esc_deriv", "esc", "eye", "heart", "hsc_and_b_cell", "ipsc", "kidney", "liver", "lung", "lymphoblastoid", "mesench", "muscle", "myosat", "neurosph", "other", "pancreas", "placenta_and_eem", "pns", "reproductive", "sm_muscle", "spleen", "stromal", "thymus", "urinary")

#Loop through the gene-annotated bivalent promoter BED for each tissue, outputting only the unique genes as new files
for (tissue in tissues) {
	file <- fread(paste(source_dir, paste0("genes_", tissue, "_bivalent_promoters.bed"), sep = "/"), header = FALSE, quote = FALSE, sep = "\t")
	df <- as.data.frame(file)
	df <- unique(df[,11])
	write.table(df, file = paste(out_dir, paste0("unique_genes_", tissue, "_bivalent_promoters.bed"), sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

