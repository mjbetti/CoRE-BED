library("data.table")

#Declare the path of the global CoRE-BED elements directory
source_dir <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/tss_5000_1000/bivalent_promoters"
out_dir <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/tss_5000_1000/bivalent_promoters"

tissues <- c("adipose", "adrenal_gland", "artery", "blood", "breast", "cultured_fibroblast", "ebv_transformed_lymphocyte", "es", "esophagus_muscularis_mucosa", "esophagus_squamous_epithelium", "heart", "intestine", "ips", "kidney", "liver", "lung", "neuron", "ovary", "pancreas", "prostate", "skeletal_muscle", "skin", "spleen", "stomach", "testis", "thyroid", "uterus", "vagina")

#Loop through the gene-annotated bivalent promoter BED for each tissue, outputting only the unique genes as new files
for (tissue in tissues) {
	file <- fread(paste(source_dir, paste0("genes_", tissue, "_bivalent_promoters_5000_1000.bed"), sep = "/"), header = FALSE, quote = FALSE, sep = "\t")
	df <- as.data.frame(file)
	df <- unique(df[,12])
	write.table(df, file = paste(out_dir, paste0("unique_genes_", tissue, "_bivalent_promoters_5000_1000.bed"), sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

