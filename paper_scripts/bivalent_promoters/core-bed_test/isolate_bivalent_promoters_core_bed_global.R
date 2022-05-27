library("data.table")

#Declare the path of the global CoRE-BED elements directory
source_dir <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/tss_5000_1000"
out_dir <- paste(source_dir, "bivalent_promoters", sep = "/")

tissues <- c("adipose", "adrenal_gland", "artery", "blood", "breast", "cultured_fibroblast", "ebv_transformed_lymphocyte", "es", "esophagus_muscularis_mucosa", "esophagus_squamous_epithelium", "heart", "intestine", "ips", "kidney", "liver", "lung", "neuron", "ovary", "pancreas", "prostate", "skeletal_muscle", "skin", "spleen", "stomach", "testis", "thyroid", "uterus", "vagina")

#Loop through the promoter BED for each tissue, isolating the bivalent promoters and outputting them to the designated output directory
for (tissue in tissues) {
	file <- fread(paste(source_dir, paste0(tissue, "_all_promoters_5000_1000.bed"), sep = "/"), header = FALSE, quote = FALSE, sep = "\t")
	df <- as.data.frame(file)
	df <- df[(df[,4] == "bivalent_promoter"),]
	write.table(df, file = paste(out_dir, paste0(tissue, "_bivalent_promoters_5000_1000.bed"), sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}