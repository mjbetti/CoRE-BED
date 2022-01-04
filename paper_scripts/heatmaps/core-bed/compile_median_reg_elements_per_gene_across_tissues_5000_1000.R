source_dir <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/rcircos/gene_density/genes_per_reg_element"
out_dir <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/rcircos/gene_density/genes_per_reg_element/reg_elements_per_gene"

tissues <- c("adipose", "adrenal_gland", "artery", "blood", "breast", "cultured_fibroblast", "ebv_transformed_lymphocyte", "es", "esophagus_muscularis_mucosa", "esophagus_squamous_epithelium", "heart", "intestine", "ips", "kidney", "liver", "lung", "neuron", "ovary", "pancreas", "prostate", "skeletal_muscle", "skin", "spleen", "stomach", "testis", "thyroid", "uterus", "vagina")

chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
new_df <- data.frame(chrs)

for (tissue in tissues) {
	file_name <- paste(source_dir, paste0(tissue, "_genes_reg_elements_per_mb.txt"), sep = "/")
	file <- read.table(file_name, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
	df <- as.data.frame(file)
	median_reg_elements <- df$median_reg_elements_per_gene
	new_df <- data.frame(new_df, median_reg_elements)
}

names(new_df) <- c("chromosome", "adipose", "adrenal_gland", "artery", "blood", "breast", "cultured_fibroblast", "ebv_transformed_lymphocyte", "es", "esophagus_muscularis_mucosa", "esophagus_squamous_epithelium", "heart", "intestine", "ips", "kidney", "liver", "lung", "neuron", "ovary", "pancreas", "prostate", "skeletal_muscle", "skin", "spleen", "stomach", "testis", "thyroid", "uterus", "vagina")

write.table(new_df, file = paste(out_dir, "median_reg_elements_per_gene_all_tissues_5000_1000.txt", sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)