is.nan.data.frame <- function(x)
do.call(cbind, lapply(x, is.nan))

#Declare the path of the directory containing the UCSC BED file counts for all genes and regulatory elements across the genome
source_dir <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/rcircos/gene_density"
out_dir <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/rcircos/gene_density/genes_per_reg_element"

#Loop through each file, converting it to UCSC BED format
tissues <- c("adipose", "adrenal_gland", "artery", "blood", "breast", "cultured_fibroblast", "ebv_transformed_lymphocyte", "es", "esophagus_muscularis_mucosa", "esophagus_squamous_epithelium", "heart", "intestine", "ips", "kidney", "liver", "lung", "neuron", "ovary", "pancreas", "prostate", "skeletal_muscle", "skin", "spleen", "stomach", "testis", "thyroid", "uterus", "vagina")
chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")

for (tissue in tissues) {
	file_name <- paste(source_dir, paste0("binned_counts_", tissue, "_all_regulatory_elements_5000_1000.sorted.bed"), sep = "/")
	file <- read.table(file_name, sep = "\t", stringsAsFactors = FALSE)
	df <- as.data.frame(file)
	stat_df <- data.frame()
	for (chr in chrs) {
		chr_df <- df[(df[,1] == chr),]
		chr_df[,6] <- chr_df[,4] / chr_df[,5]
		chr_df <- chr_df[(chr_df[,6] != Inf),]
		chr_df <- chr_df[!is.na(chr_df[,6]),]
		chr_df <- chr_df[!is.nan(chr_df[,6]),]
		
		chr_df[,7] <- chr_df[,5] / chr_df[,4]
		chr_df <- chr_df[(chr_df[,7] != Inf),]
		chr_df <- chr_df[!is.na(chr_df[,7]),]
		chr_df <- chr_df[!is.nan(chr_df[,7]),]
		
		stat_df_row <- c(chr, mean(chr_df[,6]), median(chr_df[,6]), mean(chr_df[,7]), median(chr_df[,7]))
		stat_df <- rbind(stat_df, stat_df_row)
	}
	names(stat_df) <- c("chromosome", "mean_genes_per_reg_element", "median_genes_per_reg_element", "mean_reg_elements_per_gene", "median_reg_elements_per_gene")
	write.table(stat_df, file = paste(out_dir, paste0(tissue, "_genes_reg_elements_per_mb.txt"), sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}
