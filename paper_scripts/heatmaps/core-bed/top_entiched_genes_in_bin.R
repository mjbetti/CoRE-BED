library("data.table")

is.nan.data.frame <- function(x)
do.call(cbind, lapply(x, is.nan))

#Declare the paths of the binned gene and regulatory element counts, as well as the bed file listing all genes in each bin
binned_counts_dir <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/rcircos/gene_density"
binned_genes_path <- "/home/bettimj/reference_genomes/hg38.genome.1Mb.genes.bed"

genes_per_reg_out_dir <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/rcircos/gene_density/ind_tissues_genes_per_reg_element"
regs_per_gene_out_dir <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/rcircos/gene_density/ind_tissues_reg_elements_per_gene"

#Open the binned genes file as a data frame
binned_genes_file <- fread(binned_genes_path, header = FALSE, sep = "\t", quote = "")
binned_genes_df <- as.data.frame(binned_genes_file)

#Loop through each of these files by tissue, extracting the genes in the highest enriched bin on each chromosome
tissues <- c("adipose", "adrenal_gland", "artery", "blood", "breast", "cultured_fibroblast", "ebv_transformed_lymphocyte", "es", "esophagus_muscularis_mucosa", "esophagus_squamous_epithelium", "heart", "intestine", "ips", "kidney", "liver", "lung", "neuron", "ovary", "pancreas", "prostate", "skeletal_muscle", "skin", "spleen", "stomach", "testis", "thyroid", "uterus", "vagina")

##Genes per regulatory element
for (tissue in tissues) {
	file <- fread(paste(binned_counts_dir, paste0("binned_counts_", tissue, "_all_regulatory_elements_5000_1000.sorted.bed"), sep = "/"), header = FALSE, sep = "\t", quote = "")
	df <- as.data.frame(file)
	df[,6] <- (df[,4] / df[,5])
	df <- df[!is.na(df[,6]),]
	df <- df[!is.nan(df[,6]),]
	df <- df[(df[,6] != Inf),]
	max_bin <- df[(max(df[,6])),]
	print(paste0(tissue, ":"))
	print(max_bin)
	max_genes_df <- binned_genes_df[(binned_genes_df[,1] == max_bin[,1] & binned_genes_df[,2] == max_bin[,2] & binned_genes_df[,3] == max_bin[,3]),]
	max_genes_info <- strsplit(max_genes_df[,13], ";")
	gene_name <- unlist(lapply(max_genes_info, `[`, 4))
	gene_name <- substring(gene_name, 11)
	write.table(gene_name, file <- paste(genes_per_reg_out_dir, paste0(tissue, "_top_enriched_genes_per_reg_element.txt"), sep = "/"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

print("\n")
print("\n")

##Regulatory elements per gene
for (tissue in tissues) {
	file <- fread(paste(binned_counts_dir, paste0("binned_counts_", tissue, "_all_regulatory_elements_5000_1000.sorted.bed"), sep = "/"), header = FALSE, sep = "\t", quote = "")
	df <- as.data.frame(file)
	df[,6] <- (df[,5] / df[,4])
	df <- df[!is.na(df[,6]),]
	df <- df[!is.nan(df[,6]),]
	df <- df[(df[,6] != Inf),]
	max_bin <- df[(max(df[,6])),]
	print(paste0(tissue, ":"))
	print(max_bin)
	max_genes_df <- binned_genes_df[(binned_genes_df[,1] == max_bin[,1] & binned_genes_df[,2] == max_bin[,2] & binned_genes_df[,3] == max_bin[,3]),]
	max_genes_info <- strsplit(max_genes_df[,13], ";")
	gene_name <- unlist(lapply(max_genes_info, `[`, 4))
	gene_name <- substring(gene_name, 11)
	write.table(gene_name, file <- paste(regs_per_gene_out_dir, paste0(tissue, "_top_enriched_reg_elements_per_gene.txt"), sep = "/"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}