library("pheatmap")

#Open the genes per regulatory element and regulatory elements per gene files as data frames
genes_per_reg_path <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/rcircos/gene_density/genes_per_reg_element/median_genes_per_reg_element_all_tissues_5000_1000.txt"
regs_per_gene_path <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/rcircos/gene_density/genes_per_reg_element/reg_elements_per_gene/median_reg_elements_per_gene_all_tissues_5000_1000.txt"

##genes per regulatory element
genes_per_reg_file <- read.table(genes_per_reg_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
genes_per_reg_df <- as.data.frame(genes_per_reg_file)
genes_per_reg_df <- as.matrix(genes_per_reg_df[,-1])

pdf("pheatmap_median_genes_per_reg_element_all_tissues_5000_1000.pdf")
pheatmap(genes_per_reg_df, display_numbers = T, color = colorRampPalette(c('white','red'))(100), cluster_rows = T, cluster_cols = T, fontsize_number = 5, labels_row = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"), labels_col = c("Adipose", "Adrenal gland", "Artery", "Blood", "Breast", "Cultured fibroblast", "EBV-transformed lymphocyte", "ES", "Esophagus muscularis mucosa", "Esophagus squamous epithelium", "Heart", "Intestine", "iPS", "Kidney", "Liver", "Lung", "Neuron", "Ovary", "Pancreas", "Prostate", "Skeletal muscle", "Skin", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina"), main = "Median genes per regulatory element per Mb", breaks = seq(0, 8, length.out=101))
dev.off()

##regulatory elements per gene
regs_per_gene_file <- read.table(regs_per_gene_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
regs_per_gene_df <- as.data.frame(regs_per_gene_file)
regs_per_gene_df <- as.matrix(regs_per_gene_df[,-1])

pdf("pheatmap_median_reg_elements_per_gene_all_tissues_5000_1000.pdf")
pheatmap(regs_per_gene_df, display_numbers = T, color = colorRampPalette(c('white','red'))(100), cluster_rows = T, cluster_cols = T, fontsize_number = 5, labels_row = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"), labels_col = c("Adipose", "Adrenal gland", "Artery", "Blood", "Breast", "Cultured fibroblast", "EBV-transformed lymphocyte", "ES", "Esophagus muscularis mucosa", "Esophagus squamous epithelium", "Heart", "Intestine", "iPS", "Kidney", "Liver", "Lung", "Neuron", "Ovary", "Pancreas", "Prostate", "Skeletal muscle", "Skin", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina"), main = "Median regulatory elements per gene per Mb", breaks = seq(0, 8, length.out=101))
dev.off()