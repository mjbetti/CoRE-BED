library("data.table")

#Declare the path of the global CoRE-BED elements directory
source_dir <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/tss_5000_1000/epimap/bed_format/concat"
out_dir <- paste(source_dir, "bivalent_promoters", sep = "/")

tissues <- c("adipose", "blood", "bone", "brain", "cancer", "digestive", "endocrine", "endothelial", "epithelial", "esc_deriv", "esc", "eye", "heart", "hsc_and_b_cell", "ipsc", "kidney", "liver", "lung", "lymphoblastoid", "mesench", "muscle", "myosat", "neurosph", "other", "pancreas", "placenta_and_eem", "pns", "reproductive", "sm_muscle", "spleen", "stromal", "thymus", "urinary")

#Loop through the promoter BED for each tissue, isolating the bivalent promoters and outputting them to the designated output directory
for (tissue in tissues) {
	file <- fread(paste(source_dir, paste0(tissue, "_all_regulatory_elements_epimap_5000_1000.sorted.bed"), sep = "/"), header = FALSE, quote = FALSE, sep = "\t")
	df <- as.data.frame(file)
	df <- df[(df[,4] == "bivalent_promoter"),]
	write.table(df, file = paste(out_dir, paste0(tissue, "_bivalent_promoters_all_regulatory_elements_epimap_5000_1000.sorted.bed"), sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}