library("data.table")

#Declare the paths of the CoRE-BED-training results files across the five decision tree structures reported
regs_1k5k <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/alternative_trees/epimap/1k5k/all_regulatory_elements_epimap_1000_5000.txt"
regs_2k2k <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/alternative_trees/epimap/2k2k/all_regulatory_elements_epimap_2000_2000.txt"
regs_5k1k <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/alternative_trees/epimap/5k1k/all_regulatory_elements_epimap_5000_1000.txt"
regs_notss_enh <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/alternative_trees/epimap/notss_enh_5k1k/all_regulatory_elements_epimap_5000_1000_no_tss_enh_first.txt"
regs_notss_prom <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/alternative_trees/epimap/notss_prom_5k1k/all_regulatory_elements_epimap_5000_1000_no_tss_prom_first.txt"

paths <- c(regs_1k5k, regs_2k2k, regs_5k1k, regs_notss_enh, regs_notss_prom)

out_dir <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/alternative_trees/epimap/jaccard_prep"

#Open each of these files as a data frame
for (path in paths) {
	sample_name <- strsplit(path, "/")
	sample_name <- unlist(sapply(sample_name, "[[", 10))
	file <- fread(path, header = FALSE, sep = "\t", quote = "")
	df <- as.data.frame(file)
	df <- data.frame(paste0("chr", df[,1]), df[,2:3])
	write.table(df, file = paste(out_dir, paste0(sample_name, ".all_annos.bed"), sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}