library("data.table")

#Declare the paths of the CoRE-BED-training results files across the five decision tree structures reported
regs_1k5k <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/alternative_trees/epimap/exp_only/1k5k/all_regulatory_elements_epimap_1000_5000.txt"
regs_2k2k <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/alternative_trees/epimap/exp_only/2k2k/all_regulatory_elements_epimap_2000_2000.txt"
regs_5k1k <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/alternative_trees/epimap/exp_only/5k1k/all_regulatory_elements_epimap_5000_1000.txt"
regs_notss_enh <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/alternative_trees/epimap/exp_only/notss_enh_first/all_regulatory_elements_epimap_5000_1000_no_tss_enh_first.txt"
regs_notss_prom <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/alternative_trees/epimap/exp_only/notss_prom_first/all_regulatory_elements_epimap_5000_1000_no_tss_prom_first.txt"

paths <- c(regs_1k5k, regs_2k2k, regs_5k1k, regs_notss_enh, regs_notss_prom)

#Open each of these files as a data frame
for (path in paths) {
	print(path)
	file <- fread(path, header = FALSE, sep = "\t", quote = "")
	df <- as.data.frame(file)
	total <- nrow(df)
	active_prom <- nrow(df[(df[,4] == "active_promoter"),])
	bivalent_prom <- nrow(df[(df[,4] == "bivalent_promoter"),])
	silenced_prom <- nrow(df[(df[,4] == "silenced_promoter"),])
	active_enh <- nrow(df[(df[,4] == "active_enhancer"),])
	poised_enh <- nrow(df[(df[,4] == "poised_enhancer"),])
	primed_enh <- nrow(df[(df[,4] == "primed_enhancer"),])
	shannon_ent <- -(((active_prom/total)*log2(active_prom/total)) + ((bivalent_prom/total)*log2(bivalent_prom/total)) + ((silenced_prom/total)*log2(silenced_prom/total)) + ((active_enh/total)*log2(active_enh/total)) + ((poised_enh/total)*log2(poised_enh/total)) + ((primed_enh/total)*log2(primed_enh/total)))
	print(shannon_ent)
}