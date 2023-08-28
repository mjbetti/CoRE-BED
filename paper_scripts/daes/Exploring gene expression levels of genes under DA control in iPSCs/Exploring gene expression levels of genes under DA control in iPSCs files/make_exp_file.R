library("data.table")

exp_dir <- "/data/g_gamazon_lab/abehd/Genetically_Determined/iPSC_data/ipSc/RNAseq_data/unique_individual_data"

init_file <- fread("/data/g_gamazon_lab/abehd/Genetically_Determined/iPSC_data/ipSc/RNAseq_data/unique_individual_data/HPSI0114i-bezi_1.GRCh37.75.cdna.kallisto.transcripts.abundance.rnaseq.20150415.tsv", header = TRUE, sep = "\t", quote = "")
init_df <- as.data.frame(init_file)
init_df <- init_df[,1]
init_df <- as.data.frame(init_df)
names(init_df) <- "target_id"

for (file in list.files(exp_dir, pattern = "*.tsv")) {
	print(paste(file, "..."))
	sample_name <- strsplit(file, "_")
	sample_name <- unlist(lapply(sample_name, `[[`, 1))
	exp_file <- fread(paste(exp_dir, file, sep = "/"), header = TRUE, sep = "\t", quote = "")
	exp_df <- as.data.frame(exp_file)
	init_df[[sample_name]] <- exp_df$tpm
}

write.table(init_df, file = "/home/bettimj/gamazon_rotation/mod_core-bed/primed_promoter_ipsc/hipsci_exp_all_samples.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
