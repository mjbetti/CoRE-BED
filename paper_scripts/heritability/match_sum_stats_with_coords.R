library("data.table")

path_1kg <- "/home/bettimj/gamazon_rotation/core-bed_analysis/s2g_snps/1kg_scores/allsnps.txt.gz"
path_ukbb <- "/home/bettimj/gamazon_rotation/core-bed_analysis/s2g_snps/ukbb_scores/allsnps.txt.gz"
# path_dbsnp <- "/home/bettimj/gamazon_rotation/00-All.bed"
sum_stats_dir <- "/home/bettimj/gamazon_rotation/core-bed_analysis/s2g_snps/h2/sum_stats"
out_dir <- paste(sum_stats_dir, "reformatted", sep = "/")

#Open the files as data frames
file_1kg <- fread(path_1kg, header = FALSE, sep = " ", quote = "")
file_1kg <- as.data.frame(file_1kg)

file_ukbb <- fread(path_ukbb, header = FALSE, sep = " ", quote = "")
file_ukbb <- as.data.frame(file_ukbb)

cat_df <- rbind(file_1kg, file_ukbb)
cat_df <- unique(cat_df[order(cat_df$V2, cat_df$V3),])

# file_dbsnp <- fread(path_dbsnp, header = FALSE, quote = "", sep = "\t")
# file_dbsnp <- as.data.frame(file_dbsnp)

for (i in list.files(sum_stats_dir, pattern = "*.sumstats.gz")) {
	print(paste0("Processing ", i, "..."))
	gwas_path <- paste(sum_stats_dir, i, sep = "/")
	gwas_file <- fread(gwas_path, header = TRUE, sep = "\t", quote = "")
	gwas_df <- as.data.frame(gwas_file)
	merged_df <- merge(gwas_df, cat_df, by.x = "SNP", by.y = "V1")
	print(nrow(gwas_df))
	print(nrow(merged_df))
	merged_df <- data.frame(merged_df[,(ncol(merged_df) - 1):ncol(merged_df)], merged_df[,1:(ncol(merged_df) - 2)])
	names(merged_df)[1:2] <- c("CHR", "POS")
	write.table(merged_df, file = paste(out_dir, paste0(i, ".reformat.txt"), sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}