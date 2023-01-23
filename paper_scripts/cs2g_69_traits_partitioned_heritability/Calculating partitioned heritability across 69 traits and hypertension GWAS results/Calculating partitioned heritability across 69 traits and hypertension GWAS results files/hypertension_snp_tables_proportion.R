library("data.table")
library("dplyr")

gwas_array <- c("hypertension")

#Declare the path of the working directory
in_dir <- "/home/bettimj/gamazon_rotation/mod_core-bed/uk_biobank/h2/anno_1000G/hypertension"

cat_df <- data.frame()

for (gwas in gwas_array) {
	print(gwas)
	results_file <- fread(paste0(gwas, "_h2.results"), header = TRUE, sep = "\t", quote = "")
	results_df <- as.data.frame(results_file)
	results_df <- data.frame(results_df[,3])
	results_df <- transpose(results_df)
	cat_df <- rbind(cat_df, results_df)
}

cat_df <- data.frame(gwas_array, cat_df)
names(cat_df) <- c("GWAS", "Class1", "Class2", "Class3", "Class4", "Class5", "Class6", "Class7", "Class8", "any_corebed")

write.table(cat_df, file = "all_core_bed_h2_proportions.csv", sep = ",", quote = FALSE, row.names = FALSE, col.names = TRUE)