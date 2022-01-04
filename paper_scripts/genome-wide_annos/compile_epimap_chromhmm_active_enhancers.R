library("data.table")

#Declare the path of the source directory
source_dir <- "/home/bettimj/gamazon_rotation/epimap_chromhmm/hg38/personal.broadinstitute.org/cboix/epimap/ChromHMM/observed_aux_18_hg38/CALLS"

cat_df <- data.frame()

#Loop through each of the original files, retaining only the coordinates with one of the two classes of active enhancer
for (i in list.files(source_dir, pattern = "BSS*")) {
	print(i)
	file <- fread(paste(source_dir, i, sep = "/"), header = FALSE, sep = "\t", quote = FALSE)
	df <- as.data.frame(file)
	df <- df[(df[,4] == "EnhA1" | df[,4] == "EnhA2"),]
	df <- df[,1:3]
	cat_df <- rbind(cat_df, df)
}

cat_df <- unique(cat_df)

write.table(cat_df, file = "all_epimap_active_enh_18_CALLS_segments_hg38.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

