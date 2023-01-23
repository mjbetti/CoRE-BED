library("data.table")

#Declare the path of the directory containing all of the CoRE-BED annotations across the 816 biosamples
anno_dir <- "/home/bettimj/gamazon_rotation/mod_core-bed/all_anno"

tissue_names <- c()

for (file in list.files(anno_dir, pattern = "*_hg19.genome.1kb.bed")) {
	tissue_names <- c(tissue_names, file)
}

#Open the first file as a data frame
init_file <- fread(paste(anno_dir, tissue_names[1], sep = "/"), header = FALSE, sep = "\t", quote = "")
init_df <- as.data.frame(init_file)
init_df$is_functional <- 0
init_df$is_functional[init_df[,4] == 2 | init_df[,4] == 3 | init_df[,4] == 4 | init_df[,4] == 5 | init_df[,4] == 6 | init_df[,4] == 7 | init_df[,4] == 8] <- 1

for (tissue in tissue_names[2:length(tissue_names)]) {
	print(paste0(tissue, "..."))
	new_file <- fread(paste(anno_dir, tissue, sep = "/"), header = FALSE, sep = "\t", quote = "")
	new_df <- as.data.frame(new_file)
	init_df[,4] <- new_df[,4]
	init_df$is_functional[init_df[,4] == 2 | init_df[,4] == 3 | init_df[,4] == 4 | init_df[,4] == 5 | init_df[,4] == 6 | init_df[,4] == 7 | init_df[,4] == 8] <- 1
}

#Write out the results as a new file
init_df_out <- data.frame(init_df[,1:3], init_df[,5])

write.table(init_df_out, file = "/home/bettimj/gamazon_rotation/mod_core-bed/all_anno/all_tissues_hg19_anno.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)