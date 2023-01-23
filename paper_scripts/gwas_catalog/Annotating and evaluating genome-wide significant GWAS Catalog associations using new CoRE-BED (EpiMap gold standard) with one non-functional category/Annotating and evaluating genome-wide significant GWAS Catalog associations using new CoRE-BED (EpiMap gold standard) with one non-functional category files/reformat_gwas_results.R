library("data.table")
library("stringr")

#Declare the path of the input file
in_path <- "/home/bettimj/gamazon_rotation/gwas_catalog/all_associations_v1.0.txt"

#Open the file as a data frame
in_file <- fread("/home/bettimj/gamazon_rotation/gwas_catalog/all_associations_v1.0.txt", header = TRUE, sep = "\t", quote = "")
in_df <- as.data.frame(in_file)

#Remove white space in column entries
for (i in seq(1:34)) {
	in_df[,i] = gsub(" ", "_", in_df[,i], fixed = TRUE)
}

#And in header (also replacing "-" and "/" with "_")
for (i in seq(1:34)) {
	names(in_df)[i] = gsub(" ", "_", names(in_df)[i], fixed = TRUE)
}

for (i in seq(1:34)) {
	names(in_df)[i] = gsub("-", "_", names(in_df)[i], fixed = TRUE)
}

for (i in seq(1:34)) {
	names(in_df)[i] = gsub("/", "_", names(in_df)[i], fixed = TRUE)
}

#Split the coordinate column and take only the first value
in_df$CHR_POS <- strsplit(in_df$CHR_POS, ";")
in_df$CHR_POS <- unlist(lapply(in_df$CHR_POS, `[`, 1))

in_df$CHR_POS <- strsplit(in_df$CHR_POS, "_")
in_df$CHR_POS <- unlist(lapply(in_df$CHR_POS, `[`, 1))

#Loop through each column and replace any spaces with underscores
for (i in seq(1:(ncol(in_df)))) {
	in_df[,i] <- gsub(" ", "_", in_df[,i], fixed = TRUE)
}

in_df[in_df == ""] <- NA

out_df <- in_df[!is.na(in_df$CHR_POS),]

#Output the data frame back out as a .txt file
write.table(out_df, file = "all_associations_v1.0.reformat.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)