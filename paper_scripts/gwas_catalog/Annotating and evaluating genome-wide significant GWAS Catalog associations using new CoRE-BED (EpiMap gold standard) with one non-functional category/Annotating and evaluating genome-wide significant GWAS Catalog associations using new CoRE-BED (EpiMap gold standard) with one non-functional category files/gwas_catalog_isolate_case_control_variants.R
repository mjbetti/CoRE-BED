library("data.table")

#Declare the path of the (already-reformatted) GWAS catalog variant set
in_path <- "/home/bettimj/gamazon_rotation/gwas_catalog/all_associations_v1.0.reformat.txt"
out_dir <- "/home/bettimj/gamazon_rotation/gwas_catalog"

#Open the file as a data frame
file <- fread(in_path, header = TRUE, sep = "\t", quote = "")
df <- as.data.frame(file)

df_case_control <- df[(grepl("case",df$INITIAL_SAMPLE_SIZE) | grepl("control",df$INITIAL_SAMPLE_SIZE)),]

#Write out the case-control data frame as a new file
write.table(df_case_control, file = paste(out_dir, "all_associations_v1.0.reformat.case_control.txt", sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)