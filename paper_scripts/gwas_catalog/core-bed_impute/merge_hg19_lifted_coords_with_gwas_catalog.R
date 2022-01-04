library("data.table")

#Declare the path of the pre-filtered GWAS Catalog set based on hg38, containing only case-control variants that are not missing GWAS effect size, p-value, or MAF; as well as the UCSC BED file of coordinates lifted over to hg19
hg38_path <- "/home/bettimj/gamazon_rotation/gwas_catalog/all_associations_v1.0.reformat.case_control.recal.txt"
hg19_bed_path <- "/home/bettimj/gamazon_rotation/gwas_catalog/all_associations_v1.0.reformat.case_control.recal.lifted.hg19.bed"

#Open both of the files as data frames
hg38_file <- fread(hg38_path, header = TRUE, sep = "\t", quote = FALSE)
hg38_df <- as.data.frame(hg38_file)

hg19_bed_file <- fread(hg19_bed_path, header = FALSE, sep = "\t", quote = FALSE)
hg19_bed_df <- as.data.frame(hg19_bed_file)

#Merge the two files based on common rsID
merged_df <- merge(hg38_df, hg19_bed_df, by.x = "SNPS", by.y = "V4")
merged_df$CHR_ID <- merged_df$V1
merged_df$CHR_ID <- sub("chr", "", merged_df$CHR_ID)
merged_df$CHR_POS <- merged_df$V2
merged_df <- merged_df[,1:35]

#Write the new lifted set of GWAS Catalog variants out to a .txt file
write.table(merged_df, file = "/home/bettimj/gamazon_rotation/gwas_catalog/all_associations_v1.0.reformat.case_control.recal.lifted.hg19.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
