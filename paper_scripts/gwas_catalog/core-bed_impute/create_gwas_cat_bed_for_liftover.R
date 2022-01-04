library("data.table")

#Declare the path of the pre-filtered GWAS Catalog set based on hg38, containing only case-control variants that are not missing GWAS effect size, p-value, or MAF
hg38_path <- "/home/bettimj/gamazon_rotation/gwas_catalog/all_associations_v1.0.reformat.case_control.recal.txt"

#Open the file as a data frame
hg38_file <- fread(hg38_path, header = TRUE, sep = "\t", quote = "")
hg38_df <- as.data.frame(hg38_file)

#Create a UCSC BED file containing the SNP coordinates for each variant, along with the corresponding rsID
bed_df <- data.frame(paste0("chr", hg38_df$CHR_ID), hg38_df$CHR_POS, (hg38_df$CHR_POS + 1), hg38_df$SNPS)

#Write out as a new UCSC BED file
write.table(bed_df, file = "/home/bettimj/gamazon_rotation/gwas_catalog/all_associations_v1.0.reformat.case_control.recal.hg38.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)