library("data.table")

#Disable scientific notation
#options(scipen = 999)

#Declare the path of the set of annotated case-control variants
var_path <- "/home/bettimj/gamazon_rotation/gwas_catalog/all_associations_v1.0.reformat.case_control.txt"

#Open the file as a data frame
var_file <- fread(var_path, header = TRUE, quote = "", sep = "\t")
var_df <- as.data.frame(var_file)

#First remove all variants with either a p-value, MAF, or odds ratio of NA
var_df <- var_df[!is.na(var_df$P_VALUE),]
var_df <- var_df[!is.na(var_df$RISK_ALLELE_FREQUENCY),]
var_df <- var_df[!is.na(var_df$OR_or_BETA),]

#First filter variants with a RISK_ALLELE_FREQUENCY greater than 1
var_df <- var_df[!(as.numeric(var_df$RISK_ALLELE_FREQUENCY > 1)),]

#For all remaining variants with a RISK_ALLELE_FREQUENCY greater than 0.5, recalculate MAF
mafs_needmod <- var_df[(as.numeric(var_df$RISK_ALLELE_FREQUENCY) > 0.5),]
mafs_nomod <- var_df[(as.numeric(var_df$RISK_ALLELE_FREQUENCY) <= 0.5),]

mafs_needmod["MAF"] <- (1 - as.numeric(mafs_needmod$RISK_ALLELE_FREQUENCY))
mafs_nomod["MAF"] <- as.numeric(mafs_nomod$RISK_ALLELE_FREQUENCY)

var_df <- rbind(mafs_needmod, mafs_nomod)
var_df <- var_df[!is.na(var_df$MAF),]
var_df <- var_df[order(var_df$CHR_ID, var_df$CHR_POS),]
#var_df <- var_df[(var_df$MAF >= 0.01),]

var_df <- var_df[order(var_df$OR_or_BETA, decreasing = TRUE),]

#Write out the resulting data frame as a new file
write.table(var_df, file = "all_associations_v1.0.reformat.case_control.recal.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)