library("data.table")
library("stringr")
library("dplyr")

#Declare the path of the input file (CoRE-BED annotated GWAS summary statistics)
in_path <- "/home/bettimj/gamazon_rotation/core-bed_analysis/uk_biobank/hypertension/hg19/tss_5000_1000/phecode-401-both_sexes.func_anno.all_tissues.hg19.5000_1000.tsv"

file <- fread(in_path, header = TRUE, sep = "\t", quote = "")
df <- as.data.frame(file)

#Split the input into separate sets of variants with a regulatory annotation and those without
reg_df <- df[((rowSums(df[,(ncol(df) - 27):ncol(df)] == "active_promoter") > 0) | (rowSums(df[,(ncol(df) - 27):ncol(df)] == "bivalent_promoter") > 0) | (rowSums(df[,(ncol(df) - 27):ncol(df)] == "silenced_promoter") > 0) | (rowSums(df[,(ncol(df) - 27):ncol(df)] == "active_enhancer") > 0) | (rowSums(df[,(ncol(df) - 27):ncol(df)] == "poised_enhancer") > 0) | (rowSums(df[,(ncol(df) - 27):ncol(df)] == "primed_enhancer") > 0)),]
nrow(reg_df)

noreg_df <- df[((rowSums(df[,(ncol(df) - 27):ncol(df)] == "active_promoter") == 0) & (rowSums(df[,(ncol(df) - 27):ncol(df)] == "bivalent_promoter") == 0) & (rowSums(df[,(ncol(df) - 27):ncol(df)] == "silenced_promoter") == 0) & (rowSums(df[,(ncol(df) - 27):ncol(df)] == "active_enhancer") == 0) & (rowSums(df[,(ncol(df) - 27):ncol(df)] == "poised_enhancer") == 0) & (rowSums(df[,(ncol(df) - 27):ncol(df)] == "primed_enhancer") == 0)),]
nrow(noreg_df)

#Write the split data frames out to two new files
write.table(reg_df, file = "phecode-401-both_sexes.func_anno.all_tissues.hg19.5000_1000.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

write.table(noreg_df, file = "phecode-401-both_sexes.func_anno.all_tissues.hg19.5000_1000.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)