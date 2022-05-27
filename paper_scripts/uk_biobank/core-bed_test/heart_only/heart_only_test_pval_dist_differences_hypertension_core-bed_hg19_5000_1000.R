library("data.table")
library("gap")
library("qqman")

#Declare the path of the annotated summary statistics files
sum_stats_path <- "/home/bettimj/gamazon_rotation/core-bed_analysis/uk_biobank/hypertension/hg19/tss_5000_1000/phecode-401-both_sexes.func_anno.heart.hg19.5000_1000.tsv"

#Open the summary statistics file as a data frame
sum_stats_file <- fread(sum_stats_path, header = TRUE, sep = "\t", quote = "")
df <- as.data.frame(sum_stats_file)

#Split the input into separate sets of variants with a regulatory annotation and those without
reg_anno_df <- df[(df[,47] == "active_promoter" | df[,47] == "bivalent_promoter" | df[,47] == "silenced_promoter" | df[,47] == "active_enhancer" | df[,47] == "poised_enhancer" | df[,47] == "primed_enhancer"),]
nrow(reg_anno_df)

noreg_anno_df <- df[(df[,47] != "active_promoter" & df[,47] != "bivalent_promoter" & df[,47] != "silenced_promoter" & df[,47] != "active_enhancer" & df[,47] != "poised_enhancer" & df[,47] != "primed_enhancer"),]
nrow(noreg_anno_df)
  
#Open the split annotated summary statistics as data frames
reg_pvals_meta <- reg_anno_df$pval_meta[!is.na(reg_anno_df$pval_meta)]

noreg_pvals_meta <- noreg_anno_df$pval_meta[!is.na(noreg_anno_df$pval_meta)]

print("Median p-value with annotation:")
median(reg_pvals_meta)
print("Median p-value without annotation:")
median(noreg_pvals_meta)

print("Mean p-value with annotation:")
mean(reg_pvals_meta)
print("Mean p-value without annotation:")
mean(noreg_pvals_meta)

wilcox.test(reg_pvals_meta, noreg_pvals_meta)