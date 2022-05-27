library("data.table")
library("gap")
library("qqman")

#Declare the paths of the split annotated summary statistics files
reg_anno_path <- "/home/bettimj/gamazon_rotation/core-bed_analysis/uk_biobank/hypertension/hg19/tss_5000_1000/reg_anno_phecode-401-both_sexes.func_anno.all_tissues.hg19.5000_1000.tsv"
noreg_anno_path <- "/home/bettimj/gamazon_rotation/core-bed_analysis/uk_biobank/hypertension/hg19/tss_5000_1000/noreg_anno_phecode-401-both_sexes.func_anno.all_tissues.hg19.5000_1000.tsv"
  
#Open the split annotated summary statistics as data frames
reg_anno_file <- fread(reg_anno_path, header = TRUE, quote = "", sep = "\t")
reg_anno_df <- as.data.frame(reg_anno_file)
reg_pvals_eur <- reg_anno_df$pval_EUR[!is.na(reg_anno_df$pval_EUR)]
reg_pvals_meta <- reg_anno_df$pval_meta[!is.na(reg_anno_df$pval_meta)]

noreg_anno_file <- fread(noreg_anno_path, header = TRUE, quote = "", sep = "\t")
noreg_anno_df <- as.data.frame(noreg_anno_file)
noreg_pvals_eur <- noreg_anno_df$pval_EUR[!is.na(noreg_anno_df$pval_EUR)]
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