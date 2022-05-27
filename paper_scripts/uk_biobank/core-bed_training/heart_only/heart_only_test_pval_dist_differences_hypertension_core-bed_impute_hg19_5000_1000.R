library("data.table")

#Declare the paths of the original summary statistics file and binary annotation status
sum_stats_path <- "/home/bettimj/gamazon_rotation/uk_biobank/hypertension/phecode-401-both_sexes.tsv"
anno_path <- "/home/bettimj/gamazon_rotation/core-bed_analysis/core_bed_impute/uk_biobank/tss_5000_1000/sliced_annos/heart_only_encoded_anno_status_all.epimap_phecode-401-both_sexes.func_anno.5000_1000.txt"
  
#Open the split annotated summary statistics as a data frame
sum_stats_file <- fread(sum_stats_path, header = TRUE, sep = "\t", quote = "")
sum_stats_df <- as.data.frame(sum_stats_file)

anno_file <- fread(anno_path, header = FALSE, sep = "\t", quote = "")
anno_df <- as.data.frame(anno_file)

#Add the annotations to the sumstats
sum_stats_df <- data.frame(sum_stats_df, anno_df)
names(sum_stats_df)[length(names(sum_stats_df))] <- "CORE_BED_ANNO"

#Split the summary stats based on whether or not each coordinate has a CoRE-BED Impute annotation in any of the 827 CoRE-BED Impute cell/tissue types
reg_sum_stats_df <- sum_stats_df[(sum_stats_df$CORE_BED_ANNO == 1),]
noreg_sum_stats_df <- sum_stats_df[(sum_stats_df$CORE_BED_ANNO == 0),]

nrow(reg_sum_stats_df)
nrow(noreg_sum_stats_df)

reg_pvals_eur <- reg_anno_df$pval_EUR[!is.na(reg_anno_df$pval_EUR)]
reg_pvals_meta <- reg_anno_df$pval_meta[!is.na(reg_anno_df$pval_meta)]
reg_betas_eur <- reg_anno_df$beta_EUR[!is.na(reg_anno_df$beta_EUR)]
reg_betas_meta <- reg_anno_df$beta_meta[!is.na(reg_anno_df$beta_meta)]

noreg_pvals_eur <- noreg_anno_df$pval_EUR[!is.na(noreg_anno_df$pval_EUR)]
noreg_pvals_meta <- noreg_anno_df$pval_meta[!is.na(noreg_anno_df$pval_meta)]
noreg_betas_eur <- noreg_anno_df$beta_EUR[!is.na(noreg_anno_df$beta_EUR)]
noreg_betas_meta <- noreg_anno_df$beta_meta[!is.na(noreg_anno_df$beta_meta)]

#Comparing p-values
print("Comparing p-values of variants with and without a functional annotation...")

print("Median p-value of variants with a promoter/enhancer annotation:")
print(median(reg_pvals_eur))

print("Mean p-value of variants with a promoter/enhancer annotation:")
print(mean(reg_pvals_eur))

print("Range of promoter/enhancer variant p-values:")
print(range(reg_pvals_eur))

print("Median p-value of variants without a functional annotation:")
print(median(noreg_pvals_eur))

print("Mean p-value of variants without a functional annotation:")
print(mean(noreg_pvals_eur))

print("Range of p-values of variants without a functional annotation:")
print(range(noreg_pvals_eur))

wilcox.test(reg_pvals_eur, noreg_pvals_eur)

#Comparing betas
print("Comparing betas of variants with and without a functional annotation...")

print("Median beta of variants with a promoter/enhancer annotation:")
print(median(reg_betas_eur))

print("Mean beta of variants with a promoter/enhancer annotation:")
print(mean(reg_betas_eur))

print("Range of promoter/enhancer variant betas:")
print(range(reg_betas_eur))

print("Median beta of variants without a functional annotation:")
print(median(noreg_betas_eur))

print("Mean beta of variants without a functional annotation:")
print(mean(noreg_betas_eur))

print("Range of betas of variants without a functional annotation:")
print(range(noreg_betas_eur))

wilcox.test(reg_betas_eur, noreg_betas_eur)

##Meta analysis
#Comparing p-values
print("Comparing p-values of meta-analysis variants with and without a functional annotation...")

print("Median p-value of variants with a promoter/enhancer annotation:")
print(median(reg_pvals_meta))

print("Mean p-value of variants with a promoter/enhancer annotation:")
print(mean(reg_pvals_meta))

print("Range of promoter/enhancer variant p-values:")
print(range(reg_pvals_meta))

print("Median p-value of variants without a functional annotation:")
print(median(noreg_pvals_meta))

print("Mean p-value of variants without a functional annotation:")
print(mean(noreg_pvals_meta))

print("Range of p-values of variants without a functional annotation:")
print(range(noreg_pvals_meta))

wilcox.test(reg_pvals_meta, noreg_pvals_meta)

#Comparing betas
print("Comparing betas of variants with and without a functional annotation...")

print("Median beta of variants with a promoter/enhancer annotation:")
print(median(reg_betas_meta))

print("Mean beta of variants with a promoter/enhancer annotation:")
print(mean(reg_betas_meta))

print("Range of promoter/enhancer variant betas:")
print(range(reg_betas_meta))

print("Median beta of variants without a functional annotation:")
print(median(noreg_betas_meta))


print("Range of betas of variants without a functional annotation:")
print(range(noreg_betas_meta))

wilcox.test(reg_betas_meta, noreg_betas_meta)