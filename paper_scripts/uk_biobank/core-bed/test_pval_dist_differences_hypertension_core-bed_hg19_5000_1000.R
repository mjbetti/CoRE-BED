library("data.table")

#Declare the paths of the split annotated summary statistics files
reg_anno_path <- "/home/bettimj/gamazon_rotation/core-bed_analysis/uk_biobank/hypertension/hg19/tss_5000_1000/reg_anno_phecode-401-both_sexes.func_anno.all_tissues.hg19.5000_1000.tsv"
noreg_anno_path <- "/home/bettimj/gamazon_rotation/core-bed_analysis/uk_biobank/hypertension/hg19/tss_5000_1000/noreg_anno_phecode-401-both_sexes.func_anno.all_tissues.hg19.5000_1000.tsv"
  
#Open the split annotated summary statistics as data frames
reg_anno_file <- fread(reg_anno_path, header = TRUE, quote = "", sep = "\t")
reg_anno_df <- as.data.frame(reg_anno_file)
reg_pvals_eur <- reg_anno_df$pval_EUR[!is.na(reg_anno_df$pval_EUR)]
reg_pvals_meta <- reg_anno_df$pval_meta[!is.na(reg_anno_df$pval_meta)]
reg_betas_eur <- reg_anno_df$beta_EUR[!is.na(reg_anno_df$beta_EUR)]
reg_betas_meta <- reg_anno_df$beta_meta[!is.na(reg_anno_df$beta_meta)]

noreg_anno_file <- fread(noreg_anno_path, header = TRUE, quote = "", sep = "\t")
noreg_anno_df <- as.data.frame(noreg_anno_file)
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