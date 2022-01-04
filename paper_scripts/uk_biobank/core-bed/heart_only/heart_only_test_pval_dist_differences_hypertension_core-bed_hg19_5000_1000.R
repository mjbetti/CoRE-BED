library("data.table")

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