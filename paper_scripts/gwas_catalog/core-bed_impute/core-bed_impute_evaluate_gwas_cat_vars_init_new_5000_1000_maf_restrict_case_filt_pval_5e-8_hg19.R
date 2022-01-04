library("data.table")
library("stringr")
library("dplyr")

file <- fread("/home/bettimj/gamazon_rotation/core-bed_analysis/core_bed_impute/gwas_catalog/hg19/tss_5000_1000/all_epimap_tissues_all_associations_v1.0.reformat.case_control.recal.lifted.hg19.func_anno_5000_1000.txt", header = TRUE, sep = "\t", quote = "")
df <- as.data.frame(file)

#Retain only variants with a genome-wide significant p-value of 5e-8
df <- df[(df$PVALUE_MLOG >= 7.301029995663981),]

#Retain only variants with a MAF >= 0.01
df <- df[(df$MAF >= 0.01),]

#Filter out variants from a study with < 5,000 cases
sample_size <- df$INITIAL_SAMPLE_SIZE
sample_size <- strsplit(sample_size, "_cases,")
cases <- unlist(sapply(sample_size, "[[", 1))
cases <- strsplit(cases, "_")
cases <- unlist(sapply(cases, "[[", 1))
cases <- str_replace(cases, ",", "")
cases <- as.numeric(cases)
df$CASES <- cases
df <- df[!is.na(df$CASES),]
df <- df[df$CASES >= 5000,]

print("Splitting file into those with and without promoter/enhancer annotation...")
reg_df <- df[((rowSums(df[,36:862] == "active_promoter") > 0) | (rowSums(df[,36:862] == "bivalent_promoter") > 0) | (rowSums(df[,36:862] == "silenced_promoter") > 0) | (rowSums(df[,36:862] == "active_enhancer") > 0) | (rowSums(df[,36:862] == "poised_enhancer") > 0) | (rowSums(df[,36:862] == "primed_enhancer") > 0)),]
print("Number of variants with a functional (promoter/enhancer) annotation:")
print(nrow(reg_df))

#Comparing effect sizes
noreg_df <- df[!((rowSums(df[,36:862] == "active_promoter") > 0) | (rowSums(df[,36:862] == "bivalent_promoter") > 0) | (rowSums(df[,36:862] == "silenced_promoter") > 0) | (rowSums(df[,36:862] == "active_enhancer") > 0) | (rowSums(df[,36:862] == "poised_enhancer") > 0) | (rowSums(df[,36:862] == "primed_enhancer") > 0)),]
print("Number of variants without a functional (promoter/enhancer) annotation:")
print(nrow(noreg_df))

#Comparing effect sizes
print("Comparing effect sizes of variants with and without a functional annotation...")
effect_anno <- reg_df$OR_or_BETA

print("Median odds ratio of promoter/enhancer variants:")
print(median(effect_anno))

print("Range of promoter/enhancer variant odds ratios:")
print(range(effect_anno))

effect_noanno <- noreg_df$OR_or_BETA

print("Median odds ratio of variants without a functional annotation:")
print(median(effect_noanno))

print("Range of odds ratios of variants without a functional annotation:")
print(range(effect_noanno))

wilcox.test(effect_anno, effect_noanno)

#Comparing p-values
print("Comparing p-values of variants with and without a functional annotation...")
p_anno <- reg_df$PVALUE_MLOG

print("Median p-value of promoter/enhancer variants:")
print(median(p_anno))

print("Range of promoter/enhancer variant p-values:")
print(range(p_anno))

p_noanno <- noreg_df$PVALUE_MLOG

print("Median p-value of variants without a functional annotation:")
print(median(p_noanno))

print("Range of p-values of variants without a functional annotation:")
print(range(p_noanno))

wilcox.test(p_anno, p_noanno)

#Comparing allele frequency
print("Comparing MAFs of variants with and without a functional annotation...")
maf_anno <- reg_df$RISK_ALLELE_FREQUENCY
maf_anno <- as.numeric(maf_anno)

print("Median MAF of promoter/enhancer variants:")
print(median(maf_anno))

print("Range of promoter/enhancer variant MAFs:")
print(range(maf_anno))

maf_noanno <- noreg_df$RISK_ALLELE_FREQUENCY
maf_noanno <- as.numeric(maf_noanno)

print("Median MAF of variants without a functional annotation:")
print(median(maf_noanno))

print("Range of MAFs of variants without a functional annotation:")
print(range(maf_noanno))

wilcox.test(maf_anno, maf_noanno)