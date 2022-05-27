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
reg_pvals_eur <- reg_anno_df$pval_EUR[!is.na(reg_anno_df$pval_EUR)]
reg_pvals_meta <- reg_anno_df$pval_meta[!is.na(reg_anno_df$pval_meta)]

noreg_pvals_eur <- noreg_anno_df$pval_EUR[!is.na(noreg_anno_df$pval_EUR)]
noreg_pvals_meta <- noreg_anno_df$pval_meta[!is.na(noreg_anno_df$pval_meta)]

#Plot a Q-Q plot using the meta-analysis p-values from each file
##Meta analysis p-values
par(mar=c(1,1,1,1))
tiff("heart_superimposed_qq_ht_5000_1000_meta_x8_y50.hg19.test.tiff", res = 300, width = 5, height = 5, units = 'in', compression = c( "lzw"))
qqunif(reg_pvals_meta, col = "red", lcol = "black", xlim = c(0, 8), ylim = c(0, 50))
par(new=T)
qqunif(noreg_pvals_meta, col = "blue", lcol = "black", xlim = c(0, 8), ylim = c(0, 50), axes = F)
legend('topleft', c("CoRE-BED-test annotation", "no CoRE-BED-test annotation"),
             fill = c("red", "blue"))
dev.off()