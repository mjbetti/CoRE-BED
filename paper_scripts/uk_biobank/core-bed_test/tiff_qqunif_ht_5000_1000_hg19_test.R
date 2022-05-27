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

#Plot a Q-Q plot using the meta-analysis p-values from each file
##Meta analysis p-values
par(mar=c(1,1,1,1))
tiff("superimposed_qq_ht_5000_1000_meta_x8_y50.hg19_test.tiff", res = 300, width = 5, height = 5, units = 'in', compression = c( "lzw"))
qqunif(reg_pvals_meta, col = "red", lcol = "black", xlim = c(0, 8), ylim = c(0, 50))
par(new=T)
qqunif(noreg_pvals_meta, col = "blue", lcol = "black", xlim = c(0, 8), ylim = c(0, 50), axes = F)
legend('topleft', c("CoRE-BED-test annotation", "no CoRE-BED-test annotation"),
             fill = c("red", "blue"))
dev.off()