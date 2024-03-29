library("data.table")
library("gap")
library("qqman")

#Declare the paths of the annotated summary statistics file
sum_stats_path <- "/home/bettimj/gamazon_rotation/mod_core-bed/uk_biobank/att2_bin_annos_epimap_phecode-401-both_sexes.func_anno.txt"
  
#Open the split annotated summary statistics as a data frame
sum_stats_file <- fread(sum_stats_path, header = TRUE, sep = "\t", quote = "")
sum_stats_df <- as.data.frame(sum_stats_file)

#Split the summary stats based on whether or not each coordinate has a CoRE-BED Impute annotation in any of the 816 cell/tissue types
reg_sum_stats_df <- sum_stats_df[(sum_stats_df$CORE_BED_ANNO == 1),]
noreg_sum_stats_df <- sum_stats_df[(sum_stats_df$CORE_BED_ANNO == 0),]

reg_pvals_eur <- reg_sum_stats_df$pval_EUR[!is.na(reg_sum_stats_df$pval_EUR)]
reg_pvals_meta <- reg_sum_stats_df$pval_meta[!is.na(reg_sum_stats_df$pval_meta)]

noreg_pvals_eur <- noreg_sum_stats_df$pval_EUR[!is.na(noreg_sum_stats_df$pval_EUR)]
noreg_pvals_meta <- noreg_sum_stats_df$pval_meta[!is.na(noreg_sum_stats_df$pval_meta)]

#Plot a Q-Q plot using the meta-analysis p-values from each file
##European p-values
par(mar=c(1,1,1,1))
tiff("mod_superimposed_qq_ht_eur_x8_y50.tiff", res = 300, width = 5, height = 5, units = 'in', compression = c( "lzw"))
qqunif(reg_pvals_eur, col = "red", lcol = "black", xlim = c(0, 8), ylim = c(0, 50))
par(new=T)
qqunif(noreg_pvals_eur, col = "blue", lcol = "black", xlim = c(0, 8), ylim = c(0, 50), axes = F)
legend('topleft', c("CoRE-BED annotation", "no CoRE-BED annotation"),
             fill = c("red", "blue"))
dev.off()

##Meta analysis p-values
par(mar=c(1,1,1,1))
tiff("mod_superimposed_qq_ht_meta_x8_y50.tiff", res = 300, width = 5, height = 5, units = 'in', compression = c( "lzw"))
qqunif(reg_pvals_meta, col = "red", lcol = "black", xlim = c(0, 8), ylim = c(0, 50))
par(new=T)
qqunif(noreg_pvals_meta, col = "blue", lcol = "black", xlim = c(0, 8), ylim = c(0, 50), axes = F)
legend('topleft', c("CoRE-BED annotation", "no CoRE-BED annotation"),
             fill = c("red", "blue"))
dev.off()

##Using y max of 25
###European p-values
par(mar=c(1,1,1,1))
tiff("att2_mod_superimposed_qq_ht_eur_x5_y25.tiff", res = 300, width = 5, height = 5, units = 'in', compression = c( "lzw"))
qqunif(reg_pvals_eur, col = "red", lcol = "black", xlim = c(0, 5), ylim = c(0, 25))
par(new=T)
qqunif(noreg_pvals_eur, col = "blue", lcol = "black", xlim = c(0, 5), ylim = c(0, 25), axes = F)
legend('topleft', c("CoRE-BED annotation", "no CoRE-BED annotation"),
             fill = c("red", "blue"))
dev.off()

##Meta analysis p-values
par(mar=c(1,1,1,1))
tiff("att2_mod_superimposed_qq_ht_meta_x5_y25.tiff", res = 300, width = 5, height = 5, units = 'in', compression = c( "lzw"))
qqunif(reg_pvals_meta, col = "red", lcol = "black", xlim = c(0, 5), ylim = c(0, 25))
par(new=T)
qqunif(noreg_pvals_meta, col = "blue", lcol = "black", xlim = c(0, 5), ylim = c(0, 25), axes = F)
legend('topleft', c("CoRE-BED annotation", "no CoRE-BED annotation"),
             fill = c("red", "blue"))
dev.off()

##Using x range of 1.5 to 5, y max of 25
###European p-values
par(mar=c(1,1,1,1))
tiff("att2_mod_superimposed_qq_ht_eur_x1.5_5_y25.tiff", res = 300, width = 5, height = 5, units = 'in', compression = c( "lzw"))
qqunif(reg_pvals_eur, col = "red", lcol = "black", xlim = c(1.5, 5), ylim = c(0, 25))
par(new=T)
qqunif(noreg_pvals_eur, col = "blue", lcol = "black", xlim = c(1.5, 5), ylim = c(0, 25), axes = F)
legend('topleft', c("CoRE-BED annotation", "no CoRE-BED annotation"),
             fill = c("red", "blue"))
dev.off()

##Meta analysis p-values
par(mar=c(1,1,1,1))
tiff("att2_mod_superimposed_qq_ht_meta_x1.5_5_y25.tiff", res = 300, width = 5, height = 5, units = 'in', compression = c( "lzw"))
qqunif(reg_pvals_meta, col = "red", lcol = "black", xlim = c(1.5, 5), ylim = c(0, 25))
par(new=T)
qqunif(noreg_pvals_meta, col = "blue", lcol = "black", xlim = c(1.5, 5), ylim = c(0, 25), axes = F)
legend('topleft', c("CoRE-BED annotation", "no CoRE-BED annotation"),
             fill = c("red", "blue"))
dev.off()