library("data.table")

orig_path <- "/home/bettimj/gamazon_rotation/core-bed_analysis/s2g_snps/exp_val_snps/exp_val_snps.txt"
anno_path <- "/home/bettimj/gamazon_rotation/core-bed_analysis/s2g_snps/exp_val_snps/sliced_annos/encoded_anno_status_all.exp_val_snps.func_anno.5000_1000.txt"

orig_file <- fread(orig_path, header = FALSE, sep = " ", quote = "")
orig_df <- as.data.frame(orig_file)

anno_file <- fread(anno_path, header = TRUE, sep = "\t", quote = "")
anno_df <- as.data.frame(anno_file)

cat_df <- data.frame(orig_df, anno_df)
names(cat_df) <- c("chr", "pos", "rsid", "validated_gene", "pheno", "cs2g_pred_gene", "cs2g_score", "pred_method", "active_promoter", "bivalent_promoter", "silenced_promoter", "active_enhancer", "poised_enhancer", "primed_enhancer", "any_promoter", "any_enhancer", "any_anno")

write.table(cat_df, file = "all_annos_exp_val_snps.5000_1000.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)