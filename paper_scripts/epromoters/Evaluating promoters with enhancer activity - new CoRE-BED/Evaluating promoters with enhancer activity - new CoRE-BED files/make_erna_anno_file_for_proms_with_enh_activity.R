library("data.table")

tissues <- c("Adipose_sub", "Adipose_vis", "Adrenal", "Artery_aor", "Artery_cor", "Artery_tib", "Bladder", "Brain_Amy", "Brain_Ant", "Brain_Cau", "Brain_Cerebellum", "Brain_Cer", "Brain_Cortex", "Brain_Fro", "Brain_Hippo", "Brain_Hypo", "Brain_Nuc", "Brain_Put", "Brain_Spi", "Brain_Sub", "Breast", "Cells_ebv", "Cells_fibroblast", "Cells_leukemia", "Cervix_ecto", "Cervix_endo", "Colon_sigmoid", "Colon_transverse", "Esophagus_gast", "Esophagus_muc", "Esophagus_mus", "Fallopian_tube", "Heart_Atr", "Heart_Left", "Kidney", "Liver", "Lung", "Minor_salivary_gland", "Muscle_skeletal", "Nerve", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_no_sun", "Skin_sun", "Small_intestine", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_blood")

anno_path <- "/home/bettimj/gamazon_rotation/eRNA_predixcan_models/attempt3/ensembl87_fantom5_roadmap_hg19_reg_annos_hg19.txt"
anno_file <- fread(anno_path, header = TRUE, sep = "\t", quote = "")
anno_df <- as.data.frame(anno_file)

erna_names <- c()
for (tissue in tissues) {
	path <- paste0("/home/bettimj/gamazon_rotation/herna/reformatted/", tissue, ".rds.erna_exp.txt")
	file <- fread(path, header = TRUE, sep = "\t", quote = "")
	df <- as.data.frame(file)
	erna_name_vec <- df[,1]
	erna_names <- append(erna_names, erna_name_vec)
}

erna_names <- unique(erna_names)

coord_df <- anno_df[(anno_df$geneid %in% erna_names),]

#Write out as a bed file
write.table(coord_df, file = "/home/bettimj/gamazon_rotation/herna/reformatted/erna_all_tissues.hg19.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)