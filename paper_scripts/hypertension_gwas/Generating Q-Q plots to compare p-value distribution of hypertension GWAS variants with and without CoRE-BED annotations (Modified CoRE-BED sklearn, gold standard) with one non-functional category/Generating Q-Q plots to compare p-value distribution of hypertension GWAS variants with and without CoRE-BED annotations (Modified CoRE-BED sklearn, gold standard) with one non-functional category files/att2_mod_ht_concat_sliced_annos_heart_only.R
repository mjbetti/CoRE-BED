library("data.table")

#Loop through each sliced file, concatenating them to a growing data frame and then outputting this to a final file
tissues <- c("heart_aorta_30f", "heart_aorta_34m", "heart_ascending_aorta_51f", "heart_ascending_aorta_53f", "heart_coronary_artery_51f", "heart_coronary_artery_53f", "heart_101_days", "heart_59_days_76_days_f", "heart_80_days", "heart_96_days", "heart_103_days_f", "heart_105_days_f", "heart_110_days_f", "heart_116_days_98_days_f", "heart_116_days_117_days_f", "heart_147_days_f", "heart_91_days_f", "heart_left_ventricle_53f", "heart_left_ventricle_101_days_103_days_f", "heart_left_ventricle_136_days_f", "heart_left_ventricle_34m", "heart_left_ventricle_3m", "heart_27m_35m", "heart_3m", "heart_105_days_m", "heart_110_days_m", "heart_120_days_m", "heart_72_days_76_days_m", "heart_91_days_m", "heart_96_days_m", "heart_right_ventricle_101_days_103_days_f", "heart_right_ventricle_34m", "heart_right_ventricle_3m", "heart_left_atrium_101_days_f", "heart_right_atrium_51f", "heart_right_atrium_53f", "heart_right_atrium_34m", "heart_thoracic_aorta_37m", "heart_thoracic_aorta_54m", "heart_tibial_artery_53f", "heart_tibial_artery_37m")

concat_file <- fread("epimap_phecode-401-both_sexes.func_anno.heart_101_days.tsv", header = TRUE, sep = "\t", quote = FALSE)
concat_file <- data.frame(concat_file)
concat_file[(concat_file == 1 | concat_file == 2 | concat_file == 3 | concat_file == 4 | concat_file == 5 | concat_file == 6 | concat_file == 7 | concat_file == 8)] <- "functional"
concat_file[(concat_file == 0 | concat_file == 1)] <- "nonfunctional"

concat_file[(concat_file == "functional")] <- 1
concat_file[(concat_file == "nonfunctional")] <- 0

for (tissue in tissues) {
	print(paste0("Processing ", tissue, "..."))
	name <- paste0("epimap_phecode-401-both_sexes.func_anno.", tissue, ".tsv")
	file <- fread(name, header = TRUE, sep = "\t", quote = "")
	file <- as.data.frame(file)
	concat_file[,2] <- file
	concat_file[,2][concat_file[,2] == 1 | concat_file[,2] == 2 | concat_file[,2] == 3 | concat_file[,2] == 4 | concat_file[,2] == 5 | concat_file[,2] == 6 | concat_file[,2] == 7 | concat_file[,2] == 8] <- "functional"
concat_file[,2][(concat_file[,2] == 0 | concat_file[,2] == 1)] <- "nonfunctional"

concat_file[,2][(concat_file[,2] == "functional")] <- 1
concat_file[,2][(concat_file[,2] == "nonfunctional")] <- 0

concat_file[,1][concat_file[,1] == 0 & concat_file[,2] == 1] <- 1
}

concat_file <- concat_file[,1]

write.table(concat_file, file = "att2_epimap_phecode-401-both_sexes.func_anno.all_heart.tsv", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)