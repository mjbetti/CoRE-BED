library("data.table")
library("optparse")

option_list = list(
	make_option(c("-c", "--chr"), type = "integer", default = NULL, help = "chromosome number", metavar = "integer"))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#Loop through each sliced file, concatenating them to a growing data frame and then outputting this to a final file
tissues <- c("heart_aorta_34m", "heart_ascending_aorta_51f", "heart_ascending_aorta_53f", "heart_coronary_artery_51f", "heart_coronary_artery_53f", "heart_101_days", "heart_59_days_76_days_f", "heart_80_days", "heart_96_days", "heart_103_days_f", "heart_105_days_f", "heart_110_days_f", "heart_116_days_98_days_f", "heart_116_days_117_days_f", "heart_147_days_f", "heart_91_days_f", "heart_left_ventricle_53f", "heart_left_ventricle_101_days_103_days_f", "heart_left_ventricle_136_days_f", "heart_left_ventricle_34m", "heart_left_ventricle_3m", "heart_27m_35m", "heart_3m", "heart_105_days_m", "heart_110_days_m", "heart_120_days_m", "heart_72_days_76_days_m", "heart_91_days_m", "heart_96_days_m", "heart_right_ventricle_101_days_103_days_f", "heart_right_ventricle_34m", "heart_right_ventricle_3m", "heart_left_atrium_101_days_f", "heart_right_atrium_51f", "heart_right_atrium_53f", "heart_right_atrium_34m", "heart_thoracic_aorta_37m", "heart_thoracic_aorta_54m", "heart_tibial_artery_53f", "heart_tibial_artery_37m")

concat_file <- fread(paste0("/home/bettimj/gamazon_rotation/mod_core-bed/uk_biobank/h2/anno_1000G/chr", opt$chr, "/1000G.EUR.QC.", opt$chr, ".heart_aorta_30f.bim"), header = FALSE, sep = "\t", quote = FALSE)
concat_file <- concat_file[,1]
concat_file <- as.data.frame(concat_file)
print(head(concat_file))

class1 <- concat_file
class1[(class1 == 1)] <- "func"
class1[!(class1 == "func")] <- 0
class1[(class1 == "func")] <- 1

class2 <- concat_file
class2[(class2 == 2)] <- "func"
class2[!(class2 == "func")] <- 0
class2[(class2 == "func")] <- 1

class3 <- concat_file
class3[(class3 == 3)] <- "func"
class3[!(class3 == "func")] <- 0
class3[(class3 == "func")] <- 1

class4 <- concat_file
class4[(class4 == 4)] <- "func"
class4[!(class4 == "func")] <- 0
class4[(class4 == "func")] <- 1

class5 <- concat_file
class5[(class5 == 5)] <- "func"
class5[!(class5 == "func")] <- 0
class5[(class5 == "func")] <- 1

class6 <- concat_file
class6[(class6 == 6)] <- "func"
class6[!(class6 == "func")] <- 0
class6[(class6 == "func")] <- 1

class7 <- concat_file
class7[(class7 == 7)] <- "func"
class7[!(class7 == "func")] <- 0
class7[(class7 == "func")] <- 1

class8 <- concat_file
class8[(class8 == 8)] <- "func"
class8[!(class8 == "func")] <- 0
class8[(class8 == "func")] <- 1

any_anno <- concat_file
any_anno[(any_anno == 1 | any_anno == 2 | any_anno == 3 | any_anno == 4 | any_anno == 5 | any_anno == 6 | any_anno == 7 | any_anno == 8)] <- 1
any_anno[!(any_anno == 1)] <- 0

for (tissue in tissues) {
	print(paste0("Processing ", tissue, "..."))
	name <- paste0("/home/bettimj/gamazon_rotation/mod_core-bed/uk_biobank/h2/anno_1000G/chr", opt$chr, "/1000G.EUR.QC.", opt$chr, ".", tissue, ".bim")
	file <- fread(name, header = FALSE, sep = "\t", quote = "")
	file <- file[,1]
	file <- as.data.frame(file)
	
	class1[,2] <- file
	class1[,2][class1[,2] == 1] <- "func"
	class1[,2][!(class1[,2] == "func")] <- 0
	class1[,2][(class1[,2] == "func")] <- 1
	class1[,1][class1[,1] == 0 & class1[,2] == 1] <- 1

	class2[,2] <- file
	class2[,2][class2[,2] == 2] <- "func"
	class2[,2][!(class2[,2] == "func")] <- 0
	class2[,2][(class2[,2] == "func")] <- 1
	class2[,1][class2[,1] == 0 & class2[,2] == 1] <- 1
	
	class3[,2] <- file
	class3[,2][class3[,2] == 3] <- "func"
	class3[,2][!(class3[,2] == "func")] <- 0
	class3[,2][(class3[,2] == "func")] <- 1
	class3[,1][class3[,1] == 0 & class3[,2] == 1] <- 1
	
	class4[,2] <- file
	class4[,2][class4[,2] == 4] <- "func"
	class4[,2][!(class4[,2] == "func")] <- 0
	class4[,2][(class4[,2] == "func")] <- 1
	class4[,1][class4[,1] == 0 & class4[,2] == 1] <- 1
	
	class5[,2] <- file
	class5[,2][class5[,2] == 5] <- "func"
	class5[,2][!(class5[,2] == "func")] <- 0
	class5[,2][(class5[,2] == "func")] <- 1
	class5[,1][class5[,1] == 0 & class5[,2] == 1] <- 1
	
	class6[,2] <- file
	class6[,2][class6[,2] == 6] <- "func"
	class6[,2][!(class6[,2] == "func")] <- 0
	class6[,2][(class6[,2] == "func")] <- 1
	class6[,1][class6[,1] == 0 & class6[,2] == 1] <- 1
	
	class7[,2] <- file
	class7[,2][class7[,2] == 7] <- "func"
	class7[,2][!(class7[,2] == "func")] <- 0
	class7[,2][(class7[,2] == "func")] <- 1
	class7[,1][class7[,1] == 0 & class7[,2] == 1] <- 1

	class8[,2] <- file
	class8[,2][class8[,2] == 8] <- "func"
	class8[,2][!(class8[,2] == "func")] <- 0
	class8[,2][(class8[,2] == "func")] <- 1
	class8[,1][class8[,1] == 0 & class8[,2] == 1] <- 1
	
	any_anno[,2] <- file
	any_anno[,2][any_anno[,2] == 1 | any_anno[,2] == 2 | any_anno[,2] == 3 | any_anno[,2] == 4 | any_anno[,2] == 5 | any_anno[,2] == 6 | any_anno[,2] == 7 | any_anno[,2] == 8] <- 1
	any_anno[,2][!(any_anno[,2] == 1)] <- 0
	any_anno[,1][any_anno[,1] == 0 & any_anno[,2] == 1] <- 1
}

concat_file <- data.frame(class1[,1], class2[,1], class3[,1], class4[,1], class5[,1], class6[,1], class7[,1], class8[,1], any_anno[,1])
names(concat_file) <- c("class1", "class2", "class3", "class4", "class5", "class6", "class7", "class8", "any_anno")

print(head(concat_file))

write.table(concat_file, file = paste0("/home/bettimj/gamazon_rotation/mod_core-bed/uk_biobank/h2/anno_1000G/hypertension/1000G.EUR.QC.", opt$chr, ".annot"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)