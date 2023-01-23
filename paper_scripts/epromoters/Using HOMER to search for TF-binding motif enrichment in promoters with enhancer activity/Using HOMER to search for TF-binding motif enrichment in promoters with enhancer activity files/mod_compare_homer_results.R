library("data.table")

proms_with_enh_path <- "/home/bettimj/gamazon_rotation/mod_core-bed/homer/proms_with_enh_activity/knownResults.txt"
other_shared_path <- "/home/bettimj/gamazon_rotation/mod_core-bed/homer/shared_no_prom_with_enh_activity/knownResults.txt"

#Open each of the paths as a data frame
proms_with_enh_file <- fread(proms_with_enh_path, header = TRUE, sep = "\t", quote = "")
proms_with_enh_df <- as.data.frame(proms_with_enh_file)
nrow(proms_with_enh_df)
#This p-value threshold was tested because it corrects for the 440 tests we ran (#motifs)
proms_with_enh_df <- proms_with_enh_df[(proms_with_enh_df[,3] <= 1.13e-4),]

other_shared_file <- fread(other_shared_path, header = TRUE, sep = "\t", quote = "")
other_shared_df <- as.data.frame(other_shared_file)
nrow(other_shared_df)
#This p-value threshold was tested because it corrects for the 440 tests we ran (#motifs)
other_shared_df <- other_shared_df[(other_shared_df[,3] <= 1.13e-4),]

proms_with_enh_df <- proms_with_enh_df[,1]
#other_shared_df <- other_shared_df[1:23,1]
other_shared_df <- other_shared_df[,1]

sig_uniq_prom_enh <- proms_with_enh_df[!(proms_with_enh_df %in% other_shared_df)]

length(sig_uniq_prom_enh)
#[1] 130

write.table(sig_uniq_prom_enh, file = "/home/bettimj/gamazon_rotation/mod_core-bed/homer/uniq_sig_to_enhlike_promoters.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)