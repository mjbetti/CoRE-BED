library("data.table")

#Declare the paths of the CaVEMaN, CAVIAR, and DAPG high-confidence, fine-mapped GTEx v8 eQTLs
caveman_path <- "/home/bettimj/gamazon_rotation/gtex8_eqtl_sqtl/eqtl/fine_mapped/GTEx_v8_finemapping_CaVEMaN/GTEx_v8_finemapping_CaVEMaN.txt.gz"
caviar_path <- "/home/bettimj/gamazon_rotation/gtex8_eqtl_sqtl/eqtl/fine_mapped/GTEx_v8_finemapping_CAVIAR/CAVIAR_Results_v8_GTEx_LD_HighConfidentVariants.gz"
dapg_path <- "/home/bettimj/gamazon_rotation/gtex8_eqtl_sqtl/eqtl/fine_mapped/GTEx_v8_finemapping_DAPG/GTEx_v8_finemapping_DAPG.CS95.txt.gz"

#Open each of the files as data frames and convert to bed format (chr, start, end)
caveman_file <- fread(caveman_path, header = TRUE, sep = "\t", quote = "")
caveman_df <- as.data.frame(caveman_file)
caveman_df <- data.frame(caveman_df[,4], caveman_df[,5], caveman_df[,5])
names(caveman_df) <- c("chr", "start", "end")

caviar_file <- fread(caviar_path, header = TRUE, sep = "\t", quote = "")
caviar_df <- as.data.frame(caviar_file)
caviar_df <- data.frame(caviar_df[,4], caviar_df[,5], caviar_df[,5])
caviar_df[,1] <- paste0("chr", caviar_df[,1])
names(caviar_df) <- c("chr", "start", "end")

dapg_file <- fread(dapg_path, header = TRUE, sep = "\t", quote = "")
dapg_df <- as.data.frame(dapg_file)
dapg_df <- data.frame(dapg_df[,1], dapg_df[,2], dapg_df[,2])
dapg_df[,1] <- paste0("chr", dapg_df[,1])
names(dapg_df) <- c("chr", "start", "end")

#Combine the three data frames together, and filter out duplicates
combined_df <- rbind(caveman_df, caviar_df)
combined_df <- rbind(combined_df, dapg_df)

combined_df <- unique(combined_df)

#Write the final data frame out to a new bed file
write.table(combined_df, file = "gtexv8_fine-mapped_hc_eqtls.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)