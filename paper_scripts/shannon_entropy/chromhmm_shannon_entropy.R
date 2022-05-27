library("data.table")

#Declare the paths of the CoRE-BED-training results files across the five decision tree structures reported
path_all <- "/home/bettimj/gamazon_rotation/epimap_chromhmm/hg38/personal.broadinstitute.org/cboix/epimap/ChromHMM/observed_aux_18_hg38/CALLS/all_epimap_reg_elements_18_CALLS_segments_hg38.bed"
path_enh <- "/home/bettimj/gamazon_rotation/epimap_chromhmm/hg38/personal.broadinstitute.org/cboix/epimap/ChromHMM/observed_aux_18_hg38/CALLS/all_epimap_enh_18_CALLS_segments_hg38.bed"

#Calculate Shannon entropy for ChromHMM (all regulatory states and only enhancers)
##All states
file_all <- fread(path_all, header = FALSE, sep = "\t", quote = "")
df_all <- as.data.frame(file_all)
total_all <- nrow(df_all)
quies <- nrow(df_all[(df_all[,4] == "Quies"),])
enha2 <- nrow(df_all[(df_all[,4] == "EnhA2"),])
tssa <- nrow(df_all[(df_all[,4] == "TssA"),])
tssflnk <- nrow(df_all[(df_all[,4] == "TssFlnk"),])
enhwk <- nrow(df_all[(df_all[,4] == "EnhWk"),])
reprpcwk <- nrow(df_all[(df_all[,4] == "ReprPCWk"),])
enha1 <- nrow(df_all[(df_all[,4] == "EnhA1"),])
tssflnku <- nrow(df_all[(df_all[,4] == "TssFlnkU"),])
tssbiv <- nrow(df_all[(df_all[,4] == "TssBiv"),])
enhbiv <- nrow(df_all[(df_all[,4] == "EnhBiv"),])
reprpc <- nrow(df_all[(df_all[,4] == "ReprPC"),])
enhg2 <- nrow(df_all[(df_all[,4] == "EnhG2"),])
enhg1 <- nrow(df_all[(df_all[,4] == "EnhG1"),])
tx <- nrow(df_all[(df_all[,4] == "Tx"),])
txwk <- nrow(df_all[(df_all[,4] == "TxWk"),])
tssflnkd <- nrow(df_all[(df_all[,4] == "TssFlnkD"),])
het <- nrow(df_all[(df_all[,4] == "Het"),])
znf_rpts <- nrow(df_all[(df_all[,4] == "ZNF/Rpts"),])
shannon_ent <- -(((quies/total_all)*log2(quies/total_all)) + ((enha2/total_all)*log2(enha2/total_all)) + ((tssa/total_all)*log2(tssa/total_all)) + ((tssflnk/total_all)*log2(tssflnk/total_all)) + ((enhwk/total_all)*log2(enhwk/total_all)) + ((reprpcwk/total_all)*log2(reprpcwk/total_all)) + ((enha1/total_all)*log2(enha1/total_all)) + ((tssflnku/total_all)*log2(tssflnku/total_all)) + ((tssbiv/total_all)*log2(tssbiv/total_all)) + ((enhbiv/total_all)*log2(enhbiv/total_all)) + ((reprpc/total_all)*log2(reprpc/total_all)) + ((enhg2/total_all)*log2(enhg2/total_all)) + ((enhg1/total_all)*log2(enhg1/total_all)) + ((tx/total_all)*log2(tx/total_all)) + ((txwk/total_all)*log2(txwk/total_all)) + ((tssflnkd/total_all)*log2(tssflnkd/total_all)) + ((het/total_all)*log2(het/total_all)) + ((znf_rpts/total_all)*log2(znf_rpts/total_all)))
print("All states:")
print(shannon_ent)

##Enhancers
total_enh <- (enha2 + enhwk + enha1 + enhbiv + enhg2 + enhg1)
shannon_ent <- -(((enha2/total_enh)*log2(enha2/total_enh)) + ((enhwk/total_enh)*log2(enhwk/total_enh)) + ((enha1/total_enh)*log2(enha1/total_enh)) + ((enhbiv/total_enh)*log2(enhbiv/total_enh)) + ((enhg2/total_enh)*log2(enhg2/total_enh)) + ((enhg1/total_enh)*log2(enhg1/total_enh)))
print("Enhancer states:")
print(shannon_ent)