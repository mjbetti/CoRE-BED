library("data.table")

#Declare the path of the global CoRE-BED elements directory
lost_path <- "/home/bettimj/gamazon_rotation/mod_core-bed/resolve_primed_promoters/genes_lost_h1_to_neuron.bed"
gained_path <- "/home/bettimj/gamazon_rotation/mod_core-bed/resolve_primed_promoters/genes_gained_h9_to_sm_muscle.bed"
out_dir <- "/home/bettimj/gamazon_rotation/mod_core-bed/resolve_primed_promoters"


#Lost
file <- fread(lost_path, header = FALSE, quote = FALSE, sep = "\t")
df <- as.data.frame(file)
df <- unique(df[,11])
write.table(df, file = paste(out_dir, paste0("unique_genes_lost_h1_to_neuron.bed"), sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

#Gained
file <- fread(gained_path, header = FALSE, quote = FALSE, sep = "\t")
df <- as.data.frame(file)
df <- unique(df[,11])
write.table(df, file = paste(out_dir, paste0("unique_genes_gained_h9_to_sm_muscle.bed"), sep = "/"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
