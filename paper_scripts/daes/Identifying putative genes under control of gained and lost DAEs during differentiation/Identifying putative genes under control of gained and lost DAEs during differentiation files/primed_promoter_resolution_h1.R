library("data.table")

#Declare the paths
path_esc <- "/home/bettimj/gamazon_rotation/mod_core-bed/all_anno/esc_h1-hesc_hg19.genome.1kb.bed.gz"
path_npc <- "/home/bettimj/gamazon_rotation/mod_core-bed/all_anno/esc_deriv_neural_progenitor_h1-hesc_origin_hg19.genome.1kb.bed.gz"
path_nc <- "/home/bettimj/gamazon_rotation/mod_core-bed/all_anno/esc_deriv_neural_crest_h1-hesc_origin_hg19.genome.1kb.bed.gz"

path_esc_raw <- "/home/bettimj/gamazon_rotation/mod_core-bed/epigenomes_expanded/esc_h1-hesc_epigenome_hg19.txt.gz"
path_npc_raw <- "/home/bettimj/gamazon_rotation/mod_core-bed/epigenomes_expanded/esc_deriv_neural_progenitor_h1-hesc_origin_epigenome_hg19.txt.gz"
path_nc_raw <- "/home/bettimj/gamazon_rotation/mod_core-bed/epigenomes_expanded/esc_deriv_neural_crest_h1-hesc_origin_epigenome_hg19.txt.gz"

#Open as data frames
file_esc <- fread(path_esc, header = FALSE, sep = "\t", quote = "")
df_esc <- as.data.frame(file_esc)

file_npc <- fread(path_npc, header = FALSE, sep = "\t", quote = "")
df_npc <- as.data.frame(file_npc)

file_nc <- fread(path_nc, header = FALSE, sep = "\t", quote = "")
df_nc <- as.data.frame(file_nc)

file_esc_raw <- fread(path_esc_raw, header = TRUE, sep = "\t", quote = "")
df_esc_raw <- as.data.frame(file_esc_raw)
df_esc_raw$coord <- paste(df_esc_raw$chrom, df_esc_raw$start, df_esc_raw$end, sep = "_")

file_npc_raw <- fread(path_npc_raw, header = TRUE, sep = "\t", quote = "")
df_npc_raw <- as.data.frame(file_npc_raw)
df_npc_raw$coord <- paste(df_npc_raw$chrom, df_npc_raw$start, df_npc_raw$end, sep = "_")

file_nc_raw <- fread(path_nc_raw, header = TRUE, sep = "\t", quote = "")
df_nc_raw <- as.data.frame(file_nc_raw)
df_nc_raw$coord <- paste(df_nc_raw$chrom, df_nc_raw$start, df_nc_raw$end, sep = "_")

#ESC to NC
#For the initial two files, retain only those coordinates classified as as Class 4
new_df <- data.frame(df_esc[,c(1:3)], df_esc[,4], df_nc[,4])
new_df <- new_df[(new_df[,4] == 4 | new_df[,5] == 4),]

print("Total number in ESC:")
print(nrow(new_df[new_df[,4] == 4,]))
#911

print("Total number in ESC-deriv:")
print(nrow(new_df[new_df[,5] == 4,]))
#65

print("Total number same:")
print(nrow(new_df[(new_df[,4] == 4 & new_df[,5] == 4),]))
#41

print("Total number lost from ESC to ESC-deriv:")
print(nrow(new_df[(new_df[,4] == 4 & !(new_df[,5] == 4)),]))
#870

print("Total number gained from ESC to ESC-deriv:")
print(nrow(new_df[(!(new_df[,4] == 4) & new_df[,5] == 4),]))
#24

#For those lost, what do the lost marks revolve to?
print("Class 0:")
print(nrow(new_df[(new_df[,4] == 4 & new_df[,5] == 0),]))
#847

print("Class 1:")
print(nrow(new_df[(new_df[,4] == 4 & new_df[,5] == 1),]))
#10

print("Class 2:")
print(nrow(new_df[(new_df[,4] == 4 & new_df[,5] == 2),]))
#0

print("Class 3:")
print(nrow(new_df[(new_df[,4] == 4 & new_df[,5] == 3),]))
#3

print("Class 5:")
print(nrow(new_df[(new_df[,4] == 4 & new_df[,5] == 5),]))
#7

print("Class 6:")
print(nrow(new_df[(new_df[,4] == 4 & new_df[,5] == 6),]))
#0

print("Class 7:")
print(nrow(new_df[(new_df[,4] == 4 & new_df[,5] == 7),]))
#3

print("Class 8:")
print(nrow(new_df[(new_df[,4] == 4 & new_df[,5] == 8),]))
#0

#For those resolving to Class 0, find the changes in epigenomic marks
#new_df_lost <- new_df[(new_df[,4] == 4),c(1:3)]
new_df_lost <- new_df[(new_df[,4] == 4 & !(new_df[,5] == 4)),c(1:3)]

write.table(new_df_lost, file = "lost_h1_to_neuron.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)