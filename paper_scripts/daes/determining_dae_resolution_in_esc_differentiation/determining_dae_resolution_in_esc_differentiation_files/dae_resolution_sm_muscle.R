library("data.table")

#Declare the paths
path_esc <- "/home/bettimj/gamazon_rotation/mod_core-bed/all_anno/esc_h9_hg19.genome.1kb.bed.gz"
path_npc <- "/home/bettimj/gamazon_rotation/mod_core-bed/all_anno/esc_deriv_smooth_muscle_h9_origin_hg19.genome.1kb.bed.gz"

path_esc_raw <- "/home/bettimj/gamazon_rotation/mod_core-bed/epigenomes_expanded/esc_h9_epigenome_hg19.txt.gz"
path_npc_raw <- "/home/bettimj/gamazon_rotation/mod_core-bed/epigenomes_expanded/esc_deriv_smooth_muscle_h9_origin_epigenome_hg19.txt.gz"

#Open as data frames
file_esc <- fread(path_esc, header = FALSE, sep = "\t", quote = "")
df_esc <- as.data.frame(file_esc)

file_npc <- fread(path_npc, header = FALSE, sep = "\t", quote = "")
df_npc <- as.data.frame(file_npc)

file_esc_raw <- fread(path_esc_raw, header = TRUE, sep = "\t", quote = "")
df_esc_raw <- as.data.frame(file_esc_raw)
df_esc_raw$coord <- paste(df_esc_raw$chrom, df_esc_raw$start, df_esc_raw$end, sep = "_")

file_npc_raw <- fread(path_npc_raw, header = TRUE, sep = "\t", quote = "")
df_npc_raw <- as.data.frame(file_npc_raw)
df_npc_raw$coord <- paste(df_npc_raw$chrom, df_npc_raw$start, df_npc_raw$end, sep = "_")

#ESC to NPC
#For the initial two files, retain only those coordinates classified as as Class 4
new_df <- data.frame(df_esc[,c(1:3)], df_esc[,4], df_npc[,4])
new_df <- new_df[(new_df[,4] == 4 | new_df[,5] == 4),]

print("Total number in ESC:")
print(nrow(new_df[new_df[,4] == 4,]))
#176

print("Total number in ESC-deriv:")
print(nrow(new_df[new_df[,5] == 4,]))
#3092

print("Total number same:")
print(nrow(new_df[(new_df[,4] == 4 & new_df[,5] == 4),]))
#116

print("Total number lost from ESC to ESC-deriv:")
print(nrow(new_df[(new_df[,4] == 4 & !(new_df[,5] == 4)),]))
#60

print("Total number gained from ESC to ESC-deriv:")
print(nrow(new_df[(!(new_df[,4] == 4) & new_df[,5] == 4),]))
#2976

#For those lost, what do the lost marks revolve to?
print("Class 0:")
print(nrow(new_df[(new_df[,4] == 4 & new_df[,5] == 0),]))
#32

print("Class 1:")
print(nrow(new_df[(new_df[,4] == 4 & new_df[,5] == 1),]))
#2

print("Class 2:")
print(nrow(new_df[(new_df[,4] == 4 & new_df[,5] == 2),]))
#0

print("Class 3:")
print(nrow(new_df[(new_df[,4] == 4 & new_df[,5] == 3),]))
#11

print("Class 5:")
print(nrow(new_df[(new_df[,4] == 4 & new_df[,5] == 5),]))
#4

print("Class 6:")
print(nrow(new_df[(new_df[,4] == 4 & new_df[,5] == 6),]))
#0

print("Class 7:")
print(nrow(new_df[(new_df[,4] == 4 & new_df[,5] == 7),]))
#1

print("Class 8:")
print(nrow(new_df[(new_df[,4] == 4 & new_df[,5] == 8),]))
#10

#For those resolving to Class 0, find the changes in epigenomic marks
new_df_lost <- new_df[(new_df[,4] == 4 & !(new_df[,5] == 4)),c(1:3)]
names(new_df_lost) <- c("chrom", "start", "end")
new_df_lost$coord <- paste(new_df_lost$chrom, new_df_lost$start, new_df_lost$end, sep = "_")

raw_esc <- df_esc_raw[(df_esc_raw$coord %in% new_df_lost$coord),]
raw_npc <- df_npc_raw[(df_npc_raw$coord %in% new_df_lost$coord),]

#Number gained
#For the initial two files, retain only those coordinates classified as as Class 4
new_df <- data.frame(df_esc[,c(1:3)], df_esc[,4], df_npc[,4])
new_df <- new_df[(new_df[,4] == 4 | new_df[,5] == 4),]

print("Total number in ESC:")
print(nrow(new_df[new_df[,4] == 4,]))
#911

print("Total number in neurons:")
print(nrow(new_df[new_df[,5] == 4,]))
#405

print("Total number same:")
print(nrow(new_df[(new_df[,4] == 4 & new_df[,5] == 4),]))
#227

print("Total number lost from ESC to neuron:")
print(nrow(new_df[(new_df[,4] == 4 & !(new_df[,5] == 4)),]))
#684

print("Total number gained from ESC to neuron:")
print(nrow(new_df[(!(new_df[,4] == 4) & new_df[,5] == 4),]))
#178

#For those gained, where did the gained marks come from?
print("Class 0:")
print(nrow(new_df[(new_df[,4] == 0 & new_df[,5] == 4),]))
#2707

print("Class 1:")
print(nrow(new_df[(new_df[,4] == 1 & new_df[,5] == 4),]))
#25

print("Class 2:")
print(nrow(new_df[(new_df[,4] == 2 & new_df[,5] == 4),]))
#0

print("Class 3:")
print(nrow(new_df[(new_df[,4] == 3 & new_df[,5] == 4),]))
#11

print("Class 5:")
print(nrow(new_df[(new_df[,4] == 5 & new_df[,5] == 4),]))
#22

print("Class 6:")
print(nrow(new_df[(new_df[,4] == 6 & new_df[,5] == 4),]))
#1

print("Class 7:")
print(nrow(new_df[(new_df[,4] == 7 & new_df[,5] == 4),]))
#1

print("Class 8:")
print(nrow(new_df[(new_df[,4] == 8 & new_df[,5] == 4),]))
#9

#new_df_lost <- new_df[(new_df[,4] == 4),c(1:3)]
new_df_lost <- new_df[(!(new_df[,4] == 4) & new_df[,5] == 4),c(1:3)]
names(new_df_lost) <- c("chrom", "start", "end")
new_df_lost$coord <- paste(new_df_lost$chrom, new_df_lost$start, new_df_lost$end, sep = "_")

raw_esc <- df_esc_raw[(df_esc_raw$coord %in% new_df_lost$coord),]
raw_npc <- df_npc_raw[(df_npc_raw$coord %in% new_df_lost$coord),]
