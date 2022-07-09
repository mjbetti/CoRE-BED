library("data.table")

path <- "/home/bettimj/gamazon_rotation/core-bed_analysis/s2g_snps/exp_val_snps/sliced_annos/all_annos_exp_val_snps.5000_1000.txt"
file <- fread(path, header = TRUE, sep = "\t", quote = "")
df <- as.data.frame(file)

#Split all summary stats based on whether or not each coordinate has a CoRE-BED annotation in any of the CoRE-BED-training or test cell/tissue types
any_anno_df <- df[(df[,17] == 1),]

print("In a promoter:")
print(nrow(any_anno_df[(any_anno_df$any_promoter == 1),]))
print("In an enhancer:")
print(nrow(any_anno_df[(any_anno_df$any_enhancer == 1),]))

print("In an active promoter:")
print(nrow(any_anno_df[(any_anno_df$active_promoter == 1),]))
print("In a bivalent promoter:")
print(nrow(any_anno_df[(any_anno_df$bivalent_promoter == 1),]))
print("In a silenced promoter:")
print(nrow(any_anno_df[(any_anno_df$silenced_promoter == 1),]))

print("In an active enhancer:")
print(nrow(any_anno_df[(any_anno_df$active_enhancer == 1),]))
print("In a poised enhancer:")
print(nrow(any_anno_df[(any_anno_df$poised_enhancer == 1),]))
print("In a primed enhancer:")
print(nrow(any_anno_df[(any_anno_df$primed_enhancer == 1),]))