library("data.table")

path <- "/home/bettimj/gamazon_rotation/core-bed_analysis/s2g_snps/exp_val_snps/sliced_annos/all_annos_exp_val_snps.5000_1000.txt"
file <- fread(path, header = TRUE, sep = "\t", quote = "")
df <- as.data.frame(file)

#Split all summary stats based on whether or not each coordinate has a CoRE-BED annotation in any of the CoRE-BED-training or test cell/tissue types
any_anno_df <- df[(df[,17] == 1),]
no_any_anno_df <- df[(df[,17] == 0),]

print("With annotation:")
print(nrow(any_anno_df))
print("No anno:")
print(nrow(no_any_anno_df))
print("Percent any anno:")
print(nrow(any_anno_df)/nrow(df))