library("data.table")
library("stringr")

file <- fread("/home/bettimj/gamazon_rotation/gencode.v38lift37.annotation.gtf.gz", header = FALSE, quote = "", sep = "\t")
df <- as.data.frame(file)
new_df <- df[df[,3] == "gene",]

info_cols <- strsplit(new_df[,9], "[;]")

geneid_full <- unlist(sapply(info_cols, "[[", 1))
geneid_full <- strsplit(geneid_full, " ")
geneid_full <- unlist(sapply(geneid_full, "[[", 2))
geneid_full <- noquote(geneid_full)
geneid_full <- str_replace_all(geneid_full, '"', '')

geneid_full <- unlist(sapply(info_cols, "[[", 1))
geneid_full <- strsplit(geneid_full, " ")
geneid_full <- unlist(sapply(geneid_full, "[[", 2))
geneid_full <- noquote(geneid_full)
geneid_full <- str_replace_all(geneid_full, '"', '')

geneid <- strsplit(geneid_full, "[.]")
geneid <- unlist(sapply(geneid, "[[", 1))

genetype <- unlist(sapply(info_cols, "[[", 2))
genetype <- strsplit(genetype, " ")
genetype <- unlist(sapply(genetype, "[[", 3))
genetype <- noquote(genetype)
genetype <- str_replace_all(genetype, '"', '')

genename <- unlist(sapply(info_cols, "[[", 3))
genename <- strsplit(genename, " ")
genename <- unlist(sapply(genename, "[[", 3))
genename <- noquote(genename)
genename <- str_replace_all(genename, '"', '')

#Compile the final output
out_df <- data.frame(new_df[,1], new_df[,4], new_df[,5], new_df[,7], geneid_full, geneid, genetype, genename)
names(out_df) <- c("chr", "left", "right", "strand", "geneid_full", "geneid", "genetype", "genename")

write.table(out_df, file = "gencode.v38.GRCh37.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)