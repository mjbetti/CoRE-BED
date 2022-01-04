library("data.table")

#Declare the path of the file containing all ABC enhancers
abc_path <- "/home/bettimj/gamazon_rotation/abc_enhancers/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"

#Open the file as a data frame
abc_file <- fread(abc_path, header = TRUE, quote = "", sep = "\t")
abc_df <- as.data.frame(abc_file)

abc_df <- data.frame(abc_df$chr, abc_df$start, abc_df$end)
abc_df <- unique(abc_df)

#Write out the data frame as a bed file
write.table(abc_df, file = "unique_AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)