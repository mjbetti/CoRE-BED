#Declare the path of the directory containing CoRE-BED enhancer annotations
source_dir <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/tss_5000_1000"

#Loop through each file retain only active ehnancers, and convert each resulting output to UCSC BED format
tissues <- c("adipose", "adrenal_gland", "artery", "blood", "breast", "cultured_fibroblast", "ebv_transformed_lymphocyte", "es", "esophagus_muscularis_mucosa", "esophagus_squamous_epithelium", "heart", "intestine", "ips", "kidney", "liver", "lung", "neuron", "ovary", "pancreas", "prostate", "skeletal_muscle", "skin", "spleen", "stomach", "testis", "thyroid", "uterus", "vagina")

for (tissue in tissues) {
        file_name <- paste0(tissue, "_all_enhancers_5000_1000.bed")
        file <- read.table(paste(source_dir, file_name, sep = "/"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        df <- as.data.frame(file)
        df <- df[(df[,4] == "active_enhancer"),]
        write.table(df, file = paste0(tissue, "_active_enhancers_5000_1000.bed"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

