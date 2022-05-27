library("data.table")
library("dplyr")

#Declare the path of the conservation scores file, as well as the paths of the top and bottom enriched bins in neuron
cons_path <- "/home/bettimj/gamazon_rotation/EvolutionaryRate_RegulatoryGenome/data/EvoStats_final2.txt"
binned_genes_dir <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/rcircos/gene_density/ind_tissues_reg_elements_per_gene/ensid"
top_bin_path <- paste(binned_genes_dir, "neuron_top_enriched_regs_per_gene.txt", sep = "/")
bottom_bin_path <- paste(binned_genes_dir, "neuron_bottom_enriched_regs_per_gene.txt", sep = "/")

#Open the files as data frames
cons_file <- fread(cons_path, header = TRUE, sep = "\t", quote = "")
cons_df <- as.data.frame(cons_file)

top_bin_file <- fread(top_bin_path, header = FALSE, sep = "\t", quote = "")
top_bin_df <- as.data.frame(top_bin_file)

bottom_bin_file <- fread(bottom_bin_path, header = FALSE, sep = "\t", quote = "")
bottom_bin_df <- as.data.frame(bottom_bin_file)

cons_scores_top <- merge(top_bin_df, cons_df, by.x = "V1", by.y = "GeneID")
cons_scores_top <- data.frame(cons_scores_top[,1], cons_scores_top$GeneName, cons_scores_top$Chimp_dN.dS, cons_scores_top$Type)

cons_scores_bottom <- merge(bottom_bin_df, cons_df, by.x = "V1", by.y = "GeneID")
cons_scores_bottom <- data.frame(cons_scores_bottom[,1], cons_scores_bottom$GeneName, cons_scores_bottom$Chimp_dN.dS, cons_scores_bottom$Type)

mean(na.omit(cons_scores_top[,3]))
mean(na.omit(cons_scores_bottom[,3]))

wilcox.test(na.omit(cons_scores_top[,3]), na.omit(cons_scores_bottom[,3]))

print(cons_scores_top[,2])
print(cons_scores_bottom[,2])

print(summary(cons_scores_top[,4]))

cons_scores_top %>% group_by(cons_scores_top.Type) %>% summarise(count=n())
cons_scores_bottom %>% group_by(cons_scores_bottom.Type) %>% summarise(count=n())