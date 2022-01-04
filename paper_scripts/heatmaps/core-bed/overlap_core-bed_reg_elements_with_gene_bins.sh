bedtools sort -i 
/home/bettimj/reference_genomes/hg38.genome.1Mb.genecounts.bed > /home/bettimj/reference_genomes/hg38.genome.1Mb.genecounts.sorted.bed

TISSUES=(adipose adrenal_gland artery blood breast cultured_fibroblast ebv_transformed_lymphocyte es esophagus_muscularis_mucosa esophagus_squamous_epithelium heart intestine ips kidney liver lung neuron ovary pancreas prostate skeletal_muscle skin spleen stomach testis thyroid uterus vagina)

for TISSUE in "${TISSUES[@]}"; do
bedtools sort -i /home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/rcircos/gene_density/$TISSUE\_all_regulatory_elements_5000_1000.bed > $TISSUE\_all_regulatory_elements_5000_1000.sorted.bed
done

for TISSUE in "${TISSUES[@]}"; do
bedtools intersect -a /home/bettimj/reference_genomes/hg38.genome.1Mb.genecounts.sorted.bed -b $TISSUE\_all_regulatory_elements_5000_1000.sorted.bed -c -sorted > binned_counts_$TISSUE\_all_regulatory_elements_5000_1000.sorted.bed
done