bedtools sort -i /home/bettimj/gamazon_rotation/gencode.v38.GRCh38.txt > /home/bettimj/gamazon_rotation/gencode.v38.GRCh38.sorted.bed

#Annotated with gene using bedtools closest
TISSUES=("adipose" "adrenal_gland" "artery" "blood" "breast" "cultured_fibroblast" "ebv_transformed_lymphocyte" "es" "esophagus_muscularis_mucosa" "esophagus_squamous_epithelium" "heart" "intestine" "ips" "kidney" "liver" "lung" "neuron" "ovary" "pancreas" "prostate" "skeletal_muscle" "skin" "spleen" "stomach" "testis" "thyroid" "uterus" "vagina")

for TISSUE in ${TISSUES[@]}; do
    bedtools closest \
    -a /home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/tss_5000_1000/bivalent_promoters/$TISSUE\_bivalent_promoters_5000_1000.bed \
    -b /home/bettimj/gamazon_rotation/gencode.v38.GRCh38.sorted.bed \
    > /home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/tss_5000_1000/bivalent_promoters/genes_$TISSUE\_bivalent_promoters_5000_1000.bed
done