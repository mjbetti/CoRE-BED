bedtools sort -i /home/bettimj/gamazon_rotation/gencode.v38.GRCh37.txt > /home/bettimj/gamazon_rotation/gencode.v38.GRCh38.sorted.bed

#Annotated with gene using bedtools closest:
TISSUES=("adipose" "blood" "bone" "brain" "cancer" "digestive" "endocrine" "endothelial" "epithelial" "esc_deriv" "esc" "eye" "heart" "hsc_and_b_cell" "ipsc" "kidney" "liver" "lung" "lymphoblastoid" "mesench" "muscle" "myosat" "neurosph" "other" "pancreas" "placenta_and_eem" "pns" "reproductive" "sm_muscle" "spleen" "stromal" "thymus" "urinary")

for TISSUE in ${TISSUES[@]}; do
    bedtools sort -i /home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/tss_5000_1000/epimap/bed_format/concat/bivalent_promoters/$TISSUE\_bivalent_promoters_all_regulatory_elements_epimap_5000_1000.sorted.bed > /home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/tss_5000_1000/epimap/bed_format/concat/bivalent_promoters/$TISSUE\_bivalent_promoters_all_regulatory_elements_epimap_5000_1000.sorted.sorted.bed
    bedtools closest \
    -a /home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/tss_5000_1000/epimap/bed_format/concat/bivalent_promoters/$TISSUE\_bivalent_promoters_all_regulatory_elements_epimap_5000_1000.sorted.sorted.bed \
-b /home/bettimj/gamazon_rotation/gencode.v38.GRCh37.sorted.bed \
> /home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/tss_5000_1000/epimap/bed_format/concat/bivalent_promoters/genes_$TISSUE\_bivalent_promoters_5000_1000.bed
done