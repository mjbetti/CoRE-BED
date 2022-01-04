scp phecode-401-both_sexes.func_anno.adipose.hg19.5000_1000.tsv phecode-401-both_sexes.func_anno.all_tissues.hg19.5000_1000.tsv

tissues=("adrenal_gland" "artery" "blood" "breast" "cultured_fibroblast" "ebv_transformed_lymphocyte" "es" "esophagus_muscularis_mucosa" "esophagus_squamous_epithelium" "heart" "intestine" "ips" "kidney" "liver" "lung" "neuron" "ovary" "pancreas" "prostate" "skeletal_muscle" "skin" "spleen" "stomach" "testis" "thyroid" "uterus" "vagina")

for tissue in ${tissues[@]}; do
cut -f47 phecode-401-both_sexes.func_anno.$tissue\.hg19.5000_1000.tsv | paste phecode-401-both_sexes.func_anno.all_tissues.hg19.5000_1000.tsv - | sponge phecode-401-both_sexes.func_anno.all_tissues.hg19.5000_1000.tsv
done