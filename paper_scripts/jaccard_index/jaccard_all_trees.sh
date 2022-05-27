#The coordinates in each were next sorted and merged:
for i in $(ls /home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/alternative_trees/epimap/jaccard_prep/*.all_annos.bed); do
echo $i
bedtools sort -i $i | mergeBed | uniq > $i.sorted.merged.uniq.bed
done


#...and then computed Jaccard statistic:
SOURCE_DIR=/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/alternative_trees/epimap/jaccard_prep

intervene pairwise \
-i $SOURCE_DIR\/all_regulatory_elements_epimap_1000_5000.txt.all_annos.bed.sorted.merged.uniq.bed \
$SOURCE_DIR\/all_regulatory_elements_epimap_2000_2000.txt.all_annos.bed.sorted.merged.uniq.bed \
$SOURCE_DIR\/all_regulatory_elements_epimap_5000_1000.txt.all_annos.bed.sorted.merged.uniq.bed \
$SOURCE_DIR\/all_regulatory_elements_epimap_5000_1000_no_tss_prom_first.txt.all_annos.bed.sorted.merged.uniq.bed \
$SOURCE_DIR\/all_regulatory_elements_epimap_5000_1000_no_tss_enh_first.txt.all_annos.bed \
/home/bettimj/gamazon_rotation/abc_enhancers/intervene_core_bed/core-bed_impute/merged_inputs/merged_gtexv8_fine-mapped_hc_eqtls.sorted.bed
--names="1k5k","2k2k","5k1k","notss_prom","notss_enh","GTEx_v8_fine-mapped_eQTLs" \
--compute=jaccard \
--htype=number

#1k5k
SOURCE_DIR=/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/alternative_trees/epimap/jaccard_prep

intervene pairwise \
-i $SOURCE_DIR\/all_regulatory_elements_epimap_1000_5000.txt.all_annos.bed.sorted.merged.uniq.bed \
/home/bettimj/gamazon_rotation/abc_enhancers/intervene_core_bed/core-bed_impute/merged_inputs/merged_gtexv8_fine-mapped_hc_eqtls.sorted.bed \
--names="1k5k","GTEx_v8_fine-mapped_eQTLs" \
--compute=jaccard \
--htype=number


# ...results:
#
# 1k5k    GTEx_v8_fine-mapped_eQTLs
# 1k5k    1.0     0.00172321
# GTEx_v8_fine-mapped_eQTLs       0.00172321      1.0


#2k2k
SOURCE_DIR=/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/alternative_trees/epimap/jaccard_prep

intervene pairwise \
-i $SOURCE_DIR\/all_regulatory_elements_epimap_2000_2000.txt.all_annos.bed.sorted.merged.uniq.bed \
/home/bettimj/gamazon_rotation/abc_enhancers/intervene_core_bed/core-bed_impute/merged_inputs/merged_gtexv8_fine-mapped_hc_eqtls.sorted.bed \
--names="2k2k","GTEx_v8_fine-mapped_eQTLs" \
--compute=jaccard \
--htype=number


# ...results:
# 
# 2k2k    GTEx_v8_fine-mapped_eQTLs
# 2k2k    1.0     0.00173808
# GTEx_v8_fine-mapped_eQTLs       0.00173808      1.0


#5k1k
SOURCE_DIR=/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/alternative_trees/epimap/jaccard_prep

intervene pairwise \
-i $SOURCE_DIR\/all_regulatory_elements_epimap_5000_1000.txt.all_annos.bed.sorted.merged.uniq.bed \
/home/bettimj/gamazon_rotation/abc_enhancers/intervene_core_bed/core-bed_impute/merged_inputs/merged_gtexv8_fine-mapped_hc_eqtls.sorted.bed \
--names="5k1k","GTEx_v8_fine-mapped_eQTLs" \
--compute=jaccard \
--htype=number


# ...results:
# 
#         5k1k    GTEx_v8_fine-mapped_eQTLs
# 5k1k    1.0     0.00172423
# GTEx_v8_fine-mapped_eQTLs       0.00172423      1.0


#No TSS promoter
SOURCE_DIR=/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/alternative_trees/epimap/jaccard_prep

intervene pairwise \
-i $SOURCE_DIR\/all_regulatory_elements_epimap_5000_1000_no_tss_prom_first.txt.all_annos.bed.sorted.merged.uniq.bed \
/home/bettimj/gamazon_rotation/abc_enhancers/intervene_core_bed/core-bed_impute/merged_inputs/merged_gtexv8_fine-mapped_hc_eqtls.sorted.bed \
--names="notss_prom","GTEx_v8_fine-mapped_eQTLs" \
--compute=jaccard \
--htype=number


# ...results:
# 
# notss_prom      GTEx_v8_fine-mapped_eQTLs
# notss_prom      1.0     0.00175
# GTEx_v8_fine-mapped_eQTLs       0.00175 1.0


#No TSS enhancer
SOURCE_DIR=/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/alternative_trees/epimap/jaccard_prep

intervene pairwise \
-i $SOURCE_DIR\/all_regulatory_elements_epimap_5000_1000_no_tss_enh_first.txt.all_annos.bed.sorted.merged.uniq.bed \
/home/bettimj/gamazon_rotation/abc_enhancers/intervene_core_bed/core-bed_impute/merged_inputs/merged_gtexv8_fine-mapped_hc_eqtls.sorted.bed \
--names="notss_enh","GTEx_v8_fine-mapped_eQTLs" \
--compute=jaccard \
--htype=number


# ...results:
# 
# notss_enh       GTEx_v8_fine-mapped_eQTLs
# notss_enh       1.0     0.00175
# GTEx_v8_fine-mapped_eQTLs       0.00175 1.0
