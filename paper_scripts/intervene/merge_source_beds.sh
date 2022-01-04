bedtools sort -i /home/bettimj/gamazon_rotation/abc_enhancers/unique_AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.lifted_hg38.bed | bedtools merge > merged_unique_AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.lifted_hg38.bed

bedtools sort -i /home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/tss_5000_1000/all_promoters_enhancers_5000_1000.bed | bedtools merge > merged_all_promoters_enhancers_5000_1000.bed

bedtools sort -i /home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/tss_5000_1000/epimap/all_regulatory_elements_epimap_5000_1000.sorted.uniq.lifted_hg38.bed | bedtools merge > merged_all_regulatory_elements_epimap_5000_1000.sorted.uniq.lifted_hg38.bed

bedtools sort -i /home/bettimj/gamazon_rotation/epimap_chromhmm/hg38/personal.broadinstitute.org/cboix/epimap/ChromHMM/observed_aux_18_hg38/CALLS/all_epimap_enh_18_CALLS_segments_hg38.sorted.bed | bedtools merge > merged_all_epimap_enh_18_CALLS_segments_hg38.sorted.bed

bedtools sort -i /home/bettimj/gamazon_rotation/gtex8_eqtl_sqtl/eqtl/fine_mapped/gtexv8_fine-mapped_hc_eqtls.sorted.bed | bedtools merge > merged_gtexv8_fine-mapped_hc_eqtls.sorted.bed

bedtools sort -i /home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/tss_5000_1000/all_enhancers_5000_1000.bed | bedtools merge > merged_all_enhancers_5000_1000.bed

bedtools sort -i /home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/tss_5000_1000/epimap/all_enhancers_epimap_5000_1000.sorted.uniq.lifted_hg38.bed | bedtools merge > merged_all_enhancers_epimap_5000_1000.sorted.uniq.lifted_hg38.bed