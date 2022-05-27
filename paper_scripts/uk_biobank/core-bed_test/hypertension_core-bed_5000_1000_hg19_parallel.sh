#Declare the path of the GWAS results downloaded from the GWAS catalog, as well as each of the 28 tissues
tissue_array=("adipose" "adrenal_gland" "artery" "blood" "breast" "cultured_fibroblast" "ebv_transformed_lymphocyte" "es" "esophagus_muscularis_mucosa" "esophagus_squamous_epithelium" "heart" "intestine" "ips" "kidney" "liver" "lung" "neuron" "ovary" "pancreas" "prostate" "skeletal_muscle" "skin" "spleen" "stomach" "testis" "thyroid" "uterus" "vagina")
SCRIPT_DIR=/home/bettimj/CoRE-BED
IN_FILE=/home/bettimj/gamazon_rotation/uk_biobank/hypertension/phecode-401-both_sexes.tsv

for tissue in ${tissue_array[@]}; do
	sbatch \
	--job-name=$tissue\_ht_5k1k \
	--account=aldrich_lab \
	--nodes=1 \
	--ntasks=1 \
	--cpus-per-task=1 \
	--mem=128G \
	--time=5-00:00:00 \
	--wrap="python $SCRIPT_DIR\/core-bed.py \
	-i $IN_FILE \
	-g hg19 \
	-t $tissue \
	-o phecode-401-both_sexes.func_anno.$tissue\.hg19.5000_1000.tsv  \
	-ud 5000 \
	-dd 1000 \
	-v \
	--no_multianno \
	--bed_cols 1,2,2 \
	--input_header"
done