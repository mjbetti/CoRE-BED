EPIMAP_SAMPLES=("epithelial_breast_mcf_10a")

for SAMPLE in "${EPIMAP_SAMPLES[@]}"; do
python /gpfs52/data/g_gamazon_lab/bettimj/mod_core-bed/download_bedgraph/epimap_download_bigwig.py \
-g hg19 \
-t $SAMPLE \
-v \
-r /home/bettimj/gamazon_rotation/mod_core-bed/cs2g_snps/trackplot
done