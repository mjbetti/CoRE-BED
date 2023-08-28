EPIMAP_SAMPLES=("brain_astrocyte")

for SAMPLE in "${EPIMAP_SAMPLES[@]}"; do
python /home/bettimj/gamazon_rotation/mod_core-bed/astrocyte/trackplot/epimap_download_bigwig.py \
-g hg19 \
-t $SAMPLE \
-v \
-r /home/bettimj/gamazon_rotation/mod_core-bed/astrocyte/trackplot
done