LIFTOVER_BIN=/home/bettimj/liftOver
CHAIN=/home/bettimj/reference_genomes/hg19ToHg38.over.chain.gz
WORKING_DIR=/home/bettimj/gamazon_rotation/core-bed_analysis/core_bed_impute/playground

for FILE in $(ls $WORKING_DIR\/ref_files_lifted_hg38\/); do
echo $FILE
ROOT_NAME=${FILE::-11}
$LIFTOVER_BIN \
$WORKING_DIR\/ref_files_lifted_hg38/$FILE \
$CHAIN \
$WORKING_DIR\/ref_files_lifted_hg38/$ROOT_NAME\.hg38.bed \
$WORKING_DIR\/ref_files_lifted_hg38/$ROOT_NAME\.unlifted.bed
bgzip $WORKING_DIR\/ref_files_lifted_hg38/$ROOT_NAME\.hg38.bed
done