WORKING_DIR=/home/bettimj/gamazon_rotation/core-bed_analysis/core_bed_impute/playground

for FILE in $(ls $WORKING_DIR\/ref_files); do
echo $FILE
ROOT_NAME=${FILE::-11}
gunzip -c $WORKING_DIR\/ref_files/$FILE | tail -n +2 > $WORKING_DIR\/ref_files_lifted_hg38/$ROOT_NAME\hg19.nohead.bed
done

for FILE in $(ls); do
echo $FILE
bgzip $FILE
done