#Declare the path of the GWAS results downloaded from the GWAS catalog
SCRIPT_DIR=/home/bettimj/CoRE-BED
IN_FILE=/home/bettimj/gamazon_rotation/gwas_catalog/all_associations_v1.0.reformat.case_control.recal.lifted.hg19.txt

python $SCRIPT_DIR\/core-bed_impute.py \
-i $IN_FILE \
-g hg19 \
-t all \
-o all_epimap_tissues_all_associations_v1.0.reformat.case_control.recal.lifted.hg19.func_anno_5000_1000.txt  \
-ud 5000 \
-dd 1000 \
-r /home/bettimj/gamazon_rotation/core-bed_analysis/core_bed_impute/playground/ref_files \
-v \
--no_multianno \
--bed_cols 13,14,14 \
--input_header