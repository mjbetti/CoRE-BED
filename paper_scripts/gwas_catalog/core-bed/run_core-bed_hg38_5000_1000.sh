#Declare the path of the GWAS results downloaded from the GWAS catalog
SCRIPT_DIR=/home/bettimj/CoRE-BED
IN_FILE=/home/bettimj/gamazon_rotation/gwas_catalog/all_associations_v1.0.reformat.case_control.recal.txt

python $SCRIPT_DIR\/core-bed.py \
-i $IN_FILE \
-g hg38 \
-t all \
-o all_28_tissues_all_associations_v1.0.reformat.case_control.recal.func_anno_5000_1000.txt  \
-ud 5000 \
-dd 1000 \
-v \
--no_multianno \
--bed_cols 12,13,13 \
--input_header