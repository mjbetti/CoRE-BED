#Declare the paths of the liftOver directory, as well as the working directory
LIFTOVER_BIN=/home/bettimj/liftOver
CHAIN=/home/bettimj/reference_genomes/hg38ToHg19.over.chain.gz
WORKING_DIR=/home/bettimj/gamazon_rotation/gwas_catalog

#Lift over the bed coordinates from hg38 to hg19
$LIFTOVER_BIN \
$WORKING_DIR\/all_associations_v1.0.reformat.case_control.recal.hg38.bed \
$CHAIN \
$WORKING_DIR\/all_associations_v1.0.reformat.case_control.recal.lifted.hg19.bed \
$WORKING_DIR\/all_associations_v1.0.reformat.case_control.recal.unlifted.hg38.bed