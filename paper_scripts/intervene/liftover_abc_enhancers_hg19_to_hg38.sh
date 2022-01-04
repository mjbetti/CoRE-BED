#Declare the paths of the liftOver directory, as well as the working directory
LIFTOVER_BIN=/home/bettimj/liftOver
CHAIN=/home/bettimj/reference_genomes/hg19ToHg38.over.chain.gz
WORKING_DIR=/home/bettimj/gamazon_rotation/abc_enhancers

#Lift over the bed coordinates from hg19 to hg38
$LIFTOVER_BIN \
$WORKING_DIR\/unique_AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.bed \
$CHAIN \
$WORKING_DIR\/unique_AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.lifted_hg38.bed \
$WORKING_DIR\/unique_AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.unlifted.bed