chr_array=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

for chr in ${chr_array[@]}; do
	sbatch \
    --job-name=chr$chr \
    --account=aldrich_lab \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem=16G \
    --time=2-00:00:00 \
    --wrap="Rscript mod_argin_heart_cS2G_concat_sliced_annos_1kg.R \
    -c $chr"
done