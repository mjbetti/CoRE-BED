chr_array=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

for chr in ${chr_array[@]}; do
	sbatch \
    --job-name=chr$chr \
    --account=g_gamazon_lab \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem=16G \
    --time=2-00:00:00 \
    --wrap="python /home/bettimj/ldsc/ldsc.py \
			--l2 \
			--bfile /home/bettimj/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr \
			--ld-wind-cm 1 \
			--annot /home/bettimj/gamazon_rotation/mod_core-bed/uk_biobank/h2/anno_1000G/chr$chr\/1000G.EUR.QC.$chr\.annot \
			--thin-annot \
			--out /home/bettimj/gamazon_rotation/mod_core-bed/uk_biobank/h2/anno_1000G/core_bed_EUR.QC.$chr \
			--print-snps /home/bettimj/ldsc/hapmap3_snps/hm.$chr\.snp"
done