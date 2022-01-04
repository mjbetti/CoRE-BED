bedtools makewindows -g hg38.chrom.sizes -w 1000000 > hg38.genome.1Mb.bed

bedtools sort -i hg38.genome.1Mb.bed > hg38.genome.1Mb.sorted.bed

bedtools intersect -a hg38.genome.1Mb.sorted.bed -b gencode.v38.annotation.bed -c -sorted > hg38.genome.1Mb.genecounts.bed