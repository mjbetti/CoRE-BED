bedtools makewindows -g hg19.chrom.sizes -w 1000000 > hg19.genome.1Mb.bed

bedtools sort -i hg19.genome.1Mb.bed > hg19.genome.1Mb.sorted.bed

bedtools intersect -a hg19.genome.1Mb.sorted.bed -b gencode.v38.hg19.annotation.bed -c -sorted > hg19.genome.1Mb.genecounts.bed