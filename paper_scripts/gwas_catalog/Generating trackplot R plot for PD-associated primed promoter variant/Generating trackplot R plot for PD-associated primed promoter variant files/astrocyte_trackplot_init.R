library("trackplot")

#Path to bigWig files
bigWigs = c("brain_astrocyte_H3K4me2_hg19.bigWig", "brain_astrocyte_H3K9ac_hg19.bigWig", "brain_astrocyte_ATAC-seq_hg19.bigWig", "brain_astrocyte_DNase-seq_hg19.bigWig")

#Step-1. Extract the siganl for your loci of interst
track_data = track_extract(bigWigs = bigWigs, loci = "chr1:205,688,378-205,738,378")

#Step-1a (optional). Summarize trcks by condition
track_data = track_summarize(summary_list = track_data, condition = c("H3K4me2", "H3K9ac", "ATAC-seq", "DNase-seq"), stat = "mean")

#Step-2. 
#Basic Plot 
#track_plot(summary_list = track_data)

#With gene models (by default autoamtically queries UCSC genome browser for hg19 transcripts)
#track_plot(summary_list = track_data, draw_gene_track = TRUE, build = "hg19")

#Heighlight regions of interest
#chr1:205713378
markregions = data.frame(
    chr = "chr1",
    start = 205713328,
    end = 205713428,
    name = c("rs823128")
  )

pdf(file = "brain_astrocyte_trackplot.pdf")
track_plot(
  summary_list = track_data,
  draw_gene_track = TRUE,
  show_ideogram = TRUE,
  build = "hg19",
  regions = markregions
)
dev.off()