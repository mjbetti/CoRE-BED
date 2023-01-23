library("trackplot")

#Path to bigWig files
bigWigs = c("epithelial_breast_mcf_10a_H3K27ac_hg19.bigWig", "epithelial_breast_mcf_10a_H3K27me3_hg19.bigWig", "epithelial_breast_mcf_10a_H3K4me1_hg19.bigWig", "epithelial_breast_mcf_10a_H3K4me2_hg19.bigWig", "epithelial_breast_mcf_10a_H3K4me3_hg19.bigWig", "epithelial_breast_mcf_10a_ATAC-seq_hg19.bigWig", "epithelial_breast_mcf_10a_DNase-seq_hg19.bigWig", "epithelial_breast_mcf_10a_CTCF_hg19.bigWig", "epithelial_breast_mcf_10a_EP300_hg19.bigWig", "epithelial_breast_mcf_10a_H2AFZ_hg19.bigWig", "epithelial_breast_mcf_10a_H3K36me3_hg19.bigWig", "epithelial_breast_mcf_10a_H3K79me2_hg19.bigWig", "epithelial_breast_mcf_10a_H3K9ac_hg19.bigWig", "epithelial_breast_mcf_10a_H3K9me3_hg19.bigWig", "epithelial_breast_mcf_10a_H4K20me1_hg19.bigWig", "epithelial_breast_mcf_10a_POLR2A_hg19.bigWig", "epithelial_breast_mcf_10a_RAD21_hg19.bigWig", "epithelial_breast_mcf_10a_SMC3_hg19.bigWig")

#Step-1. Extract the siganl for your loci of interst
track_data = track_extract(bigWigs = bigWigs, loci = "chr5:56,109,276-56,159,276")

#Step-1a (optional). Summarize trcks by condition
track_data = track_summarize(summary_list = track_data, condition = c("H3K27ac", "H3K27me3", "H3K4me1", "H3K4me2", "H3K4me3", "ATAC-seq", "DNase-seq", "CTCF", "EP300", "H2AFZ", "H3K36me3", "H3K79me2", "H3K9ac", "H3K9me3", "H4K20me1", "POLR2A", "RAD21", "SMC3"), stat = "mean")

#Step-2. 
#Basic Plot 
#track_plot(summary_list = track_data)

#With gene models (by default autoamtically queries UCSC genome browser for hg19 transcripts)
#track_plot(summary_list = track_data, draw_gene_track = TRUE, build = "hg19")

#Heighlight regions of interest
#chr3:56134276
markregions = data.frame(
    chr = "chr5",
    start = 56130776,
    end = 56137776,
    name = c("rs16886397")
  )

pdf(file = "cs2g_breast_cancer_trackplot.pdf")
track_plot(
  summary_list = track_data,
  draw_gene_track = TRUE,
  show_ideogram = TRUE,
  build = "hg19",
  regions = markregions
)
dev.off()

##########################

library("trackplot")

#Path to bigWig files
bigWigs = c("epithelial_breast_mcf_10a_ATAC-seq_hg19.bigWig", "epithelial_breast_mcf_10a_DNase-seq_hg19.bigWig", "epithelial_breast_mcf_10a_CTCF_hg19.bigWig", "epithelial_breast_mcf_10a_RAD21_hg19.bigWig", "epithelial_breast_mcf_10a_SMC3_hg19.bigWig")

#Step-1. Extract the siganl for your loci of interst
track_data = track_extract(bigWigs = bigWigs, loci = "chr5:56,109,276-56,159,276")

#Step-1a (optional). Summarize trcks by condition
track_data = track_summarize(summary_list = track_data, condition = c("ATAC-seq", "DNase-seq", "CTCF", "RAD21", "SMC3"), stat = "mean")

#Step-2. 
#Basic Plot 
#track_plot(summary_list = track_data)

#With gene models (by default autoamtically queries UCSC genome browser for hg19 transcripts)
#track_plot(summary_list = track_data, draw_gene_track = TRUE, build = "hg19")

#Heighlight regions of interest
#chr3:56134276
markregions = data.frame(
    chr = "chr5",
    start = 56130726,
    end = 56137826,
    name = c("rs16886397")
  )

pdf(file = "cs2g_breast_cancer_trackplot.pdf")
track_plot(
  summary_list = track_data,
  draw_gene_track = TRUE,
  show_ideogram = TRUE,
  build = "hg19",
  regions = markregions
)
dev.off()