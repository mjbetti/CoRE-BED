library("RCircos")
library("data.table")

#Declare the paths of each of the priors files containing the coordinates of all regulatory elements (promoters and enhancers) in each of the 28 CoRE-BED supported tissue types
priors_dir <- "/home/bettimj/gamazon_rotation/core-bed_analysis/bayesian_priors/tss_5000_1000"

#Open each of the files as a data frame (compatible as a BED recognized by RCircos)
tissues <- c("adipose", "adrenal_gland", "artery", "blood", "breast", "cultured_fibroblast", "ebv_transformed_lymphocyte", "es", "esophagus_muscularis_mucosa", "esophagus_squamous_epithelium", "heart", "intestine", "ips", "kidney", "liver", "lung", "neuron", "ovary", "pancreas", "prostate", "skeletal_muscle", "skin", "spleen", "stomach", "testis", "thyroid", "uterus", "vagina")

for (tissue in tissues) {
	path <- paste(priors_dir, paste0(tissue, "_all_regulatory_elements_5000_1000.txt"), sep = "/")
	file <- fread(path, header = FALSE, quote = "", sep = "\t")
	file <- as.data.frame(file)
	valid_chrs <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X", "Y")
	file <- file[(file[,1] %in% valid_chrs),]
	file[,1] <- paste0("chr", file[,1])
	file[,4] <- "red"
	names(file) <- c("Chromosome", "chromStart", "chromEnd", "PlotColor")
	#Initialize the RCircos core components
	##hg38 chromosome ideogram
	data(UCSC.HG38.Human.CytoBandIdeogram)
	cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
	tracks.inside <- 1
	tracks.outside <- 0

	RCircos.Set.Core.Components(cyto.info, chr.exclude = NULL, tracks.inside, tracks.outside)
	rcircos.params <- RCircos.Get.Plot.Parameters()
	#rcircos.params$base.per.unit <- 10

	#Initialize plotting with core components
	out.file <- paste0(tissue, "_core_bed_all_reg_elements_5000_1000_hg38.pdf")
	pdf(file=out.file, height=8, width=8, compress=TRUE)
	RCircos.Set.Plot.Area()
	RCircos.Chromosome.Ideogram.Plot()

	#Plot each of the CoRE-BED tissues as a new tile track
	track.num <- 1
	side <- "in"
	tile.colors <- file[,4];
	file["PlotColor"] <- tile.colors;
	RCircos.Tile.Plot(file, track.num, side)

	dev.off()
}