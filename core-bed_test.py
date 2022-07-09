<<<<<<< HEAD
###Developed by Michael J Betti, April 2021, updated 8 July 2022
__author__ = "Michael J Betti"
__license__ = "MIT"
=======
###Developed by Michael J Betti, April 2021, updated 25 February 2022
__author__ = "Michael J Betti"
__copyright__ = "Copyright 2021, Michael J Betti"
__license__ = "BSD"
>>>>>>> 3283f061327aed512684aab1c4b30a7a1c01e0ad
__maintainer__ = "Michael J Betti"
__email__ = "mjbetti3@gmail.com"
__status__ = "Development"

import os, sys, requests, argparse, pybedtools, pandas as pd
from itertools import chain

parser = argparse.ArgumentParser(add_help = True)
parser.add_argument("-i", "--input", type = str, required = True, help = "the input bed file (required)")
<<<<<<< HEAD
parser.add_argument("-s", "--separator", type = str, required = False, help = "the upstream boundary distance from a TSS (default: 5000 bp)", default = "\t")
=======
>>>>>>> 3283f061327aed512684aab1c4b30a7a1c01e0ad
parser.add_argument("-g", "--ref_genome", type = str, required = True, help = "the human or mouse reference genome build on which the input coordinates are based (required) (valid options: GRCh38/hg38, GRCh37/hg19, GRCm39/mm39, GRCm38/mm10, or GRCm37/mm9)")
parser.add_argument("-t", "--tissue", type = str, required = True, help = "the tissue of interest (required) (valid human options: Adipose, Adrenal_gland, Artery, Blood, Breast, Cultured_fibroblast, EBV_transformed_lymphocyte, ES, Esophagus_muscularis_mucosa, Esophagus_squamous_epithelium, Heart, Intestine, iPS, Kidney, Liver, Lung, Neuron, Ovary, Pancreas, Prostate, Skeletal_muscle, Skin, Spleen, Stomach, Testis, Thyroid, Uterus, Vagina, All, User_provided_files, User_provided_urls; valid mouse options: User_provided_files, User_provided_urls)")
parser.add_argument("-ud", "--tss_distance_upstream", type = int, required = False, help = "the upstream boundary distance from a TSS (default: 5000 bp)", default = 5000)
parser.add_argument("-dd", "--tss_distance_downstream", type = int, required = False, help = "the downstream boundary distance from a TSS (default: 1000 bp)", default = 1000)
parser.add_argument("-o", "--output", type = str, required = False, help = "the name of the output file", default = "out.bed")
parser.add_argument("-r", "--ref_dir", type = str, required = False, help = "the path of the reference file directory", default = "ref_files")
parser.add_argument("--no_multianno", required = False, help = "if a coordinate overlaps with multiple regions, keep the most significant occurance", action = "store_true")
<<<<<<< HEAD
parser.add_argument("--write_summary", required = False, help = "Write out a summary of regulatory counts as a .txt file", action = "store_true")
parser.add_argument("--write_anno_only", required = False, help = "Instead of the input file appended with an annotation column, write out only the annotation column", action = "store_true")
=======
>>>>>>> 3283f061327aed512684aab1c4b30a7a1c01e0ad
parser.add_argument("--bed_cols", type = str, required = False, help = "if the input is not in traditional UCSC BED format, specify the column numbers of chr, start, and end separated by commas", default = "1,2,3")
parser.add_argument("--input_header", required = False, help = "use if the input file has a header", action = "store_true")
parser.add_argument("--user_4me1", type = str, required = False, help = "if the User_provided_files or User_provided_urls tissue option is specified, specify either the path or URL of the user-provided H3K4me1 ChIP-seq peaks")
parser.add_argument("--user_4me3", type = str, required = False, help = "if the User_provided_files or User_provided_urls tissue option is specified, specify either the path or URL of the user-provided H3K4me3 ChIP-seq peaks")
parser.add_argument("--user_27ac", type = str, required = False, help = "if the User_provided_files or User_provided_urls tissue option is specified, specify either the path or URL of the user-provided H3K27ac ChIP-seq peaks")
parser.add_argument("--user_27me3", type = str, required = False, help = "if the User_provided_files or User_provided_urls tissue option is specified, specify either the path or URL of the user-provided H3K27me3 ChIP-seq peaks")
#parser.add_argument("--user_36me3", type = str, required = False, help = "if the User_provided_files or User_provided_urls tissue option is specified, specify either the path or URL of the user-provided H3K36me3 ChIP-seq peaks")
parser.add_argument("--user_ctcf", type = str, required = False, help = "if the User_provided_files or User_provided_urls tissue option is specified, specify either the path or URL of the user-provided CTCF ChIP-seq peaks")
parser.add_argument("--user_dnase", type = str, required = False, help = "if the User_provided_files or User_provided_urls tissue option is specified, specify either the path or URL of the user-provided DNase-seq peaks")
parser.add_argument("--user_tissue_names", type = str, required = False, help = "if the User_provided_files or User_provided_urls tissue option is specified, specify the names of each corresponding tissue type")
parser.add_argument("-v", "--verbose", required = False, help = "return logging as terminal output", action = "store_true")
<<<<<<< HEAD
parser.add_argument("-h3k4me1", "--only_h3k4me1", required = False, help = "Instead of full CoRE-BED regulatory annotations, only evaluate input for H3K4me1 overlap in the target tissue", action = "store_true")
parser.add_argument("-h3k4me3", "--only_h3k4me3", required = False, help = "Instead of full CoRE-BED regulatory annotations, only evaluate input for H3K4me3 overlap in the target tissue", action = "store_true")
parser.add_argument("-h3k27ac", "--only_h3k27ac", required = False, help = "Instead of full CoRE-BED regulatory annotations, only evaluate input for H3K27ac overlap in the target tissue", action = "store_true")
parser.add_argument("-h3k27me3", "--only_h3k27me3", required = False, help = "Instead of full CoRE-BED regulatory annotations, only evaluate input for H3K27me3 overlap in the target tissue", action = "store_true")
parser.add_argument("-dhs", "--only_dhs", required = False, help = "Instead of full CoRE-BED regulatory annotations, only evaluate input for DHS overlap in the target tissue", action = "store_true")
parser.add_argument("-ctcf", "--only_ctcf", required = False, help = "Instead of full CoRE-BED regulatory annotations, only evaluate input for CTCF overlap in the target tissue", action = "store_true")
=======
>>>>>>> 3283f061327aed512684aab1c4b30a7a1c01e0ad
args = parser.parse_args()

#Check that required arguments are specified
assert args.input, "Must specify input file (-i, --input)"
assert args.ref_genome, "Must specify reference genome build (-g, --ref_genome)"
assert args.tissue, "Must specify tissue type (-t, --tissue)"
if args.tissue.lower() == "user_provided_files" or args.tissue.lower() == "user_provided_urls":
	assert args.user_4me1 and args.user_4me3 and args.user_27ac and args.user_27me3 and args.user_ctcf and args.user_dnase and args.user_tissue_names, "Must provide histone ChIP-seq, CTCF ChIP-seq, and DNase-seq files when using the User_provided_files tissue option or URLs when using the User_provided_urls option. Must also provide the tissue/cell types names of each specified set of files."

#If the reference genome is human, check that specified tissue type is one of the 30 valid options
###Create an array with all of the parsed in tissue types
tissue_array = args.tissue.lower().split(",")
all_array = ["adipose", "adrenal_gland", "artery", "blood", "breast", "cultured_fibroblast", "ebv_transformed_lymphocyte", "es", "esophagus_muscularis_mucosa", "esophagus_squamous_epithelium", "heart", "intestine", "ips", "kidney", "liver", "lung", "neuron", "ovary", "pancreas", "prostate", "skeletal_muscle", "skin", "spleen", "stomach", "testis", "thyroid", "uterus", "vagina"]
tissue_array = list(chain.from_iterable(all_array if item == "all" else [item] for item in tissue_array))

###Create an array with all of the parsed in epigenetic arguments
if ("user_provided_files" in tissue_array or "user_provided_urls" in tissue_array):
	user_4me1_array = args.user_4me1.split(",")
	user_4me3_array = args.user_4me3.split(",")
	user_27ac_array = args.user_27ac.split(",")
	user_27me3_array = args.user_27me3.split(",")
	#user_36me3_array = args.user_36me3.split(",")
	user_ctcf_array = args.user_ctcf.split(",")
	user_dnase_array = args.user_dnase.split(",")
	user_tissue_names_array = args.user_tissue_names.lower().split(",")
	user_files_index = 0

if args.ref_genome.lower() == "hg38" or args.ref_genome.lower() == "grch38" or args.ref_genome.lower() == "hg19" or args.ref_genome.lower() == "grch37":
	for tissue in tissue_array:
		assert tissue == "adipose" or tissue == "adrenal_gland" or tissue == "artery" or tissue == "blood" or tissue == "breast" or tissue == "cultured_fibroblast" or tissue == "ebv_transformed_lymphocyte" or tissue == "es" or tissue == "esophagus_muscularis_mucosa" or tissue == "esophagus_squamous_epithelium" or tissue == "heart" or tissue == "intestine" or tissue == "ips" or tissue == "kidney" or tissue == "liver" or tissue == "lung" or tissue == "neuron" or tissue == "ovary" or tissue == "pancreas" or tissue == "prostate" or tissue == "skeletal_muscle" or tissue == "skin" or tissue == "spleen" or tissue == "stomach" or tissue == "testis" or tissue == "thyroid" or tissue == "uterus" or tissue == "vagina" or tissue == "all" or tissue == "user_provided_files" or tissue == "user_provided_urls", "Tissue type must be one of the 31 valid human options (Adipose, Adrenal_gland, Artery, Blood, Breast, Cultured_fibroblast, EBV_transformed_lymphocyte, ES, Esophagus_muscularis_mucosa, Esophagus_squamous_epithelium, Heart, Intestine, iPS, Kidney, Liver, Lung, Neuron, Ovary, Pancreas, Prostate, Skeletal_muscle, Skin, Spleen, Stomach, Testis, Thyroid, Uterus, Vagina, All, User_provided_files, or User_provided_urls)"
	
elif args.ref_genome.lower() == "mm39" or args.ref_genome.lower() == "grcm39" or args.ref_genome.lower() == "mm10" or args.ref_genome.lower() == "grcm38" or args.ref_genome.lower() == "mm9" or args.ref_genome.lower() == "grcm37":
	assert args.tissue.lower() == "user_provided_files" or args.tissue.lower() == "user_provided_urls", "Tissue type must be one of the 2 valid mouse options (User_provided_files or User_provided_urls)"

#Download the appropriate reference files based on the specified genome build and tissue arguments
<<<<<<< HEAD
if not os.path.exists(args.ref_dir):
	os.mkdir(args.ref_dir)

=======
>>>>>>> 3283f061327aed512684aab1c4b30a7a1c01e0ad
def download_ref(url_ref, out_name):
	out_path = os.path.join(args.ref_dir, out_name)
	r = requests.get(url_ref, allow_redirects = True)
	open(out_path, 'wb').write(r.content)

if args.verbose:
	print("Classifying {tissue_type} cis-regulatory regions in {input_file}...".format(input_file = args.input, tissue_type = args.tissue))
	print("\n")
if not os.path.exists(args.ref_dir):
	os.mkdir(args.ref_dir)

if args.verbose:
	print("Downloading {ref_genome} TSS coordinates and {tissue_type} histone ChIP-seq bed files...".format(ref_genome = args.ref_genome, tissue_type = args.tissue))
	print("\n")

###Human###
#GRCh38/hg38
if args.ref_genome.lower() == "hg38" or args.ref_genome.lower() == "grch38":
	if not os.path.exists(args.ref_dir + "/refTSS_v3.1_human_coordinate.hg38.bed.gz"):
		download_ref("http://reftss.clst.riken.jp/datafiles/3.1/human/refTSS_v3.1_human_coordinate.hg38.bed.gz", "refTSS_v3.1_human_coordinate.hg38.bed.gz")
	
	for tissue in tissue_array:
		#Adipose (Homo sapiens subcutaneous abdominal adipose tissue tissue nuclear fraction female adult (49 years), Homo sapiens omental fat pad tissue female adult (51 years) (DNase-seq and CTCF ChIP-seq only))
		if tissue == "adipose":
			if not os.path.exists(args.ref_dir + "/adipose_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF534VNL/@@download/ENCFF534VNL.bed.gz", "adipose_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/adipose_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF591SLF/@@download/ENCFF591SLF.bed.gz", "adipose_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/adipose_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF658NHX/@@download/ENCFF658NHX.bed.gz", "adipose_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/adipose_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF435MSC/@@download/ENCFF435MSC.bed.gz", "adipose_27me3_hg38.bed.gz")
			#if not os.path.exists("ref_files/adipose_36me3_hg38.bed.gz"):	
				#download_ref("https://www.encodeproject.org/files/ENCFF075XPG/@@download/ENCFF075XPG.bed.gz", "adipose_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/adipose_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF464DDN/@@download/ENCFF464DDN.bed.gz", "adipose_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/adipose_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF436YWU/@@download/ENCFF436YWU.bed.gz", "adipose_dnase_hg38.bed.gz")
		
		#Adrenal gland (Homo sapiens adrenal gland tissue male embryo (97 days), Homo sapiens adrenal gland tissue embryo (96 days) (DNase-seq only), and Homo sapiens adrenal gland tissue male adult (37 years) (CTCF ChIP-seq only))
		if tissue == "adrenal_gland":
			if not os.path.exists(args.ref_dir + "/adrenal_gland_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF953KZN/@@download/ENCFF953KZN.bed.gz", "adrenal_gland_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/adrenal_gland_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF959FFU/@@download/ENCFF959FFU.bed.gz", "adrenal_gland_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/adrenal_gland_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF231NZU/@@download/ENCFF231NZU.bed.gz", "adrenal_gland_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/adrenal_gland_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF852DFJ/@@download/ENCFF852DFJ.bed.gz", "adrenal_gland_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/adrenal_gland_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF852DFJ/@@download/ENCFF852DFJ.bed.gz", "adrenal_gland_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/adrenal_gland_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF297BQF/@@download/ENCFF297BQF.bed.gz", "adrenal_gland_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/adrenal_gland_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF949ZRW/@@download/ENCFF949ZRW.bed.gz", "adrenal_gland_dnase_hg38.bed.gz")
			
		#Artery (Homo sapiens aorta tissue male adult (34 years), Homo sapiens aorta tissue female adult (41 years) (DNase-seq only), and Homo sapiens ascending aorta tissue female adult (51 years) (CTCF ChIP-seq only))
		if tissue == "artery":
			if not os.path.exists(args.ref_dir + "/artery_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF824OFK/@@download/ENCFF824OFK.bed.gz", "artery_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/artery_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF876LPD/@@download/ENCFF876LPD.bed.gz", "artery_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/artery_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF872XZM/@@download/ENCFF872XZM.bed.gz", "artery_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/artery_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF093PYF/@@download/ENCFF093PYF.bed.gz", "artery_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/artery_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF587NCS/@@download/ENCFF587NCS.bed.gz", "artery_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/artery_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF218JYZ/@@download/ENCFF218JYZ.bed.gz", "artery_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/artery_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF829CJM/@@download/ENCFF829CJM.bed.gz", "artery_dnase_hg38.bed.gz")

		#Blood (Homo sapiens peripheral blood mononuclear cell male adult (39 years) and Homo sapiens K562 (DNase-seq and CTCF ChIP-seq only))	
		if tissue == "blood":
			if not os.path.exists(args.ref_dir + "/blood_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF734KLT/@@download/ENCFF734KLT.bed.gz", "blood_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/blood_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF165VDC/@@download/ENCFF165VDC.bed.gz", "blood_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/blood_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF844OYX/@@download/ENCFF844OYX.bed.gz", "blood_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/blood_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF636JEP/@@download/ENCFF636JEP.bed.gz", "blood_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/blood_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF422UFO/@@download/ENCFF422UFO.bed.gz", "blood_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/blood_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF582SNT/@@download/ENCFF582SNT.bed.gz", "blood_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/blood_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF274YGF/@@download/ENCFF274YGF.bed.gz", "blood_dnase_hg38.bed.gz")
	
		#Breast (Homo sapiens breast epithelium tissue female adult (53 years) and Homo sapiens breast epithelium tissue female adult (51 years) (DNase-seq only))
		elif tissue == "breast":
			if not os.path.exists(args.ref_dir + "/breast_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF698LQG/@@download/ENCFF698LQG.bed.gz", "breast_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/breast_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF590DGH/@@download/ENCFF590DGH.bed.gz", "breast_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/breast_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF507WEL/@@download/ENCFF507WEL.bed.gz", "breast_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/breast_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF962YZN/@@download/ENCFF962YZN.bed.gz", "breast_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/breast_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF540RWM/@@download/ENCFF540RWM.bed.gz", "breast_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/breast_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF910ACN/@@download/ENCFF910ACN.bed.gz", "breast_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/breast_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF594NFE/@@download/ENCFF594NFE.bed.gz", "breast_dnase_hg38.bed.gz")
		
		#Cultured fibroblast (Homo sapiens IMR-90)
		elif tissue == "cultured_fibroblast":
			if not os.path.exists(args.ref_dir + "/cultured_fibroblast_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF611UWF/@@download/ENCFF611UWF.bed.gz", "cultured_fibroblast_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/cultured_fibroblast_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF093NQC/@@download/ENCFF093NQC.bed.gz", "cultured_fibroblast_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/cultured_fibroblast_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF805GNH/@@download/ENCFF805GNH.bed.gz", "cultured_fibroblast_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/cultured_fibroblast_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF336IXL/@@download/ENCFF336IXL.bed.gz", "cultured_fibroblast_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/cultured_fibroblast_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF449ADN/@@download/ENCFF449ADN.bed.gz", "cultured_fibroblast_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/cultured_fibroblast_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF203SRF/@@download/ENCFF203SRF.bed.gz", "cultured_fibroblast_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/cultured_fibroblast_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF800DVI/@@download/ENCFF800DVI.bed.gz", "cultured_fibroblast_dnase_hg38.bed.gz")
		
		#EBV-transformed lymphocyte (Homo sapiens GM12878)
		elif tissue == "ebv_transformed_lymphocyte":
			if not os.path.exists(args.ref_dir + "/ebv_transformed_lymphocyte_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF321BVG/@@download/ENCFF321BVG.bed.gz", "ebv_transformed_lymphocyte_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/ebv_transformed_lymphocyte_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF998CEU/@@download/ENCFF998CEU.bed.gz", "ebv_transformed_lymphocyte_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/ebv_transformed_lymphocyte_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF023LTU/@@download/ENCFF023LTU.bed.gz", "ebv_transformed_lymphocyte_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/ebv_transformed_lymphocyte_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF291DHI/@@download/ENCFF291DHI.bed.gz", "ebv_transformed_lymphocyte_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/ebv_transformed_lymphocyte_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF432EMI/@@download/ENCFF432EMI.bed.gz", "ebv_transformed_lymphocyte_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/ebv_transformed_lymphocyte_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF796WRU/@@download/ENCFF796WRU.bed.gz", "ebv_transformed_lymphocyte_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/ebv_transformed_lymphocyte_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF759OLD/@@download/ENCFF759OLD.bed.gz", "ebv_transformed_lymphocyte_dnase_hg38.bed.gz")
		
		#Embryonic stem cell (Homo sapiens H1)
		elif tissue == "es":
			if not os.path.exists(args.ref_dir + "/es_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF613QAB/@@download/ENCFF613QAB.bed.gz", "es_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/es_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF356MCC/@@download/ENCFF356MCC.bed.gz", "es_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/es_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF689CJG/@@download/ENCFF689CJG.bed.gz", "es_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/es_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF050CUG/@@download/ENCFF050CUG.bed.gz", "es_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/es_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF813VFV/@@download/ENCFF813VFV.bed.gz", "es_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/es_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF692RPA/@@download/ENCFF692RPA.bed.gz", "es_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/es_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF983UCL/@@download/ENCFF983UCL.bed.gz", "es_dnase_hg38.bed.gz")
	
		#Esophagus muscularis mucosa (Homo sapiens esophagus muscularis mucosa tissue female adult (51 years) and Homo sapiens esophagus muscularis mucosa tissue male adult (37 years) (DNase-seq only))
		elif tissue == "esophagus_muscularis_mucosa":
			if not os.path.exists(args.ref_dir + "/esophagus_muscularis_mucosa_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF473ELV/@@download/ENCFF473ELV.bed.gz", "esophagus_muscularis_mucosa_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/esophagus_muscularis_mucosa_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF154TIT/@@download/ENCFF154TIT.bed.gz", "esophagus_muscularis_mucosa_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/esophagus_muscularis_mucosa_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF114SOB/@@download/ENCFF114SOB.bed.gz", "esophagus_muscularis_mucosa_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/esophagus_muscularis_mucosa_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF450WVA/@@download/ENCFF450WVA.bed.gz", "esophagus_muscularis_mucosa_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/esophagus_muscularis_mucosa_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF101WAE/@@download/ENCFF101WAE.bed.gz", "esophagus_muscularis_mucosa_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/esophagus_muscularis_mucosa_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF851XPN/@@download/ENCFF851XPN.bed.gz", "esophagus_muscularis_mucosa_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/esophagus_muscularis_mucosa_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF726TPF/@@download/ENCFF726TPF.bed.gz", "esophagus_muscularis_mucosa_dnase_hg38.bed.gz")
		
		#Esophagus squamous epithelium (Homo sapiens esophagus squamous epithelium tissue female adult (51 years) and Homo sapiens esophagus squamous epithelium tissue male adult (37 years) (DNase-seq only))
		elif tissue == "esophagus_squamous_epithelium":
			if not os.path.exists(args.ref_dir + "/esophagus_squamous_epithelium_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF828YQS/@@download/ENCFF828YQS.bed.gz", "esophagus_squamous_epithelium_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/esophagus_squamous_epithelium_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF470ELM/@@download/ENCFF470ELM.bed.gz", "esophagus_squamous_epithelium_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/esophagus_squamous_epithelium_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF673PCP/@@download/ENCFF673PCP.bed.gz", "esophagus_squamous_epithelium_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/esophagus_squamous_epithelium_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF905MKP/@@download/ENCFF905MKP.bed.gz", "esophagus_squamous_epithelium_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/esophagus_squamous_epithelium_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF885BKA/@@download/ENCFF885BKA.bed.gz", "esophagus_squamous_epithelium_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/esophagus_squamous_epithelium_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF384CCC/@@download/ENCFF384CCC.bed.gz", "esophagus_squamous_epithelium_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/esophagus_squamous_epithelium_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF839REP/@@download/ENCFF839REP.bed.gz", "esophagus_squamous_epithelium_dnase_hg38.bed.gz")
		
		#Heart (Homo sapiens heart left ventricle tissue female adult (53 years))
		elif tissue == "heart":
			if not os.path.exists(args.ref_dir + "/heart_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF194KSV/@@download/ENCFF194KSV.bed.gz", "heart_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/heart_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF008ODN/@@download/ENCFF008ODN.bed.gz", "heart_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/heart_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF635VBV/@@download/ENCFF635VBV.bed.gz", "heart_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/heart_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF242JDO/@@download/ENCFF242JDO.bed.gz", "heart_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/heart_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF767VOC/@@download/ENCFF767VOC.bed.gz", "heart_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/heart_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF165RCX/@@download/ENCFF165RCX.bed.gz", "heart_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/heart_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF846JLR/@@download/ENCFF846JLR.bed.gz", "heart_dnase_hg38.bed.gz")
		
		#Intestine (Homo sapiens sigmoid colon tissue male adult (37 years) and Homo sapiens sigmoid colon tissue female adult (53 years) (DNase-seq only))
		elif tissue == "intestine":
			if not os.path.exists(args.ref_dir + "/intestine_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF052XGA/@@download/ENCFF052XGA.bed.gz", "intestine_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/intestine_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF250ZHO/@@download/ENCFF250ZHO.bed.gz", "intestine_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/intestine_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF111SYA/@@download/ENCFF111SYA.bed.gz", "intestine_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/intestine_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF310IST/@@download/ENCFF310IST.bed.gz", "intestine_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/intestine_27me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF406OXE/@@download/ENCFF406OXE.bed.gz", "intestine_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/intestine_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF480QYI/@@download/ENCFF480QYI.bed.gz", "intestine_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/intestine_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF274XSP/@@download/ENCFF274XSP.bed.gz", "intestine_dnase_hg38.bed.gz")
	
		#iPS cell (Homo sapiens iPS DF 19.11 and Homo sapiens GM23338 originated from GM23248 (CTCF ChIP-seq only))
		elif tissue == "ips":
			if not os.path.exists(args.ref_dir + "/ips_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF314LSR/@@download/ENCFF314LSR.bed.gz", "ips_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/ips_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF450VVJ/@@download/ENCFF450VVJ.bed.gz", "ips_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/ips_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF037RXP/@@download/ENCFF037RXP.bed.gz", "ips_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/ips_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF844FEB/@@download/ENCFF844FEB.bed.gz", "ips_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/ips_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF855HXT/@@download/ENCFF855HXT.bed.gz", "ips_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/ips_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF651DQL/@@download/ENCFF651DQL.bed.gz", "ips_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/ips_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF172GXY/@@download/ENCFF172GXY.bed.gz", "ips_dnase_hg38.bed.gz")
		
		#Kidney (Homo sapiens kidney tissue male adult (67 years), Homo sapiens kidney tissue female embryo (120 days) for H3K27me3 - only dataset for kidney and DNase-seq, and Homo sapiens kidney epithelial cell (CTCF ChIP-seq))
		elif tissue == "kidney":
			if not os.path.exists(args.ref_dir + "/kidney_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF249OJK/@@download/ENCFF249OJK.bed.gz", "kidney_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/kidney_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF430ZMN/@@download/ENCFF430ZMN.bed.gz", "kidney_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/kidney_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF512JNT/@@download/ENCFF512JNT.bed.gz", "kidney_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/kidney_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF775YUI/@@download/ENCFF775YUI.bed.gz", "kidney_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/kidney_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF310IRC/@@download/ENCFF310IRC.bed.gz", "kidney_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/kidney_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF059SDM/@@download/ENCFF059SDM.bed.gz", "kidney_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/kidney_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF388NNQ/@@download/ENCFF388NNQ.bed.gz", "kidney_dnase_hg38.bed.gz")

		#Liver (Homo sapiens liver tissue male adult (32 years), Homo sapiens liver tissue female embryo (101 days) and female embryo (113 days) (DNase-seq only), and Homo sapiens right lobe of liver tissue female adult (53 years) (CTCF ChIP-seq only))
		elif tissue == "liver":
			if not os.path.exists(args.ref_dir + "/liver_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF625TWQ/@@download/ENCFF625TWQ.bed.gz", "liver_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/liver_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF830ODB/@@download/ENCFF830ODB.bed.gz", "liver_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/liver_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF287VIA/@@download/ENCFF287VIA.bed.gz", "liver_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/liver_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF287VIA/@@download/ENCFF287VIA.bed.gz", "liver_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/liver_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF666OGH/@@download/ENCFF666OGH.bed.gz", "liver_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/liver_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF046XEZ/@@download/ENCFF046XEZ.bed.gz", "liver_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/liver_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF933OGA/@@download/ENCFF933OGA.bed.gz", "liver_dnase_hg38.bed.gz")
		
		#Lung (Homo sapiens upper lobe of left lung tissue female adult (51 years))
		elif tissue == "lung":
			if not os.path.exists(args.ref_dir + "/lung_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF926YKP/@@download/ENCFF926YKP.bed.gz", "lung_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/lung_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF843JPR/@@download/ENCFF843JPR.bed.gz", "lung_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/lung_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF496RDA/@@download/ENCFF496RDA.bed.gz", "lung_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/lung_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF248IPH/@@download/ENCFF248IPH.bed.gz", "lung_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/lung_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF459IYH/@@download/ENCFF459IYH.bed.gz", "lung_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/lung_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF629UDP/@@download/ENCFF629UDP.bed.gz", "lung_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/lung_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF041XPT/@@download/ENCFF041XPT.bed.gz", "lung_dnase_hg38.bed.gz")
		
		#Neuron (Homo sapiens SK-N-SH)
		elif tissue == "neuron":
			if not os.path.exists(args.ref_dir + "/neuron_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF632FAM/@@download/ENCFF632FAM.bed.gz", "neuron_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/neuron_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF682JYE/@@download/ENCFF682JYE.bed.gz", "neuron_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/neuron_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF138VUT/@@download/ENCFF138VUT.bed.gz", "neuron_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/neuron_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF277NRX/@@download/ENCFF277NRX.bed.gz", "neuron_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/neuron_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF140RVM/@@download/ENCFF140RVM.bed.gz", "neuron_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/neuron_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF608JKY/@@download/ENCFF608JKY.bed.gz", "neuron_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/neuron_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF752OZB/@@download/ENCFF752OZB.bed.gz", "neuron_dnase_hg38.bed.gz")
		
		#Ovary (Homo sapiens ovary tissue female adult (30 years) and Homo sapiens ovary tissue female adult (53 years) (CTCF ChIP-seq only))
		elif tissue == "ovary":
			if not os.path.exists(args.ref_dir + "/ovary_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF458XZT/@@download/ENCFF458XZT.bed.gz", "ovary_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/ovary_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF008MLA/@@download/ENCFF008MLA.bed.gz", "ovary_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/ovary_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF328KKO/@@download/ENCFF328KKO.bed.gz", "ovary_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/ovary_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF495DAT/@@download/ENCFF495DAT.bed.gz", "ovary_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/ovary_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF558KMA/@@download/ENCFF558KMA.bed.gz", "ovary_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/ovary_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF782GIY/@@download/ENCFF782GIY.bed.gz", "ovary_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/ovary_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF418PRB/@@download/ENCFF418PRB.bed.gz", "ovary_dnase_hg38.bed.gz")
		
		#Pancreas (Homo sapiens pancreas tissue female adult (41 years))
		elif tissue == "pancreas":
			if not os.path.exists(args.ref_dir + "/pancreas_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF453GYD/@@download/ENCFF453GYD.bed.gz", "pancreas_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/pancreas_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF200ITL/@@download/ENCFF200ITL.bed.gz", "pancreas_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/pancreas_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF498SSN/@@download/ENCFF498SSN.bed.gz", "pancreas_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/pancreas_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF663DVO/@@download/ENCFF663DVO.bed.gz", "pancreas_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/pancreas_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF456HUQ/@@download/ENCFF456HUQ.bed.gz", "pancreas_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/pancreas_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF682RNL/@@download/ENCFF682RNL.bed.gz", "pancreas_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/pancreas_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF949DIO/@@download/ENCFF949DIO.bed.gz", "pancreas_dnase_hg38.bed.gz")
		
		#Prostate (Homo sapiens prostate gland tissue male adult (37 years))
		elif tissue == "prostate":
			if not os.path.exists(args.ref_dir + "/prostate_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF275KSR/@@download/ENCFF275KSR.bed.gz", "prostate_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/prostate_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF155ZYZ/@@download/ENCFF155ZYZ.bed.gz", "prostate_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/prostate_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF201VZW/@@download/ENCFF201VZW.bed.gz", "prostate_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/prostate_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF198QRW/@@download/ENCFF198QRW.bed.gz", "prostate_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/prostate_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF028CPL/@@download/ENCFF028CPL.bed.gz", "prostate_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/prostate_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF351NKT/@@download/ENCFF351NKT.bed.gz", "prostate_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/prostate_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF233XNV/@@download/ENCFF233XNV.bed.gz", "prostate_dnase_hg38.bed.gz")
		
		#Skeletal muscle (Homo sapiens muscle of leg tissue female embryo (110 days), Homo sapiens muscle of leg tissue female embryo (115 days) (DNase-seq only), and Homo sapiens gastrocnemius medialis tissue male adult (37 years) (CTCF ChIP-seq only))
		elif tissue == "skeletal_muscle":
			if not os.path.exists(args.ref_dir + "/skeletal_muscle_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF170ZFK/@@download/ENCFF170ZFK.bed.gz", "skeletal_muscle_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/skeletal_muscle_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF942MOM/@@download/ENCFF942MOM.bed.gz", "skeletal_muscle_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/skeletal_muscle_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF154ZCF/@@download/ENCFF154ZCF.bed.gz", "skeletal_muscle_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/skeletal_muscle_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF672MZK/@@download/ENCFF672MZK.bed.gz", "skeletal_muscle_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/skeletal_muscle_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF495BIY/@@download/ENCFF495BIY.bed.gz", "skeletal_muscle_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/skeletal_muscle_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF143NGY/@@download/ENCFF143NGY.bed.gz", "skeletal_muscle_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/skeletal_muscle_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF055FYP/@@download/ENCFF055FYP.bed.gz", "skeletal_muscle_dnase_hg38.bed.gz")

		#Skin (Homo sapiens suprapubic skin tissue male adult (37 years) and Homo sapiens foreskin fibroblast male newborn (DNase-seq only))
		elif tissue == "skin":
			if not os.path.exists(args.ref_dir + "/skin_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF800MDH/@@download/ENCFF800MDH.bed.gz", "skin_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/skin_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF240CFZ/@@download/ENCFF240CFZ.bed.gz", "skin_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/skin_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF398EEO/@@download/ENCFF398EEO.bed.gz", "skin_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/skin_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF686GBZ/@@download/ENCFF686GBZ.bed.gz", "skin_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/skin_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF804VWO/@@download/ENCFF804VWO.bed.gz", "skin_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/skin_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF503GIR/@@download/ENCFF503GIR.bed.gz", "skin_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/skin_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF119TMJ/@@download/ENCFF119TMJ.bed.gz", "skin_dnase_hg38.bed.gz")
		
		#Spleen (Homo sapiens spleen tissue male adult (37 years) and Homo sapiens spleen tissue embryo (112 days) (DNase-seq only))
		elif tissue == "spleen":
			if not os.path.exists(args.ref_dir + "/spleen_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF652QYY/@@download/ENCFF652QYY.bed.gz", "spleen_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/spleen_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF401XEH/@@download/ENCFF401XEH.bed.gz", "spleen_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/spleen_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF322YIK/@@download/ENCFF322YIK.bed.gz", "spleen_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/spleen_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF159IED/@@download/ENCFF159IED.bed.gz", "spleen_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/spleen_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF380BZO/@@download/ENCFF380BZO.bed.gz", "spleen_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/spleen_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF455BVG/@@download/ENCFF455BVG.bed.gz", "spleen_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/spleen_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF686XTC/@@download/ENCFF686XTC.bed.gz", "spleen_dnase_hg38.bed.gz")
		
		#Stomach (Homo sapiens stomach tissue male adult (37 years) and Homo sapiens stomach tissue male adult (34 years) (DNase-seq only))
		elif tissue == "stomach":
			if not os.path.exists(args.ref_dir + "/stomach_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF475EJU/@@download/ENCFF475EJU.bed.gz", "stomach_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/stomach_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF494BWU/@@download/ENCFF494BWU.bed.gz", "stomach_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/stomach_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF978AOD/@@download/ENCFF978AOD.bed.gz", "stomach_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/stomach_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF692HXY/@@download/ENCFF692HXY.bed.gz", "stomach_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/stomach_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF243EAW/@@download/ENCFF243EAW.bed.gz", "stomach_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/stomach_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF533OIH/@@download/ENCFF533OIH.bed.gz", "stomach_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/stomach_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF033LEB/@@download/ENCFF033LEB.bed.gz", "stomach_dnase_hg38.bed.gz")
		
		#Testis (Homo sapiens testis tissue male adult (37 years) and Homo sapiens testis tissue male adult (54 years) (DNase-seq only))
		elif tissue == "testis":
			if not os.path.exists(args.ref_dir + "/testis_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF600MJO/@@download/ENCFF600MJO.bed.gz", "testis_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/testis_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF612RSR/@@download/ENCFF612RSR.bed.gz", "testis_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/testis_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF540YUQ/@@download/ENCFF540YUQ.bed.gz", "testis_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/testis_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF389XVY/@@download/ENCFF389XVY.bed.gz", "testis_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/testis_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF469ASG/@@download/ENCFF469ASG.bed.gz", "testis_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/testis_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF387GJA/@@download/ENCFF387GJA.bed.gz", "testis_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/testis_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF914YUP/@@download/ENCFF914YUP.bed.gz", "testis_dnase_hg38.bed.gz")
		
		#Thyroid (Homo sapiens thyroid gland tissue female adult (51 years))
		elif tissue == "thyroid":
			if not os.path.exists(args.ref_dir + "/thyroid_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF346IPY/@@download/ENCFF346IPY.bed.gz", "thyroid_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/thyroid_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF925FSR/@@download/ENCFF925FSR.bed.gz", "thyroid_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/thyroid_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF138WNV/@@download/ENCFF138WNV.bed.gz", "thyroid_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/thyroid_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF607NRZ/@@download/ENCFF607NRZ.bed.gz", "thyroid_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/thyroid_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF922QTV/@@download/ENCFF922QTV.bed.gz", "thyroid_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/thyroid_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF728ZGC/@@download/ENCFF728ZGC.bed.gz", "thyroid_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/thyroid_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF336FSX/@@download/ENCFF336FSX.bed.gz", "thyroid_dnase_hg38.bed.gz")
		
		#Uterus (Homo sapiens uterus tissue female adult (51 years) and Homo sapiens uterus tissue female adult (53 years) (DNase-seq only))
		elif tissue == "uterus":
			if not os.path.exists(args.ref_dir + "/uterus_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF150BRF/@@download/ENCFF150BRF.bed.gz", "uterus_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/uterus_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF022DOY/@@download/ENCFF022DOY.bed.gz", "uterus_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/uterus_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF108SAD/@@download/ENCFF108SAD.bed.gz", "uterus_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/uterus_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF254GYM/@@download/ENCFF254GYM.bed.gz", "uterus_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/uterus_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF374OMB/@@download/ENCFF374OMB.bed.gz", "uterus_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/uterus_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF447GBU/@@download/ENCFF447GBU.bed.gz", "uterus_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/uterus_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF661ERT/@@download/ENCFF661ERT.bed.gz", "uterus_dnase_hg38.bed.gz")
		
		#Vagina (Homo sapiens vagina tissue female adult (51 years))
		elif tissue == "vagina":
			if not os.path.exists(args.ref_dir + "/vagina_4me1_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF005TOV/@@download/ENCFF005TOV.bed.gz", "vagina_4me1_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/vagina_4me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF653IVY/@@download/ENCFF653IVY.bed.gz", "vagina_4me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/vagina_27ac_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF165SPP/@@download/ENCFF165SPP.bed.gz", "vagina_27ac_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/vagina_27me3_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF866QFX/@@download/ENCFF866QFX.bed.gz", "vagina_27me3_hg38.bed.gz")
			#if not os.path.exists(args.ref_dir + "/vagina_36me3_hg38.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF823VSG/@@download/ENCFF823VSG.bed.gz", "vagina_36me3_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/vagina_ctcf_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF626XNO/@@download/ENCFF626XNO.bed.gz", "vagina_ctcf_hg38.bed.gz")
			if not os.path.exists(args.ref_dir + "/vagina_dnase_hg38.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF872AJU/@@download/ENCFF872AJU.bed.gz", "vagina_dnase_hg38.bed.gz")
			
		#User-provided files
		elif tissue == "user_provided_files":
			user_files_index += 1
		
		#User-provided URLs
		elif tissue == "user_provided_urls":
			download_ref(user_4me1_array[user_files_index], (user_tissue_names_array[user_files_index] + "_4me1_hg38.bed.gz"))
			download_ref(user_4me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_4me3_hg38.bed.gz"))
			download_ref(user_27ac_array[user_files_index], (user_tissue_names_array[user_files_index] + "_27ac_hg38.bed.gz"))
			download_ref(user_27me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_27me3_hg38.bed.gz"))
			#download_ref(user_36me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_36me3_hg38.bed.gz"))
			download_ref(user_ctcf_array[user_files_index], (user_tissue_names_array[user_files_index] + "_ctcf_hg38.bed.gz"))
			download_ref(user_dnase_array[user_files_index], (user_tissue_names_array[user_files_index] + "_dnase_hg38.bed.gz"))
			user_files_index += 1
			

###################################################################################
#GRCh37/hg19
if args.ref_genome.lower() == "hg19" or args.ref_genome.lower() == "grch37":
	#This TSS dataset was manually lifted over from the original hg38 refTSS dataset and is imported from the CoRE-BED repository
	for tissue in tissue_array:
		#Adipose (Homo sapiens subcutaneous abdominal adipose tissue tissue nuclear fraction female adult (49 years), Homo sapiens omental fat pad tissue female adult (53 years) (DNase-seq only), and Homo sapiens omental fat pad tissue female adult (51 years) (CTCF ChIP-seq only))
		if tissue == "adipose":
			if not os.path.exists(args.ref_dir + "/adipose_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF258ONX/@@download/ENCFF258ONX.bed.gz", "adipose_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/adipose_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF665CSL/@@download/ENCFF665CSL.bed.gz", "adipose_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/adipose_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF621LOW/@@download/ENCFF621LOW.bed.gz", "adipose_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/adipose_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF464CMC/@@download/ENCFF464CMC.bed.gz", "adipose_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/adipose_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF553JCE/@@download/ENCFF553JCE.bed.gz", "adipose_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/adipose_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF517HGL/@@download/ENCFF517HGL.bed.gz", "adipose_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/adipose_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF925YKV/@@download/ENCFF925YKV.bed.gz", "adipose_dnase_hg19.bed.gz")
		
		#Adrenal gland (Homo sapiens adrenal gland tissue male embryo (97 days), Homo sapiens adrenal gland tissue embryo (96 days) (DNase-seq only), and Homo sapiens adrenal gland tissue male adult (37 years) (CTCF ChIP-seq only))
		if tissue == "adrenal_gland":
			if not os.path.exists(args.ref_dir + "/adrenal_gland_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF100QRG/@@download/ENCFF100QRG.bed.gz", "adrenal_gland_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/adrenal_gland_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF825DIG/@@download/ENCFF825DIG.bed.gz", "adrenal_gland_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/adrenal_gland_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF433MRW/@@download/ENCFF433MRW.bed.gz", "adrenal_gland_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/adrenal_gland_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF810FEE/@@download/ENCFF810FEE.bed.gz", "adrenal_gland_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/adrenal_gland_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF185XYP/@@download/ENCFF185XYP.bed.gz", "adrenal_gland_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/adrenal_gland_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF622LJD/@@download/ENCFF622LJD.bed.gz", "adrenal_gland_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/adrenal_gland_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF525FRH/@@download/ENCFF525FRH.bed.gz", "adrenal_gland_dnase_hg19.bed.gz")
		
		#Artery (Homo sapiens aorta tissue male adult (34 years), Homo sapiens coronary artery tissue female adult (53 years) (DNase-seq only), and Homo sapiens ascending aorta tissue female adult (51 years) (CTCF ChIP-seq only))
		if tissue == "artery":
			if not os.path.exists(args.ref_dir + "/artery_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF028YXN/@@download/ENCFF028YXN.bed.gz", "artery_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/artery_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF870QVN/@@download/ENCFF870QVN.bed.gz", "artery_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/artery_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF195VBI/@@download/ENCFF195VBI.bed.gz", "artery_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/artery_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF251MOP/@@download/ENCFF251MOP.bed.gz", "artery_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/artery_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF498YBT/@@download/ENCFF498YBT.bed.gz", "artery_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/artery_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF171RBT/@@download/ENCFF171RBT.bed.gz", "artery_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/artery_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF021ZDC/@@download/ENCFF021ZDC.bed.gz", "artery_dnase_hg19.bed.gz")
		
		#Blood (Homo sapiens peripheral blood mononuclear cell male adult (39 years) and Homo sapiens K562 (DNase-seq and CTCF ChIP-seq only))	
		if tissue == "blood":
			if not os.path.exists(args.ref_dir + "/blood_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF553SEK/@@download/ENCFF553SEK.bed.gz", "blood_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/blood_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF477HMF/@@download/ENCFF477HMF.bed.gz", "blood_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/blood_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF442KRP/@@download/ENCFF442KRP.bed.gz", "blood_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/blood_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF987VSV/@@download/ENCFF987VSV.bed.gz", "blood_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/blood_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF068OGE/@@download/ENCFF068OGE.bed.gz", "blood_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/blood_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF085HTY/@@download/ENCFF085HTY.bed.gz", "blood_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/blood_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF621ZJY/@@download/ENCFF621ZJY.bed.gz", "blood_dnase_hg19.bed.gz")
		
		#Breast (Homo sapiens breast epithelium tissue female adult (53 years) and Homo sapiens MCF-7 (DNase-seq only))
		elif tissue == "breast":
			if not os.path.exists(args.ref_dir + "/breast_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF336DDM/@@download/ENCFF336DDM.bed.gz", "breast_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/breast_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF065TIH/@@download/ENCFF065TIH.bed.gz", "breast_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/breast_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF154XFN/@@download/ENCFF154XFN.bed.gz", "breast_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/breast_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF291WFP/@@download/ENCFF291WFP.bed.gz", "breast_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/breast_27me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF906MJM/@@download/ENCFF906MJM.bed.gz", "breast_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/breast_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF722THT/@@download/ENCFF722THT.bed.gz", "breast_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/breast_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF846DFL/@@download/ENCFF846DFL.bed.gz", "breast_dnase_hg19.bed.gz")
		
		#Cultured fibroblast (Homo sapiens IMR-90)
		elif tissue == "cultured_fibroblast":
			if not os.path.exists(args.ref_dir + "/cultured_fibroblast_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF830JTY/@@download/ENCFF830JTY.bed.gz", "cultured_fibroblast_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/cultured_fibroblast_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF154CUR/@@download/ENCFF154CUR.bed.gz", "cultured_fibroblast_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/cultured_fibroblast_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF678QLP/@@download/ENCFF678QLP.bed.gz", "cultured_fibroblast_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/cultured_fibroblast_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF041RLH/@@download/ENCFF041RLH.bed.gz", "cultured_fibroblast_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/cultured_fibroblast_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF118DFR/@@download/ENCFF118DFR.bed.gz", "cultured_fibroblast_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/cultured_fibroblast_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF453XKM/@@download/ENCFF453XKM.bed.gz", "cultured_fibroblast_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/cultured_fibroblast_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF895BJS/@@download/ENCFF895BJS.bed.gz", "cultured_fibroblast_dnase_hg19.bed.gz")
	
		#EBV-transformed lymphocyte (Homo sapiens GM12878)
		elif tissue == "ebv_transformed_lymphocyte":
			if not os.path.exists(args.ref_dir + "/ebv_transformed_lymphocyte_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF921LKB/@@download/ENCFF921LKB.bed.gz", "ebv_transformed_lymphocyte_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/ebv_transformed_lymphocyte_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF295GNH/@@download/ENCFF295GNH.bed.gz", "ebv_transformed_lymphocyte_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/ebv_transformed_lymphocyte_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF816AHV/@@download/ENCFF816AHV.bed.gz", "ebv_transformed_lymphocyte_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/ebv_transformed_lymphocyte_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF247VUO/@@download/ENCFF247VUO.bed.gz", "ebv_transformed_lymphocyte_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/ebv_transformed_lymphocyte_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF233NNT/@@download/ENCFF233NNT.bed.gz", "ebv_transformed_lymphocyte_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/ebv_transformed_lymphocyte_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF710VEH/@@download/ENCFF710VEH.bed.gz", "ebv_transformed_lymphocyte_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/ebv_transformed_lymphocyte_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF273MVV/@@download/ENCFF273MVV.bed.gz", "ebv_transformed_lymphocyte_dnase_hg19.bed.gz")
		
		#Embryonic stem cell (Homo sapiens H1)
		elif tissue == "es":
			if not os.path.exists(args.ref_dir + "/es_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF834NWA/@@download/ENCFF834NWA.bed.gz", "es_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/es_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF865TYO/@@download/ENCFF865TYO.bed.gz", "es_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/es_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF335JSA/@@download/ENCFF335JSA.bed.gz", "es_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/es_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF900RIU/@@download/ENCFF900RIU.bed.gz", "es_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/es_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF452OKA/@@download/ENCFF452OKA.bed.gz", "es_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/es_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF093VEE/@@download/ENCFF093VEE.bed.gz", "es_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/es_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF810YLL/@@download/ENCFF810YLL.bed.gz", "es_dnase_hg19.bed.gz")
		
		#Esophagus muscularis mucosa (Homo sapiens esophagus muscularis mucosa tissue female adult (51 years) and Homo sapiens esophagus muscularis mucosa tissue male adult (37 years) (DNase-seq only))
		elif tissue == "esophagus_muscularis_mucosa":
			if not os.path.exists(args.ref_dir + "/esophagus_muscularis_mucosa_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF278OGS/@@download/ENCFF278OGS.bed.gz", "esophagus_muscularis_mucosa_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/esophagus_muscularis_mucosa_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF821CCT/@@download/ENCFF821CCT.bed.gz", "esophagus_muscularis_mucosa_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/esophagus_muscularis_mucosa_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF333FTU/@@download/ENCFF333FTU.bed.gz", "esophagus_muscularis_mucosa_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/esophagus_muscularis_mucosa_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF683CXO/@@download/ENCFF683CXO.bed.gz", "esophagus_muscularis_mucosa_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/esophagus_muscularis_mucosa_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF040XZW/@@download/ENCFF040XZW.bed.gz", "esophagus_muscularis_mucosa_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/esophagus_muscularis_mucosa_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF081PSJ/@@download/ENCFF081PSJ.bed.gz", "esophagus_muscularis_mucosa_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/esophagus_muscularis_mucosa_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF747AGO/@@download/ENCFF747AGO.bed.gz", "esophagus_muscularis_mucosa_dnase_hg19.bed.gz")
		
		#Esophagus squamous epithelium (Homo sapiens esophagus squamous epithelium tissue female adult (51 years) and Homo sapiens esophagus squamous epithelium tissue male adult (37 years) (DNase-seq only))
		elif tissue == "esophagus_squamous_epithelium":
			if not os.path.exists(args.ref_dir + "/esophagus_squamous_epithelium_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF693JWR/@@download/ENCFF693JWR.bed.gz", "esophagus_squamous_epithelium_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/esophagus_squamous_epithelium_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF659VMQ/@@download/ENCFF659VMQ.bed.gz", "esophagus_squamous_epithelium_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/esophagus_squamous_epithelium_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF081YZL/@@download/ENCFF081YZL.bed.gz", "esophagus_squamous_epithelium_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/esophagus_squamous_epithelium_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF639ADZ/@@download/ENCFF639ADZ.bed.gz", "esophagus_squamous_epithelium_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/esophagus_squamous_epithelium_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF641FAA/@@download/ENCFF641FAA.bed.gz", "esophagus_squamous_epithelium_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/esophagus_squamous_epithelium_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF401XPE/@@download/ENCFF401XPE.bed.gz", "esophagus_squamous_epithelium_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/esophagus_squamous_epithelium_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF197ZNE/@@download/ENCFF197ZNE.bed.gz", "esophagus_squamous_epithelium_dnase_hg19.bed.gz")
		
		#Heart (Homo sapiens heart left ventricle tissue female adult (53 years))
		elif tissue == "heart":
			if not os.path.exists(args.ref_dir + "/heart_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF266FCG/@@download/ENCFF266FCG.bed.gz", "heart_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/heart_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF086QMV/@@download/ENCFF086QMV.bed.gz", "heart_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/heart_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF546TJN/@@download/ENCFF546TJN.bed.gz", "heart_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/heart_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF614LBC/@@download/ENCFF614LBC.bed.gz", "heart_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/heart_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF079NXF/@@download/ENCFF079NXF.bed.gz", "heart_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/heart_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF240UFV/@@download/ENCFF240UFV.bed.gz", "heart_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/heart_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF939RCS/@@download/ENCFF939RCS.bed.gz", "heart_dnase_hg19.bed.gz")
		
		#Intestine (Homo sapiens sigmoid colon tissue male adult (37 years) and Homo sapiens sigmoid colon tissue female adult (53 years) (DNase-seq only))
		elif tissue == "intestine":
			if not os.path.exists(args.ref_dir + "/intestine_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF409HTF/@@download/ENCFF409HTF.bed.gz", "intestine_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/intestine_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF878SLR/@@download/ENCFF878SLR.bed.gz", "intestine_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/intestine_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF314GGO/@@download/ENCFF314GGO.bed.gz", "intestine_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/intestine_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF099PRP/@@download/ENCFF099PRP.bed.gz", "intestine_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/intestine_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF419JUG/@@download/ENCFF419JUG.bed.gz", "intestine_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/intestine_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF637HAC/@@download/ENCFF637HAC.bed.gz", "intestine_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/intestine_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF882CST/@@download/ENCFF882CST.bed.gz", "intestine_dnase_hg19.bed.gz")
		
		#iPS cell (Homo sapiens iPS DF 19.11 and Homo sapiens GM23338 originated from GM23248 (CTCF ChIP-seq only))
		elif tissue == "ips":
			if not os.path.exists(args.ref_dir + "/ips_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF485HSK/@@download/ENCFF485HSK.bed.gz", "ips_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/ips_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF538JKF/@@download/ENCFF538JKF.bed.gz", "ips_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/ips_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF631JHU/@@download/ENCFF631JHU.bed.gz", "ips_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/ips_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF396DQE/@@download/ENCFF396DQE.bed.gz", "ips_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/ips_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF510FJC/@@download/ENCFF510FJC.bed.gz", "ips_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/ips_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF760HUV/@@download/ENCFF760HUV.bed.gz", "ips_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/ips_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF779LNR/@@download/ENCFF779LNR.bed.gz", "ips_dnase_hg19.bed.gz")
		
		#Kidney (Homo sapiens kidney tissue male adult (67 years), Homo sapiens kidney tissue female embryo (120 days) for H3K27me3 - only dataset for kidney and DNase-seq, and Homo sapiens kidney epithelial cell (CTCF ChIP-seq))
		elif tissue == "kidney":
			if not os.path.exists(args.ref_dir + "/kidney_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF388YQF/@@download/ENCFF388YQF.bed.gz", "kidney_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/kidney_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF913FVX/@@download/ENCFF913FVX.bed.gz", "kidney_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/kidney_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF953SWO/@@download/ENCFF953SWO.bed.gz", "kidney_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/kidney_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF754PGB/@@download/ENCFF754PGB.bed.gz", "kidney_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/kidney_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF216OKG/@@download/ENCFF216OKG.bed.gz", "kidney_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/kidney_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF168VBK/@@download/ENCFF168VBK.bed.gz", "kidney_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/kidney_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF330CNE/@@download/ENCFF330CNE.bed.gz", "kidney_dnase_hg19.bed.gz")
	
		#Liver (Homo sapiens liver tissue male adult (32 years), Homo sapiens liver tissue female embryo (101 days) and female embryo (113 days) (DNase-seq), and Homo sapiens right lobe of liver tissue female adult (53 years) (CTCF ChIP-seq only))
		elif tissue == "liver":
			if not os.path.exists(args.ref_dir + "/liver_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF327BBS/@@download/ENCFF327BBS.bed.gz", "liver_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/liver_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF065NRN/@@download/ENCFF065NRN.bed.gz", "liver_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/liver_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF752QSK/@@download/ENCFF752QSK.bed.gz", "liver_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/liver_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF933NEK/@@download/ENCFF933NEK.bed.gz", "liver_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/liver_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF744HEP/@@download/ENCFF744HEP.bed.gz", "liver_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/liver_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF665OBP/@@download/ENCFF665OBP.bed.gz", "liver_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/liver_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF286LYP/@@download/ENCFF286LYP.bed.gz", "liver_dnase_hg19.bed.gz")
		
		#Lung (Homo sapiens upper lobe of left lung tissue female adult (51 years) and Homo sapiens upper lobe of left lung tissue male adult (37 years) (DNase-seq only))
		elif tissue == "lung":
			if not os.path.exists(args.ref_dir + "/lung_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF433RBT/@@download/ENCFF433RBT.bed.gz", "lung_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/lung_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF142ZXD/@@download/ENCFF142ZXD.bed.gz", "lung_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/lung_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF307QJK/@@download/ENCFF307QJK.bed.gz", "lung_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/lung_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF957XLF/@@download/ENCFF957XLF.bed.gz", "lung_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/lung_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF344NKJ/@@download/ENCFF344NKJ.bed.gz", "lung_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/lung_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF232DOG/@@download/ENCFF232DOG.bed.gz", "lung_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/lung_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF889ADL/@@download/ENCFF889ADL.bed.gz", "lung_dnase_hg19.bed.gz")
		
		#Neuron (Homo sapiens SK-N-SH and Homo sapiens tibial nerve tissue female adult (51 years) (DNase-seq only))
		elif tissue == "neuron":
			if not os.path.exists(args.ref_dir + "/neuron_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF580GTZ/@@download/ENCFF580GTZ.bed.gz", "neuron_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/neuron_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF363ZFM/@@download/ENCFF363ZFM.bed.gz", "neuron_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/neuron_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF362OBM/@@download/ENCFF362OBM.bed.gz", "neuron_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/neuron_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF102YDR/@@download/ENCFF102YDR.bed.gz", "neuron_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/neuron_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF712RJX/@@download/ENCFF712RJX.bed.gz", "neuron_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/neuron_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF861DPF/@@download/ENCFF861DPF.bed.gz", "neuron_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/neuron_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF424BYY/@@download/ENCFF424BYY.bed.gz", "neuron_dnase_hg19.bed.gz")
		
		#Ovary (Homo sapiens ovary tissue female adult (30 years) and Homo sapiens ovary tissue female adult (53 years) (CTCF ChIP-seq only))
		elif tissue == "ovary":
			if not os.path.exists(args.ref_dir + "/ovary_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF917PWI/@@download/ENCFF917PWI.bed.gz", "ovary_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/ovary_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF320JHG/@@download/ENCFF320JHG.bed.gz", "ovary_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/ovary_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF657AUA/@@download/ENCFF657AUA.bed.gz", "ovary_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/ovary_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF712UCB/@@download/ENCFF712UCB.bed.gz", "ovary_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/ovary_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF302DXB/@@download/ENCFF302DXB.bed.gz", "ovary_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/ovary_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF261BWI/@@download/ENCFF261BWI.bed.gz", "ovary_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/ovary_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF883WWT/@@download/ENCFF883WWT.bed.gz", "ovary_dnase_hg19.bed.gz")
		
		#Pancreas (Homo sapiens pancreas tissue male adult (34 years) and Homo sapiens body of pancreas tissue male adult (37 years) (CTCF ChIP-seq only))
		elif tissue == "pancreas":
			if not os.path.exists(args.ref_dir + "/pancreas_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF668HLF/@@download/ENCFF668HLF.bed.gz", "pancreas_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/pancreas_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF340YEE/@@download/ENCFF340YEE.bed.gz", "pancreas_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/pancreas_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF583QFI/@@download/ENCFF583QFI.bed.gz", "pancreas_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/pancreas_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF973GJC/@@download/ENCFF973GJC.bed.gz", "pancreas_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/pancreas_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF544VYY/@@download/ENCFF544VYY.bed.gz", "pancreas_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/pancreas_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF025MZL/@@download/ENCFF025MZL.bed.gz", "pancreas_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/pancreas_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF897PRD/@@download/ENCFF897PRD.bed.gz", "pancreas_dnase_hg19.bed.gz")
		
		#Prostate (Homo sapiens prostate gland tissue male adult (37 years) and Homo sapiens epithelial cell of prostate (DNase-seq only))
		elif tissue == "prostate":
			if not os.path.exists(args.ref_dir + "/prostate_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF942RGZ/@@download/ENCFF942RGZ.bed.gz", "prostate_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/prostate_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF495DBS/@@download/ENCFF495DBS.bed.gz", "prostate_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/prostate_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF002NFG/@@download/ENCFF002NFG.bed.gz", "prostate_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/prostate_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF318QJC/@@download/ENCFF318QJC.bed.gz", "prostate_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/prostate_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF122XKF/@@download/ENCFF122XKF.bed.gz", "prostate_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/prostate_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF899MQP/@@download/ENCFF899MQP.bed.gz", "prostate_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/prostate_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF608WCU/@@download/ENCFF608WCU.bed.gz", "prostate_dnase_hg19.bed.gz")
		
		#Skeletal muscle (Homo sapiens muscle of leg tissue female embryo (110 days), Homo sapiens muscle of leg tissue female embryo (115 days) (DNase-seq only), and Homo sapiens gastrocnemius medialis tissue male adult (37 years) (CTCF ChIP-seq only))
		elif tissue == "skeletal_muscle":
			if not os.path.exists(args.ref_dir + "/skeletal_muscle_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF746HHI/@@download/ENCFF746HHI.bed.gz", "skeletal_muscle_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/skeletal_muscle_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF027HBT/@@download/ENCFF027HBT.bed.gz", "skeletal_muscle_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/skeletal_muscle_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF563PBB/@@download/ENCFF563PBB.bed.gz", "skeletal_muscle_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/skeletal_muscle_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF191PRP/@@download/ENCFF191PRP.bed.gz", "skeletal_muscle_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/skeletal_muscle_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF678XMX/@@download/ENCFF678XMX.bed.gz", "skeletal_muscle_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/skeletal_muscle_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF588XCW/@@download/ENCFF588XCW.bed.gz", "skeletal_muscle_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/skeletal_muscle_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF470QNA/@@download/ENCFF470QNA.bed.gz", "skeletal_muscle_dnase_hg19.bed.gz")
		
		#Skin (Homo sapiens suprapubic skin tissue male adult (37 years) and Homo sapiens foreskin fibroblast male newborn (DNase-seq only))
		elif tissue == "skin":
			if not os.path.exists(args.ref_dir + "/skin_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF904DTJ/@@download/ENCFF904DTJ.bed.gz", "skin_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/skin_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF778GBQ/@@download/ENCFF778GBQ.bed.gz", "skin_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/skin_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF390JYY/@@download/ENCFF390JYY.bed.gz", "skin_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/skin_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF752AID/@@download/ENCFF752AID.bed.gz", "skin_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/skin_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF106JDV/@@download/ENCFF106JDV.bed.gz", "skin_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/skin_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF129ZVT/@@download/ENCFF129ZVT.bed.gz", "skin_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/skin_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF169THH/@@download/ENCFF169THH.bed.gz", "skin_dnase_hg19.bed.gz")
		
		#Spleen (Homo sapiens spleen tissue male adult (37 years) and Homo sapiens spleen tissue embryo (112 days) (DNase-seq only))
		elif tissue == "spleen":
			if not os.path.exists(args.ref_dir + "/spleen_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF585JQZ/@@download/ENCFF585JQZ.bed.gz", "spleen_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/spleen_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF736XRS/@@download/ENCFF736XRS.bed.gz", "spleen_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/spleen_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF538SIJ/@@download/ENCFF538SIJ.bed.gz", "spleen_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/spleen_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF065IBC/@@download/ENCFF065IBC.bed.gz", "spleen_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/spleen_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF382JFO/@@download/ENCFF382JFO.bed.gz", "spleen_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/spleen_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF459AHK/@@download/ENCFF459AHK.bed.gz", "spleen_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/spleen_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF587YNA/@@download/ENCFF587YNA.bed.gz", "spleen_dnase_hg19.bed.gz")
		
		#Stomach (Homo sapiens stomach tissue male adult (37 years) and Homo sapiens stomach tissue male adult (34 years) (DNase-seq only))
		elif tissue == "stomach":
			if not os.path.exists(args.ref_dir + "/stomach_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF324NLX/@@download/ENCFF324NLX.bed.gz", "stomach_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/stomach_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF161ERT/@@download/ENCFF161ERT.bed.gz", "stomach_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/stomach_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF506AKR/@@download/ENCFF506AKR.bed.gz", "stomach_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/stomach_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF435CWM/@@download/ENCFF435CWM.bed.gz", "stomach_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/stomach_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF984QYQ/@@download/ENCFF984QYQ.bed.gz", "stomach_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/stomach_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF703CWY/@@download/ENCFF703CWY.bed.gz", "stomach_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/stomach_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF033LEB/@@download/ENCFF033LEB.bed.gz", "stomach_dnase_hg19.bed.gz")
		
		#Testis (Homo sapiens testis tissue male adult (37 years) and Homo sapiens testis tissue male adult (54 years) (DNase-seq only))
		elif tissue == "testis":
			if not os.path.exists(args.ref_dir + "/testis_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF620AJW/@@download/ENCFF620AJW.bed.gz", "testis_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/testis_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF047XWN/@@download/ENCFF047XWN.bed.gz", "testis_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/testis_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF567XUE/@@download/ENCFF567XUE.bed.gz", "testis_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/testis_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF563BJS/@@download/ENCFF563BJS.bed.gz", "testis_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/testis_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF074JIX/@@download/ENCFF074JIX.bed.gz", "testis_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/testis_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF432XLE/@@download/ENCFF432XLE.bed.gz", "testis_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/testis_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF914YUP/@@download/ENCFF914YUP.bed.gz", "testis_dnase_hg19.bed.gz")
		
		#Thyroid (Homo sapiens thyroid gland tissue female adult (51 years) and Homo sapiens thymus tissue female embryo (147 days) (DNase-seq only))
		elif tissue == "thyroid":
			if not os.path.exists(args.ref_dir + "/thyroid_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF346LTA/@@download/ENCFF346LTA.bed.gz", "thyroid_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/thyroid_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF228FJJ/@@download/ENCFF228FJJ.bed.gz", "thyroid_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/thyroid_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF086VAY/@@download/ENCFF086VAY.bed.gz", "thyroid_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/thyroid_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF957XLF/@@download/ENCFF957XLF.bed.gz", "thyroid_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/thyroid_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF839YKD/@@download/ENCFF839YKD.bed.gz", "thyroid_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/thyroid_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF520JZM/@@download/ENCFF520JZM.bed.gz", "thyroid_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/thyroid_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF494PNZ/@@download/ENCFF494PNZ.bed.gz", "thyroid_dnase_hg19.bed.gz")
		
		#Uterus (Homo sapiens uterus tissue female adult (51 years) and Homo sapiens HeLa-S3 (DNase-seq only))
		elif tissue == "uterus":
			if not os.path.exists(args.ref_dir + "/uterus_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF995XPF/@@download/ENCFF995XPF.bed.gz", "uterus_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/uterus_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF155UOS/@@download/ENCFF155UOS.bed.gz", "uterus_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/uterus_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF471JES/@@download/ENCFF471JES.bed.gz", "uterus_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/uterus_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF119JPA/@@download/ENCFF119JPA.bed.gz", "uterus_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/uterus_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF551UHB/@@download/ENCFF551UHB.bed.gz", "uterus_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/uterus_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF450VON/@@download/ENCFF450VON.bed.gz", "uterus_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/uterus_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF775GBM/@@download/ENCFF775GBM.bed.gz", "uterus_dnase_hg19.bed.gz")
		
		#Vagina (Homo sapiens vagina tissue female adult (51 years))
		elif tissue == "vagina":
			if not os.path.exists(args.ref_dir + "/vagina_4me1_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF324OUK/@@download/ENCFF324OUK.bed.gz", "vagina_4me1_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/vagina_4me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF367PZY/@@download/ENCFF367PZY.bed.gz", "vagina_4me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/vagina_27ac_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF435LEQ/@@download/ENCFF435LEQ.bed.gz", "vagina_27ac_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/vagina_27me3_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF993WTL/@@download/ENCFF993WTL.bed.gz", "vagina_27me3_hg19.bed.gz")
			#if not os.path.exists(args.ref_dir + "/vagina_36me3_hg19.bed.gz"):
				#download_ref("https://www.encodeproject.org/files/ENCFF643BMR/@@download/ENCFF643BMR.bed.gz", "vagina_36me3_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/vagina_ctcf_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF498YGJ/@@download/ENCFF498YGJ.bed.gz", "vagina_ctcf_hg19.bed.gz")
			if not os.path.exists(args.ref_dir + "/vagina_dnase_hg19.bed.gz"):
				download_ref("https://www.encodeproject.org/files/ENCFF600WRF/@@download/ENCFF600WRF.bed.gz", "vagina_dnase_hg19.bed.gz")
		
		#User-provided files
		elif tissue == "user_provided_files":
			user_files_index += 1
		
		#User-provided URLs
		elif tissue == "user_provided_urls":
			download_ref(user_4me1_array[user_files_index], (user_tissue_names_array[user_files_index] + "_4me1_hg19.bed.gz"))
			download_ref(user_4me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_4me3_hg19.bed.gz"))
			download_ref(user_27ac_array[user_files_index], (user_tissue_names_array[user_files_index] + "_27ac_hg19.bed.gz"))
			download_ref(user_27me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_27me3_hg19.bed.gz"))
			#download_ref(user_36me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_36me3_hg19.bed.gz"))
			download_ref(user_ctcf_array[user_files_index], (user_tissue_names_array[user_files_index] + "_ctcf_hg19.bed.gz"))
			download_ref(user_dnase_array[user_files_index], (user_tissue_names_array[user_files_index] + "_dnase_hg19.bed.gz"))
			user_files_index += 1
		
###Mouse###
#GRCm39/mm39
elif args.ref_genome.lower() == "mm39" or args.ref_genome.lower() == "grcm39":
	#This TSS dataset was manually lifted over from the original mm10 refTSS dataset and is imported from the CoRE-BED repository
	for tissue in tissue_array:
		#User-provided files
		if tissue == "user_provided_files":
			user_files_index += 1
		
		#User-provided URLs
		elif tissue == "user_provided_urls":
			download_ref(user_4me1_array[user_files_index], (user_tissue_names_array[user_files_index] + "_4me1_mm39.bed.gz"))
			download_ref(user_4me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_4me3_mm39.bed.gz"))
			download_ref(user_27ac_array[user_files_index], (user_tissue_names_array[user_files_index] + "_27ac_mm39.bed.gz"))
			download_ref(user_27me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_27me3_mm39.bed.gz"))
			#download_ref(user_36me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_36me3_mm39.bed.gz"))
			download_ref(user_ctcf_array[user_files_index], (user_tissue_names_array[user_files_index] + "_ctcf_mm39.bed.gz"))
			download_ref(user_dnase_array[user_files_index], (user_tissue_names_array[user_files_index] + "_dnase_mm39.bed.gz"))
			user_files_index += 1

#GRCm38/mm10
elif args.ref_genome.lower() == "mm10" or args.ref_genome.lower() == "grcm38":
	if not os.path.exists(args.ref_dir + "/refTSS_v3.1_mouse_coordinate.mm10.bed.gz"):
		download_ref("http://reftss.clst.riken.jp/datafiles/current/mouse/refTSS_v3.1_mouse_coordinate.mm10.bed.gz", "refTSS_v3.1_mouse_coordinate.mm10.bed.gz")
	for tissue in tissue_array:
		#User-provided files
		if tissue == "user_provided_files":
			user_files_index += 1
		
		#User-provided URLs
		elif tissue == "user_provided_urls":
			download_ref(user_4me1_array[user_files_index], (user_tissue_names_array[user_files_index] + "_4me1_mm10.bed.gz"))
			download_ref(user_4me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_4me3_mm10.bed.gz"))
			download_ref(user_27ac_array[user_files_index], (user_tissue_names_array[user_files_index] + "_27ac_mm10.bed.gz"))
			download_ref(user_27me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_27me3_mm10.bed.gz"))
			#download_ref(user_36me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_36me3_mm10.bed.gz"))
			download_ref(user_ctcf_array[user_files_index], (user_tissue_names_array[user_files_index] + "_ctcf_mm10.bed.gz"))
			download_ref(user_dnase_array[user_files_index], (user_tissue_names_array[user_files_index] + "_dnase_mm10.bed.gz"))
			user_files_index += 1
		
	
#GRCm37/mm9
elif args.ref_genome.lower() == "mm9" or args.ref_genome.lower() == "grcm37":
	#This TSS dataset was manually lifted over from the original mm10 refTSS dataset and is imported from the CoRE-BED repository
	for tissue in tissue_array:
		#User-provided files
		if tissue == "user_provided_files":
			user_files_index += 1
		
		#User-provided URLs
		elif tissue == "user_provided_urls":
			download_ref(user_4me1_array[user_files_index], (user_tissue_names_array[user_files_index] + "_4me1_mm9.bed.gz"))
			download_ref(user_4me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_4me3_mm9.bed.gz"))
			download_ref(user_27ac_array[user_files_index], (user_tissue_names_array[user_files_index] + "_27ac_mm9.bed.gz"))
			download_ref(user_27me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_27me3_mm9.bed.gz"))
			#download_ref(user_36me3_array[user_files_index], (user_tissue_names_array[user_files_index] + "_36me3_mm9.bed.gz"))
			download_ref(user_ctcf_array[user_files_index], (user_tissue_names_array[user_files_index] + "_ctcf_mm9.bed.gz"))
			download_ref(user_dnase_array[user_files_index], (user_tissue_names_array[user_files_index] + "_dnase_mm9.bed.gz"))
			user_files_index += 1
	
###Import all of the bed files to be worked on (input and histone mark ChIP-seq) using pybedtools
##Human
if args.ref_genome.lower() == "hg38" or args.ref_genome.lower() == "grch38":
	bed_ref = "hg38"
	bed_tss_path = args.ref_dir + "/refTSS_v3.1_human_coordinate.hg38.bed.gz"
elif args.ref_genome.lower() == "hg19" or args.ref_genome.lower() == "grch37":
	bed_ref = "hg19"
	script_dir = os.path.dirname(os.path.realpath(__file__))
	bed_tss_path = script_dir + "/" + "lifted_reftss_files/" + "refTSS_v3.1_human_coordinate.lifted.sorted.hg19.bed.gz"

##Mouse
if args.ref_genome.lower() == "mm39" or args.ref_genome.lower() == "grcm39":
	bed_ref = "mm39"
	script_dir = os.path.dirname(os.path.realpath(__file__))
	bed_tss_path = script_dir + "/" + "lifted_reftss_files/" + "refTSS_v3.1_mouse_coordinate.lifted.sorted.mm39.bed.gz"
elif args.ref_genome.lower() == "mm10" or args.ref_genome.lower() == "grcm38":
	bed_ref = "mm10"
	bed_tss_path = args.ref_dir + "/refTSS_v3.1_mouse_coordinate.mm10.bed.gz"
elif args.ref_genome.lower() == "mm9" or args.ref_genome.lower() == "mm9":
	bed_ref = "mm9"
	script_dir = os.path.dirname(os.path.realpath(__file__))
	bed_tss_path = script_dir + "/" + "lifted_reftss_files/" + "refTSS_v3.1_mouse_coordinate.lifted.sorted.mm9.bed.gz"

#Manipulate the TSS bed file so that the coordinates are expanded +-2kb (based on the enrichment we see in Mas et al. 2018 - https://www.nature.com/articles/s41588-018-0218-5?proof=t)
if args.verbose:
	print("Importing TSS and input files...")
	print("\n")
bed_tss_orig = pd.read_csv(bed_tss_path, sep = "\t", header = None)
bed_tss_orig.iloc[:,1] = bed_tss_orig.iloc[:,1] - args.tss_distance_upstream
bed_tss_orig.iloc[:,2] = bed_tss_orig.iloc[:,2] + args.tss_distance_downstream
bed_tss_orig = bed_tss_orig[bed_tss_orig.iloc[:,1] >= 0]
bed_tss = pybedtools.BedTool().from_dataframe(bed_tss_orig)

#Import the input file, along with the histone ChIP-seq bed files
##Because the raw input file will not necessarily be in UCSC BED format, first import it as a pandas data frame, sort by coordinates, and then input coordinates as a pybedtools object
<<<<<<< HEAD
if args.separator == " " or args.separator == "\t":
	if args.input_header:
		file_input = pd.read_csv(args.input, delim_whitespace = True, low_memory = False)
	else:
		file_input = pd.read_csv(args.input, delim_whitespace = True, header = None, low_memory = False)
else:
	if args.input_header:
		file_input = pd.read_csv(args.input, delimiter = args.separator, low_memory = False)
	else:
		file_input = pd.read_csv(args.input, delimiter = args.separator, header = None, low_memory = False)
=======
if args.input_header:
	file_input = pd.read_csv(args.input, delim_whitespace = True, low_memory = False)
else:
	file_input = pd.read_csv(args.input, delim_whitespace = True, header = None, low_memory = False)
>>>>>>> 3283f061327aed512684aab1c4b30a7a1c01e0ad
	
cols_array = args.bed_cols.split(",")
chr_col = int(cols_array[0]) - 1
start_col = int(cols_array[1]) - 1
end_col = int(cols_array[2]) - 1

col_names = file_input.columns
file_input = file_input.sort_values(by = [col_names[chr_col], col_names[start_col], col_names[end_col]])

##Sort the input file by coordinates
bed_input = file_input[[col_names[chr_col], col_names[start_col], col_names[end_col]]]
if str(bed_input.iloc[0,0]).startswith("chr") == False:
	bed_input.iloc[:,0] = "chr" + bed_input.iloc[:,0].astype(str)
bed_input = pybedtools.BedTool.from_dataframe(bed_input)

#Loop through each reference tissue parsed in from the input arguments, generating contextual functional annotations and appending each set as a new column onto the original input file
user_files_index = 0
n_cols_input_file = len(file_input.columns)
anno_df = pd.DataFrame()
for tissue in tissue_array:
	if args.verbose:
		if tissue == "user_provided_files" or tissue == "user_provided_urls":
			print("Finding overlap of {input_file} coordinates with {tissue_type} histone ChIP-seq peaks based on {ref_genome}...".format(input_file = args.input, tissue_type = user_tissue_names_array[user_files_index], ref_genome = args.ref_genome))
			print("\n")
		else:
			print("Finding overlap of {input_file} coordinates with {tissue_type} histone ChIP-seq peaks based on {ref_genome}...".format(input_file = args.input, tissue_type = tissue, ref_genome = args.ref_genome))
			print("\n")	
	if tissue == "user_provided_files":
		bed_4me1 = pybedtools.BedTool(user_4me1_array[user_files_index])
		bed_4me3 = pybedtools.BedTool(user_4me3_array[user_files_index])
		bed_27ac = pybedtools.BedTool(user_27ac_array[user_files_index])
		bed_27me3 = pybedtools.BedTool(user_27me3_array[user_files_index])
		#bed_36me3 = pybedtools.BedTool(user_36me3_array[user_files_index])
		bed_ctcf = pybedtools.BedTool(user_ctcf_array[user_files_index])
		bed_dnase = pybedtools.BedTool(user_dnase_array[user_files_index])
	elif tissue == "user_provided_urls":
		bed_4me1 = pybedtools.BedTool(args.ref_dir + "/" + user_tissue_names_array[user_files_index] + "_4me1_" + bed_ref + ".bed.gz")
		bed_4me3 = pybedtools.BedTool(args.ref_dir + "/" + user_tissue_names_array[user_files_index] + "_4me3_" + bed_ref + ".bed.gz")
		bed_27ac = pybedtools.BedTool(args.ref_dir + "/" + user_tissue_names_array[user_files_index] + "_27ac_" + bed_ref + ".bed.gz")
		bed_27me3 = pybedtools.BedTool(args.ref_dir + "/" + user_tissue_names_array[user_files_index] + "_27me3_" + bed_ref + ".bed.gz")
		#bed_36me3 = pybedtools.BedTool(args.ref_dir + "/" + user_tissue_names_array[user_files_index] + "_36me3_" + bed_ref + ".bed.gz")
		bed_ctcf = pybedtools.BedTool(args.ref_dir + "/" + user_tissue_names_array[user_files_index] + "_ctcf_" + bed_ref + ".bed.gz")
		bed_dnase = pybedtools.BedTool(args.ref_dir + "/" + user_tissue_names_array[user_files_index] + "_dnase_" + bed_ref + ".bed.gz")
	else:
		bed_4me1 = pybedtools.BedTool(args.ref_dir + "/" + tissue + "_4me1_" + bed_ref + ".bed.gz")
		bed_4me3 = pybedtools.BedTool(args.ref_dir + "/" + tissue + "_4me3_" + bed_ref + ".bed.gz")
		bed_27ac = pybedtools.BedTool(args.ref_dir + "/" + tissue + "_27ac_" + bed_ref + ".bed.gz")
		bed_27me3 = pybedtools.BedTool(args.ref_dir + "/" + tissue + "_27me3_" + bed_ref + ".bed.gz")
		#bed_36me3 = pybedtools.BedTool(args.ref_dir + "/" + tissue + "_36me3_" + bed_ref + ".bed.gz")
		bed_ctcf = pybedtools.BedTool(args.ref_dir + "/" + tissue + "_ctcf_" + bed_ref + ".bed.gz")
		bed_dnase = pybedtools.BedTool(args.ref_dir + "/" + tissue + "_dnase_" + bed_ref + ".bed.gz")
<<<<<<< HEAD
	
	#If one of the special flags specifying only one epigenomic mark is enabled, then evaluate this overlap. Otherwise, move on to the normal workflow.
	if args.only_h3k4me1:
		overlaps_4me1 = bed_input.intersect(bed_4me1, u = True)
		overlaps_4me1 = overlaps_4me1.sort()
		
		no_4me1 = bed_input.intersect(bed_4me1, v = True)
		no_4me1 = no_4me1.sort()
		
		if overlaps_4me1.count() > 0:
			overlaps_4me1_df = pybedtools.BedTool.to_dataframe(overlaps_4me1)
			overlaps_4me1_df["region_classification"] = "True"
		else:
			overlaps_4me1_df = pd.DataFrame()
	
		if no_4me1.count() > 0:
			no_4me1_df = pybedtools.BedTool.to_dataframe(no_4me1)
			no_4me1_df["region_classification"] = "False"
		else:
			no_4me1_df = pd.DataFrame()
		
		all_data_frames = [overlaps_4me1_df, no_4me1_df]

		concatenated_df = pd.concat(all_data_frames)
		concatenated_df = concatenated_df.drop_duplicates()
		concatenated_df.iloc[:,0] = concatenated_df.iloc[:,0].str[3:]
		col_names = ["chr", "start", "end", "func_anno"]
		concatenated_df.columns = col_names
		concatenated_df = concatenated_df.sort_values(by = [col_names[0], col_names[1], col_names[2]])

		#For rows that are completely identical (all entries except for the regulatory annotation), the rows are collapsed into one, with the multiple annotations concatenated in a comma-separated string)
		n_cols = len(concatenated_df.columns)
		names_array = []
		counter = 1
		for i in range(0, n_cols):
			colname = "col_" + str(counter)
			names_array.append(colname)
			counter += 1
		concatenated_df.columns = names_array
		last_col = names_array[-1]
		concatenated_df[last_col] = concatenated_df.groupby(names_array[0:-1])[last_col].transform(lambda x: ','.join(x))
		concatenated_df = concatenated_df.drop_duplicates()

		#If a coordinate overlaps with multiple regulatory elements and the user specifically wishes to keep only one, then we will keep only the most significant (i.e. active_promoter over unclassified_within_2kb_of_tss)
		if args.no_multianno:
			from operator import itemgetter
			annotations = concatenated_df.iloc[:,-1]
			annotations = annotations.str.split(',')
			firsts = list(map(itemgetter(0), annotations))
			concatenated_df.drop(concatenated_df.columns[-1], axis = 1, inplace = True)
			concatenated_df[names_array[-1]] = firsts
	
		#Merge the original input file with the functional annotations from the concatenated data frame
		concatenated_df["col_1"] = concatenated_df["col_1"].astype(str)
		concatenated_df["col_2"] = concatenated_df["col_2"].astype(str)
		concatenated_df["col_3"] = concatenated_df["col_3"].astype(str)
	
		file_input_colnames = file_input.columns
		chr_file_input = file_input_colnames[chr_col]
		start_file_input = file_input_colnames[start_col]
		end_file_input = file_input_colnames[end_col]
	
		file_input[chr_file_input] = file_input[chr_file_input].astype(str)
		file_input[start_file_input] = file_input[start_file_input].astype(str)
		file_input[end_file_input] = file_input[end_file_input].astype(str)
	
		merged_df = file_input.merge(concatenated_df, how = "left", left_on = [chr_file_input, start_file_input, end_file_input], right_on = ["col_1", "col_2", "col_3"])
		merged_df = merged_df.drop(columns = ["col_1", "col_2", "col_3"])
	
		#Rename the functional annotation column with the tissue of interest
		n_cols_merged_df = len(merged_df.columns)
		names_array_merged_df = []
		for i in range(0, n_cols_merged_df):
			colname = merged_df.columns[i]
			names_array_merged_df.append(colname)
		
		if tissue == "user_provided_files" or tissue == "user_provided_urls":
			func_names = user_tissue_names_array[user_files_index].upper() + "_FUNCTIONAL_ANNOTATION"
			names_array_merged_df[-1] = func_names
			user_files_index += 1
		else:
			func_names = tissue.upper() + "_FUNCTIONAL_ANNOTATION"
			names_array_merged_df[-1] = func_names
		
		merged_df.columns = names_array_merged_df
	
		last_merged_df_column = names_array_merged_df[-1]
		anno_df[func_names] = merged_df[func_names]
	
	elif args.only_h3k4me3:
		overlaps_4me3 = bed_input.intersect(bed_4me3, u = True)
		overlaps_4me3 = overlaps_4me3.sort()
		
		no_4me3 = bed_input.intersect(bed_4me3, v = True)
		no_4me3 = no_4me3.sort()
		
		if overlaps_4me3.count() > 0:
			overlaps_4me3_df = pybedtools.BedTool.to_dataframe(overlaps_4me3)
			overlaps_4me3_df["region_classification"] = "True"
		else:
			overlaps_4me3_df = pd.DataFrame()
	
		if no_4me3.count() > 0:
			no_4me3_df = pybedtools.BedTool.to_dataframe(no_4me3)
			no_4me3_df["region_classification"] = "False"
		else:
			no_4me3_df = pd.DataFrame()
		
		all_data_frames = [overlaps_4me3_df, no_4me3_df]

		concatenated_df = pd.concat(all_data_frames)
		concatenated_df = concatenated_df.drop_duplicates()
		concatenated_df.iloc[:,0] = concatenated_df.iloc[:,0].str[3:]
		col_names = ["chr", "start", "end", "func_anno"]
		concatenated_df.columns = col_names
		concatenated_df = concatenated_df.sort_values(by = [col_names[0], col_names[1], col_names[2]])

		#For rows that are completely identical (all entries except for the regulatory annotation), the rows are collapsed into one, with the multiple annotations concatenated in a comma-separated string)
		n_cols = len(concatenated_df.columns)
		names_array = []
		counter = 1
		for i in range(0, n_cols):
			colname = "col_" + str(counter)
			names_array.append(colname)
			counter += 1
		concatenated_df.columns = names_array
		last_col = names_array[-1]
		concatenated_df[last_col] = concatenated_df.groupby(names_array[0:-1])[last_col].transform(lambda x: ','.join(x))
		concatenated_df = concatenated_df.drop_duplicates()

		#If a coordinate overlaps with multiple regulatory elements and the user specifically wishes to keep only one, then we will keep only the most significant (i.e. active_promoter over unclassified_within_2kb_of_tss)
		if args.no_multianno:
			from operator import itemgetter
			annotations = concatenated_df.iloc[:,-1]
			annotations = annotations.str.split(',')
			firsts = list(map(itemgetter(0), annotations))
			concatenated_df.drop(concatenated_df.columns[-1], axis = 1, inplace = True)
			concatenated_df[names_array[-1]] = firsts
	
		#Merge the original input file with the functional annotations from the concatenated data frame
		concatenated_df["col_1"] = concatenated_df["col_1"].astype(str)
		concatenated_df["col_2"] = concatenated_df["col_2"].astype(str)
		concatenated_df["col_3"] = concatenated_df["col_3"].astype(str)
	
		file_input_colnames = file_input.columns
		chr_file_input = file_input_colnames[chr_col]
		start_file_input = file_input_colnames[start_col]
		end_file_input = file_input_colnames[end_col]
	
		file_input[chr_file_input] = file_input[chr_file_input].astype(str)
		file_input[start_file_input] = file_input[start_file_input].astype(str)
		file_input[end_file_input] = file_input[end_file_input].astype(str)
	
		merged_df = file_input.merge(concatenated_df, how = "left", left_on = [chr_file_input, start_file_input, end_file_input], right_on = ["col_1", "col_2", "col_3"])
		merged_df = merged_df.drop(columns = ["col_1", "col_2", "col_3"])
	
		#Rename the functional annotation column with the tissue of interest
		n_cols_merged_df = len(merged_df.columns)
		names_array_merged_df = []
		for i in range(0, n_cols_merged_df):
			colname = merged_df.columns[i]
			names_array_merged_df.append(colname)
		
		if tissue == "user_provided_files" or tissue == "user_provided_urls":
			func_names = user_tissue_names_array[user_files_index].upper() + "_FUNCTIONAL_ANNOTATION"
			names_array_merged_df[-1] = func_names
			user_files_index += 1
		else:
			func_names = tissue.upper() + "_FUNCTIONAL_ANNOTATION"
			names_array_merged_df[-1] = func_names
		
		merged_df.columns = names_array_merged_df
	
		last_merged_df_column = names_array_merged_df[-1]
		anno_df[func_names] = merged_df[func_names]
		
	elif args.only_h3k27ac:
		overlaps_27ac = bed_input.intersect(bed_27ac, u = True)
		overlaps_27ac = overlaps_27ac.sort()
		
		no_27ac = bed_input.intersect(bed_27ac, v = True)
		no_27ac = no_27ac.sort()
		
		if overlaps_27ac.count() > 0:
			overlaps_27ac_df = pybedtools.BedTool.to_dataframe(overlaps_27ac)
			overlaps_27ac_df["region_classification"] = "True"
		else:
			overlaps_27ac_df = pd.DataFrame()
	
		if no_27ac.count() > 0:
			no_27ac_df = pybedtools.BedTool.to_dataframe(no_27ac)
			no_27ac_df["region_classification"] = "False"
		else:
			no_27ac_df = pd.DataFrame()
		
		all_data_frames = [overlaps_27ac_df, no_27ac_df]

		concatenated_df = pd.concat(all_data_frames)
		concatenated_df = concatenated_df.drop_duplicates()
		concatenated_df.iloc[:,0] = concatenated_df.iloc[:,0].str[3:]
		col_names = ["chr", "start", "end", "func_anno"]
		concatenated_df.columns = col_names
		concatenated_df = concatenated_df.sort_values(by = [col_names[0], col_names[1], col_names[2]])

		#For rows that are completely identical (all entries except for the regulatory annotation), the rows are collapsed into one, with the multiple annotations concatenated in a comma-separated string)
		n_cols = len(concatenated_df.columns)
		names_array = []
		counter = 1
		for i in range(0, n_cols):
			colname = "col_" + str(counter)
			names_array.append(colname)
			counter += 1
		concatenated_df.columns = names_array
		last_col = names_array[-1]
		concatenated_df[last_col] = concatenated_df.groupby(names_array[0:-1])[last_col].transform(lambda x: ','.join(x))
		concatenated_df = concatenated_df.drop_duplicates()

		#If a coordinate overlaps with multiple regulatory elements and the user specifically wishes to keep only one, then we will keep only the most significant (i.e. active_promoter over unclassified_within_2kb_of_tss)
		if args.no_multianno:
			from operator import itemgetter
			annotations = concatenated_df.iloc[:,-1]
			annotations = annotations.str.split(',')
			firsts = list(map(itemgetter(0), annotations))
			concatenated_df.drop(concatenated_df.columns[-1], axis = 1, inplace = True)
			concatenated_df[names_array[-1]] = firsts
	
		#Merge the original input file with the functional annotations from the concatenated data frame
		concatenated_df["col_1"] = concatenated_df["col_1"].astype(str)
		concatenated_df["col_2"] = concatenated_df["col_2"].astype(str)
		concatenated_df["col_3"] = concatenated_df["col_3"].astype(str)
	
		file_input_colnames = file_input.columns
		chr_file_input = file_input_colnames[chr_col]
		start_file_input = file_input_colnames[start_col]
		end_file_input = file_input_colnames[end_col]
	
		file_input[chr_file_input] = file_input[chr_file_input].astype(str)
		file_input[start_file_input] = file_input[start_file_input].astype(str)
		file_input[end_file_input] = file_input[end_file_input].astype(str)
	
		merged_df = file_input.merge(concatenated_df, how = "left", left_on = [chr_file_input, start_file_input, end_file_input], right_on = ["col_1", "col_2", "col_3"])
		merged_df = merged_df.drop(columns = ["col_1", "col_2", "col_3"])
	
		#Rename the functional annotation column with the tissue of interest
		n_cols_merged_df = len(merged_df.columns)
		names_array_merged_df = []
		for i in range(0, n_cols_merged_df):
			colname = merged_df.columns[i]
			names_array_merged_df.append(colname)
		
		if tissue == "user_provided_files" or tissue == "user_provided_urls":
			func_names = user_tissue_names_array[user_files_index].upper() + "_FUNCTIONAL_ANNOTATION"
			names_array_merged_df[-1] = func_names
			user_files_index += 1
		else:
			func_names = tissue.upper() + "_FUNCTIONAL_ANNOTATION"
			names_array_merged_df[-1] = func_names
		
		merged_df.columns = names_array_merged_df
	
		last_merged_df_column = names_array_merged_df[-1]
		anno_df[func_names] = merged_df[func_names]
		
	elif args.only_h3k27me3:
		overlaps_27me3 = bed_input.intersect(bed_27me3, u = True)
		overlaps_27me3 = overlaps_27me3.sort()
		
		no_27me3 = bed_input.intersect(bed_27me3, v = True)
		no_27me3 = no_27me3.sort()
		
		if overlaps_27me3.count() > 0:
			overlaps_27me3_df = pybedtools.BedTool.to_dataframe(overlaps_27me3)
			overlaps_27me3_df["region_classification"] = "True"
		else:
			overlaps_27me3_df = pd.DataFrame()
	
		if no_27me3.count() > 0:
			no_27me3_df = pybedtools.BedTool.to_dataframe(no_27me3)
			no_27me3_df["region_classification"] = "False"
		else:
			no_27me3_df = pd.DataFrame()
		
		all_data_frames = [overlaps_27me3_df, no_27me3_df]

		concatenated_df = pd.concat(all_data_frames)
		concatenated_df = concatenated_df.drop_duplicates()
		concatenated_df.iloc[:,0] = concatenated_df.iloc[:,0].str[3:]
		col_names = ["chr", "start", "end", "func_anno"]
		concatenated_df.columns = col_names
		concatenated_df = concatenated_df.sort_values(by = [col_names[0], col_names[1], col_names[2]])

		#For rows that are completely identical (all entries except for the regulatory annotation), the rows are collapsed into one, with the multiple annotations concatenated in a comma-separated string)
		n_cols = len(concatenated_df.columns)
		names_array = []
		counter = 1
		for i in range(0, n_cols):
			colname = "col_" + str(counter)
			names_array.append(colname)
			counter += 1
		concatenated_df.columns = names_array
		last_col = names_array[-1]
		concatenated_df[last_col] = concatenated_df.groupby(names_array[0:-1])[last_col].transform(lambda x: ','.join(x))
		concatenated_df = concatenated_df.drop_duplicates()

		#If a coordinate overlaps with multiple regulatory elements and the user specifically wishes to keep only one, then we will keep only the most significant (i.e. active_promoter over unclassified_within_2kb_of_tss)
		if args.no_multianno:
			from operator import itemgetter
			annotations = concatenated_df.iloc[:,-1]
			annotations = annotations.str.split(',')
			firsts = list(map(itemgetter(0), annotations))
			concatenated_df.drop(concatenated_df.columns[-1], axis = 1, inplace = True)
			concatenated_df[names_array[-1]] = firsts
	
		#Merge the original input file with the functional annotations from the concatenated data frame
		concatenated_df["col_1"] = concatenated_df["col_1"].astype(str)
		concatenated_df["col_2"] = concatenated_df["col_2"].astype(str)
		concatenated_df["col_3"] = concatenated_df["col_3"].astype(str)
	
		file_input_colnames = file_input.columns
		chr_file_input = file_input_colnames[chr_col]
		start_file_input = file_input_colnames[start_col]
		end_file_input = file_input_colnames[end_col]
	
		file_input[chr_file_input] = file_input[chr_file_input].astype(str)
		file_input[start_file_input] = file_input[start_file_input].astype(str)
		file_input[end_file_input] = file_input[end_file_input].astype(str)
	
		merged_df = file_input.merge(concatenated_df, how = "left", left_on = [chr_file_input, start_file_input, end_file_input], right_on = ["col_1", "col_2", "col_3"])
		merged_df = merged_df.drop(columns = ["col_1", "col_2", "col_3"])
	
		#Rename the functional annotation column with the tissue of interest
		n_cols_merged_df = len(merged_df.columns)
		names_array_merged_df = []
		for i in range(0, n_cols_merged_df):
			colname = merged_df.columns[i]
			names_array_merged_df.append(colname)
		
		if tissue == "user_provided_files" or tissue == "user_provided_urls":
			func_names = user_tissue_names_array[user_files_index].upper() + "_FUNCTIONAL_ANNOTATION"
			names_array_merged_df[-1] = func_names
			user_files_index += 1
		else:
			func_names = tissue.upper() + "_FUNCTIONAL_ANNOTATION"
			names_array_merged_df[-1] = func_names
		
		merged_df.columns = names_array_merged_df
	
		last_merged_df_column = names_array_merged_df[-1]
		anno_df[func_names] = merged_df[func_names]
		
	elif args.only_dhs:
		overlaps_dhs = bed_input.intersect(bed_dhs, u = True)
		overlaps_dhs = overlaps_dhs.sort()
		
		no_dhs = bed_input.intersect(bed_dhs, v = True)
		no_dhs = no_dhs.sort()
		
		if overlaps_dhs.count() > 0:
			overlaps_dhs_df = pybedtools.BedTool.to_dataframe(overlaps_dhs)
			overlaps_dhs_df["region_classification"] = "True"
		else:
			overlaps_dhs_df = pd.DataFrame()
	
		if no_dhs.count() > 0:
			no_dhs_df = pybedtools.BedTool.to_dataframe(no_dhs)
			no_dhs_df["region_classification"] = "False"
		else:
			no_dhs_df = pd.DataFrame()
		
		all_data_frames = [overlaps_dhs_df, no_dhs_df]

		concatenated_df = pd.concat(all_data_frames)
		concatenated_df = concatenated_df.drop_duplicates()
		concatenated_df.iloc[:,0] = concatenated_df.iloc[:,0].str[3:]
		col_names = ["chr", "start", "end", "func_anno"]
		concatenated_df.columns = col_names
		concatenated_df = concatenated_df.sort_values(by = [col_names[0], col_names[1], col_names[2]])

		#For rows that are completely identical (all entries except for the regulatory annotation), the rows are collapsed into one, with the multiple annotations concatenated in a comma-separated string)
		n_cols = len(concatenated_df.columns)
		names_array = []
		counter = 1
		for i in range(0, n_cols):
			colname = "col_" + str(counter)
			names_array.append(colname)
			counter += 1
		concatenated_df.columns = names_array
		last_col = names_array[-1]
		concatenated_df[last_col] = concatenated_df.groupby(names_array[0:-1])[last_col].transform(lambda x: ','.join(x))
		concatenated_df = concatenated_df.drop_duplicates()

		#If a coordinate overlaps with multiple regulatory elements and the user specifically wishes to keep only one, then we will keep only the most significant (i.e. active_promoter over unclassified_within_2kb_of_tss)
		if args.no_multianno:
			from operator import itemgetter
			annotations = concatenated_df.iloc[:,-1]
			annotations = annotations.str.split(',')
			firsts = list(map(itemgetter(0), annotations))
			concatenated_df.drop(concatenated_df.columns[-1], axis = 1, inplace = True)
			concatenated_df[names_array[-1]] = firsts
	
		#Merge the original input file with the functional annotations from the concatenated data frame
		concatenated_df["col_1"] = concatenated_df["col_1"].astype(str)
		concatenated_df["col_2"] = concatenated_df["col_2"].astype(str)
		concatenated_df["col_3"] = concatenated_df["col_3"].astype(str)
	
		file_input_colnames = file_input.columns
		chr_file_input = file_input_colnames[chr_col]
		start_file_input = file_input_colnames[start_col]
		end_file_input = file_input_colnames[end_col]
	
		file_input[chr_file_input] = file_input[chr_file_input].astype(str)
		file_input[start_file_input] = file_input[start_file_input].astype(str)
		file_input[end_file_input] = file_input[end_file_input].astype(str)
	
		merged_df = file_input.merge(concatenated_df, how = "left", left_on = [chr_file_input, start_file_input, end_file_input], right_on = ["col_1", "col_2", "col_3"])
		merged_df = merged_df.drop(columns = ["col_1", "col_2", "col_3"])
	
		#Rename the functional annotation column with the tissue of interest
		n_cols_merged_df = len(merged_df.columns)
		names_array_merged_df = []
		for i in range(0, n_cols_merged_df):
			colname = merged_df.columns[i]
			names_array_merged_df.append(colname)
		
		if tissue == "user_provided_files" or tissue == "user_provided_urls":
			func_names = user_tissue_names_array[user_files_index].upper() + "_FUNCTIONAL_ANNOTATION"
			names_array_merged_df[-1] = func_names
			user_files_index += 1
		else:
			func_names = tissue.upper() + "_FUNCTIONAL_ANNOTATION"
			names_array_merged_df[-1] = func_names
		
		merged_df.columns = names_array_merged_df
	
		last_merged_df_column = names_array_merged_df[-1]
		anno_df[func_names] = merged_df[func_names]
		
	elif args.only_ctcf:
		overlaps_ctcf = bed_input.intersect(bed_ctcf, u = True)
		overlaps_ctcf = overlaps_ctcf.sort()
	
		no_ctcf = bed_input.intersect(bed_ctcf, v = True)
		no_ctcf = no_ctcf.sort()
	
		if overlaps_ctcf.count() > 0:
			overlaps_ctcf_df = pybedtools.BedTool.to_dataframe(overlaps_ctcf)
			overlaps_ctcf_df["region_classification"] = "True"
		else:
			overlaps_ctcf_df = pd.DataFrame()

		if no_ctcf.count() > 0:
			no_ctcf_df = pybedtools.BedTool.to_dataframe(no_ctcf)
			no_ctcf_df["region_classification"] = "False"
		else:
			no_ctcf_df = pd.DataFrame()
	
		all_data_frames = [overlaps_ctcf_df, no_ctcf_df]

		concatenated_df = pd.concat(all_data_frames)
		concatenated_df = concatenated_df.drop_duplicates()
		concatenated_df.iloc[:,0] = concatenated_df.iloc[:,0].str[3:]
		col_names = ["chr", "start", "end", "func_anno"]
		concatenated_df.columns = col_names
		concatenated_df = concatenated_df.sort_values(by = [col_names[0], col_names[1], col_names[2]])

		#For rows that are completely identical (all entries except for the regulatory annotation), the rows are collapsed into one, with the multiple annotations concatenated in a comma-separated string)
		n_cols = len(concatenated_df.columns)
		names_array = []
		counter = 1
		for i in range(0, n_cols):
			colname = "col_" + str(counter)
			names_array.append(colname)
			counter += 1
		concatenated_df.columns = names_array
		last_col = names_array[-1]
		concatenated_df[last_col] = concatenated_df.groupby(names_array[0:-1])[last_col].transform(lambda x: ','.join(x))
		concatenated_df = concatenated_df.drop_duplicates()

		#If a coordinate overlaps with multiple regulatory elements and the user specifically wishes to keep only one, then we will keep only the most significant (i.e. active_promoter over unclassified_within_2kb_of_tss)
		if args.no_multianno:
			from operator import itemgetter
			annotations = concatenated_df.iloc[:,-1]
			annotations = annotations.str.split(',')
			firsts = list(map(itemgetter(0), annotations))
			concatenated_df.drop(concatenated_df.columns[-1], axis = 1, inplace = True)
			concatenated_df[names_array[-1]] = firsts

		#Merge the original input file with the functional annotations from the concatenated data frame
		concatenated_df["col_1"] = concatenated_df["col_1"].astype(str)
		concatenated_df["col_2"] = concatenated_df["col_2"].astype(str)
		concatenated_df["col_3"] = concatenated_df["col_3"].astype(str)

		file_input_colnames = file_input.columns
		chr_file_input = file_input_colnames[chr_col]
		start_file_input = file_input_colnames[start_col]
		end_file_input = file_input_colnames[end_col]

		file_input[chr_file_input] = file_input[chr_file_input].astype(str)
		file_input[start_file_input] = file_input[start_file_input].astype(str)
		file_input[end_file_input] = file_input[end_file_input].astype(str)

		merged_df = file_input.merge(concatenated_df, how = "left", left_on = [chr_file_input, start_file_input, end_file_input], right_on = ["col_1", "col_2", "col_3"])
		merged_df = merged_df.drop(columns = ["col_1", "col_2", "col_3"])

		#Rename the functional annotation column with the tissue of interest
		n_cols_merged_df = len(merged_df.columns)
		names_array_merged_df = []
		for i in range(0, n_cols_merged_df):
			colname = merged_df.columns[i]
			names_array_merged_df.append(colname)
	
		if tissue == "user_provided_files" or tissue == "user_provided_urls":
			func_names = user_tissue_names_array[user_files_index].upper() + "_FUNCTIONAL_ANNOTATION"
			names_array_merged_df[-1] = func_names
			user_files_index += 1
		else:
			func_names = tissue.upper() + "_FUNCTIONAL_ANNOTATION"
			names_array_merged_df[-1] = func_names
			
		merged_df.columns = names_array_merged_df
	
		last_merged_df_column = names_array_merged_df[-1]
		anno_df[func_names] = merged_df[func_names]

	###Compare the input file to the reference bed files, starting with the modified TSS bed
	##Do the peaks fall within boundary of a TSS?
	else:
		overlaps_tss = bed_input.intersect(bed_tss, u = True)
		overlaps_tss = overlaps_tss.sort()
		#print(overlaps_tss.count())
		no_tss = bed_input.intersect(bed_tss, v = True)
		overlaps_tss = overlaps_tss.sort()
		#print(no_tss.count())

		#Classifying putative promoters
		overlaps_4me3 = overlaps_tss.intersect(bed_4me3, u = True)
		overlaps_4me3 = overlaps_4me3.sort()
		#print(overlaps_4me3.count())
		no_4me3 = overlaps_tss.intersect(bed_4me3, v = True)
		no_4me3 = no_4me3.sort()
		#print(no_4me3.count())
	
		active_promoter = overlaps_4me3.intersect(bed_27me3, v = True)
		active_promoter = active_promoter.sort()
		#print(active_promoter.count())
 
		bivalent_promoter = overlaps_4me3.intersect(bed_27me3, u = True)
		bivalent_promoter = bivalent_promoter.sort()
		#print(bivalent_promoter.count())
	
		silenced_promoter = no_4me3.intersect(bed_27me3, u = True)
		silenced_promoter = silenced_promoter.sort()
		#print(silenced_promoter.count())

		###Further classify unclassified promoter-like regions
		unclassified_overlaps_tss = overlaps_tss.intersect(active_promoter, v = True)
		unclassified_overlaps_tss = unclassified_overlaps_tss.intersect(bivalent_promoter, v = True)
		unclassified_overlaps_tss = unclassified_overlaps_tss.intersect(silenced_promoter, v = True)
		#print(unclassified_overlaps_tss.count())
		tss_overlaps_ctcf = unclassified_overlaps_tss.intersect(bed_ctcf, u = True)
		tss_no_ctcf = unclassified_overlaps_tss.intersect(bed_ctcf, v = True)

		ctcf_open_within_tss = tss_overlaps_ctcf.intersect(bed_dnase, u = True)
		ctcf_open_within_tss = ctcf_open_within_tss.sort()

		ctcf_closed_within_tss = tss_overlaps_ctcf.intersect(bed_dnase, v = True)
		ctcf_closed_within_tss = ctcf_closed_within_tss.sort()
	
		no_ctcf_open_within_tss = tss_no_ctcf.intersect(bed_dnase, u = True)
		no_ctcf_open_within_tss = no_ctcf_open_within_tss.sort()
	
		no_ctcf_closed_within_tss = tss_no_ctcf.intersect(bed_dnase, v = True)
		no_ctcf_closed_within_tss = no_ctcf_closed_within_tss.sort()
	
		#Classifying putative enhancers
		overlaps_4me1 = no_tss.intersect(bed_4me1, u = True)

		active_enhancer = overlaps_4me1.intersect(bed_27ac, u = True)
		active_enhancer = active_enhancer.sort()

		poised_enhancer = overlaps_4me1.intersect(bed_27me3, u = True)
		poised_enhancer = poised_enhancer.sort()

		primed_enhancer = overlaps_4me1.intersect(active_enhancer, v = True)
		primed_enhancer = primed_enhancer.intersect(poised_enhancer, v = True)
		primed_enhancer = primed_enhancer.sort()

		###Further classify unclassified enhancer-like regions
		unclassified_no_tss = no_tss.intersect(active_enhancer, v = True)
		unclassified_no_tss = no_tss.intersect(poised_enhancer, v = True)
		unclassified_no_tss = no_tss.intersect(primed_enhancer, v = True)
	
		no_tss_overlaps_ctcf = unclassified_no_tss.intersect(bed_ctcf, u = True)
		no_tss_no_ctcf = unclassified_no_tss.intersect(bed_ctcf, v = True)

		ctcf_open_no_tss = no_tss_overlaps_ctcf.intersect(bed_dnase, u = True)
		ctcf_open_no_tss = ctcf_open_no_tss.sort()

		ctcf_closed_no_tss = no_tss_overlaps_ctcf.intersect(bed_dnase, v = True)
		ctcf_closed_no_tss = ctcf_closed_no_tss.sort()

		no_ctcf_open_no_tss = no_tss_no_ctcf.intersect(bed_dnase, u = True)
		no_ctcf_open_no_tss = no_ctcf_open_no_tss.sort()

		no_ctcf_closed_no_tss = no_tss_no_ctcf.intersect(bed_dnase, v = True)
		no_ctcf_closed_no_tss = no_ctcf_closed_no_tss.sort()

		#If they have one or more coordinate entries, convert each of the classifier variables to a pandas DataFrame object so that we can append the classifications as an additional column
		if active_promoter.count() > 0:
			active_promoter_df = pybedtools.BedTool.to_dataframe(active_promoter)
			active_promoter_df["region_classification"] = "active_promoter"
		else:
			active_promoter_df = pd.DataFrame()
	
		if bivalent_promoter.count() > 0:
			bivalent_promoter_df = pybedtools.BedTool.to_dataframe(bivalent_promoter)
			bivalent_promoter_df["region_classification"] = "bivalent_promoter"
		else:
			bivalent_promoter_df = pd.DataFrame()

		if silenced_promoter.count() > 0:
			silenced_promoter_df = pybedtools.BedTool.to_dataframe(silenced_promoter)
			silenced_promoter_df["region_classification"] = "silenced_promoter"
		else:
			silenced_promoter_df = pd.DataFrame()

		if ctcf_open_within_tss.count() > 0:
			ctcf_open_within_tss_df = pybedtools.BedTool.to_dataframe(ctcf_open_within_tss)
			ctcf_open_within_tss_df["region_classification"] = ("unclassified_open_chromatin;ctcf_binding_site;within_" + str(args.tss_distance_upstream) + "_bp_upstream_" + str(args.tss_distance_downstream) + "_bp_downstream_of_tss")
		else:
			ctcf_open_within_tss_df = pd.DataFrame()

		if ctcf_closed_within_tss.count() > 0:
			ctcf_closed_within_tss_df = pybedtools.BedTool.to_dataframe(ctcf_closed_within_tss)
			ctcf_closed_within_tss_df["region_classification"] = ("unclassified_closed_chromatin;ctcf_binding_site;within_" + str(args.tss_distance_upstream) + "_bp_upstream_" + str(args.tss_distance_downstream) + "_bp_downstream_of_tss")
		else:
			ctcf_closed_within_tss_df = pd.DataFrame()
	
		if no_ctcf_open_within_tss.count() > 0:
			no_ctcf_open_within_tss_df = pybedtools.BedTool.to_dataframe(no_ctcf_open_within_tss)
			no_ctcf_open_within_tss_df["region_classification"] = ("unclassified_open_chromatin_within_" + str(args.tss_distance_upstream) + "_bp_upstream_" + str(args.tss_distance_downstream) + "_bp_downstream_of_tss")
		else:
			no_ctcf_open_within_tss_df = pd.DataFrame()
	
		if no_ctcf_closed_within_tss.count() > 0:
			no_ctcf_closed_within_tss_df = pybedtools.BedTool.to_dataframe(no_ctcf_closed_within_tss)
			no_ctcf_closed_within_tss_df["region_classification"] = ("unclassified_closed_chromatin_within_" + str(args.tss_distance_upstream) + "_bp_upstream_" + str(args.tss_distance_downstream) + "_bp_downstream_of_tss")
		else:
			no_ctcf_closed_within_tss_df = pd.DataFrame()

		if active_enhancer.count() > 0:
			active_enhancer_df = pybedtools.BedTool.to_dataframe(active_enhancer)
			active_enhancer_df["region_classification"] = "active_enhancer"
		else:
			active_enhancer_df = pd.DataFrame()

		if poised_enhancer.count() > 0:
			poised_enhancer_df = pybedtools.BedTool.to_dataframe(poised_enhancer)
			poised_enhancer_df["region_classification"] = "poised_enhancer"
		else:
			poised_enhancer_df = pd.DataFrame()

		if primed_enhancer.count() > 0:
			primed_enhancer_df = pybedtools.BedTool.to_dataframe(primed_enhancer)
			primed_enhancer_df["region_classification"] = "primed_enhancer"
		else:
			primed_enhancer_df = pd.DataFrame()

		if ctcf_open_no_tss.count() > 0:
			ctcf_open_no_tss_df = pybedtools.BedTool.to_dataframe(ctcf_open_no_tss)
			ctcf_open_no_tss_df["region_classification"] = ("unclassified_open_chromatin;ctcf_binding_site;beyond_" + str(args.tss_distance_upstream) + "_bp_upstream_" + str(args.tss_distance_downstream) + "_bp_downstream_of_tss")
		else:
			ctcf_open_no_tss_df = pd.DataFrame()
	
		if ctcf_closed_no_tss.count() > 0:
			ctcf_closed_no_tss_df = pybedtools.BedTool.to_dataframe(ctcf_closed_no_tss)
			ctcf_closed_no_tss_df["region_classification"] = ("unclassified_closed_chromatin;ctcf_binding_site;beyond_" + str(args.tss_distance_upstream) + "_bp_upstream_" + str(args.tss_distance_downstream) + "_bp_downstream_of_tss")
		else:
			ctcf_closed_no_tss_df = pd.DataFrame()
	
		if no_ctcf_open_no_tss.count() > 0:
			no_ctcf_open_no_tss_df = pybedtools.BedTool.to_dataframe(no_ctcf_open_no_tss)
			no_ctcf_open_no_tss_df["region_classification"] = ("unclassified_open_chromatin_beyond_" + str(args.tss_distance_upstream) + "_bp_upstream_" + str(args.tss_distance_downstream) + "_bp_downstream_of_tss")
		else:
			no_ctcf_open_no_tss_df = pd.DataFrame()
	
		if no_ctcf_closed_no_tss.count() > 0:
			no_ctcf_closed_no_tss_df = pybedtools.BedTool.to_dataframe(no_ctcf_closed_no_tss)
			no_ctcf_closed_no_tss_df["region_classification"] = ("unclassified_closed_chromatin_beyond_" + str(args.tss_distance_upstream) + "_bp_upstream_" + str(args.tss_distance_downstream) + "_bp_downstream_of_tss")
		else:
			no_ctcf_closed_no_tss_df = pd.DataFrame()
	
		#Concatenate all of the pandas DataFrame objects into one, remove all duplicate lines, and sort by coordinate
		all_data_frames = [active_promoter_df, bivalent_promoter_df, silenced_promoter_df, ctcf_open_within_tss_df, ctcf_closed_within_tss_df, no_ctcf_open_within_tss_df, no_ctcf_closed_within_tss_df, active_enhancer_df, poised_enhancer_df, primed_enhancer_df, ctcf_open_no_tss_df, ctcf_closed_no_tss_df, no_ctcf_open_no_tss_df, no_ctcf_closed_no_tss_df]

		concatenated_df = pd.concat(all_data_frames)
		concatenated_df = concatenated_df.drop_duplicates()
		concatenated_df.iloc[:,0] = concatenated_df.iloc[:,0].str[3:]
		col_names = ["chr", "start", "end", "func_anno"]
		concatenated_df.columns = col_names
		concatenated_df = concatenated_df.sort_values(by = [col_names[0], col_names[1], col_names[2]])

		#For rows that are completely identical (all entries except for the regulatory annotation), the rows are collapsed into one, with the multiple annotations concatenated in a comma-separated string)
		n_cols = len(concatenated_df.columns)
		names_array = []
		counter = 1
		for i in range(0, n_cols):
			colname = "col_" + str(counter)
			names_array.append(colname)
			counter += 1
		concatenated_df.columns = names_array
		last_col = names_array[-1]
		concatenated_df[last_col] = concatenated_df.groupby(names_array[0:-1])[last_col].transform(lambda x: ','.join(x))
		concatenated_df = concatenated_df.drop_duplicates()

		#If a coordinate overlaps with multiple regulatory elements and the user specifically wishes to keep only one, then we will keep only the most significant (i.e. active_promoter over unclassified_within_2kb_of_tss)
		if args.no_multianno:
			from operator import itemgetter
			annotations = concatenated_df.iloc[:,-1]
			annotations = annotations.str.split(',')
			firsts = list(map(itemgetter(0), annotations))
			concatenated_df.drop(concatenated_df.columns[-1], axis = 1, inplace = True)
			concatenated_df[names_array[-1]] = firsts
	
		#Merge the original input file with the functional annotations from the concatenated data frame
		concatenated_df["col_1"] = concatenated_df["col_1"].astype(str)
		concatenated_df["col_2"] = concatenated_df["col_2"].astype(str)
		concatenated_df["col_3"] = concatenated_df["col_3"].astype(str)
	
		file_input_colnames = file_input.columns
		chr_file_input = file_input_colnames[chr_col]
		start_file_input = file_input_colnames[start_col]
		end_file_input = file_input_colnames[end_col]
	
		file_input[chr_file_input] = file_input[chr_file_input].astype(str)
		file_input[start_file_input] = file_input[start_file_input].astype(str)
		file_input[end_file_input] = file_input[end_file_input].astype(str)
	
		merged_df = file_input.merge(concatenated_df, how = "left", left_on = [chr_file_input, start_file_input, end_file_input], right_on = ["col_1", "col_2", "col_3"])
		merged_df = merged_df.drop(columns = ["col_1", "col_2", "col_3"])
	
		#Rename the functional annotation column with the tissue of interest
		n_cols_merged_df = len(merged_df.columns)
		names_array_merged_df = []
		for i in range(0, n_cols_merged_df):
			colname = merged_df.columns[i]
			names_array_merged_df.append(colname)
		
		if tissue == "user_provided_files" or tissue == "user_provided_urls":
			func_names = user_tissue_names_array[user_files_index].upper() + "_FUNCTIONAL_ANNOTATION"
			names_array_merged_df[-1] = func_names
			user_files_index += 1
		else:
			func_names = tissue.upper() + "_FUNCTIONAL_ANNOTATION"
			names_array_merged_df[-1] = func_names
		
		merged_df.columns = names_array_merged_df
	
		last_merged_df_column = names_array_merged_df[-1]
		anno_df[func_names] = merged_df[func_names]
	
		#Write out the number of peaks/coordinates in each category if the --write_summary flag is specified. Also print this statement to the console if verbose mode is enabled
		if tissue == "user_provided_files" or tissue == "user_provided_urls":
			reg_element_counts_str = "Identified regions in {specified_tissue_type}:\nPutative Promoters\n{active_promoter_count} regions in an active promoter\n{bivalent_promoter_count} regions in a bivalent promoter\n{silenced_promoter_count} regions in a silenced promoter\n{unclassified_open_chromatin_ctcf_within_range_of_tss_count} regions in a CTCF binding site within an open region of chromatin, within {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n{unclassified_closed_chromatin_ctcf_within_range_of_tss_count} regions in a CTCF binding site within a closed region of chromatin, within {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n{unclassified_open_chromatin_within_range_of_tss_count} regions within an open region of chromatin (non-CTCF binding site), within {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n{unclassified_closed_chromatin_within_range_of_tss_count} regions within a closed region of chromatin (non-CTCF binding site), within {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n\nPutative Enhancers\n{active_enhancer_count} regions in an active enhancer\n{poised_enhancer_count} regions in a poised enhancer\n{primed_enhancer_count} in a primed enhancer\n{unclassified_open_chromatin_ctcf_beyond_range_of_tss_count} regions in a CTCF binding site within an open region of chromatin, beyond {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n{unclassified_closed_chromatin_ctcf_beyond_range_of_tss_count} regions in a CTCF binding site within a closed region of chromatin, beyond {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n{unclassified_closed_chromatin_beyond_range_of_tss_count} regions within a closed region of chromatin (non-CTCF binding site), beyond {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n".format(specified_tissue_type = user_tissue_names_array[user_files_index], active_promoter_count = list(map(lambda x: str(x).startswith("active_promoter"), merged_df[last_merged_df_column])).count(True), bivalent_promoter_count = list(map(lambda x: str(x).startswith("bivalent_promoter"), merged_df[last_merged_df_column])).count(True), silenced_promoter_count = list(map(lambda x: str(x).startswith("silenced_promoter"), merged_df[last_merged_df_column])).count(True), unclassified_open_chromatin_ctcf_within_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_open_chromatin;ctcf_binding_site;within_"), merged_df[last_merged_df_column])).count(True), unclassified_closed_chromatin_ctcf_within_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_closed_chromatin;ctcf_binding_site;within_"), merged_df[last_merged_df_column])).count(True), unclassified_open_chromatin_within_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_open_chromatin_within_"), merged_df[last_merged_df_column])).count(True), unclassified_closed_chromatin_within_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_closed_chromatin_within_"), merged_df[last_merged_df_column])).count(True), active_enhancer_count = list(map(lambda x: str(x).startswith("active_enhancer"), merged_df[last_merged_df_column])).count(True), poised_enhancer_count = list(map(lambda x: str(x).startswith("poised_enhancer"), merged_df[last_merged_df_column])).count(True), primed_enhancer_count = list(map(lambda x: str(x).startswith("primed_enhancer"), merged_df[last_merged_df_column])).count(True), unclassified_open_chromatin_ctcf_beyond_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_open_chromatin;ctcf_binding_site;beyond_"), merged_df[last_merged_df_column])).count(True), unclassified_closed_chromatin_ctcf_beyond_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_closed_chromatin;ctcf_binding_site;beyond_"), merged_df[last_merged_df_column])).count(True), unclassified_open_chromatin_beyond_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_closed_chromatin_beyond_"), merged_df[last_merged_df_column])).count(True), unclassified_closed_chromatin_beyond_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_closed_chromatin_beyond_"), merged_df[last_merged_df_column])).count(True), tss_dist_upstream = str(args.tss_distance_upstream), tss_dist_downstream = str(args.tss_distance_downstream))
			if args.write_summary:
				with open((user_tissue_names_array[user_files_index] + ".func_anno_summary.txt"), "w") as sum_file:
					sum_file.write(reg_element_counts_str)
			if args.verbose:
				print(reg_element_counts_str)
		else:
			reg_element_counts_str = "Identified regions in {specified_tissue_type}:\nPutative Promoters\n{active_promoter_count} regions in an active promoter\n{bivalent_promoter_count} regions in a bivalent promoter\n{silenced_promoter_count} regions in a silenced promoter\n{unclassified_open_chromatin_ctcf_within_range_of_tss_count} regions in a CTCF binding site within an open region of chromatin, within {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n{unclassified_closed_chromatin_ctcf_within_range_of_tss_count} regions in a CTCF binding site within a closed region of chromatin, within {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n{unclassified_open_chromatin_within_range_of_tss_count} regions within an open region of chromatin (non-CTCF binding site), within {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n{unclassified_closed_chromatin_within_range_of_tss_count} regions within a closed region of chromatin (non-CTCF binding site), within {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n\nPutative Enhancers\n{active_enhancer_count} regions in an active enhancer\n{poised_enhancer_count} regions in a poised enhancer\n{primed_enhancer_count} in a primed enhancer\n{unclassified_open_chromatin_ctcf_beyond_range_of_tss_count} regions in a CTCF binding site within an open region of chromatin, beyond {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n{unclassified_closed_chromatin_ctcf_beyond_range_of_tss_count} regions in a CTCF binding site within a closed region of chromatin, beyond {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n{unclassified_open_chromatin_beyond_range_of_tss_count} regions within an open region of chromatin (non-CTCF binding site), beyond {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n{unclassified_closed_chromatin_beyond_range_of_tss_count} regions within a closed region of chromatin (non-CTCF binding site), beyond {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n".format(specified_tissue_type = tissue, active_promoter_count = list(map(lambda x: str(x).startswith("active_promoter"), merged_df[last_merged_df_column])).count(True), bivalent_promoter_count = list(map(lambda x: str(x).startswith("bivalent_promoter"), merged_df[last_merged_df_column])).count(True), silenced_promoter_count = list(map(lambda x: str(x).startswith("silenced_promoter"), merged_df[last_merged_df_column])).count(True), unclassified_open_chromatin_ctcf_within_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_open_chromatin;ctcf_binding_site;within_"), merged_df[last_merged_df_column])).count(True), unclassified_closed_chromatin_ctcf_within_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_closed_chromatin;ctcf_binding_site;within_"), merged_df[last_merged_df_column])).count(True), unclassified_open_chromatin_within_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_open_chromatin_within_"), merged_df[last_merged_df_column])).count(True), unclassified_closed_chromatin_within_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_closed_chromatin_within_"), merged_df[last_merged_df_column])).count(True), active_enhancer_count = list(map(lambda x: str(x).startswith("active_enhancer"), merged_df[last_merged_df_column])).count(True), poised_enhancer_count = list(map(lambda x: str(x).startswith("poised_enhancer"), merged_df[last_merged_df_column])).count(True), primed_enhancer_count = list(map(lambda x: str(x).startswith("primed_enhancer"), merged_df[last_merged_df_column])).count(True), unclassified_open_chromatin_ctcf_beyond_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_open_chromatin;ctcf_binding_site;beyond_"), merged_df[last_merged_df_column])).count(True), unclassified_closed_chromatin_ctcf_beyond_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_closed_chromatin;ctcf_binding_site;beyond_"), merged_df[last_merged_df_column])).count(True), unclassified_open_chromatin_beyond_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_open_chromatin_beyond_"), merged_df[last_merged_df_column])).count(True), unclassified_closed_chromatin_beyond_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_closed_chromatin_beyond_"), merged_df[last_merged_df_column])).count(True), tss_dist_upstream = str(args.tss_distance_upstream), tss_dist_downstream = str(args.tss_distance_downstream))
			if args.write_summary:
				with open((tissue + ".func_anno_summary.txt"), "w") as sum_file:
					sum_file.write(reg_element_counts_str)
			if args.verbose:
					print(reg_element_counts_str)
=======

	###Compare the input file to the reference bed files, starting with the modified TSS bed
	##Do the peaks fall within boundary of a TSS?
	overlaps_tss = bed_input.intersect(bed_tss, u = True)
	overlaps_tss = overlaps_tss.sort()
	#print(overlaps_tss.count())
	no_tss = bed_input.intersect(bed_tss, v = True)
	overlaps_tss = overlaps_tss.sort()
	#print(no_tss.count())

	#Classifying putative promoters
	overlaps_4me3 = overlaps_tss.intersect(bed_4me3, u = True)
	overlaps_4me3 = overlaps_4me3.sort()
	#print(overlaps_4me3.count())
	no_4me3 = overlaps_tss.intersect(bed_4me3, v = True)
	no_4me3 = no_4me3.sort()
	#print(no_4me3.count())
	
	active_promoter = overlaps_4me3.intersect(bed_27me3, v = True)
	active_promoter = active_promoter.sort()
	#print(active_promoter.count())
 
	bivalent_promoter = overlaps_4me3.intersect(bed_27me3, u = True)
	bivalent_promoter = bivalent_promoter.sort()
	#print(bivalent_promoter.count())
 	
	silenced_promoter = no_4me3.intersect(bed_27me3, u = True)
	silenced_promoter = silenced_promoter.sort()
	#print(silenced_promoter.count())

	###Further classify unclassified promoter-like regions
	unclassified_overlaps_tss = overlaps_tss.intersect(active_promoter, v = True)
	unclassified_overlaps_tss = unclassified_overlaps_tss.intersect(bivalent_promoter, v = True)
	unclassified_overlaps_tss = unclassified_overlaps_tss.intersect(silenced_promoter, v = True)
	#print(unclassified_overlaps_tss.count())
	tss_overlaps_ctcf = unclassified_overlaps_tss.intersect(bed_ctcf, u = True)
	tss_no_ctcf = unclassified_overlaps_tss.intersect(bed_ctcf, v = True)

	ctcf_open_within_tss = tss_overlaps_ctcf.intersect(bed_dnase, u = True)
	ctcf_open_within_tss = ctcf_open_within_tss.sort()

	ctcf_closed_within_tss = tss_overlaps_ctcf.intersect(bed_dnase, v = True)
	ctcf_closed_within_tss = ctcf_closed_within_tss.sort()
	
	no_ctcf_open_within_tss = tss_no_ctcf.intersect(bed_dnase, u = True)
	no_ctcf_open_within_tss = no_ctcf_open_within_tss.sort()
	
	no_ctcf_closed_within_tss = tss_no_ctcf.intersect(bed_dnase, v = True)
	no_ctcf_closed_within_tss = no_ctcf_closed_within_tss.sort()
	
	#Classifying putative enhancers
	overlaps_4me1 = no_tss.intersect(bed_4me1, u = True)

	active_enhancer = overlaps_4me1.intersect(bed_27ac, u = True)
	active_enhancer = active_enhancer.sort()

	poised_enhancer = overlaps_4me1.intersect(bed_27me3, u = True)
	poised_enhancer = poised_enhancer.sort()

	primed_enhancer = overlaps_4me1.intersect(active_enhancer, v = True)
	primed_enhancer = primed_enhancer.intersect(poised_enhancer, v = True)
	primed_enhancer = primed_enhancer.sort()

	###Further classify unclassified enhancer-like regions
	unclassified_no_tss = no_tss.intersect(active_enhancer, v = True)
	unclassified_no_tss = no_tss.intersect(poised_enhancer, v = True)
	unclassified_no_tss = no_tss.intersect(primed_enhancer, v = True)
	
	no_tss_overlaps_ctcf = unclassified_no_tss.intersect(bed_ctcf, u = True)
	no_tss_no_ctcf = unclassified_no_tss.intersect(bed_ctcf, v = True)

	ctcf_open_no_tss = no_tss_overlaps_ctcf.intersect(bed_dnase, u = True)
	ctcf_open_no_tss = ctcf_open_no_tss.sort()

	ctcf_closed_no_tss = no_tss_overlaps_ctcf.intersect(bed_dnase, v = True)
	ctcf_closed_no_tss = ctcf_closed_no_tss.sort()

	no_ctcf_open_no_tss = no_tss_no_ctcf.intersect(bed_dnase, u = True)
	no_ctcf_open_no_tss = no_ctcf_open_no_tss.sort()

	no_ctcf_closed_no_tss = no_tss_no_ctcf.intersect(bed_dnase, v = True)
	no_ctcf_closed_no_tss = no_ctcf_closed_no_tss.sort()

	#If they have one or more coordinate entries, convert each of the classifier variables to a pandas DataFrame object so that we can append the classifications as an additional column
	if active_promoter.count() > 0:
		active_promoter_df = pybedtools.BedTool.to_dataframe(active_promoter)
		active_promoter_df["region_classification"] = "active_promoter"
	else:
		active_promoter_df = pd.DataFrame()
	
	if bivalent_promoter.count() > 0:
		bivalent_promoter_df = pybedtools.BedTool.to_dataframe(bivalent_promoter)
		bivalent_promoter_df["region_classification"] = "bivalent_promoter"
	else:
		bivalent_promoter_df = pd.DataFrame()

	if silenced_promoter.count() > 0:
		silenced_promoter_df = pybedtools.BedTool.to_dataframe(silenced_promoter)
		silenced_promoter_df["region_classification"] = "silenced_promoter"
	else:
		silenced_promoter_df = pd.DataFrame()

	if ctcf_open_within_tss.count() > 0:
		ctcf_open_within_tss_df = pybedtools.BedTool.to_dataframe(ctcf_open_within_tss)
		ctcf_open_within_tss_df["region_classification"] = ("unclassified_open_chromatin;ctcf_binding_site;within_" + str(args.tss_distance_upstream) + "_bp_upstream_" + str(args.tss_distance_downstream) + "_bp_downstream_of_tss")
	else:
		ctcf_open_within_tss_df = pd.DataFrame()

	if ctcf_closed_within_tss.count() > 0:
		ctcf_closed_within_tss_df = pybedtools.BedTool.to_dataframe(ctcf_closed_within_tss)
		ctcf_closed_within_tss_df["region_classification"] = ("unclassified_closed_chromatin;ctcf_binding_site;within_" + str(args.tss_distance_upstream) + "_bp_upstream_" + str(args.tss_distance_downstream) + "_bp_downstream_of_tss")
	else:
		ctcf_closed_within_tss_df = pd.DataFrame()
	
	if no_ctcf_open_within_tss.count() > 0:
		no_ctcf_open_within_tss_df = pybedtools.BedTool.to_dataframe(no_ctcf_open_within_tss)
		no_ctcf_open_within_tss_df["region_classification"] = ("unclassified_open_chromatin_within_" + str(args.tss_distance_upstream) + "_bp_upstream_" + str(args.tss_distance_downstream) + "_bp_downstream_of_tss")
	else:
		no_ctcf_open_within_tss_df = pd.DataFrame()
	
	if no_ctcf_closed_within_tss.count() > 0:
		no_ctcf_closed_within_tss_df = pybedtools.BedTool.to_dataframe(no_ctcf_closed_within_tss)
		no_ctcf_closed_within_tss_df["region_classification"] = ("unclassified_closed_chromatin_within_" + str(args.tss_distance_upstream) + "_bp_upstream_" + str(args.tss_distance_downstream) + "_bp_downstream_of_tss")
	else:
		no_ctcf_closed_within_tss_df = pd.DataFrame()

	if active_enhancer.count() > 0:
		active_enhancer_df = pybedtools.BedTool.to_dataframe(active_enhancer)
		active_enhancer_df["region_classification"] = "active_enhancer"
	else:
		active_enhancer_df = pd.DataFrame()

	if poised_enhancer.count() > 0:
		poised_enhancer_df = pybedtools.BedTool.to_dataframe(poised_enhancer)
		poised_enhancer_df["region_classification"] = "poised_enhancer"
	else:
		poised_enhancer_df = pd.DataFrame()

	if primed_enhancer.count() > 0:
		primed_enhancer_df = pybedtools.BedTool.to_dataframe(primed_enhancer)
		primed_enhancer_df["region_classification"] = "primed_enhancer"
	else:
		primed_enhancer_df = pd.DataFrame()

	if ctcf_open_no_tss.count() > 0:
		ctcf_open_no_tss_df = pybedtools.BedTool.to_dataframe(ctcf_open_no_tss)
		ctcf_open_no_tss_df["region_classification"] = ("unclassified_open_chromatin;ctcf_binding_site;beyond_" + str(args.tss_distance_upstream) + "_bp_upstream_" + str(args.tss_distance_downstream) + "_bp_downstream_of_tss")
	else:
		ctcf_open_no_tss_df = pd.DataFrame()
	
	if ctcf_closed_no_tss.count() > 0:
		ctcf_closed_no_tss_df = pybedtools.BedTool.to_dataframe(ctcf_closed_no_tss)
		ctcf_closed_no_tss_df["region_classification"] = ("unclassified_closed_chromatin;ctcf_binding_site;beyond_" + str(args.tss_distance_upstream) + "_bp_upstream_" + str(args.tss_distance_downstream) + "_bp_downstream_of_tss")
	else:
		ctcf_closed_no_tss_df = pd.DataFrame()
	
	if no_ctcf_open_no_tss.count() > 0:
		no_ctcf_open_no_tss_df = pybedtools.BedTool.to_dataframe(no_ctcf_open_no_tss)
		no_ctcf_open_no_tss_df["region_classification"] = ("unclassified_open_chromatin_beyond_" + str(args.tss_distance_upstream) + "_bp_upstream_" + str(args.tss_distance_downstream) + "_bp_downstream_of_tss")
	else:
		no_ctcf_open_no_tss_df = pd.DataFrame()
	
	if no_ctcf_closed_no_tss.count() > 0:
		no_ctcf_closed_no_tss_df = pybedtools.BedTool.to_dataframe(no_ctcf_closed_no_tss)
		no_ctcf_closed_no_tss_df["region_classification"] = ("unclassified_closed_chromatin_beyond_" + str(args.tss_distance_upstream) + "_bp_upstream_" + str(args.tss_distance_downstream) + "_bp_downstream_of_tss")
	else:
		no_ctcf_closed_no_tss_df = pd.DataFrame()
	
	#Concatenate all of the pandas DataFrame objects into one, remove all duplicate lines, and sort by coordinate
	all_data_frames = [active_promoter_df, bivalent_promoter_df, silenced_promoter_df, ctcf_open_within_tss_df, ctcf_closed_within_tss_df, no_ctcf_open_within_tss_df, no_ctcf_closed_within_tss_df, active_enhancer_df, poised_enhancer_df, primed_enhancer_df, ctcf_open_no_tss_df, ctcf_closed_no_tss_df, no_ctcf_open_no_tss_df, no_ctcf_closed_no_tss_df]

	concatenated_df = pd.concat(all_data_frames)
	concatenated_df = concatenated_df.drop_duplicates()
	concatenated_df.iloc[:,0] = concatenated_df.iloc[:,0].str[3:]
	col_names = ["chr", "start", "end", "func_anno"]
	concatenated_df.columns = col_names
	concatenated_df = concatenated_df.sort_values(by = [col_names[0], col_names[1], col_names[2]])

	#For rows that are completely identical (all entries except for the regulatory annotation), the rows are collapsed into one, with the multiple annotations concatenated in a comma-separated string)
	n_cols = len(concatenated_df.columns)
	names_array = []
	counter = 1
	for i in range(0, n_cols):
		colname = "col_" + str(counter)
		names_array.append(colname)
		counter += 1
	concatenated_df.columns = names_array
	last_col = names_array[-1]
	concatenated_df[last_col] = concatenated_df.groupby(names_array[0:-1])[last_col].transform(lambda x: ','.join(x))
	concatenated_df = concatenated_df.drop_duplicates()

	#If a coordinate overlaps with multiple regulatory elements and the user specifically wishes to keep only one, then we will keep only the most significant (i.e. active_promoter over unclassified_within_2kb_of_tss)
	if args.no_multianno:
		from operator import itemgetter
		annotations = concatenated_df.iloc[:,-1]
		annotations = annotations.str.split(',')
		firsts = list(map(itemgetter(0), annotations))
		concatenated_df.drop(concatenated_df.columns[-1], axis = 1, inplace = True)
		concatenated_df[names_array[-1]] = firsts
	
	#Merge the original input file with the functional annotations from the concatenated data frame
	concatenated_df["col_1"] = concatenated_df["col_1"].astype(str)
	concatenated_df["col_2"] = concatenated_df["col_2"].astype(str)
	concatenated_df["col_3"] = concatenated_df["col_3"].astype(str)
	
	file_input_colnames = file_input.columns
	chr_file_input = file_input_colnames[chr_col]
	start_file_input = file_input_colnames[start_col]
	end_file_input = file_input_colnames[end_col]
	
	file_input[chr_file_input] = file_input[chr_file_input].astype(str)
	file_input[start_file_input] = file_input[start_file_input].astype(str)
	file_input[end_file_input] = file_input[end_file_input].astype(str)
	
	merged_df = file_input.merge(concatenated_df, how = "left", left_on = [chr_file_input, start_file_input, end_file_input], right_on = ["col_1", "col_2", "col_3"])
	merged_df = merged_df.drop(columns = ["col_1", "col_2", "col_3"])
	
	#Rename the functional annotation column with the tissue of interest
	n_cols_merged_df = len(merged_df.columns)
	names_array_merged_df = []
	for i in range(0, n_cols_merged_df):
		colname = merged_df.columns[i]
		names_array_merged_df.append(colname)
		
	if tissue == "user_provided_files" or tissue == "user_provided_urls":
		func_names = user_tissue_names_array[user_files_index].upper() + "_FUNCTIONAL_ANNOTATION"
		names_array_merged_df[-1] = func_names
		user_files_index += 1
	else:
		func_names = tissue.upper() + "_FUNCTIONAL_ANNOTATION"
		names_array_merged_df[-1] = func_names
		
	merged_df.columns = names_array_merged_df
	
	last_merged_df_column = names_array_merged_df[-1]
	anno_df[func_names] = merged_df[func_names]
	
	#Write out the number of peaks/coordinates in each category. Also print this statement to the console if verbose mode is enabled
	if tissue == "user_provided_files" or tissue == "user_provided_urls":
		reg_element_counts_str = "Identified regions in {specified_tissue_type}:\nPutative Promoters\n{active_promoter_count} regions in an active promoter\n{bivalent_promoter_count} regions in a bivalent promoter\n{silenced_promoter_count} regions in a silenced promoter\n{unclassified_open_chromatin_ctcf_within_range_of_tss_count} regions in a CTCF binding site within an open region of chromatin, within {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n{unclassified_closed_chromatin_ctcf_within_range_of_tss_count} regions in a CTCF binding site within a closed region of chromatin, within {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n{unclassified_open_chromatin_within_range_of_tss_count} regions within an open region of chromatin (non-CTCF binding site), within {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n{unclassified_closed_chromatin_within_range_of_tss_count} regions within a closed region of chromatin (non-CTCF binding site), within {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n\nPutative Enhancers\n{active_enhancer_count} regions in an active enhancer\n{poised_enhancer_count} regions in a poised enhancer\n{primed_enhancer_count} in a primed enhancer\n{unclassified_open_chromatin_ctcf_beyond_range_of_tss_count} regions in a CTCF binding site within an open region of chromatin, beyond {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n{unclassified_closed_chromatin_ctcf_beyond_range_of_tss_count} regions in a CTCF binding site within a closed region of chromatin, beyond {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n{unclassified_closed_chromatin_beyond_range_of_tss_count} regions within a closed region of chromatin (non-CTCF binding site), beyond {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n".format(specified_tissue_type = user_tissue_names_array[user_files_index], active_promoter_count = list(map(lambda x: str(x).startswith("active_promoter"), merged_df[last_merged_df_column])).count(True), bivalent_promoter_count = list(map(lambda x: str(x).startswith("bivalent_promoter"), merged_df[last_merged_df_column])).count(True), silenced_promoter_count = list(map(lambda x: str(x).startswith("silenced_promoter"), merged_df[last_merged_df_column])).count(True), unclassified_open_chromatin_ctcf_within_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_open_chromatin;ctcf_binding_site;within_"), merged_df[last_merged_df_column])).count(True), unclassified_closed_chromatin_ctcf_within_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_closed_chromatin;ctcf_binding_site;within_"), merged_df[last_merged_df_column])).count(True), unclassified_open_chromatin_within_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_open_chromatin_within_"), merged_df[last_merged_df_column])).count(True), unclassified_closed_chromatin_within_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_closed_chromatin_within_"), merged_df[last_merged_df_column])).count(True), active_enhancer_count = list(map(lambda x: str(x).startswith("active_enhancer"), merged_df[last_merged_df_column])).count(True), poised_enhancer_count = list(map(lambda x: str(x).startswith("poised_enhancer"), merged_df[last_merged_df_column])).count(True), primed_enhancer_count = list(map(lambda x: str(x).startswith("primed_enhancer"), merged_df[last_merged_df_column])).count(True), unclassified_open_chromatin_ctcf_beyond_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_open_chromatin;ctcf_binding_site;beyond_"), merged_df[last_merged_df_column])).count(True), unclassified_closed_chromatin_ctcf_beyond_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_closed_chromatin;ctcf_binding_site;beyond_"), merged_df[last_merged_df_column])).count(True), unclassified_open_chromatin_beyond_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_closed_chromatin_beyond_"), merged_df[last_merged_df_column])).count(True), unclassified_closed_chromatin_beyond_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_closed_chromatin_beyond_"), merged_df[last_merged_df_column])).count(True), tss_dist_upstream = str(args.tss_distance_upstream), tss_dist_downstream = str(args.tss_distance_downstream))
		with open((user_tissue_names_array[user_files_index] + ".func_anno_summary.txt"), "w") as sum_file:
			sum_file.write(reg_element_counts_str)
		if args.verbose:
			print(reg_element_counts_str)
	else:
		reg_element_counts_str = "Identified regions in {specified_tissue_type}:\nPutative Promoters\n{active_promoter_count} regions in an active promoter\n{bivalent_promoter_count} regions in a bivalent promoter\n{silenced_promoter_count} regions in a silenced promoter\n{unclassified_open_chromatin_ctcf_within_range_of_tss_count} regions in a CTCF binding site within an open region of chromatin, within {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n{unclassified_closed_chromatin_ctcf_within_range_of_tss_count} regions in a CTCF binding site within a closed region of chromatin, within {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n{unclassified_open_chromatin_within_range_of_tss_count} regions within an open region of chromatin (non-CTCF binding site), within {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n{unclassified_closed_chromatin_within_range_of_tss_count} regions within a closed region of chromatin (non-CTCF binding site), within {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n\nPutative Enhancers\n{active_enhancer_count} regions in an active enhancer\n{poised_enhancer_count} regions in a poised enhancer\n{primed_enhancer_count} in a primed enhancer\n{unclassified_open_chromatin_ctcf_beyond_range_of_tss_count} regions in a CTCF binding site within an open region of chromatin, beyond {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n{unclassified_closed_chromatin_ctcf_beyond_range_of_tss_count} regions in a CTCF binding site within a closed region of chromatin, beyond {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n{unclassified_open_chromatin_beyond_range_of_tss_count} regions within an open region of chromatin (non-CTCF binding site), beyond {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n{unclassified_closed_chromatin_beyond_range_of_tss_count} regions within a closed region of chromatin (non-CTCF binding site), beyond {tss_dist_upstream} bp upstream or {tss_dist_downstream} bp downstream of a TSS\n".format(specified_tissue_type = tissue, active_promoter_count = list(map(lambda x: str(x).startswith("active_promoter"), merged_df[last_merged_df_column])).count(True), bivalent_promoter_count = list(map(lambda x: str(x).startswith("bivalent_promoter"), merged_df[last_merged_df_column])).count(True), silenced_promoter_count = list(map(lambda x: str(x).startswith("silenced_promoter"), merged_df[last_merged_df_column])).count(True), unclassified_open_chromatin_ctcf_within_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_open_chromatin;ctcf_binding_site;within_"), merged_df[last_merged_df_column])).count(True), unclassified_closed_chromatin_ctcf_within_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_closed_chromatin;ctcf_binding_site;within_"), merged_df[last_merged_df_column])).count(True), unclassified_open_chromatin_within_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_open_chromatin_within_"), merged_df[last_merged_df_column])).count(True), unclassified_closed_chromatin_within_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_closed_chromatin_within_"), merged_df[last_merged_df_column])).count(True), active_enhancer_count = list(map(lambda x: str(x).startswith("active_enhancer"), merged_df[last_merged_df_column])).count(True), poised_enhancer_count = list(map(lambda x: str(x).startswith("poised_enhancer"), merged_df[last_merged_df_column])).count(True), primed_enhancer_count = list(map(lambda x: str(x).startswith("primed_enhancer"), merged_df[last_merged_df_column])).count(True), unclassified_open_chromatin_ctcf_beyond_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_open_chromatin;ctcf_binding_site;beyond_"), merged_df[last_merged_df_column])).count(True), unclassified_closed_chromatin_ctcf_beyond_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_closed_chromatin;ctcf_binding_site;beyond_"), merged_df[last_merged_df_column])).count(True), unclassified_open_chromatin_beyond_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_open_chromatin_beyond_"), merged_df[last_merged_df_column])).count(True), unclassified_closed_chromatin_beyond_range_of_tss_count = list(map(lambda x: str(x).startswith("unclassified_closed_chromatin_beyond_"), merged_df[last_merged_df_column])).count(True), tss_dist_upstream = str(args.tss_distance_upstream), tss_dist_downstream = str(args.tss_distance_downstream))
		with open((tissue + ".func_anno_summary.txt"), "w") as sum_file:
			sum_file.write(reg_element_counts_str)
		if args.verbose:
				print(reg_element_counts_str)
>>>>>>> 3283f061327aed512684aab1c4b30a7a1c01e0ad

#Output the final annotated file
merged_df = merged_df.rename({last_merged_df_column: "dropcol"}, axis=1)
merged_df = merged_df.drop(columns = ["dropcol"])
anno_df_names = anno_df.columns
merged_df[anno_df_names] = anno_df

if args.verbose:
		print("Preparing output file ({output_file})...".format(output_file = args.output))
		print("\n")
<<<<<<< HEAD
if args.write_anno_only:
	merged_df = merged_df.iloc[:,-1]
	merged_df.to_csv(args.output, sep = args.separator, index = False)
elif args.input_header:
	merged_df.to_csv(args.output, sep = args.separator, index = False)
else:
	merged_df.to_csv(args.output, sep = args.separator, header = False, index = False)
=======
if args.input_header:
	merged_df.to_csv(args.output, sep = "\t", index = False)
else:
	merged_df.to_csv(args.output, sep = "\t", header = False, index = False)
>>>>>>> 3283f061327aed512684aab1c4b30a7a1c01e0ad
