###Developed by Michael J Betti, April 2021, updated 3 June 2021
__author__ = "Michael J Betti"
__copyright__ = "Copyright 2021, Michael J Betti"
__license__ = "BSD"
__maintainer__ = "Michael J Betti"
__email__ = "mjbetti3@gmail.com"
__status__ = "Development"

import os, sys, requests, argparse, pybedtools, pandas as pd

parser = argparse.ArgumentParser(add_help = True)
parser.add_argument("-i", "--input", type = str, required = True, help = "the input bed file (required)")
parser.add_argument("-g", "--ref_genome", type = str, required = True, help = "the human reference genome build on which the input coordinates are based (required) (valid options: GRCh38/hg38 and GRCh37/hg19)")
parser.add_argument("-t", "--tissue", type = str, required = True, help = "the tissue of interest (required) (valid options: Adipose, Adrenal_gland, Artery, Blood, Breast, Cultured_fibroblast, EBV_transformed_lymphocyte, ES, Esophagus_muscularis_mucosa, Esophagus_squamous_epithelium, Heart, Intestine, iPS, Kidney, Liver, Lung, Neuron, Ovary, Pancreas, Prostate, Skeletal_muscle, Skin, Spleen, Stomach, Testis, Thyroid, Uterus, Vagina)")
parser.add_argument("-ud", "--tss_distance_upstream", type = int, required = False, help = "the upstream boundary distance from a TSS", default = 2000)
parser.add_argument("-dd", "--tss_distance_downstream", type = int, required = False, help = "the downstream boundary distance from a TSS", default = 2000)
parser.add_argument("-o", "--output", type = str, required = False, help = "the name of the output file", default = "out.bed")
parser.add_argument("--no_multianno", required = False, help = "if a coordinate overlaps with multiple regions, keep the most significant occurance", action = "store_true")
parser.add_argument("--bed_cols", type = str, required = False, help = "if the input is not in traditional UCSC BED format, specify the column numbers of chr, start, and end separated by commas", default = "1,2,3")
parser.add_argument("--input_header", required = False, help = "indicate whether the input file has a header (indicate true or false)", action = "store_true")
parser.add_argument("-v", "--verbose", required = False, help = "return logging as terminal output", action = "store_true")
args = parser.parse_args()

#Check that required arguments are specified
assert args.input, "Must specify input file (-i, --input)"
assert args.ref_genome, "Must specify reference genome build (-g, --ref_genome)"
assert args.tissue, "Must specify tissue type (-t, --tissue)"

#Check that specified tissue type is one of the 28 valid options
assert args.tissue.lower() == "adipose" or args.tissue.lower() == "adrenal_gland" or args.tissue.lower() == "artery" or args.tissue.lower() == "blood" or args.tissue.lower() == "breast" or args.tissue.lower() == "cultured_fibroblast" or args.tissue.lower() == "ebv_transformed_lymphocyte" or args.tissue.lower() == "es" or args.tissue.lower() == "esophagus_muscularis_mucosa" or args.tissue.lower() == "esophagus_squamous_epithelium" or args.tissue.lower() == "heart" or args.tissue.lower() == "intestine" or args.tissue.lower() == "ips" or args.tissue.lower() == "kidney" or args.tissue.lower() == "liver" or args.tissue.lower() == "lung" or args.tissue.lower() == "neuron" or args.tissue.lower() == "ovary" or args.tissue.lower() == "pancreas" or args.tissue.lower() == "prostate" or args.tissue.lower() == "skeletal_muscle" or args.tissue.lower() == "skin" or args.tissue.lower() == "spleen" or args.tissue.lower() == "stomach" or args.tissue.lower() == "testis" or args.tissue.lower() == "thyroid" or args.tissue.lower() == "uterus" or args.tissue.lower() == "vagina", "Tissue type must be one of the 28 valid options (Adipose, Adrenal_gland, Artery, Blood, Breast, Cultured_fibroblast, EBV_transformed_lymphocyte, ES, Esophagus_muscularis_mucosa, Esophagus_squamous_epithelium, Heart, Intestine, iPS, Kidney, Liver, Lung, Neuron, Ovary, Pancreas, Prostate, Skeletal_muscle, Skin, Spleen, Stomach, Testis, Thyroid, Uterus, Vagina)"

#Download the appropriate reference files based on the specified genome build and tissue arguments
def download_ref(url_ref, out_name):
	out_path = os.path.join("ref_files", out_name)
	r = requests.get(url_ref, allow_redirects = True)
	open(out_path, 'wb').write(r.content)

if args.verbose:
	print("Classifying {tissue_type} cis-regulatory regions in {input_file}...".format(input_file = args.input, tissue_type = args.tissue))
	print("\n")
if not os.path.exists("ref_files"):
	os.mkdir("ref_files")

if args.verbose:
	print("Downloading {ref_genome} TSS coordinates and {tissue_type} histone ChIP-seq bed files...".format(ref_genome = args.ref_genome, tissue_type = args.tissue))
	print("\n")

#GRCh38/hg38
if args.ref_genome.lower() == "hg38" or args.ref_genome.lower() == "grch38":
	download_ref("http://reftss.clst.riken.jp/datafiles/3.1/human/refTSS_v3.1_human_coordinate.hg38.bed.gz", "refTSS_v3.1_human_coordinate.hg38.bed.gz")
	
	#Adipose (Homo sapiens subcutaneous abdominal adipose tissue tissue nuclear fraction female adult (49 years))
	if args.tissue.lower() == "adipose":
		download_ref("https://www.encodeproject.org/files/ENCFF534VNL/@@download/ENCFF534VNL.bed.gz", "adipose_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF591SLF/@@download/ENCFF591SLF.bed.gz", "adipose_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF658NHX/@@download/ENCFF658NHX.bed.gz", "adipose_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF435MSC/@@download/ENCFF435MSC.bed.gz", "adipose_27me3_hg38.bed.gz")	
		download_ref("https://www.encodeproject.org/files/ENCFF075XPG/@@download/ENCFF075XPG.bed.gz", "adipose_36me3_hg38.bed.gz")
		
	#Adrenal gland (Homo sapiens adrenal gland tissue male embryo (97 days))
	if args.tissue.lower() == "adrenal_gland":
		download_ref("https://www.encodeproject.org/files/ENCFF953KZN/@@download/ENCFF953KZN.bed.gz", "adrenal_gland_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF959FFU/@@download/ENCFF959FFU.bed.gz", "adrenal_gland_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF231NZU/@@download/ENCFF231NZU.bed.gz", "adrenal_gland_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF852DFJ/@@download/ENCFF852DFJ.bed.gz", "adrenal_gland_27me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF852DFJ/@@download/ENCFF852DFJ.bed.gz", "adrenal_gland_36me3_hg38.bed.gz")
		
	#Artery (Homo sapiens aorta tissue male adult (34 years))
	if args.tissue.lower() == "artery":
		download_ref("https://www.encodeproject.org/files/ENCFF824OFK/@@download/ENCFF824OFK.bed.gz", "artery_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF876LPD/@@download/ENCFF876LPD.bed.gz", "artery_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF872XZM/@@download/ENCFF872XZM.bed.gz", "artery_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF093PYF/@@download/ENCFF093PYF.bed.gz", "artery_27me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF587NCS/@@download/ENCFF587NCS.bed.gz", "artery_36me3_hg38.bed.gz")

	#Blood (Homo sapiens peripheral blood mononuclear cell male adult (39 years))	
	if args.tissue.lower() == "blood":
		download_ref("https://www.encodeproject.org/files/ENCFF734KLT/@@download/ENCFF734KLT.bed.gz", "blood_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF165VDC/@@download/ENCFF165VDC.bed.gz", "blood_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF844OYX/@@download/ENCFF844OYX.bed.gz", "blood_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF636JEP/@@download/ENCFF636JEP.bed.gz", "blood_27me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF422UFO/@@download/ENCFF422UFO.bed.gz", "blood_36me3_hg38.bed.gz")
	
	#Breast (Homo sapiens breast epithelium tissue female adult (53 years))
	elif args.tissue.lower() == "breast":
		download_ref("https://www.encodeproject.org/files/ENCFF698LQG/@@download/ENCFF698LQG.bed.gz", "breast_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF590DGH/@@download/ENCFF590DGH.bed.gz", "breast_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF507WEL/@@download/ENCFF507WEL.bed.gz", "breast_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF962YZN/@@download/ENCFF962YZN.bed.gz", "breast_27me3_hg38.bed.gz")
		
	#Cultured fibroblast (Homo sapiens IMR-90)
	elif args.tissue.lower() == "cultured_fibroblast":
		download_ref("https://www.encodeproject.org/files/ENCFF611UWF/@@download/ENCFF611UWF.bed.gz", "cultured_fibroblast_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF093NQC/@@download/ENCFF093NQC.bed.gz", "cultured_fibroblast_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF805GNH/@@download/ENCFF805GNH.bed.gz", "cultured_fibroblast_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF336IXL/@@download/ENCFF336IXL.bed.gz", "cultured_fibroblast_27me3_hg38.bed.gz")
		
	#EBV-transformed lymphocyte (Homo sapiens GM12878)
	elif args.tissue.lower() == "ebv_transformed_lymphocyte":
		download_ref("https://www.encodeproject.org/files/ENCFF321BVG/@@download/ENCFF321BVG.bed.gz", "ebv_transformed_lymphocyte_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF998CEU/@@download/ENCFF998CEU.bed.gz", "ebv_transformed_lymphocyte_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF023LTU/@@download/ENCFF023LTU.bed.gz", "ebv_transformed_lymphocyte_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF291DHI/@@download/ENCFF291DHI.bed.gz", "ebv_transformed_lymphocyte_27me3_hg38.bed.gz")
		
	#Embryonic stem cell (Homo sapiens H1)
	elif args.tissue.lower() == "es":
		download_ref("https://www.encodeproject.org/files/ENCFF613QAB/@@download/ENCFF613QAB.bed.gz", "es_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF356MCC/@@download/ENCFF356MCC.bed.gz", "es_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF689CJG/@@download/ENCFF689CJG.bed.gz", "es_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF050CUG/@@download/ENCFF050CUG.bed.gz", "es_27me3_hg38.bed.gz")
	
	#Esophagus muscularis mucosa (Homo sapiens esophagus muscularis mucosa tissue female adult (51 years))
	elif args.tissue.lower() == "esophagus_muscularis_mucosa":
		download_ref("https://www.encodeproject.org/files/ENCFF473ELV/@@download/ENCFF473ELV.bed.gz", "esophagus_muscularis_mucosa_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF154TIT/@@download/ENCFF154TIT.bed.gz", "esophagus_muscularis_mucosa_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF114SOB/@@download/ENCFF114SOB.bed.gz", "esophagus_muscularis_mucosa_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF450WVA/@@download/ENCFF450WVA.bed.gz", "esophagus_muscularis_mucosa_27me3_hg38.bed.gz")
		
	#Esophagus squamous epithelium (Homo sapiens esophagus squamous epithelium tissue female adult (51 years))
	elif args.tissue.lower() == "esophagus_squamous_epithelium":
		download_ref("https://www.encodeproject.org/files/ENCFF828YQS/@@download/ENCFF828YQS.bed.gz", "esophagus_squamous_epithelium_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF470ELM/@@download/ENCFF470ELM.bed.gz", "esophagus_squamous_epithelium_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF673PCP/@@download/ENCFF673PCP.bed.gz", "esophagus_squamous_epithelium_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF905MKP/@@download/ENCFF905MKP.bed.gz", "esophagus_squamous_epithelium_27me3_hg38.bed.gz")
		
	#Heart (Homo sapiens heart left ventricle tissue female adult (53 years))
	elif args.tissue.lower() == "heart":
		download_ref("https://www.encodeproject.org/files/ENCFF194KSV/@@download/ENCFF194KSV.bed.gz", "heart_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF008ODN/@@download/ENCFF008ODN.bed.gz", "heart_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF635VBV/@@download/ENCFF635VBV.bed.gz", "heart_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF242JDO/@@download/ENCFF242JDO.bed.gz", "heart_27me3_hg38.bed.gz")
		
	#Intestine (Homo sapiens sigmoid colon tissue male adult (37 years))
	elif args.tissue.lower() == "intestine":
		download_ref("https://www.encodeproject.org/files/ENCFF052XGA/@@download/ENCFF052XGA.bed.gz", "intestine_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF250ZHO/@@download/ENCFF250ZHO.bed.gz", "intestine_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF111SYA/@@download/ENCFF111SYA.bed.gz", "intestine_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF310IST/@@download/ENCFF310IST.bed.gz", "intestine_27me3_hg38.bed.gz")
	
	#iPS cell (Homo sapiens iPS DF 19.11)
	elif args.tissue.lower() == "ips":
		download_ref("https://www.encodeproject.org/files/ENCFF314LSR/@@download/ENCFF314LSR.bed.gz", "ips_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF450VVJ/@@download/ENCFF450VVJ.bed.gz", "ips_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF037RXP/@@download/ENCFF037RXP.bed.gz", "ips_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF844FEB/@@download/ENCFF844FEB.bed.gz", "ips_27me3_hg38.bed.gz")
		
	#Kidney (Homo sapiens kidney tissue male adult (67 years), Homo sapiens kidney tissue female embryo (120 days) for H3K27me3 - only dataset for kidney)
	elif args.tissue.lower() == "kidney":
		download_ref("https://www.encodeproject.org/files/ENCFF249OJK/@@download/ENCFF249OJK.bed.gz", "kidney_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF430ZMN/@@download/ENCFF430ZMN.bed.gz", "kidney_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF512JNT/@@download/ENCFF512JNT.bed.gz", "kidney_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF775YUI/@@download/ENCFF775YUI.bed.gz", "kidney_27me3_hg38.bed.gz")

	#Liver (Homo sapiens liver tissue male adult (32 years))
	elif args.tissue.lower() == "liver":
		download_ref("https://www.encodeproject.org/files/ENCFF625TWQ/@@download/ENCFF625TWQ.bed.gz", "liver_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF830ODB/@@download/ENCFF830ODB.bed.gz", "liver_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF287VIA/@@download/ENCFF287VIA.bed.gz", "liver_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF287VIA/@@download/ENCFF287VIA.bed.gz", "liver_27me3_hg38.bed.gz")
		
	#Lung (Homo sapiens upper lobe of left lung tissue female adult (51 years))
	elif args.tissue.lower() == "lung":
		download_ref("https://www.encodeproject.org/files/ENCFF926YKP/@@download/ENCFF926YKP.bed.gz", "lung_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF843JPR/@@download/ENCFF843JPR.bed.gz", "lung_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF496RDA/@@download/ENCFF496RDA.bed.gz", "lung_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF248IPH/@@download/ENCFF248IPH.bed.gz", "lung_27me3_hg38.bed.gz")
		
	#Neuron (Homo sapiens SK-N-SH)
	elif args.tissue.lower() == "neuron":
		download_ref("https://www.encodeproject.org/files/ENCFF632FAM/@@download/ENCFF632FAM.bed.gz", "neuron_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF682JYE/@@download/ENCFF682JYE.bed.gz", "neuron_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF138VUT/@@download/ENCFF138VUT.bed.gz", "neuron_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF277NRX/@@download/ENCFF277NRX.bed.gz", "neuron_27me3_hg38.bed.gz")
		
	#Ovary (Homo sapiens ovary tissue female adult (30 years))
	elif args.tissue.lower() == "ovary":
		download_ref("https://www.encodeproject.org/files/ENCFF458XZT/@@download/ENCFF458XZT.bed.gz", "ovary_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF008MLA/@@download/ENCFF008MLA.bed.gz", "ovary_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF328KKO/@@download/ENCFF328KKO.bed.gz", "ovary_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF495DAT/@@download/ENCFF495DAT.bed.gz", "ovary_27me3_hg38.bed.gz")
		
	#Pancreas (Homo sapiens pancreas tissue female adult (41 years))
	elif args.tissue.lower() == "pancreas":
		download_ref("https://www.encodeproject.org/files/ENCFF453GYD/@@download/ENCFF453GYD.bed.gz", "pancreas_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF200ITL/@@download/ENCFF200ITL.bed.gz", "pancreas_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF498SSN/@@download/ENCFF498SSN.bed.gz", "pancreas_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF663DVO/@@download/ENCFF663DVO.bed.gz", "pancreas_27me3_hg38.bed.gz")
		
	#Prostate (Homo sapiens prostate gland tissue male adult (37 years))
	elif args.tissue.lower() == "prostate":
		download_ref("https://www.encodeproject.org/files/ENCFF275KSR/@@download/ENCFF275KSR.bed.gz", "prostate_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF155ZYZ/@@download/ENCFF155ZYZ.bed.gz", "prostate_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF201VZW/@@download/ENCFF201VZW.bed.gz", "prostate_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF198QRW/@@download/ENCFF198QRW.bed.gz", "prostate_27me3_hg38.bed.gz")
		
	#Skeletal muscle (Homo sapiens muscle of leg tissue female embryo (110 days))
	elif args.tissue.lower() == "skeletal_muscle":
		download_ref("https://www.encodeproject.org/files/ENCFF170ZFK/@@download/ENCFF170ZFK.bed.gz", "skeletal_muscle_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF942MOM/@@download/ENCFF942MOM.bed.gz", "skeletal_muscle_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF154ZCF/@@download/ENCFF154ZCF.bed.gz", "skeletal_muscle_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF672MZK/@@download/ENCFF672MZK.bed.gz", "skeletal_muscle_27me3_hg38.bed.gz")

	#Skin (Homo sapiens suprapubic skin tissue male adult (37 years))
	elif args.tissue.lower() == "skin":
		download_ref("https://www.encodeproject.org/files/ENCFF800MDH/@@download/ENCFF800MDH.bed.gz", "skin_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF240CFZ/@@download/ENCFF240CFZ.bed.gz", "skin_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF398EEO/@@download/ENCFF398EEO.bed.gz", "skin_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF686GBZ/@@download/ENCFF686GBZ.bed.gz", "skin_27me3_hg38.bed.gz")
		
	#Spleen (Homo sapiens spleen tissue male adult (37 years))
	elif args.tissue.lower() == "spleen":
		download_ref("https://www.encodeproject.org/files/ENCFF652QYY/@@download/ENCFF652QYY.bed.gz", "spleen_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF401XEH/@@download/ENCFF401XEH.bed.gz", "spleen_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF322YIK/@@download/ENCFF322YIK.bed.gz", "spleen_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF159IED/@@download/ENCFF159IED.bed.gz", "spleen_27me3_hg38.bed.gz")
		
	#Stomach (Homo sapiens stomach tissue male adult (37 years))
	elif args.tissue.lower() == "stomach":
		download_ref("https://www.encodeproject.org/files/ENCFF475EJU/@@download/ENCFF475EJU.bed.gz", "stomach_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF494BWU/@@download/ENCFF494BWU.bed.gz", "stomach_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF978AOD/@@download/ENCFF978AOD.bed.gz", "stomach_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF692HXY/@@download/ENCFF692HXY.bed.gz", "stomach_27me3_hg38.bed.gz")
		
	#Testis (Homo sapiens testis tissue male adult (37 years))
	elif args.tissue.lower() == "testis":
		download_ref("https://www.encodeproject.org/files/ENCFF600MJO/@@download/ENCFF600MJO.bed.gz", "testis_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF612RSR/@@download/ENCFF612RSR.bed.gz", "testis_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF540YUQ/@@download/ENCFF540YUQ.bed.gz", "testis_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF389XVY/@@download/ENCFF389XVY.bed.gz", "testis_27me3_hg38.bed.gz")
		
	#Thyroid (Homo sapiens thyroid gland tissue female adult (51 years))
	elif args.tissue.lower() == "thyroid":
		download_ref("https://www.encodeproject.org/files/ENCFF346IPY/@@download/ENCFF346IPY.bed.gz", "thyroid_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF925FSR/@@download/ENCFF925FSR.bed.gz", "thyroid_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF138WNV/@@download/ENCFF138WNV.bed.gz", "thyroid_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF607NRZ/@@download/ENCFF607NRZ.bed.gz", "thyroid_27me3_hg38.bed.gz")
		
	#Uterus (Homo sapiens uterus tissue female adult (51 years))
	elif args.tissue.lower() == "uterus":
		download_ref("https://www.encodeproject.org/files/ENCFF150BRF/@@download/ENCFF150BRF.bed.gz", "uterus_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF022DOY/@@download/ENCFF022DOY.bed.gz", "uterus_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF108SAD/@@download/ENCFF108SAD.bed.gz", "uterus_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF254GYM/@@download/ENCFF254GYM.bed.gz", "uterus_27me3_hg38.bed.gz")
		
	#Vagina (Homo sapiens vagina tissue female adult (51 years))
	elif args.tissue.lower() == "vagina":
		download_ref("https://www.encodeproject.org/files/ENCFF005TOV/@@download/ENCFF005TOV.bed.gz", "vagina_4me1_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF653IVY/@@download/ENCFF653IVY.bed.gz", "vagina_4me3_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF165SPP/@@download/ENCFF165SPP.bed.gz", "vagina_27ac_hg38.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF866QFX/@@download/ENCFF866QFX.bed.gz", "vagina_27me3_hg38.bed.gz")

###################################################################################
#GRCh37/hg19
elif args.ref_genome.lower() == "hg19" or args.ref_genome.lower() == "grch37":
	download_ref("https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/hg19/refGene_hg19_TSS.bed", "refGene_hg19_TSS.bed")
	
	#Adipose (Homo sapiens subcutaneous abdominal adipose tissue tissue nuclear fraction female adult (49 years))
	if args.tissue.lower() == "adipose":
		download_ref("https://www.encodeproject.org/files/ENCFF258ONX/@@download/ENCFF258ONX.bed.gz", "adipose_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF665CSL/@@download/ENCFF665CSL.bed.gz", "adipose_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF621LOW/@@download/ENCFF621LOW.bed.gz", "adipose_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF464CMC/@@download/ENCFF464CMC.bed.gz", "adipose_27me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF553JCE/@@download/ENCFF553JCE.bed.gz", "adipose_36me3_hg19.bed.gz")
		
	#Adrenal gland (Homo sapiens adrenal gland tissue male embryo (97 days))
	if args.tissue.lower() == "adrenal_gland":
		download_ref("https://www.encodeproject.org/files/ENCFF100QRG/@@download/ENCFF100QRG.bed.gz", "adrenal_gland_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF825DIG/@@download/ENCFF825DIG.bed.gz", "adrenal_gland_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF433MRW/@@download/ENCFF433MRW.bed.gz", "adrenal_gland_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF810FEE/@@download/ENCFF810FEE.bed.gz", "adrenal_gland_27me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF185XYP/@@download/ENCFF185XYP.bed.gz", "adrenal_gland_36me3_hg19.bed.gz")
		
	#Artery (Homo sapiens aorta tissue male adult (34 years))
	if args.tissue.lower() == "artery":
		download_ref("https://www.encodeproject.org/files/ENCFF028YXN/@@download/ENCFF028YXN.bed.gz", "artery_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF870QVN/@@download/ENCFF870QVN.bed.gz", "artery_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF195VBI/@@download/ENCFF195VBI.bed.gz", "artery_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF251MOP/@@download/ENCFF251MOP.bed.gz", "artery_27me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF498YBT/@@download/ENCFF498YBT.bed.gz", "artery_36me3_hg19.bed.gz")
		
	#Blood (Homo sapiens peripheral blood mononuclear cell male adult (39 years))	
	if args.tissue.lower() == "blood":
		download_ref("https://www.encodeproject.org/files/ENCFF553SEK/@@download/ENCFF553SEK.bed.gz", "blood_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF477HMF/@@download/ENCFF477HMF.bed.gz", "blood_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF442KRP/@@download/ENCFF442KRP.bed.gz", "blood_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF987VSV/@@download/ENCFF987VSV.bed.gz", "blood_27me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF068OGE/@@download/ENCFF068OGE.bed.gz", "blood_36me3_hg19.bed.gz")
		
	#Breast (Homo sapiens breast epithelium tissue female adult (53 years))
	elif args.tissue.lower() == "breast":
		download_ref("https://www.encodeproject.org/files/ENCFF336DDM/@@download/ENCFF336DDM.bed.gz", "breast_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF065TIH/@@download/ENCFF065TIH.bed.gz", "breast_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF154XFN/@@download/ENCFF154XFN.bed.gz", "breast_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF291WFP/@@download/ENCFF291WFP.bed.gz", "breast_27me3_hg19.bed.gz")
		
	#Cultured fibroblast (Homo sapiens IMR-90)
	elif args.tissue.lower() == "cultured_fibroblast":
		download_ref("https://www.encodeproject.org/files/ENCFF830JTY/@@download/ENCFF830JTY.bed.gz", "cultured_fibroblast_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF154CUR/@@download/ENCFF154CUR.bed.gz", "cultured_fibroblast_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF678QLP/@@download/ENCFF678QLP.bed.gz", "cultured_fibroblast_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF041RLH/@@download/ENCFF041RLH.bed.gz", "cultured_fibroblasts_27me3_hg19.bed.gz")
	
	#EBV-transformed lymphocyte (Homo sapiens GM12878)
	elif args.tissue.lower() == "ebv_transformed_lymphocyte":
		download_ref("https://www.encodeproject.org/files/ENCFF921LKB/@@download/ENCFF921LKB.bed.gz", "ebv_transformed_lymphocyte_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF295GNH/@@download/ENCFF295GNH.bed.gz", "ebv_transformed_lymphocyte_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF816AHV/@@download/ENCFF816AHV.bed.gz", "ebv_transformed_lymphocyte_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF247VUO/@@download/ENCFF247VUO.bed.gz", "ebv_transformed_lymphocyte_27me3_hg19.bed.gz")
		
	#Embryonic stem cell (Homo sapiens H1)
	elif args.tissue.lower() == "es":
		download_ref("https://www.encodeproject.org/files/ENCFF834NWA/@@download/ENCFF834NWA.bed.gz", "es_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF865TYO/@@download/ENCFF865TYO.bed.gz", "es_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF335JSA/@@download/ENCFF335JSA.bed.gz", "es_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF900RIU/@@download/ENCFF900RIU.bed.gz", "es_27me3_hg19.bed.gz")
		
	#Esophagus muscularis mucosa (Homo sapiens esophagus muscularis mucosa tissue female adult (51 years))
	elif args.tissue.lower() == "esophagus_muscularis_mucosa":
		download_ref("https://www.encodeproject.org/files/ENCFF278OGS/@@download/ENCFF278OGS.bed.gz", "esophagus_muscularis_mucosa_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF821CCT/@@download/ENCFF821CCT.bed.gz", "esophagus_muscularis_mucosa_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF333FTU/@@download/ENCFF333FTU.bed.gz", "esophagus_muscularis_mucosa_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF683CXO/@@download/ENCFF683CXO.bed.gz", "esophagus_muscularis_mucosa_27me3_hg19.bed.gz")
		
	#Esophagus squamous epithelium (Homo sapiens esophagus squamous epithelium tissue female adult (51 years))
	elif args.tissue.lower() == "esophagus_squamous_epithelium":
		download_ref("https://www.encodeproject.org/files/ENCFF693JWR/@@download/ENCFF693JWR.bed.gz", "esophagus_squamous_epithelium_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF659VMQ/@@download/ENCFF659VMQ.bed.gz", "esophagus_squamous_epithelium_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF081YZL/@@download/ENCFF081YZL.bed.gz", "esophagus_squamous_epithelium_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF639ADZ/@@download/ENCFF639ADZ.bed.gz", "esophagus_squamous_epithelium_27me3_hg19.bed.gz")
		
	#Heart (Homo sapiens heart left ventricle tissue female adult (53 years))
	elif args.tissue.lower() == "heart":
		download_ref("https://www.encodeproject.org/files/ENCFF266FCG/@@download/ENCFF266FCG.bed.gz", "heart_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF086QMV/@@download/ENCFF086QMV.bed.gz", "heart_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF546TJN/@@download/ENCFF546TJN.bed.gz", "heart_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF614LBC/@@download/ENCFF614LBC.bed.gz", "heart_27me3_hg19.bed.gz")
		
	#Intestine (Homo sapiens sigmoid colon tissue male adult (37 years))
	elif args.tissue.lower() == "intestine":
		download_ref("https://www.encodeproject.org/files/ENCFF409HTF/@@download/ENCFF409HTF.bed.gz", "intestine_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF878SLR/@@download/ENCFF878SLR.bed.gz", "intestine_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF314GGO/@@download/ENCFF314GGO.bed.gz", "intestine_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF099PRP/@@download/ENCFF099PRP.bed.gz", "intestine_27me3_hg19.bed.gz")
		
	#iPS cell (Homo sapiens iPS DF 19.11)
	elif args.tissue.lower() == "ips":
		download_ref("https://www.encodeproject.org/files/ENCFF485HSK/@@download/ENCFF485HSK.bed.gz", "ips_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF538JKF/@@download/ENCFF538JKF.bed.gz", "ips_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF631JHU/@@download/ENCFF631JHU.bed.gz", "ips_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF396DQE/@@download/ENCFF396DQE.bed.gz", "ips_27me3_hg19.bed.gz")
		
	#Kidney (Homo sapiens kidney tissue male adult (67 years), Homo sapiens kidney tissue female embryo (120 days) for H3K27me3 - only dataset for kidney)
	elif args.tissue.lower() == "kidney":
		download_ref("https://www.encodeproject.org/files/ENCFF388YQF/@@download/ENCFF388YQF.bed.gz", "kidney_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF913FVX/@@download/ENCFF913FVX.bed.gz", "kidney_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF953SWO/@@download/ENCFF953SWO.bed.gz", "kidney_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF754PGB/@@download/ENCFF754PGB.bed.gz", "kidney_27me3_hg19.bed.gz")
	
	#Liver (Homo sapiens liver tissue male adult (32 years))
	elif args.tissue.lower() == "liver":
		download_ref("https://www.encodeproject.org/files/ENCFF327BBS/@@download/ENCFF327BBS.bed.gz", "liver_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF065NRN/@@download/ENCFF065NRN.bed.gz", "liver_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF752QSK/@@download/ENCFF752QSK.bed.gz", "liver_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF933NEK/@@download/ENCFF933NEK.bed.gz", "liver_27me3_hg19.bed.gz")
		
	#Lung (Homo sapiens upper lobe of left lung tissue female adult (51 years))
	elif args.tissue.lower() == "lung":
		download_ref("https://www.encodeproject.org/files/ENCFF433RBT/@@download/ENCFF433RBT.bed.gz", "lung_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF142ZXD/@@download/ENCFF142ZXD.bed.gz", "lung_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF307QJK/@@download/ENCFF307QJK.bed.gz", "lung_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF957XLF/@@download/ENCFF957XLF.bed.gz", "lung_27me3_hg19.bed.gz")
		
	#Neuron (Homo sapiens SK-N-SH)
	elif args.tissue.lower() == "neuron":
		download_ref("https://www.encodeproject.org/files/ENCFF580GTZ/@@download/ENCFF580GTZ.bed.gz", "neuron_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF363ZFM/@@download/ENCFF363ZFM.bed.gz", "neuron_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF362OBM/@@download/ENCFF362OBM.bed.gz", "neuron_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF102YDR/@@download/ENCFF102YDR.bed.gz", "neuron_27me3_hg19.bed.gz")
		
	#Ovary (Homo sapiens ovary tissue female adult (30 years))
	elif args.tissue.lower() == "ovary":
		download_ref("https://www.encodeproject.org/files/ENCFF917PWI/@@download/ENCFF917PWI.bed.gz", "ovary_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF320JHG/@@download/ENCFF320JHG.bed.gz", "ovary_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF657AUA/@@download/ENCFF657AUA.bed.gz", "ovary_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF712UCB/@@download/ENCFF712UCB.bed.gz", "ovary_27me3_hg19.bed.gz")
		
	#Pancreas (Homo sapiens pancreas tissue male adult (34 years))
	elif args.tissue.lower() == "pancreas":
		download_ref("https://www.encodeproject.org/files/ENCFF668HLF/@@download/ENCFF668HLF.bed.gz", "pancreas_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF340YEE/@@download/ENCFF340YEE.bed.gz", "pancreas_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF583QFI/@@download/ENCFF583QFI.bed.gz", "pancreas_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF973GJC/@@download/ENCFF973GJC.bed.gz", "pancreas_27me3_hg19.bed.gz")
		
	#Prostate (Homo sapiens prostate gland tissue male adult (37 years))
	elif args.tissue.lower() == "prostate":
		download_ref("https://www.encodeproject.org/files/ENCFF942RGZ/@@download/ENCFF942RGZ.bed.gz", "prostate_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF495DBS/@@download/ENCFF495DBS.bed.gz", "prostate_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF002NFG/@@download/ENCFF002NFG.bed.gz", "prostate_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF318QJC/@@download/ENCFF318QJC.bed.gz", "prostate_27me3_hg19.bed.gz")
		
	#Skeletal muscle (Homo sapiens muscle of leg tissue female embryo (110 days))
	elif args.tissue.lower() == "skeletal_muscle":
		download_ref("https://www.encodeproject.org/files/ENCFF746HHI/@@download/ENCFF746HHI.bed.gz", "skeletal_muscle_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF027HBT/@@download/ENCFF027HBT.bed.gz", "skeletal_muscle_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF563PBB/@@download/ENCFF563PBB.bed.gz", "skeletal_muscle_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF191PRP/@@download/ENCFF191PRP.bed.gz", "skeletal_muscle_27me3_hg19.bed.gz")
		
	#Skin (Homo sapiens suprapubic skin tissue male adult (37 years))
	elif args.tissue.lower() == "skin":
		download_ref("https://www.encodeproject.org/files/ENCFF904DTJ/@@download/ENCFF904DTJ.bed.gz", "skin_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF778GBQ/@@download/ENCFF778GBQ.bed.gz", "skin_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF390JYY/@@download/ENCFF390JYY.bed.gz", "skin_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF752AID/@@download/ENCFF752AID.bed.gz", "skin_27me3_hg19.bed.gz")
		
	#Spleen (Homo sapiens spleen tissue male adult (37 years))
	elif args.tissue.lower() == "spleen":
		download_ref("https://www.encodeproject.org/files/ENCFF585JQZ/@@download/ENCFF585JQZ.bed.gz", "spleen_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF736XRS/@@download/ENCFF736XRS.bed.gz", "spleen_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF538SIJ/@@download/ENCFF538SIJ.bed.gz", "spleen_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF065IBC/@@download/ENCFF065IBC.bed.gz", "spleen_27me3_hg19.bed.gz")
		
	#Stomach (Homo sapiens stomach tissue male adult (37 years))
	elif args.tissue.lower() == "stomach":
		download_ref("https://www.encodeproject.org/files/ENCFF324NLX/@@download/ENCFF324NLX.bed.gz", "stomach_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF161ERT/@@download/ENCFF161ERT.bed.gz", "stomach_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF506AKR/@@download/ENCFF506AKR.bed.gz", "stomach_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF435CWM/@@download/ENCFF435CWM.bed.gz", "stomach_27me3_hg19.bed.gz")
		
	#Testis (Homo sapiens testis tissue male adult (37 years))
	elif args.tissue.lower() == "testis":
		download_ref("https://www.encodeproject.org/files/ENCFF620AJW/@@download/ENCFF620AJW.bed.gz", "testis_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF047XWN/@@download/ENCFF047XWN.bed.gz", "testis_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF567XUE/@@download/ENCFF567XUE.bed.gz", "testis_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF563BJS/@@download/ENCFF563BJS.bed.gz", "testis_27me3_hg19.bed.gz")
		
	#Thyroid (Homo sapiens thyroid gland tissue female adult (51 years))
	elif args.tissue.lower() == "thyroid":
		download_ref("https://www.encodeproject.org/files/ENCFF346LTA/@@download/ENCFF346LTA.bed.gz", "thyroid_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF228FJJ/@@download/ENCFF228FJJ.bed.gz", "thyroid_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF086VAY/@@download/ENCFF086VAY.bed.gz", "thyroid_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF957XLF/@@download/ENCFF957XLF.bed.gz", "thyroid_27me3_hg19.bed.gz")
		
	#Uterus (Homo sapiens uterus tissue female adult (51 years))
	elif args.tissue.lower() == "uterus":
		download_ref("https://www.encodeproject.org/files/ENCFF995XPF/@@download/ENCFF995XPF.bed.gz", "uterus_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF155UOS/@@download/ENCFF155UOS.bed.gz", "uterus_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF471JES/@@download/ENCFF471JES.bed.gz", "uterus_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF119JPA/@@download/ENCFF119JPA.bed.gz", "uterus_27me3_hg19.bed.gz")
		
	#Vagina (Homo sapiens vagina tissue female adult (51 years))
	elif args.tissue.lower() == "vagina":
		download_ref("https://www.encodeproject.org/files/ENCFF324OUK/@@download/ENCFF324OUK.bed.gz", "vagina_4me1_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF367PZY/@@download/ENCFF367PZY.bed.gz", "vagina_4me3_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF435LEQ/@@download/ENCFF435LEQ.bed.gz", "vagina_27ac_hg19.bed.gz")
		download_ref("https://www.encodeproject.org/files/ENCFF993WTL/@@download/ENCFF993WTL.bed.gz", "vagina_27me3_hg19.bed.gz")
	
###Import all of the bed files to be worked on (input and histone mark ChIP-seq) using pybedtools
if args.ref_genome.lower() == "hg38" or args.ref_genome.lower() == "grch38":
	bed_ref = "hg38"
	bed_tss_path = "ref_files/refTSS_v3.1_human_coordinate.hg38.bed.gz"
elif args.ref_genome.lower() == "hg19" or args.ref_genome.lower() == "grch37":
	bed_ref = "hg19"
	bed_tss_path = "ref_files/refGene_hg19_TSS.bed"

#Manipulate the TSS bed file so that the coordinates are expanded +-2kb (based on the enrichment we see in Mas et al. 2018 - https://www.nature.com/articles/s41588-018-0218-5?proof=t)
if args.verbose:
	print("Finding overlap of {input_file} coordinates with {tissue_type} histone ChIP-seq peaks based on {ref_genome}...".format(input_file = args.input, tissue_type = args.tissue, ref_genome = args.ref_genome))
	print("\n")

bed_tss_orig = pd.read_csv(bed_tss_path, sep = "\t", header = None)
bed_tss_orig.iloc[:,1] = bed_tss_orig.iloc[:,1] - args.tss_distance_upstream
bed_tss_orig.iloc[:,2] = bed_tss_orig.iloc[:,2] + args.tss_distance_downstream
bed_tss_orig = bed_tss_orig[bed_tss_orig.iloc[:,1] >= 0]
bed_tss = pybedtools.BedTool().from_dataframe(bed_tss_orig)

#Import the input bed, along with the histone ChIP-seq bed files
##Because the raw input file will not necessarily be in UCSC BED format, first import it as a pandas data frame, sort by coordinates, and then input coordinates as a pybedtools object
if args.input_header:
	file_input = pd.read_csv(args.input, delim_whitespace = True)
else:
	file_input = pd.read_csv(args.input, delim_whitespace = True, header = None)
	
cols_array = args.bed_cols.split(",")
chr_col = int(cols_array[0]) - 1
start_col = int(cols_array[1]) - 1
end_col = int(cols_array[2]) - 1

col_names = file_input.columns
file_input.sort_values(by = [col_names[chr_col], col_names[start_col], col_names[end_col]])

##Sort the input file by coordinates
bed_input = file_input[[col_names[chr_col], col_names[start_col], col_names[end_col]]]
if str(bed_input.iloc[:,0]).startswith("chr") == False:
	bed_input.iloc[:,0] = "chr" + bed_input.iloc[:,0].astype(str)
bed_input = pybedtools.BedTool.from_dataframe(bed_input)
bed_4me1 = pybedtools.BedTool("ref_files/" + args.tissue.lower() + "_4me1_" + bed_ref + ".bed.gz")
bed_4me3 = pybedtools.BedTool("ref_files/" + args.tissue.lower() + "_4me3_" + bed_ref + ".bed.gz")
bed_27ac = pybedtools.BedTool("ref_files/" + args.tissue.lower() + "_27ac_" + bed_ref + ".bed.gz")
bed_27me3 = pybedtools.BedTool("ref_files/" + args.tissue.lower() + "_27me3_" + bed_ref + ".bed.gz")

###Compare the input file to the reference bed files, starting with the modified TSS bed
##Do the peaks fall within 2 kb of a TSS?
overlaps_tss = bed_input.intersect(bed_tss)
no_tss = bed_input.intersect(bed_tss, v = True)

#Classifying putative promoters
overlaps_4me3 = overlaps_tss.intersect(bed_4me3)
no_4me3 =overlaps_tss.intersect(bed_4me3, v = True)

active_promoter = overlaps_4me3.intersect(bed_27me3, v = True)
active_promoter = active_promoter.sort()

bivalent_promoter = overlaps_4me3.intersect(bed_27me3)
bivalent_promoter = bivalent_promoter.sort()

silenced_promoter = no_4me3.intersect(bed_27me3)
silenced_promoter = silenced_promoter.sort()

marked_promoters = active_promoter + bivalent_promoter + silenced_promoter

unclassified_overlaps_tss = overlaps_tss.intersect(marked_promoters, v = True)
unclassified_overlaps_tss = unclassified_overlaps_tss.sort()

#Classifying putative enhancers
overlaps_4me1 = no_tss.intersect(bed_4me1)

active_enhancer = overlaps_4me1.intersect(bed_27ac)
active_enhancer = active_enhancer.sort()

poised_enhancer = overlaps_4me1.intersect(bed_27me3)
poised_enhancer = poised_enhancer.sort()

active_and_poised = active_enhancer + poised_enhancer
primed_enhancer = overlaps_4me1.intersect(active_and_poised, v = True)
primed_enhancer = primed_enhancer.sort()

active_poised_primed = active_and_poised + primed_enhancer
unclassified_no_tss = no_tss.intersect(active_poised_primed, v = True)
unclassified_no_tss = unclassified_no_tss.sort()

#If they have one or more coordinate entries, convert each of the classifier variables to a pandas DataFrame object so that we can append the classifications as an additional column
if args.verbose:
	print("Preparing output file ({output_file})...".format(output_file = args.output))
	print("\n")

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

if unclassified_overlaps_tss.count() > 0:
	unclassified_overlaps_tss_df = pybedtools.BedTool.to_dataframe(unclassified_overlaps_tss)
	unclassified_overlaps_tss_df["region_classification"] = "unclassified_within_2kb_of_tss"
else:
	unclassified_overlaps_tss_df = pd.DataFrame()

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

if unclassified_no_tss.count() > 0:
	unclassified_no_tss_df = pybedtools.BedTool.to_dataframe(unclassified_no_tss)
	unclassified_no_tss_df["region_classification"] = "unclassified_beyond_2kb_of_tss"
else:
	unclassified_no_tss_df = pd.DataFrame()
	
#Concatenate all of the pandas DataFrame objects into one, remove all duplicate lines, and sort by coordinate
all_data_frames = [active_promoter_df, bivalent_promoter_df, silenced_promoter_df, unclassified_overlaps_tss_df, active_enhancer_df, poised_enhancer_df, primed_enhancer_df, unclassified_no_tss_df]

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
	
#Print out the number of peaks/coordinates in each category
if args.verbose:
	print("Identified regions:\nPutative Promoters\n{active_promoter_count} regions in an active promoter\n{bivalent_promoter_count} regions in a bivalent promoter\n{silenced_promoter_count} regions in a silenced promoter\n{unclassified_overlaps_tss_count} regions in an unclassified region within 2 kb of a TSS\n\nPutative Enhancers\n{active_enhancer_count} regions in an active enhancer\n{poised_enhancer_count} regions in a poised enhancer\n{primed_enhancer_count} in a primed enhancer\n{unclassified_no_tss_count} regions in an unclassified region greater than 2 kb from a TSS\n".format(active_promoter_count = list(map(lambda x: x.startswith("active_promoter"), concatenated_df[last_col])).count(True), bivalent_promoter_count = list(map(lambda x: x.startswith("bivalent_promoter"), concatenated_df[last_col])).count(True), silenced_promoter_count = list(map(lambda x: x.startswith("silenced_promoter"), concatenated_df[last_col])).count(True), unclassified_overlaps_tss_count = list(map(lambda x: x.startswith("unclassified_within_2kb_of_tss"), concatenated_df[last_col])).count(True), active_enhancer_count = list(map(lambda x: x.startswith("active_enhancer"), concatenated_df[last_col])).count(True), poised_enhancer_count = list(map(lambda x: x.startswith("poised_enhancer"), concatenated_df[last_col])).count(True), primed_enhancer_count = list(map(lambda x: x.startswith("primed_enhancer"), concatenated_df[last_col])).count(True), unclassified_no_tss_count = list(map(lambda x: x.startswith("unclassified_beyond_2kb_of_tss"), concatenated_df[last_col])).count(True)))

#Append the functional annotations onto the original input
annotations = concatenated_df.iloc[:,-1]
annotations = annotations.reset_index(drop = True)
file_input = file_input.reset_index(drop = True)
file_input["func_anno"] = annotations

if args.input_header:
	file_input.to_csv(args.output, sep = "\t", index = False)
else:
	file_input.to_csv(args.output, sep = "\t", header = False, index = False)