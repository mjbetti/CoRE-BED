###Developed by Michael Betti, April 2021

import os, sys, requests, argparse, logging, pybedtools, pandas as pd

parser = argparse.ArgumentParser(add_help = True)
parser.add_argument("-i", "--input", type = str, help = "the input bed file")
parser.add_argument("-g", "--ref_genome", type = str, help = "the human reference genome build on which the input coordinates are based (valid options: GRCh38/hg38 and GRCh37/hg19)")
parser.add_argument("-t", "--tissue", type = str, help = "the tissue of interest (valid options: Blood, Brain, ES, Heart, Intestine, Liver, Lung, Muscle, Skin, Thyroid)")
parser.add_argument("-o", "--output", type = str, help = "the name of the output file", default = "out.bed")
parser.add_argument("-v", "--verbose", help = "return script progress as terminal output", action = "store_true")
args = parser.parse_args()

#Check that required arguments are specified
assert args.input, "Must specify input file (-i, --input)"
assert args.ref_genome, "Must specify reference genome build (-g, --ref_genome)"
assert args.tissue, "Must specify tissue type (-t, --tissue)"

#Check that specified tissue type is one of the 10 valid options
assert args.tissue.lower() == "blood" or args.tissue.lower() == "brain" or args.tissue.lower() == "es" or args.tissue.lower() == "heart" or args.tissue.lower() == "intestine" or args.tissue.lower() == "liver" or args.tissue.lower() == "lung" or args.tissue.lower() == "muscle" or args.tissue.lower() == "skin" or args.tissue.lower() == "thyroid", "Tissue type must be one of the 10 valid options (Blood, Brain, ES, Heart, Intestine, Liver, Lung, Muscle, Skin, Thyroid)"

#Check that not too many arguments are specified
#assert len(args) < 6, "Too many arguments specified"

#Download the appropriate reference files based on the specified genome build and tissue arguments
if args.verbose:
	print("Classifying {tissue_type} cis-regulatory regions in {input_file}...".format(input_file = args.input, tissue_type = args.tissue))
	print("\n")
if not os.path.exists('ref_files'):
	os.mkdir('ref_files')

if args.verbose:
	print("Downloading {ref_genome} TSS coordinates and {tissue_type} histone ChIP-seq bed files...".format(ref_genome = args.ref_genome, tissue_type = args.tissue))
	print("\n")

#GRCh38/hg38
if args.ref_genome.lower() == "hg38" or args.ref_genome.lower() == "grch38":
	url_tss = "http://reftss.clst.riken.jp/datafiles/3.1/human/refTSS_v3.1_human_coordinate.hg38.bed.gz"
	out_path_tss = os.path.join('ref_files', "refTSS_v3.1_human_coordinate.hg38.bed.gz")
	r = requests.get(url_tss, allow_redirects=True)
	open(out_path_tss, 'wb').write(r.content)
	
	#Blood (Homo sapiens peripheral blood mononuclear cell male adult (39 years))	
	if args.tissue.lower() == "blood":
		url_4me1 = "https://www.encodeproject.org/files/ENCFF734KLT/@@download/ENCFF734KLT.bed.gz"
		out_path_4me1 = os.path.join('ref_files', "blood_4me1_hg38.bed.gz")
		r = requests.get(url_4me1, allow_redirects=True)
		open(out_path_4me1, 'wb').write(r.content)
		url_4me3 = "https://www.encodeproject.org/files/ENCFF165VDC/@@download/ENCFF165VDC.bed.gz"
		out_path_4me3 = os.path.join('ref_files', "blood_4me3_hg38.bed.gz")
		r = requests.get(url_4me3, allow_redirects=True)
		open(out_path_4me3, 'wb').write(r.content)
		url_27ac = "https://www.encodeproject.org/files/ENCFF844OYX/@@download/ENCFF844OYX.bed.gz"
		out_path_27ac = os.path.join('ref_files', "blood_27ac_hg38.bed.gz")
		r = requests.get(url_27ac, allow_redirects=True)
		open(out_path_27ac, 'wb').write(r.content)
		url_27me3 = "https://www.encodeproject.org/files/ENCFF636JEP/@@download/ENCFF636JEP.bed.gz"
		out_path_27me3 = os.path.join('ref_files', "blood_27me3_hg38.bed.gz")
		r = requests.get(url_27me3, allow_redirects=True)
		open(out_path_27me3, 'wb').write(r.content)
		
	#Brain (Homo sapiens SK-N-SH)
	elif args.tissue.lower() == "brain":
		url_4me1 = "https://www.encodeproject.org/files/ENCFF632FAM/@@download/ENCFF632FAM.bed.gz"
		out_path_4me1 = os.path.join('ref_files', "brain_4me1_hg38.bed.gz")
		r = requests.get(url_4me1, allow_redirects=True)
		open(out_path_4me1, 'wb').write(r.content)
		url_4me3 = "https://www.encodeproject.org/files/ENCFF682JYE/@@download/ENCFF682JYE.bed.gz"
		out_path_4me3 = os.path.join('ref_files', "brain_4me3_hg38.bed.gz")
		r = requests.get(url_4me3, allow_redirects=True)
		open(out_path_4me3, 'wb').write(r.content)
		url_27ac = "https://www.encodeproject.org/files/ENCFF138VUT/@@download/ENCFF138VUT.bed.gz"
		out_path_27ac = os.path.join('ref_files', "brain_27ac_hg38.bed.gz")
		r = requests.get(url_27ac, allow_redirects=True)
		open(out_path_27ac, 'wb').write(r.content)
		url_27me3 = "https://www.encodeproject.org/files/ENCFF277NRX/@@download/ENCFF277NRX.bed.gz"
		out_path_27me3 = os.path.join('ref_files', "brain_27me3_hg38.bed.gz")
		r = requests.get(url_27me3, allow_redirects=True)
		open(out_path_27me3, 'wb').write(r.content)
		
	#Embryonic stem cell (Homo sapiens H1)
	elif args.tissue.lower() == "es":
		url_4me1 = "https://www.encodeproject.org/files/ENCFF613QAB/@@download/ENCFF613QAB.bed.gz"
		out_path_4me1 = os.path.join('ref_files', "es_4me1_hg38.bed.gz")
		r = requests.get(url_4me1, allow_redirects=True)
		open(out_path_4me1, 'wb').write(r.content)
		url_4me3 = "https://www.encodeproject.org/files/ENCFF356MCC/@@download/ENCFF356MCC.bed.gz"
		out_path_4me3 = os.path.join('ref_files', "es_4me3_hg38.bed.gz")
		r = requests.get(url_4me3, allow_redirects=True)
		open(out_path_4me3, 'wb').write(r.content)
		url_27ac = "https://www.encodeproject.org/files/ENCFF689CJG/@@download/ENCFF689CJG.bed.gz"
		out_path_27ac = os.path.join('ref_files', "es_27ac_hg38.bed.gz")
		r = requests.get(url_27ac, allow_redirects=True)
		open(out_path_27ac, 'wb').write(r.content)
		url_27me3 = "https://www.encodeproject.org/files/ENCFF050CUG/@@download/ENCFF050CUG.bed.gz"
		out_path_27me3 = os.path.join('ref_files', "es_27me3_hg38.bed.gz")
		r = requests.get(url_27me3, allow_redirects=True)
		open(out_path_27me3, 'wb').write(r.content)
		
	#Heart (Homo sapiens heart left ventricle tissue female adult (53 years))
	elif args.tissue.lower() == "heart":
		url_4me1 = "https://www.encodeproject.org/files/ENCFF194KSV/@@download/ENCFF194KSV.bed.gz"
		out_path_4me1 = os.path.join('ref_files', "heart_4me1_hg38.bed.gz")
		r = requests.get(url_4me1, allow_redirects=True)
		open(out_path_4me1, 'wb').write(r.content)
		url_4me3 = "https://www.encodeproject.org/files/ENCFF008ODN/@@download/ENCFF008ODN.bed.gz"
		out_path_4me3 = os.path.join('ref_files', "heart_4me3_hg38.bed.gz")
		r = requests.get(url_4me3, allow_redirects=True)
		open(out_path_4me3, 'wb').write(r.content)
		url_27ac = "https://www.encodeproject.org/files/ENCFF635VBV/@@download/ENCFF635VBV.bed.gz"
		out_path_27ac = os.path.join('ref_files', "heart_27ac_hg38.bed.gz")
		r = requests.get(url_27ac, allow_redirects=True)
		open(out_path_27ac, 'wb').write(r.content)
		url_27me3 = "https://www.encodeproject.org/files/ENCFF242JDO/@@download/ENCFF242JDO.bed.gz"
		out_path_27me3 = os.path.join('ref_files', "heart_27me3_hg38.bed.gz")
		r = requests.get(url_27me3, allow_redirects=True)
		open(out_path_27me3, 'wb').write(r.content)
		
	#Intestine (Homo sapiens sigmoid colon tissue male adult (37 years))
	elif args.tissue.lower() == "intestine":
		url_4me1 = "https://www.encodeproject.org/files/ENCFF052XGA/@@download/ENCFF052XGA.bed.gz"
		out_path_4me1 = os.path.join('ref_files', "intestine_4me1_hg38.bed.gz")
		r = requests.get(url_4me1, allow_redirects=True)
		open(out_path_4me1, 'wb').write(r.content)
		url_4me3 = "https://www.encodeproject.org/files/ENCFF250ZHO/@@download/ENCFF250ZHO.bed.gz"
		out_path_4me3 = os.path.join('ref_files', "intestine_4me3_hg38.bed.gz")
		r = requests.get(url_4me3, allow_redirects=True)
		open(out_path_4me3, 'wb').write(r.content)
		url_27ac = "https://www.encodeproject.org/files/ENCFF111SYA/@@download/ENCFF111SYA.bed.gz"
		out_path_27ac = os.path.join('ref_files', "intestine_27ac_hg38.bed.gz")
		r = requests.get(url_27ac, allow_redirects=True)
		open(out_path_27ac, 'wb').write(r.content)
		url_27me3 = "https://www.encodeproject.org/files/ENCFF310IST/@@download/ENCFF310IST.bed.gz"
		out_path_27me3 = os.path.join('ref_files', "intestine_27me3_hg38.bed.gz")
		r = requests.get(url_27me3, allow_redirects=True)
		open(out_path_27me3, 'wb').write(r.content)

	#Liver (Homo sapiens liver tissue male adult (32 years))
	elif args.tissue.lower() == "liver":
		url_4me1 = "https://www.encodeproject.org/files/ENCFF625TWQ/@@download/ENCFF625TWQ.bed.gz"
		out_path_4me1 = os.path.join('ref_files', "liver_4me1_hg38.bed.gz")
		r = requests.get(url_4me1, allow_redirects=True)
		open(out_path_4me1, 'wb').write(r.content)
		url_4me3 = "https://www.encodeproject.org/files/ENCFF830ODB/@@download/ENCFF830ODB.bed.gz"
		out_path_4me3 = os.path.join('ref_files', "liver_4me3_hg38.bed.gz")
		r = requests.get(url_4me3, allow_redirects=True)
		open(out_path_4me3, 'wb').write(r.content)
		url_27ac = "https://www.encodeproject.org/files/ENCFF287VIA/@@download/ENCFF287VIA.bed.gz"
		out_path_27ac = os.path.join('ref_files', "liver_27ac_hg38.bed.gz")
		r = requests.get(url_27ac, allow_redirects=True)
		open(out_path_27ac, 'wb').write(r.content)
		url_27me3 = "https://www.encodeproject.org/files/ENCFF287VIA/@@download/ENCFF287VIA.bed.gz"
		out_path_27me3 = os.path.join('ref_files', "liver_27me3_hg38.bed.gz")
		r = requests.get(url_27me3, allow_redirects=True)
		open(out_path_27me3, 'wb').write(r.content)
		
	#Lung (Homo sapiens upper lobe of left lung tissue female adult (51 years))
	elif args.tissue.lower() == "lung":
		url_4me1 = "https://www.encodeproject.org/files/ENCFF926YKP/@@download/ENCFF926YKP.bed.gz"
		out_path_4me1 = os.path.join('ref_files', "lung_4me1_hg38.bed.gz")
		r = requests.get(url_4me1, allow_redirects=True)
		open(out_path_4me1, 'wb').write(r.content)
		url_4me3 = "https://www.encodeproject.org/files/ENCFF843JPR/@@download/ENCFF843JPR.bed.gz"
		out_path_4me3 = os.path.join('ref_files', "lung_4me3_hg38.bed.gz")
		r = requests.get(url_4me3, allow_redirects=True)
		open(out_path_4me3, 'wb').write(r.content)
		url_27ac = "https://www.encodeproject.org/files/ENCFF496RDA/@@download/ENCFF496RDA.bed.gz"
		out_path_27ac = os.path.join('ref_files', "lung_27ac_hg38.bed.gz")
		r = requests.get(url_27ac, allow_redirects=True)
		open(out_path_27ac, 'wb').write(r.content)
		url_27me3 = "https://www.encodeproject.org/files/ENCFF248IPH/@@download/ENCFF248IPH.bed.gz"
		out_path_27me3 = os.path.join('ref_files', "lung_27me3_hg38.bed.gz")
		r = requests.get(url_27me3, allow_redirects=True)
		open(out_path_27me3, 'wb').write(r.content)
		
	#Muscle (Homo sapiens muscle of leg tissue female embryo (110 days))
	elif args.tissue.lower() == "muscle":
		url_4me1 = "https://www.encodeproject.org/files/ENCFF170ZFK/@@download/ENCFF170ZFK.bed.gz"
		out_path_4me1 = os.path.join('ref_files', "muscle_4me1_hg38.bed.gz")
		r = requests.get(url_4me1, allow_redirects=True)
		open(out_path_4me1, 'wb').write(r.content)
		url_4me3 = "https://www.encodeproject.org/files/ENCFF942MOM/@@download/ENCFF942MOM.bed.gz"
		out_path_4me3 = os.path.join('ref_files', "muscle_4me3_hg38.bed.gz")
		r = requests.get(url_4me3, allow_redirects=True)
		open(out_path_4me3, 'wb').write(r.content)
		url_27ac = "https://www.encodeproject.org/files/ENCFF154ZCF/@@download/ENCFF154ZCF.bed.gz"
		out_path_27ac = os.path.join('ref_files', "muscle_27ac_hg38.bed.gz")
		r = requests.get(url_27ac, allow_redirects=True)
		open(out_path_27ac, 'wb').write(r.content)
		url_27me3 = "https://www.encodeproject.org/files/ENCFF672MZK/@@download/ENCFF672MZK.bed.gz"
		out_path_27me3 = os.path.join('ref_files', "muscle_27me3_hg38.bed.gz")
		r = requests.get(url_27me3, allow_redirects=True)
		open(out_path_27me3, 'wb').write(r.content)
		
	#Skin (Homo sapiens suprapubic skin tissue male adult (37 years))
	elif args.tissue.lower() == "skin":
		url_4me1 = "https://www.encodeproject.org/files/ENCFF800MDH/@@download/ENCFF800MDH.bed.gz"
		out_path_4me1 = os.path.join('ref_files', "skin_4me1_hg38.bed.gz")
		r = requests.get(url_4me1, allow_redirects=True)
		open(out_path_4me1, 'wb').write(r.content)
		url_4me3 = "https://www.encodeproject.org/files/ENCFF240CFZ/@@download/ENCFF240CFZ.bed.gz"
		out_path_4me3 = os.path.join('ref_files', "skin_4me3_hg38.bed.gz")
		r = requests.get(url_4me3, allow_redirects=True)
		open(out_path_4me3, 'wb').write(r.content)
		url_27ac = "https://www.encodeproject.org/files/ENCFF398EEO/@@download/ENCFF398EEO.bed.gz"
		out_path_27ac = os.path.join('ref_files', "skin_27ac_hg38.bed.gz")
		r = requests.get(url_27ac, allow_redirects=True)
		open(out_path_27ac, 'wb').write(r.content)
		url_27me3 = "https://www.encodeproject.org/files/ENCFF686GBZ/@@download/ENCFF686GBZ.bed.gz"
		out_path_27me3 = os.path.join('ref_files', "skin_27me3_hg38.bed.gz")
		r = requests.get(url_27me3, allow_redirects=True)
		open(out_path_27me3, 'wb').write(r.content)
		
	#Thyroid (Homo sapiens thyroid gland tissue female adult (51 years))
	elif args.tissue.lower() == "thyroid":
		url_4me1 = "https://www.encodeproject.org/files/ENCFF346IPY/@@download/ENCFF346IPY.bed.gz"
		out_path_4me1 = os.path.join('ref_files', "thyroid_4me1_hg38.bed.gz")
		r = requests.get(url_4me1, allow_redirects=True)
		open(out_path_4me1, 'wb').write(r.content)
		url_4me3 = "https://www.encodeproject.org/files/ENCFF925FSR/@@download/ENCFF925FSR.bed.gz"
		out_path_4me3 = os.path.join('ref_files', "thyroid_4me3_hg38.bed.gz")
		r = requests.get(url_4me3, allow_redirects=True)
		open(out_path_4me3, 'wb').write(r.content)
		url_27ac = "https://www.encodeproject.org/files/ENCFF138WNV/@@download/ENCFF138WNV.bed.gz"
		out_path_27ac = os.path.join('ref_files', "thyroid_27ac_hg38.bed.gz")
		r = requests.get(url_27ac, allow_redirects=True)
		open(out_path_27ac, 'wb').write(r.content)
		url_27me3 = "https://www.encodeproject.org/files/ENCFF607NRZ/@@download/ENCFF607NRZ.bed.gz"
		out_path_27me3 = os.path.join('ref_files', "thyroid_27me3_hg38.bed.gz")
		r = requests.get(url_27me3, allow_redirects=True)
		open(out_path_27me3, 'wb').write(r.content)

###################################################################################
#GRCh37/hg19
elif args.ref_genome.lower() == "hg19" or args.ref_genome.lower() == "grch37":
	url_tss = "https://raw.githubusercontent.com/Shicheng-Guo/AnnotationDatabase/master/hg19/refGene_hg19_TSS.bed"
	out_path_tss = os.path.join('ref_files', "refGene_hg19_TSS.bed")
	r = requests.get(url_tss, allow_redirects=True)
	open(out_path_tss, 'wb').write(r.content)
	
	#Blood (Homo sapiens peripheral blood mononuclear cell male adult (39 years))	
	if args.tissue.lower() == "blood":
		url_4me1 = "https://www.encodeproject.org/files/ENCFF553SEK/@@download/ENCFF553SEK.bed.gz"
		out_path_4me1 = os.path.join('ref_files', "blood_4me1_hg19.bed.gz")
		r = requests.get(url_4me1, allow_redirects=True)
		open(out_path_4me1, 'wb').write(r.content)
		url_4me3 = "https://www.encodeproject.org/files/ENCFF477HMF/@@download/ENCFF477HMF.bed.gz"
		out_path_4me3 = os.path.join('ref_files', "blood_4me3_hg19.bed.gz")
		r = requests.get(url_4me3, allow_redirects=True)
		open(out_path_4me3, 'wb').write(r.content)
		url_27ac = "https://www.encodeproject.org/files/ENCFF442KRP/@@download/ENCFF442KRP.bed.gz"
		out_path_27ac = os.path.join('ref_files', "blood_27ac_hg19.bed.gz")
		r = requests.get(url_27ac, allow_redirects=True)
		open(out_path_27ac, 'wb').write(r.content)
		url_27me3 = "https://www.encodeproject.org/files/ENCFF987VSV/@@download/ENCFF987VSV.bed.gz"
		out_path_27me3 = os.path.join('ref_files', "blood_27me3_hg19.bed.gz")
		r = requests.get(url_27me3, allow_redirects=True)
		open(out_path_27me3, 'wb').write(r.content)
		
	#Brain (Homo sapiens SK-N-SH)
	elif args.tissue.lower() == "brain":
		url_4me1 = "https://www.encodeproject.org/files/ENCFF580GTZ/@@download/ENCFF580GTZ.bed.gz"
		out_path_4me1 = os.path.join('ref_files', "brain_4me1_hg19.bed.gz")
		r = requests.get(url_4me1, allow_redirects=True)
		open(out_path_4me1, 'wb').write(r.content)
		url_4me3 = "https://www.encodeproject.org/files/ENCFF363ZFM/@@download/ENCFF363ZFM.bed.gz"
		out_path_4me3 = os.path.join('ref_files', "brain_4me3_hg19.bed.gz")
		r = requests.get(url_4me3, allow_redirects=True)
		open(out_path_4me3, 'wb').write(r.content)
		url_27ac = "https://www.encodeproject.org/files/ENCFF362OBM/@@download/ENCFF362OBM.bed.gz"
		out_path_27ac = os.path.join('ref_files', "brain_27ac_hg19.bed.gz")
		r = requests.get(url_27ac, allow_redirects=True)
		open(out_path_27ac, 'wb').write(r.content)
		url_27me3 = "https://www.encodeproject.org/files/ENCFF102YDR/@@download/ENCFF102YDR.bed.gz"
		out_path_27me3 = os.path.join('ref_files', "brain_27me3_hg19.bed.gz")
		r = requests.get(url_27me3, allow_redirects=True)
		open(out_path_27me3, 'wb').write(r.content)
		
	#Embryonic stem cell (Homo sapiens H1)
	elif args.tissue.lower() == "es":
		url_4me1 = "https://www.encodeproject.org/files/ENCFF834NWA/@@download/ENCFF834NWA.bed.gz"
		out_path_4me1 = os.path.join('ref_files', "es_4me1_hg19.bed.gz")
		r = requests.get(url_4me1, allow_redirects=True)
		open(out_path_4me1, 'wb').write(r.content)
		url_4me3 = "https://www.encodeproject.org/files/ENCFF865TYO/@@download/ENCFF865TYO.bed.gz"
		out_path_4me3 = os.path.join('ref_files', "es_4me3_hg19.bed.gz")
		r = requests.get(url_4me3, allow_redirects=True)
		open(out_path_4me3, 'wb').write(r.content)
		url_27ac = "https://www.encodeproject.org/files/ENCFF335JSA/@@download/ENCFF335JSA.bed.gz"
		out_path_27ac = os.path.join('ref_files', "es_27ac_hg19.bed.gz")
		r = requests.get(url_27ac, allow_redirects=True)
		open(out_path_27ac, 'wb').write(r.content)
		url_27me3 = "https://www.encodeproject.org/files/ENCFF900RIU/@@download/ENCFF900RIU.bed.gz"
		out_path_27me3 = os.path.join('ref_files', "es_27me3_hg19.bed.gz")
		r = requests.get(url_27me3, allow_redirects=True)
		open(out_path_27me3, 'wb').write(r.content)
		
	#Heart (Homo sapiens heart left ventricle tissue female adult (53 years))
	elif args.tissue.lower() == "heart":
		url_4me1 = "https://www.encodeproject.org/files/ENCFF266FCG/@@download/ENCFF266FCG.bed.gz"
		out_path_4me1 = os.path.join('ref_files', "heart_4me1_hg19.bed.gz")
		r = requests.get(url_4me1, allow_redirects=True)
		open(out_path_4me1, 'wb').write(r.content)
		url_4me3 = "https://www.encodeproject.org/files/ENCFF086QMV/@@download/ENCFF086QMV.bed.gz"
		out_path_4me3 = os.path.join('ref_files', "heart_4me3_hg19.bed.gz")
		r = requests.get(url_4me3, allow_redirects=True)
		open(out_path_4me3, 'wb').write(r.content)
		url_27ac = "https://www.encodeproject.org/files/ENCFF546TJN/@@download/ENCFF546TJN.bed.gz"
		out_path_27ac = os.path.join('ref_files', "heart_27ac_hg19.bed.gz")
		r = requests.get(url_27ac, allow_redirects=True)
		open(out_path_27ac, 'wb').write(r.content)
		url_27me3 = "https://www.encodeproject.org/files/ENCFF614LBC/@@download/ENCFF614LBC.bed.gz"
		out_path_27me3 = os.path.join('ref_files', "heart_27me3_hg19.bed.gz")
		r = requests.get(url_27me3, allow_redirects=True)
		open(out_path_27me3, 'wb').write(r.content)
		
	#Intestine (Homo sapiens sigmoid colon tissue male adult (37 years))
	elif args.tissue.lower() == "intestine":
		url_4me1 = "https://www.encodeproject.org/files/ENCFF409HTF/@@download/ENCFF409HTF.bed.gz"
		out_path_4me1 = os.path.join('ref_files', "intestine_4me1_hg19.bed.gz")
		r = requests.get(url_4me1, allow_redirects=True)
		open(out_path_4me1, 'wb').write(r.content)
		url_4me3 = "https://www.encodeproject.org/files/ENCFF878SLR/@@download/ENCFF878SLR.bed.gz"
		out_path_4me3 = os.path.join('ref_files', "intestine_4me3_hg19.bed.gz")
		r = requests.get(url_4me3, allow_redirects=True)
		open(out_path_4me3, 'wb').write(r.content)
		url_27ac = "https://www.encodeproject.org/files/ENCFF314GGO/@@download/ENCFF314GGO.bed.gz"
		out_path_27ac = os.path.join('ref_files', "intestine_27ac_hg19.bed.gz")
		r = requests.get(url_27ac, allow_redirects=True)
		open(out_path_27ac, 'wb').write(r.content)
		url_27me3 = "https://www.encodeproject.org/files/ENCFF099PRP/@@download/ENCFF099PRP.bed.gz"
		out_path_27me3 = os.path.join('ref_files', "intestine_27me3_hg19.bed.gz")
		r = requests.get(url_27me3, allow_redirects=True)
		open(out_path_27me3, 'wb').write(r.content)
	
	#Liver (Homo sapiens liver tissue male adult (32 years))
	elif args.tissue.lower() == "liver":
		url_4me1 = "https://www.encodeproject.org/files/ENCFF327BBS/@@download/ENCFF327BBS.bed.gz"
		out_path_4me1 = os.path.join('ref_files', "liver_4me1_hg19.bed.gz")
		r = requests.get(url_4me1, allow_redirects=True)
		open(out_path_4me1, 'wb').write(r.content)
		url_4me3 = "https://www.encodeproject.org/files/ENCFF065NRN/@@download/ENCFF065NRN.bed.gz"
		out_path_4me3 = os.path.join('ref_files', "liver_4me3_hg19.bed.gz")
		r = requests.get(url_4me3, allow_redirects=True)
		open(out_path_4me3, 'wb').write(r.content)
		url_27ac = "https://www.encodeproject.org/files/ENCFF752QSK/@@download/ENCFF752QSK.bed.gz"
		out_path_27ac = os.path.join('ref_files', "liver_27ac_hg19.bed.gz")
		r = requests.get(url_27ac, allow_redirects=True)
		open(out_path_27ac, 'wb').write(r.content)
		url_27me3 = "https://www.encodeproject.org/files/ENCFF933NEK/@@download/ENCFF933NEK.bed.gz"
		out_path_27me3 = os.path.join('ref_files', "liver_27me3_hg19.bed.gz")
		r = requests.get(url_27me3, allow_redirects=True)
		open(out_path_27me3, 'wb').write(r.content)
		
	#Lung (Homo sapiens upper lobe of left lung tissue female adult (51 years))
	elif args.tissue.lower() == "lung":
		url_4me1 = "https://www.encodeproject.org/files/ENCFF433RBT/@@download/ENCFF433RBT.bed.gz"
		out_path_4me1 = os.path.join('ref_files', "lung_4me1_hg19.bed.gz")
		r = requests.get(url_4me1, allow_redirects=True)
		open(out_path_4me1, 'wb').write(r.content)
		url_4me3 = "https://www.encodeproject.org/files/ENCFF142ZXD/@@download/ENCFF142ZXD.bed.gz"
		out_path_4me3 = os.path.join('ref_files', "lung_4me3_hg19.bed.gz")
		r = requests.get(url_4me3, allow_redirects=True)
		open(out_path_4me3, 'wb').write(r.content)
		url_27ac = "https://www.encodeproject.org/files/ENCFF307QJK/@@download/ENCFF307QJK.bed.gz"
		out_path_27ac = os.path.join('ref_files', "lung_27ac_hg19.bed.gz")
		r = requests.get(url_27ac, allow_redirects=True)
		open(out_path_27ac, 'wb').write(r.content)
		url_27me3 = "https://www.encodeproject.org/files/ENCFF957XLF/@@download/ENCFF957XLF.bed.gz"
		out_path_27me3 = os.path.join('ref_files', "lung_27me3_hg19.bed.gz")
		r = requests.get(url_27me3, allow_redirects=True)
		open(out_path_27me3, 'wb').write(r.content)
		
	#Muscle (Homo sapiens muscle of leg tissue female embryo (110 days))
	elif args.tissue.lower() == "muscle":
		url_4me1 = "https://www.encodeproject.org/files/ENCFF746HHI/@@download/ENCFF746HHI.bed.gz"
		out_path_4me1 = os.path.join('ref_files', "muscle_4me1_hg19.bed.gz")
		r = requests.get(url_4me1, allow_redirects=True)
		open(out_path_4me1, 'wb').write(r.content)
		url_4me3 = "https://www.encodeproject.org/files/ENCFF027HBT/@@download/ENCFF027HBT.bed.gz"
		out_path_4me3 = os.path.join('ref_files', "muscle_4me3_hg19.bed.gz")
		r = requests.get(url_4me3, allow_redirects=True)
		open(out_path_4me3, 'wb').write(r.content)
		url_27ac = "https://www.encodeproject.org/files/ENCFF563PBB/@@download/ENCFF563PBB.bed.gz"
		out_path_27ac = os.path.join('ref_files', "muscle_27ac_hg19.bed.gz")
		r = requests.get(url_27ac, allow_redirects=True)
		open(out_path_27ac, 'wb').write(r.content)
		url_27me3 = "https://www.encodeproject.org/files/ENCFF191PRP/@@download/ENCFF191PRP.bed.gz"
		out_path_27me3 = os.path.join('ref_files', "muscle_27me3_hg19.bed.gz")
		r = requests.get(url_27me3, allow_redirects=True)
		open(out_path_27me3, 'wb').write(r.content)
		
	#Skin (Homo sapiens suprapubic skin tissue male adult (37 years))
	elif args.tissue.lower() == "skin":
		url_4me1 = "https://www.encodeproject.org/files/ENCFF904DTJ/@@download/ENCFF904DTJ.bed.gz"
		out_path_4me1 = os.path.join('ref_files', "skin_4me1_hg19.bed.gz")
		r = requests.get(url_4me1, allow_redirects=True)
		open(out_path_4me1, 'wb').write(r.content)
		url_4me3 = "https://www.encodeproject.org/files/ENCFF778GBQ/@@download/ENCFF778GBQ.bed.gz"
		out_path_4me3 = os.path.join('ref_files', "skin_4me3_hg19.bed.gz")
		r = requests.get(url_4me3, allow_redirects=True)
		open(out_path_4me3, 'wb').write(r.content)
		url_27ac = "https://www.encodeproject.org/files/ENCFF390JYY/@@download/ENCFF390JYY.bed.gz"
		out_path_27ac = os.path.join('ref_files', "skin_27ac_hg19.bed.gz")
		r = requests.get(url_27ac, allow_redirects=True)
		open(out_path_27ac, 'wb').write(r.content)
		url_27me3 = "https://www.encodeproject.org/files/ENCFF752AID/@@download/ENCFF752AID.bed.gz"
		out_path_27me3 = os.path.join('ref_files', "skin_27me3_hg19.bed.gz")
		r = requests.get(url_27me3, allow_redirects=True)
		open(out_path_27me3, 'wb').write(r.content)
		
	#Thyroid (Homo sapiens thyroid gland tissue female adult (51 years))
	elif args.tissue.lower() == "thyroid":
		url_4me1 = "https://www.encodeproject.org/files/ENCFF346LTA/@@download/ENCFF346LTA.bed.gz"
		out_path_4me1 = os.path.join('ref_files', "thyroid_4me1_hg19.bed.gz")
		r = requests.get(url_4me1, allow_redirects=True)
		open(out_path_4me1, 'wb').write(r.content)
		url_4me3 = "https://www.encodeproject.org/files/ENCFF228FJJ/@@download/ENCFF228FJJ.bed.gz"
		out_path_4me3 = os.path.join('ref_files', "thyroid_4me3_hg19.bed.gz")
		r = requests.get(url_4me3, allow_redirects=True)
		open(out_path_4me3, 'wb').write(r.content)
		url_27ac = "https://www.encodeproject.org/files/ENCFF086VAY/@@download/ENCFF086VAY.bed.gz"
		out_path_27ac = os.path.join('ref_files', "thyroid_27ac_hg19.bed.gz")
		r = requests.get(url_27ac, allow_redirects=True)
		open(out_path_27ac, 'wb').write(r.content)
		url_27me3 = "https://www.encodeproject.org/files/ENCFF957XLF/@@download/ENCFF957XLF.bed.gz"
		out_path_27me3 = os.path.join('ref_files', "thyroid_27me3_hg19.bed.gz")
		r = requests.get(url_27me3, allow_redirects=True)
		open(out_path_27me3, 'wb').write(r.content)
	
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

bed_tss_orig = pd.read_csv(bed_tss_path, sep = "\t", header = None, usecols = [0, 1, 2], names = ["chr", "start", "end"])
bed_tss_orig["start"] = bed_tss_orig["start"] - 2000
bed_tss_orig["end"] = bed_tss_orig["end"] + 2000
bed_tss_orig = bed_tss_orig[bed_tss_orig.iloc[:,1] >= 0]
bed_tss = pybedtools.BedTool().from_dataframe(bed_tss_orig)

#Import the input bed, along with the histone ChIP-seq bed files
bed_input = pybedtools.BedTool(args.input)
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
active_promoter = active_promoter.sort().merge()

bivalent_promoter = overlaps_4me3.intersect(bed_27me3)
bivalent_promoter = bivalent_promoter.sort().merge()

silenced_promoter = no_4me3.intersect(bed_27me3)
silenced_promoter = silenced_promoter.sort().merge()

marked_promoters = active_promoter + bivalent_promoter + silenced_promoter

unclassified_overlaps_tss = overlaps_tss.intersect(marked_promoters, v = True)
unclassified_overlaps_tss = unclassified_overlaps_tss.sort().merge()

#Classifying putative enhancers
overlaps_4me1 = no_tss.intersect(bed_4me1)

active_enhancer = overlaps_4me1.intersect(bed_27ac)
active_enhancer = active_enhancer.sort().merge()

poised_enhancer = overlaps_4me1.intersect(bed_27me3)
poised_enhancer = poised_enhancer.sort().merge()

active_and_poised = active_enhancer + poised_enhancer
primed_enhancer = overlaps_4me1.intersect(active_and_poised, v = True)
primed_enhancer = primed_enhancer.sort().merge()

active_poised_primed = active_and_poised + primed_enhancer
unclassified_no_tss = no_tss.intersect(active_poised_primed, v = True)
unclassified_no_tss = unclassified_no_tss.sort().merge()

#Print out the number of peaks/coordinates in each category
if args.verbose:
	print("Identified regions:\nPutative Promoters\n{active_promoter_count} regions in an active promoter\n{bivalent_promoter_count} regions in a bivalent promoter\n{silenced_promoter_count} regions in a silenced promoter\n{unclassified_overlaps_tss_count} regions in an unclassified region within 2 kb of a TSS\n\nPutative Enhancers\n{active_enhancer_count} regions in an active enhancer\n{poised_enhancer_count} regions in a poised enhancer\n{primed_enhancer_count} in a primed enhancer\n{unclassified_no_tss_count} regions in an unclassified region greater than 2 kb from a TSS\n".format(active_promoter_count = active_promoter.count(), bivalent_promoter_count = bivalent_promoter.count(), silenced_promoter_count = silenced_promoter.count(), unclassified_overlaps_tss_count = unclassified_overlaps_tss.count(), active_enhancer_count = active_enhancer.count(), poised_enhancer_count = poised_enhancer.count(), primed_enhancer_count = primed_enhancer.count(), unclassified_no_tss_count = unclassified_no_tss.count()))
	print("\n")

#If they have one or more coordinate entries, convert each of the classifier variables to a pandas DataFrame object so that we can append the classifications as an additional column
if args.verbose:
	print("Preparing output file ({output_file})...".format(output_file = args.output))

if active_promoter.count() > 0:
	active_promoter_df = pybedtools.BedTool.to_dataframe(active_promoter)
	active_promoter_df["region_classification"] = "active_promoter"
else:
	active_promoter_df = pd.DataFrame(columns = ["chrom", "start", "end", "region_classification"])

if bivalent_promoter.count() > 0:
	bivalent_promoter_df = pybedtools.BedTool.to_dataframe(bivalent_promoter)
	bivalent_promoter_df["region_classification"] = "bivalent_promoter"
else:
	bivalent_promoter_df = pd.DataFrame(columns = ["chrom", "start", "end", "region_classification"])

if silenced_promoter.count() > 0:
	silenced_promoter_df = pybedtools.BedTool.to_dataframe(silenced_promoter)
	silenced_promoter_df["region_classification"] = "silenced_promoter"
else:
	silenced_promoter_df = pd.DataFrame(columns = ["chrom", "start", "end", "region_classification"])

if unclassified_overlaps_tss.count() > 0:
	unclassified_overlaps_tss_df = pybedtools.BedTool.to_dataframe(unclassified_overlaps_tss)
	unclassified_overlaps_tss_df["region_classification"] = "unclassified_within_2kb_of_tss"
else:
	unclassified_overlaps_tss_df = pd.DataFrame(columns = ["chrom", "start", "end", "region_classification"])

if active_enhancer.count() > 0:
	active_enhancer_df = pybedtools.BedTool.to_dataframe(active_enhancer)
	active_enhancer_df["region_classification"] = "active_enhancer"
else:
	active_enhancer_df = pd.DataFrame(columns = ["chrom", "start", "end", "region_classification"])

if poised_enhancer.count() > 0:
	poised_enhancer_df = pybedtools.BedTool.to_dataframe(poised_enhancer)
	poised_enhancer_df["region_classification"] = "poised_enhancer"
else:
	poised_enhancer_df = pd.DataFrame(columns = ["chrom", "start", "end", "region_classification"])

if primed_enhancer.count() > 0:
	primed_enhancer_df = pybedtools.BedTool.to_dataframe(primed_enhancer)
	primed_enhancer_df["region_classification"] = "primed_enhancer"
else:
	primed_enhancer_df = pd.DataFrame(columns = ["chrom", "start", "end", "region_classification"])

if unclassified_no_tss.count() > 0:
	unclassified_no_tss_df = pybedtools.BedTool.to_dataframe(unclassified_no_tss)
	unclassified_no_tss_df["region_classification"] = "unclassified_beyond_2kb_of_tss"
else:
	unclassified_no_tss_df = pd.DataFrame(columns = ["chrom", "start", "end", "region_classification"])
	
#Concatenate all of the pandas DataFrame objects into one
all_data_frames = [active_promoter_df, bivalent_promoter_df, silenced_promoter_df, unclassified_overlaps_tss_df, active_enhancer_df, poised_enhancer_df, primed_enhancer_df, unclassified_no_tss_df]

concatenated_df = pd.concat(all_data_frames)

#Convert the concatenated data frame back to pybedtools objects and sort by coordinates
concatenated_bed = pybedtools.BedTool().from_dataframe(concatenated_df).sort()
concatenated_bed.saveas(args.output)

###Developed by Michael J Betti, April 2021