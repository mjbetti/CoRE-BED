import os, sys, requests, argparse, pybedtools, pandas as pd
from itertools import chain

#Download the appropriate reference files
def download_ref(url_ref, out_name):
	out_path = os.path.join("ref_files", out_name)
	r = requests.get(url_ref, allow_redirects = True)
	open(out_path, 'wb').write(r.content)

if not os.path.exists("ref_files"):
	os.mkdir("ref_files")
	
tissue_array = ["adipose", "adrenal_gland", "artery", "blood", "breast", "cultured_fibroblast", "ebv_transformed_lymphocyte", "es", "esophagus_muscularis_mucosa", "esophagus_squamous_epithelium", "heart", "intestine", "ips", "kidney", "liver", "lung", "neuron", "ovary", "pancreas", "prostate", "skeletal_muscle", "skin", "spleen", "stomach", "testis", "thyroid", "uterus", "vagina"]

###Human###
#GRCh38/hg38
download_ref("http://reftss.clst.riken.jp/datafiles/3.1/human/refTSS_v3.1_human_coordinate.hg38.bed.gz", "refTSS_v3.1_human_coordinate.hg38.bed.gz")
	
		#Adipose (Homo sapiens subcutaneous abdominal adipose tissue tissue nuclear fraction female adult (49 years), Homo sapiens omental fat pad tissue female adult (51 years) (DNase-seq and CTCF ChIP-seq only))
download_ref("https://www.encodeproject.org/files/ENCFF534VNL/@@download/ENCFF534VNL.bed.gz", "adipose_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF591SLF/@@download/ENCFF591SLF.bed.gz", "adipose_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF658NHX/@@download/ENCFF658NHX.bed.gz", "adipose_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF435MSC/@@download/ENCFF435MSC.bed.gz", "adipose_27me3_hg38.bed.gz")	
			#download_ref("https://www.encodeproject.org/files/ENCFF075XPG/@@download/ENCFF075XPG.bed.gz", "adipose_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF464DDN/@@download/ENCFF464DDN.bed.gz", "adipose_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF436YWU/@@download/ENCFF436YWU.bed.gz", "adipose_dnase_hg38.bed.gz")
		
		#Adrenal gland (Homo sapiens adrenal gland tissue male embryo (97 days), Homo sapiens adrenal gland tissue embryo (96 days) (DNase-seq only), and Homo sapiens adrenal gland tissue male adult (37 years) (CTCF ChIP-seq only))
download_ref("https://www.encodeproject.org/files/ENCFF953KZN/@@download/ENCFF953KZN.bed.gz", "adrenal_gland_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF959FFU/@@download/ENCFF959FFU.bed.gz", "adrenal_gland_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF231NZU/@@download/ENCFF231NZU.bed.gz", "adrenal_gland_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF852DFJ/@@download/ENCFF852DFJ.bed.gz", "adrenal_gland_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF852DFJ/@@download/ENCFF852DFJ.bed.gz", "adrenal_gland_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF297BQF/@@download/ENCFF297BQF.bed.gz", "adrenal_gland_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF949ZRW/@@download/ENCFF949ZRW.bed.gz", "adrenal_gland_dnase_hg38.bed.gz")
			
		#Artery (Homo sapiens aorta tissue male adult (34 years), Homo sapiens aorta tissue female adult (41 years) (DNase-seq only), and Homo sapiens ascending aorta tissue female adult (51 years) (CTCF ChIP-seq only))
download_ref("https://www.encodeproject.org/files/ENCFF824OFK/@@download/ENCFF824OFK.bed.gz", "artery_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF876LPD/@@download/ENCFF876LPD.bed.gz", "artery_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF872XZM/@@download/ENCFF872XZM.bed.gz", "artery_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF093PYF/@@download/ENCFF093PYF.bed.gz", "artery_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF587NCS/@@download/ENCFF587NCS.bed.gz", "artery_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF218JYZ/@@download/ENCFF218JYZ.bed.gz", "artery_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF829CJM/@@download/ENCFF829CJM.bed.gz", "artery_dnase_hg38.bed.gz")

		#Blood (Homo sapiens peripheral blood mononuclear cell male adult (39 years) and Homo sapiens K562 (DNase-seq and CTCF ChIP-seq only))	
download_ref("https://www.encodeproject.org/files/ENCFF734KLT/@@download/ENCFF734KLT.bed.gz", "blood_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF165VDC/@@download/ENCFF165VDC.bed.gz", "blood_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF844OYX/@@download/ENCFF844OYX.bed.gz", "blood_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF636JEP/@@download/ENCFF636JEP.bed.gz", "blood_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF422UFO/@@download/ENCFF422UFO.bed.gz", "blood_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF582SNT/@@download/ENCFF582SNT.bed.gz", "blood_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF274YGF/@@download/ENCFF274YGF.bed.gz", "blood_dnase_hg38.bed.gz")
	
		#Breast (Homo sapiens breast epithelium tissue female adult (53 years) and Homo sapiens breast epithelium tissue female adult (51 years) (DNase-seq only))
download_ref("https://www.encodeproject.org/files/ENCFF698LQG/@@download/ENCFF698LQG.bed.gz", "breast_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF590DGH/@@download/ENCFF590DGH.bed.gz", "breast_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF507WEL/@@download/ENCFF507WEL.bed.gz", "breast_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF962YZN/@@download/ENCFF962YZN.bed.gz", "breast_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF540RWM/@@download/ENCFF540RWM.bed.gz", "breast_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF910ACN/@@download/ENCFF910ACN.bed.gz", "breast_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF594NFE/@@download/ENCFF594NFE.bed.gz", "breast_dnase_hg38.bed.gz")
		
		#Cultured fibroblast (Homo sapiens IMR-90)
download_ref("https://www.encodeproject.org/files/ENCFF611UWF/@@download/ENCFF611UWF.bed.gz", "cultured_fibroblast_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF093NQC/@@download/ENCFF093NQC.bed.gz", "cultured_fibroblast_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF805GNH/@@download/ENCFF805GNH.bed.gz", "cultured_fibroblast_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF336IXL/@@download/ENCFF336IXL.bed.gz", "cultured_fibroblast_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF449ADN/@@download/ENCFF449ADN.bed.gz", "cultured_fibroblast_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF203SRF/@@download/ENCFF203SRF.bed.gz", "cultured_fibroblast_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF800DVI/@@download/ENCFF800DVI.bed.gz", "cultured_fibroblast_dnase_hg38.bed.gz")
		
		#EBV-transformed lymphocyte (Homo sapiens GM12878)
download_ref("https://www.encodeproject.org/files/ENCFF321BVG/@@download/ENCFF321BVG.bed.gz", "ebv_transformed_lymphocyte_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF998CEU/@@download/ENCFF998CEU.bed.gz", "ebv_transformed_lymphocyte_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF023LTU/@@download/ENCFF023LTU.bed.gz", "ebv_transformed_lymphocyte_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF291DHI/@@download/ENCFF291DHI.bed.gz", "ebv_transformed_lymphocyte_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF432EMI/@@download/ENCFF432EMI.bed.gz", "ebv_transformed_lymphocyte_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF796WRU/@@download/ENCFF796WRU.bed.gz", "ebv_transformed_lymphocyte_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF759OLD/@@download/ENCFF759OLD.bed.gz", "ebv_transformed_lymphocyte_dnase_hg38.bed.gz")
		
		#Embryonic stem cell (Homo sapiens H1)
download_ref("https://www.encodeproject.org/files/ENCFF613QAB/@@download/ENCFF613QAB.bed.gz", "es_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF356MCC/@@download/ENCFF356MCC.bed.gz", "es_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF689CJG/@@download/ENCFF689CJG.bed.gz", "es_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF050CUG/@@download/ENCFF050CUG.bed.gz", "es_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF813VFV/@@download/ENCFF813VFV.bed.gz", "es_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF692RPA/@@download/ENCFF692RPA.bed.gz", "es_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF983UCL/@@download/ENCFF983UCL.bed.gz", "es_dnase_hg38.bed.gz")
	
		#Esophagus muscularis mucosa (Homo sapiens esophagus muscularis mucosa tissue female adult (51 years) and Homo sapiens esophagus muscularis mucosa tissue male adult (37 years) (DNase-seq only))
download_ref("https://www.encodeproject.org/files/ENCFF473ELV/@@download/ENCFF473ELV.bed.gz", "esophagus_muscularis_mucosa_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF154TIT/@@download/ENCFF154TIT.bed.gz", "esophagus_muscularis_mucosa_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF114SOB/@@download/ENCFF114SOB.bed.gz", "esophagus_muscularis_mucosa_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF450WVA/@@download/ENCFF450WVA.bed.gz", "esophagus_muscularis_mucosa_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF101WAE/@@download/ENCFF101WAE.bed.gz", "esophagus_muscularis_mucosa_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF851XPN/@@download/ENCFF851XPN.bed.gz", "esophagus_muscularis_mucosa_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF726TPF/@@download/ENCFF726TPF.bed.gz", "esophagus_muscularis_mucosa_dnase_hg38.bed.gz")
		
		#Esophagus squamous epithelium (Homo sapiens esophagus squamous epithelium tissue female adult (51 years) and Homo sapiens esophagus squamous epithelium tissue male adult (37 years) (DNase-seq only))
download_ref("https://www.encodeproject.org/files/ENCFF828YQS/@@download/ENCFF828YQS.bed.gz", "esophagus_squamous_epithelium_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF470ELM/@@download/ENCFF470ELM.bed.gz", "esophagus_squamous_epithelium_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF673PCP/@@download/ENCFF673PCP.bed.gz", "esophagus_squamous_epithelium_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF905MKP/@@download/ENCFF905MKP.bed.gz", "esophagus_squamous_epithelium_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF885BKA/@@download/ENCFF885BKA.bed.gz", "esophagus_squamous_epithelium_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF384CCC/@@download/ENCFF384CCC.bed.gz", "esophagus_squamous_epithelium_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF839REP/@@download/ENCFF839REP.bed.gz", "esophagus_squamous_epithelium_dnase_hg38.bed.gz")
		
		#Heart (Homo sapiens heart left ventricle tissue female adult (53 years))
download_ref("https://www.encodeproject.org/files/ENCFF194KSV/@@download/ENCFF194KSV.bed.gz", "heart_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF008ODN/@@download/ENCFF008ODN.bed.gz", "heart_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF635VBV/@@download/ENCFF635VBV.bed.gz", "heart_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF242JDO/@@download/ENCFF242JDO.bed.gz", "heart_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF767VOC/@@download/ENCFF767VOC.bed.gz", "heart_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF165RCX/@@download/ENCFF165RCX.bed.gz", "heart_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF846JLR/@@download/ENCFF846JLR.bed.gz", "heart_dnase_hg38.bed.gz")
		
		#Intestine (Homo sapiens sigmoid colon tissue male adult (37 years) and Homo sapiens sigmoid colon tissue female adult (53 years) (DNase-seq only))
download_ref("https://www.encodeproject.org/files/ENCFF052XGA/@@download/ENCFF052XGA.bed.gz", "intestine_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF250ZHO/@@download/ENCFF250ZHO.bed.gz", "intestine_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF111SYA/@@download/ENCFF111SYA.bed.gz", "intestine_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF310IST/@@download/ENCFF310IST.bed.gz", "intestine_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF406OXE/@@download/ENCFF406OXE.bed.gz", "intestine_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF480QYI/@@download/ENCFF480QYI.bed.gz", "intestine_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF274XSP/@@download/ENCFF274XSP.bed.gz", "intestine_dnase_hg38.bed.gz")
	
		#iPS cell (Homo sapiens iPS DF 19.11 and Homo sapiens GM23338 originated from GM23248 (CTCF ChIP-seq only))
download_ref("https://www.encodeproject.org/files/ENCFF314LSR/@@download/ENCFF314LSR.bed.gz", "ips_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF450VVJ/@@download/ENCFF450VVJ.bed.gz", "ips_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF037RXP/@@download/ENCFF037RXP.bed.gz", "ips_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF844FEB/@@download/ENCFF844FEB.bed.gz", "ips_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF855HXT/@@download/ENCFF855HXT.bed.gz", "ips_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF651DQL/@@download/ENCFF651DQL.bed.gz", "ips_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF172GXY/@@download/ENCFF172GXY.bed.gz", "ips_dnase_hg38.bed.gz")
		
		#Kidney (Homo sapiens kidney tissue male adult (67 years), Homo sapiens kidney tissue female embryo (120 days) for H3K27me3 - only dataset for kidney and DNase-seq, and Homo sapiens kidney epithelial cell (CTCF ChIP-seq))
download_ref("https://www.encodeproject.org/files/ENCFF249OJK/@@download/ENCFF249OJK.bed.gz", "kidney_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF430ZMN/@@download/ENCFF430ZMN.bed.gz", "kidney_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF512JNT/@@download/ENCFF512JNT.bed.gz", "kidney_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF775YUI/@@download/ENCFF775YUI.bed.gz", "kidney_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF310IRC/@@download/ENCFF310IRC.bed.gz", "kidney_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF059SDM/@@download/ENCFF059SDM.bed.gz", "kidney_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF388NNQ/@@download/ENCFF388NNQ.bed.gz", "kidney_dnase_hg38.bed.gz")

		#Liver (Homo sapiens liver tissue male adult (32 years), Homo sapiens liver tissue female embryo (101 days) and female embryo (113 days) (DNase-seq only), and Homo sapiens right lobe of liver tissue female adult (53 years) (CTCF ChIP-seq only))
download_ref("https://www.encodeproject.org/files/ENCFF625TWQ/@@download/ENCFF625TWQ.bed.gz", "liver_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF830ODB/@@download/ENCFF830ODB.bed.gz", "liver_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF287VIA/@@download/ENCFF287VIA.bed.gz", "liver_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF287VIA/@@download/ENCFF287VIA.bed.gz", "liver_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF666OGH/@@download/ENCFF666OGH.bed.gz", "liver_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF046XEZ/@@download/ENCFF046XEZ.bed.gz", "liver_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF933OGA/@@download/ENCFF933OGA.bed.gz", "liver_dnase_hg38.bed.gz")
		
		#Lung (Homo sapiens upper lobe of left lung tissue female adult (51 years))
download_ref("https://www.encodeproject.org/files/ENCFF926YKP/@@download/ENCFF926YKP.bed.gz", "lung_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF843JPR/@@download/ENCFF843JPR.bed.gz", "lung_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF496RDA/@@download/ENCFF496RDA.bed.gz", "lung_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF248IPH/@@download/ENCFF248IPH.bed.gz", "lung_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF459IYH/@@download/ENCFF459IYH.bed.gz", "lung_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF629UDP/@@download/ENCFF629UDP.bed.gz", "lung_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF041XPT/@@download/ENCFF041XPT.bed.gz", "lung_dnase_hg38.bed.gz")
		
		#Neuron (Homo sapiens SK-N-SH)
download_ref("https://www.encodeproject.org/files/ENCFF632FAM/@@download/ENCFF632FAM.bed.gz", "neuron_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF682JYE/@@download/ENCFF682JYE.bed.gz", "neuron_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF138VUT/@@download/ENCFF138VUT.bed.gz", "neuron_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF277NRX/@@download/ENCFF277NRX.bed.gz", "neuron_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF140RVM/@@download/ENCFF140RVM.bed.gz", "neuron_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF608JKY/@@download/ENCFF608JKY.bed.gz", "neuron_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF752OZB/@@download/ENCFF752OZB.bed.gz", "neuron_dnase_hg38.bed.gz")
		
		#Ovary (Homo sapiens ovary tissue female adult (30 years) and Homo sapiens ovary tissue female adult (53 years) (CTCF ChIP-seq only))
download_ref("https://www.encodeproject.org/files/ENCFF458XZT/@@download/ENCFF458XZT.bed.gz", "ovary_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF008MLA/@@download/ENCFF008MLA.bed.gz", "ovary_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF328KKO/@@download/ENCFF328KKO.bed.gz", "ovary_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF495DAT/@@download/ENCFF495DAT.bed.gz", "ovary_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF558KMA/@@download/ENCFF558KMA.bed.gz", "ovary_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF782GIY/@@download/ENCFF782GIY.bed.gz", "ovary_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF418PRB/@@download/ENCFF418PRB.bed.gz", "ovary_dnase_hg38.bed.gz")
		
		#Pancreas (Homo sapiens pancreas tissue female adult (41 years))
download_ref("https://www.encodeproject.org/files/ENCFF453GYD/@@download/ENCFF453GYD.bed.gz", "pancreas_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF200ITL/@@download/ENCFF200ITL.bed.gz", "pancreas_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF498SSN/@@download/ENCFF498SSN.bed.gz", "pancreas_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF663DVO/@@download/ENCFF663DVO.bed.gz", "pancreas_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF456HUQ/@@download/ENCFF456HUQ.bed.gz", "pancreas_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF682RNL/@@download/ENCFF682RNL.bed.gz", "pancreas_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF949DIO/@@download/ENCFF949DIO.bed.gz", "pancreas_dnase_hg38.bed.gz")
		
		#Prostate (Homo sapiens prostate gland tissue male adult (37 years))
download_ref("https://www.encodeproject.org/files/ENCFF275KSR/@@download/ENCFF275KSR.bed.gz", "prostate_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF155ZYZ/@@download/ENCFF155ZYZ.bed.gz", "prostate_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF201VZW/@@download/ENCFF201VZW.bed.gz", "prostate_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF198QRW/@@download/ENCFF198QRW.bed.gz", "prostate_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF028CPL/@@download/ENCFF028CPL.bed.gz", "prostate_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF351NKT/@@download/ENCFF351NKT.bed.gz", "prostate_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF233XNV/@@download/ENCFF233XNV.bed.gz", "prostate_dnase_hg38.bed.gz")
		
		#Skeletal muscle (Homo sapiens muscle of leg tissue female embryo (110 days), Homo sapiens muscle of leg tissue female embryo (115 days) (DNase-seq only), and Homo sapiens gastrocnemius medialis tissue male adult (37 years) (CTCF ChIP-seq only))
download_ref("https://www.encodeproject.org/files/ENCFF170ZFK/@@download/ENCFF170ZFK.bed.gz", "skeletal_muscle_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF942MOM/@@download/ENCFF942MOM.bed.gz", "skeletal_muscle_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF154ZCF/@@download/ENCFF154ZCF.bed.gz", "skeletal_muscle_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF672MZK/@@download/ENCFF672MZK.bed.gz", "skeletal_muscle_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF495BIY/@@download/ENCFF495BIY.bed.gz", "skeletal_muscle_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF143NGY/@@download/ENCFF143NGY.bed.gz", "skeletal_muscle_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF055FYP/@@download/ENCFF055FYP.bed.gz", "skeletal_muscle_dnase_hg38.bed.gz")

		#Skin (Homo sapiens suprapubic skin tissue male adult (37 years) and Homo sapiens foreskin fibroblast male newborn (DNase-seq only))
download_ref("https://www.encodeproject.org/files/ENCFF800MDH/@@download/ENCFF800MDH.bed.gz", "skin_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF240CFZ/@@download/ENCFF240CFZ.bed.gz", "skin_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF398EEO/@@download/ENCFF398EEO.bed.gz", "skin_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF686GBZ/@@download/ENCFF686GBZ.bed.gz", "skin_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF804VWO/@@download/ENCFF804VWO.bed.gz", "skin_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF503GIR/@@download/ENCFF503GIR.bed.gz", "skin_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF119TMJ/@@download/ENCFF119TMJ.bed.gz", "skin_dnase_hg38.bed.gz")
		
		#Spleen (Homo sapiens spleen tissue male adult (37 years) and Homo sapiens spleen tissue embryo (112 days) (DNase-seq only))
download_ref("https://www.encodeproject.org/files/ENCFF652QYY/@@download/ENCFF652QYY.bed.gz", "spleen_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF401XEH/@@download/ENCFF401XEH.bed.gz", "spleen_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF322YIK/@@download/ENCFF322YIK.bed.gz", "spleen_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF159IED/@@download/ENCFF159IED.bed.gz", "spleen_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF380BZO/@@download/ENCFF380BZO.bed.gz", "spleen_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF455BVG/@@download/ENCFF455BVG.bed.gz", "spleen_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF686XTC/@@download/ENCFF686XTC.bed.gz", "spleen_dnase_hg38.bed.gz")
		
		#Stomach (Homo sapiens stomach tissue male adult (37 years) and Homo sapiens stomach tissue male adult (34 years) (DNase-seq only))
download_ref("https://www.encodeproject.org/files/ENCFF475EJU/@@download/ENCFF475EJU.bed.gz", "stomach_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF494BWU/@@download/ENCFF494BWU.bed.gz", "stomach_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF978AOD/@@download/ENCFF978AOD.bed.gz", "stomach_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF692HXY/@@download/ENCFF692HXY.bed.gz", "stomach_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF243EAW/@@download/ENCFF243EAW.bed.gz", "stomach_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF533OIH/@@download/ENCFF533OIH.bed.gz", "stomach_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF033LEB/@@download/ENCFF033LEB.bed.gz", "stomach_dnase_hg38.bed.gz")
		
		#Testis (Homo sapiens testis tissue male adult (37 years) and Homo sapiens testis tissue male adult (54 years) (DNase-seq only))
download_ref("https://www.encodeproject.org/files/ENCFF600MJO/@@download/ENCFF600MJO.bed.gz", "testis_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF612RSR/@@download/ENCFF612RSR.bed.gz", "testis_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF540YUQ/@@download/ENCFF540YUQ.bed.gz", "testis_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF389XVY/@@download/ENCFF389XVY.bed.gz", "testis_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF469ASG/@@download/ENCFF469ASG.bed.gz", "testis_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF387GJA/@@download/ENCFF387GJA.bed.gz", "testis_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF914YUP/@@download/ENCFF914YUP.bed.gz", "testis_dnase_hg38.bed.gz")
		
		#Thyroid (Homo sapiens thyroid gland tissue female adult (51 years))
download_ref("https://www.encodeproject.org/files/ENCFF346IPY/@@download/ENCFF346IPY.bed.gz", "thyroid_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF925FSR/@@download/ENCFF925FSR.bed.gz", "thyroid_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF138WNV/@@download/ENCFF138WNV.bed.gz", "thyroid_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF607NRZ/@@download/ENCFF607NRZ.bed.gz", "thyroid_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF922QTV/@@download/ENCFF922QTV.bed.gz", "thyroid_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF728ZGC/@@download/ENCFF728ZGC.bed.gz", "thyroid_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF336FSX/@@download/ENCFF336FSX.bed.gz", "thyroid_dnase_hg38.bed.gz")
		
		#Uterus (Homo sapiens uterus tissue female adult (51 years) and Homo sapiens uterus tissue female adult (53 years) (DNase-seq only))
download_ref("https://www.encodeproject.org/files/ENCFF150BRF/@@download/ENCFF150BRF.bed.gz", "uterus_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF022DOY/@@download/ENCFF022DOY.bed.gz", "uterus_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF108SAD/@@download/ENCFF108SAD.bed.gz", "uterus_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF254GYM/@@download/ENCFF254GYM.bed.gz", "uterus_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF374OMB/@@download/ENCFF374OMB.bed.gz", "uterus_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF447GBU/@@download/ENCFF447GBU.bed.gz", "uterus_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF661ERT/@@download/ENCFF661ERT.bed.gz", "uterus_dnase_hg38.bed.gz")
		
		#Vagina (Homo sapiens vagina tissue female adult (51 years))
download_ref("https://www.encodeproject.org/files/ENCFF005TOV/@@download/ENCFF005TOV.bed.gz", "vagina_4me1_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF653IVY/@@download/ENCFF653IVY.bed.gz", "vagina_4me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF165SPP/@@download/ENCFF165SPP.bed.gz", "vagina_27ac_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF866QFX/@@download/ENCFF866QFX.bed.gz", "vagina_27me3_hg38.bed.gz")
			#download_ref("https://www.encodeproject.org/files/ENCFF823VSG/@@download/ENCFF823VSG.bed.gz", "vagina_36me3_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF626XNO/@@download/ENCFF626XNO.bed.gz", "vagina_ctcf_hg38.bed.gz")
download_ref("https://www.encodeproject.org/files/ENCFF872AJU/@@download/ENCFF872AJU.bed.gz", "vagina_dnase_hg38.bed.gz")
	
###Import all of the bed files to be worked on (input and histone mark ChIP-seq) using pybedtools
##Human
bed_ref = "hg38"
bed_tss_path = "ref_files/refTSS_v3.1_human_coordinate.hg38.bed.gz"

#Manipulate the TSS bed file so that the coordinates are expanded -5 kb and +1 kb
bed_tss_orig = pd.read_csv(bed_tss_path, sep = "\t", header = None)
bed_tss_orig.iloc[:,1] = bed_tss_orig.iloc[:,1] - 5000
bed_tss_orig.iloc[:,2] = bed_tss_orig.iloc[:,2] + 1000
bed_tss_orig = bed_tss_orig[bed_tss_orig.iloc[:,1] >= 0]
bed_tss = pybedtools.BedTool().from_dataframe(bed_tss_orig)

#Loop through each reference tissue parsed in from the input arguments, generating contextual functional annotations and appending each set as a new column onto the original input file
for tissue in tissue_array:	
	bed_4me1 = pybedtools.BedTool("ref_files/" + tissue + "_4me1_" + bed_ref + ".bed.gz")
	bed_4me3 = pybedtools.BedTool("ref_files/" + tissue + "_4me3_" + bed_ref + ".bed.gz")
	bed_27ac = pybedtools.BedTool("ref_files/" + tissue + "_27ac_" + bed_ref + ".bed.gz")
	bed_27me3 = pybedtools.BedTool("ref_files/" + tissue + "_27me3_" + bed_ref + ".bed.gz")

	#Classifying putative promoters
	tss_4me3 = bed_tss.intersect(bed_4me3)
	tss_no_4me3 = bed_tss.intersect(bed_4me3, v = True)

	active_promoter = tss_4me3.intersect(bed_27me3, v = True)
	active_promoter = active_promoter.sort().merge()

	bivalent_promoter = tss_4me3.intersect(bed_27me3)
	bivalent_promoter = bivalent_promoter.sort().merge()

	silenced_promoter = tss_no_4me3.intersect(bed_27me3)
	silenced_promoter = silenced_promoter.sort().merge()

	marked_promoters = active_promoter + bivalent_promoter + silenced_promoter

	#Classifying putative enhancers
	no_tss_4me1 = bed_4me1.intersect(bed_tss, v = True)

	active_enhancer = no_tss_4me1.intersect(bed_27ac)
	active_enhancer = active_enhancer.sort().merge()

	poised_enhancer = no_tss_4me1.intersect(bed_27me3)
	poised_enhancer = poised_enhancer.sort().merge()

	active_and_poised = active_enhancer + poised_enhancer
	primed_enhancer = no_tss_4me1.intersect(active_and_poised, v = True)
	primed_enhancer = primed_enhancer.sort().merge()

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
	
	#Concatenate all of the pandas DataFrame objects into one, remove all duplicate lines, and sort by coordinate
	all_data_frames = [active_promoter_df, bivalent_promoter_df, silenced_promoter_df, active_enhancer_df, poised_enhancer_df, primed_enhancer_df]

	concatenated_df = pd.concat(all_data_frames)
	col_names = ["chr", "start", "end", "func_anno"]
	concatenated_df.columns = col_names
	concatenated_df = concatenated_df.drop_duplicates()
	concatenated_df.iloc[:,0] = concatenated_df.iloc[:,0].str[3:]

	#Output the final output file of all regulatory elements in the tissue type
	out_name = tissue + "_all_regulatory_elements_5000_1000.txt"
	concatenated_df.to_csv(out_name, sep = "\t", header = False, index = False)