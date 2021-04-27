# CoRE-BED
## Introduction
CoRE-BED, or Classifier of Regulatory Elements in a BED file, is a Python-based tool developed to determine whether coordinates in a UCSC BED file fall within a putative cis-regulatory element (promoter or enhancer). Based on each coordinate's overlap (or lack thereof) with annotated transcription start sites (TSSs), as well as H3K4me1, H3K4me3, H3K27ac, and H3K27me3 histone marks, it will be classified as either an active promoter, bivalent promoter, silenced promoter, unclassified element within 2 kb of a TSS, active enhancer, poised enhancer, primed enhancer, or unclassified element not within 2 kb of a TSS. This classification method can be best visualized as a tree structure:

<img width="1023" alt="Screen Shot 2021-04-27 at 11 04 58 AM" src="https://user-images.githubusercontent.com/63562495/116274604-722d7e00-a748-11eb-8605-98c98d2e48b4.png">

The main script ```core-bed.py``` takes in a user-provided UCSC BED file and downloads all reference files (transcription start site coordinates, as well as H3K4me1, H3K4me3, H3K27ac, and H3K27me3 histone ChIP-seq peak coordinates generated by the ENCODE Consortium) based on the user-specified human genome build and tissue type. The following biological samples from ENCODE were used to represent each of the 10 supported tissue types:
* Blood - *Homo sapiens* peripheral blood mononuclear cell male adult (39 years)
* Brain - *Homo sapiens* SK-N-SH
* Embryonic stem cell - *Homo sapiens* H1
* Heart - *Homo sapiens* heart left ventricle tissue female adult (53 years)
* Intestine - *Homo sapiens* sigmoid colon tissue male adult (37 years)
* Liver - *Homo sapiens* liver tissue male adult (32 years)
* Lung - *Homo sapiens* upper lobe of left lung tissue female adult (51 years)
* Muscle - *Homo sapiens* muscle of leg tissue female embryo (110 days)
* Skin - *Homo sapiens* suprapubic skin tissue male adult (37 years)
* Thyroid - *Homo sapiens* thyroid gland tissue female adult (51 years)

Upon preparation of all reference files, the BED containing transcription start sites (TSSs) is modified by subtracting 2 kb from each start coordinate and adding 2 kb to each end coordinate. The script then uses these modified TSS coordinates to determine whether each input coordinate overlaps with a genomic interval within 2 kb of a TSS. Those input coordinates that do have overlap are considered "candidate promoters" and those that do not are considered "candidate enhancers."

For the set of candidate promoters, the script next determines 
