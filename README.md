# CoRE-BED
## Introduction
CoRE-BED, or Classifier of Regulatory Elements in a BED file, is a Python-based tool developed to determine whether coordinates in a UCSC BED file fall within a putative cis-regulatory element (promoter or enhancer). Based on each coordinate's overlap (or lack thereof) with annotated transcription start sites (TSSs), as well as H3K4me1, H3K4me3, H3K27ac, and H3K27me3 histone marks, it will be classified as either an active promoter, bivalent promoter, silenced promoter, unclassified element within 2 kb of a TSS, active enhancer, poised enhancer, primed enhancer, or unclassified element not within 2 kb of a TSS. This classification method can be best visualized as a tree structure:

<img width="1023" alt="Screen Shot 2021-04-27 at 11 04 58 AM" src="https://user-images.githubusercontent.com/63562495/116274604-722d7e00-a748-11eb-8605-98c98d2e48b4.png">
