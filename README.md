# CoRE-BED
## Introduction
CoRE-BED, or Classifier of Regulatory Elements in a BED file, is a Python-based tool developed to determine whether coordinates in a UCSC BED file fall within a putative cis-regulatory element (promoter or enhancer). Based on each coordinate's overlap (or lack thereof) with annotated transcription start sites (TSSs), as well as H3K4me1, H3K4me3, H3K27ac, and H3K27me3 histone marks, it will be classified as either an active promoter, bivalent promoter, silenced promoter, unclassified element within 2 kb of a TSS, active enhancer, poised enhancer, primed enhancer, or unclassified element not within 2 kb of a TSS. This classification method can be best visualized as a tree structure:

<img width="1023" alt="Screen Shot 2021-04-27 at 11 04 58 AM" src="https://user-images.githubusercontent.com/63562495/116274604-722d7e00-a748-11eb-8605-98c98d2e48b4.png">

The main script ```core-bed.py``` takes in a user-provided UCSC BED file and downloads all reference files (transcription start site coordinates, as well as H3K4me1, H3K4me3, H3K27ac, and H3K27me3 histone ChIP-seq peak coordinates generated by the ENCODE Consortium) based on the user-specified human genome build and tissue type. The following biological samples from ENCODE were used to represent each of the 10 supported tissue types:
* ```Blood``` - *Homo sapiens* peripheral blood mononuclear cell male adult (39 years)
* ```Brain``` - *Homo sapiens* SK-N-SH
* ```Embryonic stem cell``` - *Homo sapiens* H1
* ```Heart``` - *Homo sapiens* heart left ventricle tissue female adult (53 years)
* ```Intestine``` - *Homo sapiens* sigmoid colon tissue male adult (37 years)
* ```Liver``` - *Homo sapiens* liver tissue male adult (32 years)
* ```Lung``` - *Homo sapiens* upper lobe of left lung tissue female adult (51 years)
* ```Muscle``` - *Homo sapiens* muscle of leg tissue female embryo (110 days)
* ```Skin``` - *Homo sapiens* suprapubic skin tissue male adult (37 years)
* ```Thyroid``` - *Homo sapiens* thyroid gland tissue female adult (51 years)

Upon preparation of all reference files, the BED containing transcription start sites (TSSs) is modified by subtracting 2 kb from each start coordinate and adding 2 kb to each end coordinate. The script then uses these modified TSS coordinates to determine whether each input coordinate overlaps with a genomic interval within 2 kb of a TSS. Those input coordinates that do have overlap are considered "candidate promoters" and those that do not are considered "candidate enhancers."

For the set of candidate promoters, the script next determines whether each coordinate overlaps with the H3K4me3, a histone modification commonly associated with promoter regions in active chromatin. Those that do overlap with H3K4me3 are subsequently evaluated for overlap with H3K27me3, a repressive histone mark. If a coordinate has both H3K4me3 and H3K27me3, it is classified as a "bivalent promoter," meaning that while not presently active, it is "poised" to become either active or more permanently repressed by losing one of the two histone marks. If no H3K27me3 is present, then the element is classified as an "active promoter."

Conversely, coordinates that do not have any overlap with H3K4me3 are next evaluated for overlap with H3K27me3. Elements that are marked with only H3K27me3 are classified as "silenced promoters," while those that do not overlap with either of the two marks are classified as an "unclassified element within 2 kb of a TSS."

With our set of candidate enhancers, CoRE-BED will search for overlap with three enhancer-associated histone marks, H3K4me1, H3K27ac, and H3K27me3. The script first determines whether there is any overlap with H3K4me1, a mark commonly used to denote the presence of an enhancer. If H3K4me1 is present, then we next search for co-localization of the active mark H3K27ac. If both marks localize at a single element, then we classify this as an "active enhancer." Otherwise, CoRE-BED will search for overlap with the repressive mark H3K27me3. Elements that show a co-localization of H3K4me1 and H3K27me3 will be classified as "poised enhancers," while those marked with only H3K4me1 are classified as "primed enhancers." The remaining candidate enhancer regions overlapping with none of these histone marks are classified as "unclassified elements beyond 2 kb of a TSS."

## Dependencies
All dependencies (including ```python``` (3.6.13), ```requests``` (2.25.1), ```pybedtools``` (0.7.10), and ```pandas``` (1.1.3)) can be easily installed as an Anaconda environment using the included ```core-bed_env.yml``` file:
```
conda env create -f core-bed_env.yml
```
...and then activated using the following command:
```
conda activate core-bed_env
```

## Usage
Use the ```-h``` or ```--help``` flag to view all available options:
```
python core-bed.py -h

usage: core-bed.py [-h] -i INPUT -g REF_GENOME -t TISSUE [-o OUTPUT] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        the input bed file (required)
  -g REF_GENOME, --ref_genome REF_GENOME
                        the human reference genome build on which the input
                        coordinates are based (required) (valid options:
                        GRCh38/hg38 and GRCh37/hg19)
  -t TISSUE, --tissue TISSUE
                        the tissue of interest (required) (valid options:
                        Blood, Brain, ES, Heart, Intestine, Liver, Lung,
                        Muscle, Skin, Thyroid)
  -o OUTPUT, --output OUTPUT
                        the name of the output file
  -v, --verbose         return logging as terminal output
```

A typical run of the CoRE-BED method would look something like the following:
```
python core-bed.py \
-i input_file_based_on_hg38.bed \
-g hg38 \
-t blood \
-o annotated_input_coordinates.bed \
-v
```
...with each of the user-specified arguments specifying:
* ```-i``` or ```--input``` - A set of genomic coordinates in UCSC BED format (https://genome.ucsc.edu/FAQ/FAQformat.html)
* ```-g``` or ```--ref_genome``` - The reference genome build on which the input coordinates are based
* ```-t``` or ```--tissue``` - The tissue type to which the input coordinates will be compared
* ```-o``` or ```--output``` - The desired name of the output file (default is ```out.bed```)
* ```-v``` or ```--verbose``` - Enable verbosity (i.e. print the script progress out to the console)

The final output file will be a new BED file, in which the input coordinates are annotated based on how the genomic region in which they fall was classified. The BED will consist of three data columns: chr, start, end, and annotation.

## References
Dale RK, Pedersen BS, Quinlan AR. Pybedtools: a flexible Python library for manipulating genomic datasets and annotations. Bioinformatics. 2011;27: 3423–3424.

Mas G, Blanco E, Ballaré C, Sansó M, Spill YG, Hu D, et al. Promoter bivalency favors an open chromatin architecture in embryonic stem cells. Nat Genet. 2018;50: 1452–1462.

McKinney, Wes. 2010. “Data Structures for Statistical Computing in Python.” Proceedings of the 9th Python in Science Conference. https://doi.org/10.25080/majora-92bf1922-00a.

Sha K, Boyer LA. The chromatin signature of pluripotent cells. StemBook. Cambridge (MA): Harvard Stem Cell Institute; 2009.

Shlyueva D, Stampfel G, Stark A. Transcriptional enhancers: from properties to genome-wide predictions. Nat Rev Genet. 2014;15: 272–286.
