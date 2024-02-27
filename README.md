# marsden_obj1_preservationMethod

## 1. Introduction

This document serves as a comprehensive guide and supplementary information for the bioinformatic and statistical analysis of the manuscript entitled "Unlocking Antarctic molecular time-capsules - recovering historical environmental DNA from museum-preserved sponges" by Jeunen *et al*., 2024. This work is associated with the Marsden Fast-Start fund (MFP-UOO002116).

## 2. Environment set up

### 2.1 Starting files

The bioinformatic and statistical analysis requires access to several files, including the sequencing data and metadata files.

The two raw g-zipped sequencing files can be downloaded from the Sequence Read Archive (SRA) project number PRJNA1019816 and submission number SUB13856226. The files names are: "**8373-P1-00-01_S1_L001_R1_001.fastq.gz**" and "**8373-P2-00-01_S1_L001_R1_001.fastq.gz**". Download both sequencing data files and place them in the starting folder.

Besides the sequencing data, we also need two metadatafiles, both of which can be downloaded from this GitHub repository in the folder "**metadata_files**". The first metadata file named "**barcodeMetadata8373.fasta**" will be used for demultiplexing, while the second metadata file named "**sampleMetadata8373.csv**" will be used for the statistical analysis. Place both metadata files into the starting folder.

In the Terminal, move to the starting folder. When listing the files, we should have access to all starting files.

```{code-block} bash
cd move/to/starting/folder
ls -ltr
```
