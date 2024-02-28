# marsden_obj1_preservationMethod

## 1. Introduction

This document serves as a comprehensive guide and supplementary information for the bioinformatic and statistical analysis of the manuscript entitled "Unlocking Antarctic molecular time-capsules - recovering historical environmental DNA from museum-preserved sponges" by Jeunen *et al*., 2024. This work is associated with the Marsden Fast-Start fund (MFP-UOO002116).

## 2. Environment set up

### 2.1 Starting files

The bioinformatic and statistical analysis requires access to several files, including the sequencing data and metadata files.

The two raw g-zipped sequencing files can be downloaded from the Sequence Read Archive (SRA) project number PRJNA1019816 and submission number SUB13856226. The files names are: "**8373-P1-00-01_S1_L001_R1_001.fastq.gz**" and "**8373-P2-00-01_S1_L001_R1_001.fastq.gz**". Download both sequencing data files and place them in the starting folder.

Besides the sequencing data, two metadatafiles are necessary to conduct the bioinformatic and statistical analysis, both of which can be downloaded from this GitHub repository in the folder "**metadata_files**". The first metadata file named "**barcodeMetadata8373.fasta**" will be used for demultiplexing, while the second metadata file named "**sampleMetadata8373.csv**" will be used for the statistical analysis. Place both metadata files into the starting folder.

In the Terminal, move to the starting folder. When listing the files, all starting files should be listed in the Terminal window.

```{code-block} bash
cd move/to/starting/folder
ls -ltr
```

## 3. Bioinformatic analysis

### 3.1 Initial data check

Before starting the bioinformatic pipeline, check the structure and quality of the raw sequencing files using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). First, unzip the raw sequencing data files using `gunzip`.

```{code-block} bash
gunzip *.gz
```

Next, analyse the quality and structure of raw sequencing data files using `fastqc`. Both libraries were sequenced on a MiSeq using the sequencing kit V2 1x300 bp.

```{code-block} bash
mkdir fastqc_raw
fastqc *.fastq -t 8 -o fastqc_raw/
```

Samples were sequenced across two sequencing runs. To simplify the code, combine both files using the `cat` command.

```{code-block} bash
cat *.fastq > 8373Combined.fastq
```

### 3.2 Demultiplexing

#### 3.2.1 Adapter removal

Prior to demultiplexing, remove the Illumina sequencing adapter on the 3' end of the reads using [cutadapt](https://cutadapt.readthedocs.io/en/stable/). Removing these adapters enables the anchoring of the barcodes and primers for demultiplexing due to the library build specifications, i.e., reads will start and end with the barcode and primer sequences. Roughly 82% of reads should pass this filtering step.

```{code-block} bash
cutadapt 8373Combined.fastq -a ATCTCGTATGCCGTCTTCTGCTTG -o 8373Combinedp7removed.fastq --discard-untrimmed --no-indels --cores=0
```

#### 3.2.2 Assigning reads

Once Illumina adapters are removed, demultiplex data the data using [cutadapt](https://cutadapt.readthedocs.io/en/stable/). Note that in the metadata file, information starts with `^` and ends with `$`, indicating the anchoring of barcode sequences. Additionally, since primer regions need to be removed as well, the metadata file specifies the barcode + primer region as the adapter to be removed. Since multiple files will be generated during demultiplexing, one per sample, place newly generated files in a new folder using the code below. Roughly 65% of reads sshould pass this filtering step. Low number of demultiplexed reads is due to the library of these samples being pooled with other libraries in a single sequencing run.

```{code-block} bash
mkdir demux
cutadapt 8373Combinedp7removed.fastq -g file:'barcodeMetadata8373.fasta' -o demux/{name}.fastq --discard-untrimmed --no-indels -e 2 --cores=0
```

### 3.3 Quality filtering

Once reads have been demultiplexed, filter reads based on various quality parameters using [VSEARCH](https://github.com/torognes/vsearch). To combine all files post quality filtering, rename sequence headers to sample names.

```{code-block} bash
cd demux/
mkdir qual_fasta qual_fastq
for fq in *.fastq
do
echo "\n\n\nAnalysing: ${fq}"
fasta=${fq/.fastq/.fasta}
vsearch --fastq_filter ${fq} --fastq_maxee 1.0 --fastq_minlen 190 --fastq_maxlen 220 --fastq_maxns 0 --fastqout qual_fastq/${fq} --fastaout qual_fasta/${fasta} --relabel ${fq/.fastq/}.
done
```

Once quality filtering is completed, combine all files into a single sequence file using the `cat` command.

```{code-block} bash
cat qual_fastq/*.fastq > qual_fastq/combined.fastq
cat qual_fasta/*.fasta > qual_fasta/combined.fasta
```

Prior to dereplication, check quality of "**combined.fastq**" using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to verify quality filtering was successful.

```{code-block} bash
fastqc qual_fastq/combined.fastq -t 8 
```

### 3.4 Dereplication

Next, dereplicate the data using [VSEARCH](https://github.com/torognes/vsearch), i.e., find unique sequences in the "**combined.fasta**" file.

```{code-block} bash
cd qual_fasta
mkdir derep
vsearch --derep_fulllength combined.fasta --relabel uniq. --output derep/uniques.fasta --sizeout
```

### 3.5 Denoising

Using the unique sequences, denoise the remaining reads to find all biologically relevant sequences using [USEARCH](https://www.drive5.com/usearch/). USEARCH will automatically remove chimeric sequences during this step.

```{code-block} bash
cd derep/
usearch -unoise3 uniques.fasta -zotus asv.fasta
```

### 3.6 Count table

Once a list of ASV's is generated, match this list to the sequence data to create a count table through a global alignment search using [USEARCH](https://www.drive5.com/usearch/).

```{code-block} bash
usearch -otutab combined.fasta -zotus asv.fasta -otutabout zotutable.txt
```

## 4. Taxonomy assignment

[IDTAXA](http://www2.decipher.codes/Classification.html) was used as the taxonomic classifier in this project. First, a local reference database was generated using [CRABS](https://github.com/gjeunen/reference_database_creator).

### 4.1 CRABS - local reference database

#### 4.1.1 db_download

First, download sequencing data from online repositories, including the mitofish database and NCBI.

```{code-block} bash
crabs db_download -s mitofish -o mitofish.fasta
crabs db_download -s ncbi -db nucleotide -q '16S[All Fields] AND (animals[filter] AND mitochondrion[filter])' -o ncbi16S.fasta -e gjeunen@gmail.com
```

Additionally, the NCBI taxonomy information needs to be downloaded to complete the reference database creation with CRABS.

```{code-block} bash
crabs db_download -s taxonomy
```

#### 4.1.2 db_merge

Second, merge the sequencing data from both online repositories that were downloaded previously.

```{code-block} bash
crabs db_merge -o merged.fasta -u yes -i ncbi16S.fasta mitofish.fasta
```

#### 4.1.3 insilico_pcr

Third, extract amplicons from reference sequences using an *in silico* PCR analysis.

```{code-block} bash
crabs insilico_pcr -i merged.fasta -o insilico.fasta -f GACCCTATGGAGCTTTAGAC -r CGCTGTTATCCCTADRGTAACT
```

#### 4.1.4 pga

Fourth, recover reference sequences without primer-binding regions through a pairwise global alignment analysis.

```{code-block} bash
crabs pga -i merged.fasta -o pga.fasta -db insilico.fasta -f GACCCTATGGAGCTTTAGAC -r CGCTGTTATCCCTADRGTAACT
```

#### 4.1.5 assign_tax

Fifth, assign a taxonomic lineage to each sequence.

```{code-block} bash
crabs assign_tax -i pga.fasta -o taxdb.tsv -a nucl_gb.accession2taxid -t nodes.dmp -n names.dmp -w yes
```

#### 4.1.6 dereplicate

Sixth, to reduce file size and remove reduntant sequences, dereplicate the data.

```{code-block} bash
crabs dereplicate -i taxdb.tsv -o derep.tsv -m uniq_species
```

#### 4.1.7 seq_cleanup

Seventh, prior to exporting the reference database, use various filtering parameters to retain only high quality references in the local database.

```{code-block} bash
crabs seq_cleanup -i derep.tsv -o cleandb.tsv -e yes -s yes -n 0
```

#### 4.1.8 tax_format

Eigth, export the database to SINTAX format, as CRABS currently doesn't support the IDTAXA format.

```{code-block} bash
crabs tax_format -i cleandb.tsv -o sintax.fasta -f sintax
```

Once the database is exported to SINTAX format, use the python script below to reformat the local reference database to IDTAXA specifications.
