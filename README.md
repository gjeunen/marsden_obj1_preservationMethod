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

Once the database is exported to SINTAX format, use the python script below to reformat the local reference database to IDTAXA specifications. Note that IDTAXA cannot resolve homonyms, i.e., identical taxon names. In the resulting database "**sintax.fasta**", both a gastropod and fish family are named "Chilodontidae". To ensure IDTAXA runs without errors, the python code below uses names of two ranks pasted together to generate a unique name for each taxon. For example, the fish family name of Chilodontidae will be transformed to the order name + the family name: CharaciformesChilodontidae. This output of this "hack" will be corrected at a later stage to revert taxon names back to their appropriate format.

```{code-block} bash
nano idtaxaInitFormat.py
```

```{code-block} python
#! /usr/bin/env python3

import sys

inFile = sys.argv[1]
dbOut = sys.argv[2]
taxOut = sys.argv[3]

indexNameDict = {'Root' : 0}
parentDict = {'Root' : -1}
levelDict = {'Root' : 0}
rankDict = {'Root' : 'rootrank'}
indexNumber = 0
ranks = {0: 'rootrank', 1 : 'domain', 2 : 'phylum', 3 : 'class', 4 : 'order', 5: 'family', 6 : 'genus', 7 : 'species'}

idTaxaSeqDict = {}

with open(inFile, 'r') as infile:
    for line in infile:
        previousName = 'Root'
        previousTaxonName = 'Root'
        oldName = 'start'
        levelTax = 0
        if line.startswith('>'):
            line = line.rstrip('\n')
            seqID = line.split(';')[0]
            headerStart = f'{seqID};Root'
            lineage = line.split(':')[1:8]
            for taxon in lineage:
                levelTax += 1
                taxonName = taxon.split(',')[0]
                if taxonName == '':
                    taxonName = f'{previousName}_insertaeCedis'
                uniqueTaxonName = f'{previousTaxonName}{taxonName}'
                headerStart = f'{headerStart};{uniqueTaxonName}'
                if uniqueTaxonName not in indexNameDict:
                    indexNumber += 1
                    indexNameDict[uniqueTaxonName] = indexNumber
                    parentDict[uniqueTaxonName] = indexNameDict[previousName]
                    levelDict[uniqueTaxonName] = levelTax
                    rankDict[uniqueTaxonName] = ranks[levelTax]
                taxonNumber = indexNameDict[uniqueTaxonName]
                oldName = previousName
                previousName = uniqueTaxonName
                previousTaxonName = taxonName
        else:
            idTaxaSeqDict[headerStart] = line.rstrip('\n')                

with open(taxOut, 'w') as outfile:
    for indexName in indexNameDict:
        printLine = f'{indexNameDict[indexName]}*{indexName}*{parentDict[indexName]}*{levelDict[indexName]}*{rankDict[indexName]}\n'
        _ = outfile.write(printLine)

with open(dbOut, 'w') as seqFile:
    for key in idTaxaSeqDict:
        _ = seqFile.write(f'{key}\n{idTaxaSeqDict[key]}\n')
```

Press `ctrl + x` to exit out of the editor, followed by `y` and `return`.

```{code-block} bash
chmod +x idtaxaInitFormat.py
```

```{code-block} bash
./idtaxaInitFormat.py sintax.fasta idtaxa.fasta idtaxa.txt
```

### 4.2 IDTAXA classifier training

Once the reference database has been formatted to IDTAXA specifications, train the IDTAXA classifier on the training set, i.e., the formatted reference database. IDTAXA is incorporated in the R package [DECIPHER](http://www2.decipher.codes). Hence, the following code is run in R.

```{code-block} R
library(DECIPHER)
setwd('/path/to/starting/folder/demux/qual_fasta/derep/')

seqIDTAXA <- readDNAStringSet('idtaxa.fasta')
rankIDTAXA <- read.table('idtaxa.txt', header = FALSE, col.names = c('Index', 'Name', 'Parent', 'Level', 'Rank'), sep = '*', quote = '', stringsAsFactors = FALSE)

groups <- names(seqIDTAXA)
head(groups)
groups <- gsub("(.*)(Root;)", "\\2", groups)
head(groups)
groupCounts <- table(groups)
uniqueGroups <- names(groupCounts)
length(uniqueGroups)

# count number of sequences per group and, optionally, select only a subset of sequences if the group is deemed too large (10 seqs)
maxGroupSize <- 10 # max sequences per label (>= 1)
remove <- logical(length(seqIDTAXA))
for (i in which(groupCounts > maxGroupSize)) {
  index <- which(groups==uniqueGroups[i])
  keep <- sample(length(index),
                 maxGroupSize)
  remove[index[-keep]] <- TRUE
}
sum(remove) # number of sequences eliminated

# train the classifier. Set 'maxIterations' to number of times the classifier will be trained.
maxIterations <- 5 # must be >= 1
allowGroupRemoval <- FALSE
probSeqsPrev <- integer() # suspected problem sequences from prior iteration
for (i in seq_len(maxIterations)) {
  cat("Training iteration: ", i, "\n", sep="")
  # train the classifier
  trainingSet <- LearnTaxa(seqIDTAXA[!remove],
                           names(seqIDTAXA)[!remove],
                           rankIDTAXA)
  # look for problem sequences
  probSeqs <- trainingSet$problemSequences$Index
  if (length(probSeqs)==0) {
    cat("No problem sequences remaining.\n")
    break
  } else if (length(probSeqs)==length(probSeqsPrev) &&
             all(probSeqsPrev==probSeqs)) {
    cat("Iterations converged.\n")
    break
  }
  if (i==maxIterations)
    break
  probSeqsPrev <- probSeqs
  # remove any problem sequences
  index <- which(!remove)[probSeqs]
  remove[index] <- TRUE # remove all problem sequences
  if (!allowGroupRemoval) {
    # replace any removed groups
    missing <- !(uniqueGroups %in% groups[!remove])
    missing <- uniqueGroups[missing]
    if (length(missing) > 0) {
      index <- index[groups[index] %in% missing]
      remove[index] <- FALSE # don't remove
    }
  }
}
sum(remove)
length(probSeqs)
trainingSet
plot(trainingSet)
```

### 4.3 Assign taxonomic ID to ASV

Once the classifier is trained, use the training set to assign a taxonomic ID to each ASV using the IDTAXA algorithm. Again, the code below should be run in an R environment.

```{code-block} R
library(DECIPHER)
setwd('/path/to/starting/folder/demux/qual_fasta/derep/')

# read in the ZOTU fasta file
fastaIDTAXA <- readDNAStringSet('asv.fasta')

# classify sequences
set.seed(123)
ids <- IdTaxa(fastaIDTAXA, trainingSet, type = 'extended', strand = 'top', threshold = 60, processors = NULL)
ids

# plot results
plot(ids, trainingSet)

# export output
output <- sapply(ids,
                 function(id) {
                   paste(id$taxon,
                         " (",
                         round(id$confidence, digits = 1),
                         "%)",
                         sep = "",
                         collapse = ": ")
                 })
# Create an empty character vector to store the output
outputFile <- character(0)
# Combine the zotu ID and information and add them to the output vector
for (zotu_id in names(output)) {
  entry <- paste(zotu_id, output[zotu_id], sep = "\t")  # Use tab as a separator
  outputFile <- c(outputFile, entry)
}
# Write output to file
writeLines(outputFile, 'taxonomyIDTAXA.txt')
```

### 4.4 Format taxonomic ID output file

To easily read the IDTAXA results into R for downstream processing, parse the document using the following python script.

```{code-block} bash
nano idtaxaOutputFormat.py
```

```{code-block} python
#! /usr/bin/env python3

import sys
import collections

def removeBeforeSecondUppercase(input_string):
    uppercase_count = 0
    result = ''
    count = 0
    for char in input_string:
        count += 1
        if char.isupper():
            uppercase_count += 1
            if uppercase_count == 2:
                result = input_string[count -1:]
                break
    return result

idTaxaInputFile = sys.argv[1]
idTaxaLineageOutputFile = sys.argv[2]
idTaxaScoreOutputFile = sys.argv[3]

headerRanks = ['root', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
idLineageDict = collections.defaultdict(dict)
idScoreDict = collections.defaultdict(dict)

with open(idTaxaInputFile, 'r') as inputFile:
  for line in inputFile:
    line = line.rstrip('\n')
    seqID = line.split('\t')[0]
    lineage = line.split('\t')[1]
    lineageItems = lineage.split(': ')
    for item in range(len(headerRanks)):
      try:
        if 'unclassified' in lineageItems[item]:
          idLineageDict[seqID][headerRanks[item]] = 'NA'
          idScoreDict[seqID][headerRanks[item]] = 'NA'
        else:
          lineageString = lineageItems[item].split(' ')[0]
          if lineageString != 'Root':
            lineageString = removeBeforeSecondUppercase(lineageString)
          scoreString = lineageItems[item].split('(')[1].rstrip(')')
          idScoreDict[seqID][headerRanks[item]] = scoreString
          idLineageDict[seqID][headerRanks[item]] = lineageString
      except IndexError:
        idLineageDict[seqID][headerRanks[item]] = 'NA'
        idScoreDict[seqID][headerRanks[item]] = 'NA'

with open(idTaxaLineageOutputFile, 'w') as outfile:
  _ = outfile.write('ID' + '\t' + '\t'.join(headerRanks) + '\n')
  for i in idLineageDict:
    _ = outfile.write(i)
    for j, g in idLineageDict[i].items():
      _ = outfile.write('\t' + g)
    _ = outfile.write('\n')

with open(idTaxaScoreOutputFile, 'w') as outfile:
  _ = outfile.write('ID' + '\t' + '\t'.join(headerRanks) + '\n')
  for i in idScoreDict:
    _ = outfile.write(i)
    for j, g in idScoreDict[i].items():
      _ = outfile.write('\t' + g)
    _ = outfile.write('\n')
```

Press `ctrl + x` to exit out of the editor, followed by `y` and `return`.

```{code-block} bash
chmod +x idtaxaOutputFormat.py
```

```{code-block} bash
./idtaxaOutputFormat.py taxonomyIDTAXA.txt idtaxaLineage.txt idtaxaScore.txt
```
