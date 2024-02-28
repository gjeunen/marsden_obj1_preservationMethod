# marsden_obj1_preservationMethod

## 1. Introduction

This document serves as a comprehensive guide and supplementary information for the bioinformatic and statistical analysis of the manuscript entitled "Unlocking Antarctic molecular time-capsules - recovering historical environmental DNA from museum-preserved sponges" by Jeunen *et al*., 2024. This work is associated with the Marsden Fast-Start fund (MFP-UOO002116).

## 2. Environment set up

### 2.1 Starting files

The bioinformatic and statistical analysis requires access to several files, including the sequencing data and metadata files.

The two raw g-zipped sequencing files can be downloaded from the Sequence Read Archive (SRA) project number PRJNA1019816 and submission number SUB13856226. The files names are: "**8373-P1-00-01_S1_L001_R1_001.fastq.gz**" and "**8373-P2-00-01_S1_L001_R1_001.fastq.gz**". Download both sequencing data files and place them in the starting folder.

Besides the sequencing data, three metadatafiles are necessary to conduct the bioinformatic and statistical analysis, both of which can be downloaded from this GitHub repository in the folder "**metadata_files**". The first metadata file named "**barcodeMetadata8373.fasta**" will be used for demultiplexing, while the second metadata file named "**sampleMetadata8373.csv**" will be used for the statistical analysis. Place both metadata files into the starting folder. The third metadata file named "**Porifera_NIWA_22Sept2023.xlsx**" contains the information of the NIWA Invertebrate Collection for Figure 1.

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

## 5. Data pre-processing

Prior to the statistical analysis, data tables will be filtered and parsed according to various specifications. All code within this section should be run within an R environment.

### 5.1 Read data into R

First, all data files need to be read into R, including the sample metadata file "**sampleMetadata8373.csv**", the count table "**zotutable.txt**", the taxonomy table "**idtaxaLineage.txt**", and the ASV sequence file "**asv.fasta**".

```{code-block} R
#########################
# PREPARE R ENVIRONMENT #
#########################
library(readxl)
library(Biostrings)
library(dplyr)

# set working directory
setwd('/path/to/starting/folder/demux/qual_fasta/derep/')

##################
# READ DATA IN R #
##################
metaData <- read_csv('../../../sampleMetadata8373.csv')
freqTable <- read.table('zotutable.txt', header = TRUE, sep = '\t', row.names = 1, check.names = FALSE, comment.char = '')
taxonomyTable <- read.table('idtaxaLineage.txt', header = TRUE, sep = '\t', row.names = 1, check.names = FALSE, comment.char = '')
sequenceTable <- readDNAStringSet('asv.fasta')
```

### 5.2 Negative controls

To filter the data based on our negative controls, (1) sum all negative samples to column `freqTable$NEGSUM`, (2) identify which ASVs had a positive detection in our negative controls, (3) compare the abundance of these ASVs to samples, (4) remove all detections made up of 3 reads or less across the full count table, (5) set an abundance threshold of 10x the number of reads in the negative controls for ZOTUs with a positive detection in the negative controls. For example, ZOTU2 shows 3 reads in `freqTable$NEGSUM`. Hence, remove all detections of ZOTU2 in samples when read count is lower than 30.

```{code-block} R
#####################
# NEGATIVE CONTROLS #
#####################

# 1. find columns that contain NEG in header, 2. sum rows that contain NEG in header to column NEGSUM, 3. drop individual NEG columns
negColumns <- grep('NEG', names(freqTable), value = TRUE)
freqTable$NEGSUM <- rowSums(freqTable[negColumns])
freqTable <- freqTable[, !names(freqTable) %in% negColumns]

# 1. create list of index names (ZOTU numbers) for which freqTable$NEGSUM > 0, 2. print taxonomy of ZOTUs in list
negZOTUs <- rownames(freqTable)[which(freqTable$NEGSUM > 0)]
print(taxonomyTable[negZOTUs, ])

# 1. iterate over index names in negZOTUs list and print out the necessary values
for (x in negZOTUs) {
  if (x %in% rownames(freqTable)) {
    negValue <- freqTable[x, 'NEGSUM']
    negValuePerc <- negValue / rowSums(freqTable[x, ]) * 100
    rowMeanValue <- rowMeans(freqTable[x, ])
    positiveDetections <- sum(freqTable[x, ] > 0)
    positiveDetectionsPerc <- sum(freqTable[x, ] > 0) / nrow(freqTable) * 100
    cat(x, 'value in NEG:', negValue, '% of ZOTU reads:', negValuePerc, 'mean # reads in sample:', rowMeanValue, '# +ve detections:', positiveDetections, '% +ve detections:', positiveDetectionsPerc, '\n')
  }
}

# 1. set detections of 3 or lower to 0
freqTable.negFilter <- as.data.frame(apply(freqTable, 2, function(x) ifelse(x < 4, 0, x)))

# 1. set abundance threshold 30 for ZOTU2, 2. set abundance threshold 10 for ZOTU8
freqTable.negFilter['Zotu2', -1] <- ifelse(freqTable.negFilter['Zotu2', -1] >= 30, freqTable.negFilter['Zotu2', -1], 0)
freqTable.negFilter['Zotu8', -1] <- ifelse(freqTable.negFilter['Zotu8', -1] >= 10, freqTable.negFilter['Zotu8', -1], 0)
```

### 5.3 Low confidence tax IDs

After filtering the data based on negative control samples, remove low-confidence taxonomy assignments from the data. Low-confidence taxonomy assignments are defined as sequences for which IDTAXA could not provide a taxonomic ID at the order or higher level. Additionally, remove temperate fish species from the data.

```{code-block} R
#########################
# LOW CONFIDENCE TAX ID #
#########################

# 1. find ZOTUs which do not have a taxonomic ID at order level, 2. remove this list from the frequency table
lowConfidenceTax <- rownames(taxonomyTable)[which(is.na(taxonomyTable$order))]
freqTable.lowConfidenceTax <- freqTable.negFilter[!(rownames(freqTable.negFilter) %in% lowConfidenceTax), ]

# 1. create ZOTU list that matched to temperate fish species, 2. remove this list from the frequency table
temperateZOTUs <- c('Zotu56', 'Zotu71', 'Zotu88')
freqTable.temperateFish <- freqTable.lowConfidenceTax[!(rownames(freqTable.lowConfidenceTax) %in% temperateZOTUs), ]
```

### 5.4 Artefact sequences

Remove artefact sequences using a taxon-dependent co-occurrence pattern of similar sequences.

```{code-block} R
##################
# ARTEFACT ZOTUS #
##################

# 1. create ZOTU list for additional removal, 2. remove this list from the frequency table
removeExtraZOTUs <- c('Zotu62', 'Zotu152')
freqTable.artefact <- freqTable.temperateFish[!(rownames(freqTable.temperateFish) %in% removeExtraZOTUs), ]

# 1. merge daughters with parent ZOTUs, as identified by visual inspection of the frequency table, 2. remove daughter ZOTUs
# ZOTU32 is parent of ZOTU41, ZOTU45, and ZOTU49, all belonging to the genus bathylagus
freqTable.artefact['Zotu32', ] <- freqTable.artefact['Zotu32', ] + freqTable.artefact['Zotu41', ] + freqTable.artefact['Zotu45', ] + freqTable.artefact['Zotu49', ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu41'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu45'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu49'), ]

# ZOTU15 is parent of ZOTU17, ZOTU93, and ZOTU146, all belonging to Notolepis coatsorum
freqTable.artefact['Zotu15', ] <- freqTable.artefact['Zotu15', ] + freqTable.artefact['Zotu17', ] + freqTable.artefact['Zotu93', ] + freqTable.artefact['Zotu146', ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu17'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu93'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu146'), ]

# ZOTU2 is parent of ZOTU61, ZOTU72, ZOTU86, ZOTU100, ZOTU104, ZOTU127, ZOT133, and ZOTU142, all belonging to the genus Macrourus
freqTable.artefact['Zotu2', ] <- freqTable.artefact['Zotu2', ] + freqTable.artefact['Zotu61', ] + freqTable.artefact['Zotu72', ] + freqTable.artefact['Zotu86', ] + freqTable.artefact['Zotu100', ] + freqTable.artefact['Zotu104', ] + freqTable.artefact['Zotu127', ] + freqTable.artefact['Zotu133', ] + freqTable.artefact['Zotu142', ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu61'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu72'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu86'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu100'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu104'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu127'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu133'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu142'), ]

# ZOTU59 is the parent of ZOTU78, belonging to Antimora rostrata
freqTable.artefact['Zotu59', ] <- freqTable.artefact['Zotu59', ] + freqTable.artefact['Zotu78', ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu78'), ]

# ZOTU128 is the parent of ZOTU149, belonging to Gymnoscopelus braueri
freqTable.artefact['Zotu128', ] <- freqTable.artefact['Zotu128', ] + freqTable.artefact['Zotu149', ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu149'), ]

# ZOTU35 is the parent of ZOTU80, belonging to the genus Bathydraco
freqTable.artefact['Zotu35', ] <- freqTable.artefact['Zotu35', ] + freqTable.artefact['Zotu80', ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu80'), ]

# ZOTU58 is the parent of ZOTU64 and ZOTU65, belonging to Gerlachea australis
freqTable.artefact['Zotu58', ] <- freqTable.artefact['Zotu58', ] + freqTable.artefact['Zotu64', ] + freqTable.artefact['Zotu65', ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu64'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu65'), ]

# ZOTU40 is the parent of ZOTU43, belonging to Psilodraco breviceps
freqTable.artefact['Zotu40', ] <- freqTable.artefact['Zotu40', ] + freqTable.artefact['Zotu43', ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu43'), ]

# ZOTU95 is the parent of ZOTU124 and ZOTU151, belonging to the genus Chionobathyscus
freqTable.artefact['Zotu95', ] <- freqTable.artefact['Zotu95', ] + freqTable.artefact['Zotu124', ] + freqTable.artefact['Zotu151', ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu124'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu151'), ]

# ZOTU19 is the parent of ZOTU34, belonging to Chionobathyscus dewitti
freqTable.artefact['Zotu19', ] <- freqTable.artefact['Zotu19', ] + freqTable.artefact['Zotu34', ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu34'), ]

# ZOTU16 is parent of ZOTU25, ZOTU120, ZOTU131, and ZOTU135, all belonging to the genus Chionodraco
freqTable.artefact['Zotu16', ] <- freqTable.artefact['Zotu16', ] + freqTable.artefact['Zotu25', ] + freqTable.artefact['Zotu120', ] + freqTable.artefact['Zotu131', ] + freqTable.artefact['Zotu135', ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu25'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu120'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu131'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu135'), ]

# ZOTU20 is parent of ZOTU39, ZOTU129, and ZOTU130, all belonging to the genus Chionodraco
freqTable.artefact['Zotu20', ] <- freqTable.artefact['Zotu20', ] + freqTable.artefact['Zotu39', ] + freqTable.artefact['Zotu129', ] + freqTable.artefact['Zotu130', ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu39'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu129'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu130'), ]

# ZOTU13 is the parent of ZOTU22, belonging to Pagetopsis maculatus
freqTable.artefact['Zotu13', ] <- freqTable.artefact['Zotu13', ] + freqTable.artefact['Zotu22', ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu22'), ]

# ZOTU26 is the parent of ZOTU85, belonging to the genus Paraliparis
freqTable.artefact['Zotu26', ] <- freqTable.artefact['Zotu26', ] + freqTable.artefact['Zotu85', ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu85'), ]

# ZOTU1 is the parent of ZOTU109, ZOTU110, ZOTU136, ZOTU138, ZOTU140, and ZOTU144, all belonging to Dissostichus mawsoni
freqTable.artefact['Zotu1', ] <- freqTable.artefact['Zotu1', ] + freqTable.artefact['Zotu109', ] + freqTable.artefact['Zotu110', ] + freqTable.artefact['Zotu136', ] + freqTable.artefact['Zotu138', ] + freqTable.artefact['Zotu140', ] + freqTable.artefact['Zotu144', ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu109'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu110'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu136'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu138'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu140'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu144'), ]

# ZOTU3 is the parent of ZOTU8 and ZOTU103, belonging to Pleuragramma antarctica
freqTable.artefact['Zotu3', ] <- freqTable.artefact['Zotu3', ] + freqTable.artefact['Zotu8', ] + freqTable.artefact['Zotu103', ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu8'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu103'), ]

# ZOTU73 is the parent of ZOTU97 and ZOTU125, belonging to the genus Pagrus
freqTable.artefact['Zotu73', ] <- freqTable.artefact['Zotu73', ] + freqTable.artefact['Zotu97', ] + freqTable.artefact['Zotu125', ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu97'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu125'), ]

# ZOTU44 is the parent of ZOTU52, ZOTU63, ZOTU89, and ZOTU92, all belonging to Cyclothone pygmaea
freqTable.artefact['Zotu44', ] <- freqTable.artefact['Zotu44', ] + freqTable.artefact['Zotu52', ] + freqTable.artefact['Zotu63', ] + freqTable.artefact['Zotu89', ] + freqTable.artefact['Zotu92', ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu52'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu63'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu89'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu92'), ]

# ZOTU37 is the parent of ZOTU46 and ZOTU54, belonging to the genus Hoplostethus
freqTable.artefact['Zotu37', ] <- freqTable.artefact['Zotu37', ] + freqTable.artefact['Zotu46', ] + freqTable.artefact['Zotu54', ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu46'), ]
freqTable.artefact <- freqTable.artefact[-which(rownames(freqTable.artefact) == 'Zotu54'), ]
```

### 5.5 Sample read count

The final filtering step of the data is to remove low-abundant samples from the analysis. For this experiment, the threshold is set at 10,000 reads. Since most samples from *Cinachyra* sponges failed to amplify and contain more than 10,000 reads, remove all samples from *Cinachyra* sponges.

```{code-block} R
#####################
# SAMPLE READ COUNT #
#####################

# 1. list all samples with reads lower than threshold, 2. iterate over list and print metadata for samples, as well as read count, 3. drop columns from rawTable, 4. drop rows from metaData
lowCoverageSamples <- names(which(colSums(freqTable.artefact) < 10000))
for (x in lowCoverageSamples) {
  if (x %in% metaData$sampleID) {
    row <- paste(metaData[metaData$sampleID == x, ], collapse = ' ')
    seqCount <- sum(freqTable.artefact$x)
    cat(row, seqCount, '/n')
  } else {
    seqCount <- sum(freqTable.artefact$x)
    cat('sampleID', x, 'not found in metaData, seqCount = ', seqCount, '\n')
  }
}

# 1. list samples to remove based on abundance threshold and sponge species, 2. remove samples from the frequency table.
cinachyraSamples<- metaData %>%
  filter(spongeID == 'Cinachyra barbata') %>%
  pull(sampleID)
samplesToRemove <- unique(c(cinachyraSamples, lowCoverageSamples))
freqTable.sampleThreshold <- freqTable.artefact[, !(names(freqTable.artefact) %in% samplesToRemove)]

# 1. check if any ZOTUs have a total count of 0, 2. remove these ZOTUs from the frequency table
removeZOTU <- names(which(rowSums(freqTable.sampleThreshold) == 0))
removeZOTU
```

### 5.6 Export data tables

Once the data is filtered, export and write the updated tables to new files.

```{code-block} R
#################
# UPDATE TABLES #
#################

# 1. update metaData based on headers in freqTable.sampleThreshold
metaData.sampleThreshold <- metaData[(metaData$sampleID) %in% rownames(t(freqTable.sampleThreshold)), ]

# 1. update taxonomyTable based on rownames freqTable.sampleThreshold
taxonomyTable.sampleThreshold <- taxonomyTable[(rownames(taxonomyTable) %in% rownames(freqTable.sampleThreshold)), ]

# 1. update sequenceTable based on rownames taxonomyTable.sampleThreshold
sequenceTable.sampleThreshold <- sequenceTable[rownames(taxonomyTable.sampleThreshold)]

# 1. write freqTable.sampleThreshold to output, 2. write metaData.sampleThreshold to output, 3. write taxonomyTable.sampleThreshold to output, 4. write sequenceTable.sampleThreshold to output
write.table(freqTable.sampleThreshold, 'zotutableFiltered.txt', append = FALSE, sep = '\t', dec = '.', row.names = TRUE, col.names = NA)
write.table(metaData.sampleThreshold, '../../../sampleMetadata8373Filtered.txt', append = FALSE, sep = '\t', dec = '.', row.names = FALSE)
write.table(taxonomyTable.sampleThreshold, 'taxonomyFiltered.txt', append = FALSE, sep = '\t', dec = '.', row.names = TRUE, col.names = NA)
writeXStringSet(sequenceTable.sampleThreshold, file = 'asvFiltered.fasta')
```

## 6. Statistical analysis

All code for the statistical analysis is run in an R environment.

### 6.1 Reading data into R

To complete the statistical analysis, 6 files need to be read into R, all of which were created during the bioinformatic processing or are made available on this GitHub repository. The 6 files include: (1) the filtered count table "**zotutableFiltered.txt**", (2) the updated sample metadata file "**sampleMetadata8373Filtered.txt**", (3) the filtered ASV sequence file "**asvFiltered.fasta**", (4) the phylogenetic tree "**asvFiltered-tree.tree**", (5) the updated taxonomy table "**taxonomyFiltered.txt**", and (6) the NIWA Invertebrate Collection database "**Porifera_NIWA_22Sept2023.xlsx**".

```{code-block} R
#########################
# PREPARE R ENVIRONMENT #
#########################
library(Biostrings)
library(phyloseq)
library(ggplot2)
library(vegan)
library(iNEXT)
library(microbiome)
library(ape)
library(scales)
library(readxl)
library(tidyverse)
library(ampvis2)
library(car)
library(FSA)
library(dplyr)
library(agricolae)
library(speedyseq)
library(BiodiversityR)
library(indicspecies)

# set working directory
setwd('/path/to/starting/folder/demux/qual_fasta/derep/')

##################
# READ DATA IN R #
##################
metaData <- read.table('../../../sampleMetadata8373Filtered.txt', header = TRUE, sep = '\t', row.names = 1, check.names = FALSE, comment.char = '')
freqTable <- read.table('zotutableFiltered.txt', header = TRUE, sep = '\t', row.names = 1, check.names = FALSE, comment.char = '')
taxonomyTable <- read.table('taxonomyFiltered.txt', header = TRUE, sep = '\t', row.names = 1, check.names = FALSE, comment.char = '')
sequenceTable <- readDNAStringSet('asvFiltered.fasta')
phyloTree <- read.tree('asvFiltered-tree.tree')
metadata <- read_xlsx('../../../Porifera_NIWA_22Sept2023.xlsx', sheet = 'Porifera', col_names = TRUE, na = '')
```

### 6.2 Figure 1

Figure 1 of the manuscript provides information on the number of sponge specimens in the NIWA Invertebrate Collection (NIC) per decade and facetted by preservation method, including dry, ethanol, frozen, formalin, isopropanol, and other. Specimens included in ‘other’ include preservation methods listed as Alcohol Unknown, Ethanol – Previously Unknown, and Slide. Number above bars represent number of specimens. Y-axis reported as square root transformed to increase readability of low-abundant collection numbers. For NIC specimen data, see [Porifera_NIWA_22Sept2023.xlsx](https://nzobisipt.niwa.co.nz/resource?r=obisspecify).

```{code-block} R
#######################################
# FIGURE 1: PRESERVATION OF SPECIMENS #
#######################################
# remove entries that do not have year information
metadata_year <- metadata[!is.na(metadata$`Date (Year)`), ]

# determine groups of preservation types
unique(metadata_year[c("Pres Type")])

# change all preservation type not in list to other
shouldBecomeOther <- !(metadata_year$`Pres Type` %in% c('Frozen', 'Formalin', 'Ethanol - no formalin', 'Dry', 'Isopropanol - orig unknown', 'Ethanol - orig formalin', 'Ethanol - orig Isopropanol'))
metadata_year$`Pres Type`[shouldBecomeOther] <- 'other'

# change 'Isopropanol - orig unknown' and 'Ethanol - orig Isopropanol' to Isopropanol
shouldBecomeIsopropanol <- (metadata_year$`Pres Type` %in% c('Isopropanol - orig unknown', 'Ethanol - orig Isopropanol'))
metadata_year$`Pres Type`[shouldBecomeIsopropanol] <- 'Isopropanol'

# change 'Ethanol - orig formalin' to Formalin
shouldBecomeFormalin <- (metadata_year$`Pres Type` %in% c('Ethanol - orig formalin'))
metadata_year$`Pres Type`[shouldBecomeFormalin] <- 'Formalin'

# change 'Ethanol - no formalin' to Ethanol
shouldBecomeEthanol <- (metadata_year$`Pres Type` %in% c('Ethanol - no formalin'))
metadata_year$`Pres Type`[shouldBecomeEthanol] <- 'Ethanol'

# see if data is formatted appropriately
unique(metadata_year[c("Pres Type")])

# only keep the necessary columns
metadata_year_subset <- subset(metadata_year, select = c('Pres Type', 'Date (Year)'))

# find minimum and maximum year to group samples
max(metadata_year_subset$`Date (Year)`)
min(metadata_year_subset$`Date (Year)`)

# Create 10-year bins
metadata_year_subset <- metadata_year_subset %>%
  mutate(`Year Bin` = cut(`Date (Year)`, breaks = seq(1910, 2030, by = 10), right = FALSE))

# Create all possible combinations of 'Year Bin' and 'Pres Type'
all_combinations <- expand.grid(
  `Year Bin` = unique(metadata_year_subset$`Year Bin`),
  `Pres Type` = unique(metadata_year_subset$`Pres Type`)
)

# Calculate the number of items in each 10-year bin for each 'Pres Type'
bin_counts <- metadata_year_subset %>%
  group_by(`Year Bin`, `Pres Type`) %>%
  summarize(count = n(), .groups = "drop")

# Merge with all combinations to fill in missing counts with zeros
bin_counts <- merge(all_combinations, bin_counts, by = c("Year Bin", "Pres Type"), all.x = TRUE)
bin_counts[is.na(bin_counts$count), "count"] <- 0

# set order of plots
facet_order <- c("Dry", "Ethanol", "Frozen", "Formalin", "Isopropanol", "other")

# Create separate bar plots for each 'Pres Type'
ggplot(bin_counts, aes(x = `Year Bin`, y = count, fill = `Pres Type`)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count), vjust = -0.5, size = 3) +
  labs(title = "Number of Items in 10-Year Bins by Press Type",
       x = "Year Bin",
       y = "Number of Items") +
  #theme_classic() +
  theme_grey() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  scale_y_sqrt(expand = c(0, 0), limits = c(0, 6000)) +
  facet_wrap(~ factor(`Pres Type`, levels = facet_order), scales = "free_x")
```

### 6.3 Figure 2

Figure 2 of the manuscript displays the map of the Ross Sea, Antarctica and depicts specimen collection locations. Points are coloured by preservation method, including dry (yellow), ethanol (blue), and frozen (red). Point shape is dictated by sponge ID, including Cinachyra sp. (inverted triangle), Homaxinella sp. (circle), Inflatella belli (square), Rossella nuda (diamond), and Rossella villosa (triangle).

```{code-block} R
#read in data
if (!require("pacman")) install.packages("pacman")
pacman::p_load(terra, tidyterra, ggplot2, ggnewscale, patchwork)
ATAshp <- simplifyGeom(vect("Coastline_high_res_polygon/Coastline_high_res_polygon.shp"), tolerance=500)
ATAshp$col <- sapply(ATAshp$surface, function(x){ifelse(x == "land", "grey85", "grey98")})
sponge_dat <- read.csv("../../../sampleMetadata8373.csv")
sponge_unique <- sponge_dat[!duplicated(sponge_dat$sampleSet), ]
sponge_unique <- sponge_unique[!is.na(sponge_unique$preservationType), ]
sponge_unique <- sponge_unique[sponge_unique$spongeID != "Cinachyra barbata", ]
sponge_pts <- project(vect(sponge_unique, geom = c("longitude", "latitude"), keepgeom = TRUE, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"), ATAshp)
grat <- project(vect(sf::st_graticule(lon=seq(-175, 180, 5), lat = seq(-85, -60, 5), ndiscr = 5000)), ATAshp)
ibsco <- rast("IBCSO_v2_ice-surface.nc")
sample_colors <- c("dry" = "lightgoldenrod", "ethanol" = "steelblue", "frozen" = "firebrick")
spp_names <- c(expression(italic("Homaxinella")), expression(italic("Inflatella belli")), expression(italic("Rossella nuda")), expression(italic("Rossella villosa")))

#define xy limits in m (projected coordinate space) for plot
xmn <- -1000000
xmx <- 1500000
ymn <- -2500000
ymx <- -1000000

#plot Ross Sea region
RSR_plot <- ggplot() +
 geom_spatraster(data = ibsco, show.legend = FALSE) +
 geom_spatvector(data = ATAshp, alpha = 1, fill = ATAshp$col, col = "grey20") +
 scale_fill_gradient(low = "#528B8B", high = "#CBDCDC") +
 geom_spatvector(data = grat, col = "grey50", alpha = 0.5) +
 new_scale_fill() +
 geom_spatvector(data = sponge_pts, aes(shape = spongeID, fill = preservationType), size = 3, stroke = 0.5) +
 scale_fill_manual(values = sample_colors, labels = c("Dry", "Ethanol", "Frozen") , name = "Preservation type") +
 guides(fill = guide_legend(override.aes = list(pch=21))) +
 scale_shape_manual(values = c(21, 22, 23, 24, 1), labels = spp_names, name = "Species") +
 scale_x_continuous(breaks = seq(-180, 180, by = 10)) +
 theme(panel.background = element_blank(), panel.border = element_blank(), legend.position = c(0.12, 0.37), legend.text = element_text(hjust = 0), plot.margin = unit(c(45, 5, 5, 5), "pt")) +
 coord_sf(crs = crs(ATAshp), expand = FALSE, xlim = c(xmn, xmx), ylim = c(ymn, ymx)) +
 annotate(geom = "rect", xmin = xmn, xmax = xmx, ymin = ymn, ymax = ymx, fill = NA, col = "grey10")

#plot inset map Antarctica
inset_map <- ggplot() +
 geom_rect(aes(xmin = xmn, xmax = xmx, ymin = ymn, ymax = ymx), fill = "white", col = NA, alpha = 1) +
 geom_spatvector(data = ATAshp, alpha = 1, fill = ATAshp$col, col = "black") +
 geom_rect(aes(xmin = xmn, xmax = xmx, ymin = ymn, ymax = ymx), fill = "dark red", col = "black", alpha = 0.2) +
 theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_blank(), plot.background = element_blank(), panel.border = element_blank(), axis.line = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_blank())
 
#output to pdf
pdf("RSR_spongemap_v1.pdf", width = 8, height = 5.5)
RSR_plot + annotation_custom(ggplotGrob(inset_map), xmin = xmx - (xmx-xmn)/3, xmax = xmx, ymin = ymx - (ymx-ymn)/3, ymax = ymx + 280000)
dev.off()
```

### 6.4 Figure 3

Figure 3 of the manuscript displays a bayesian phylogenetic tree generated of all 64 ZOTU sequences. Bayesian tree generated in Beast2. Tip labels represent ZOTU number. Taxonomic ID for each ZOTU can be retrieved from SUPPLEMENT 4. Inner bar graph showing the number of detections of each ZOTU sequence within the nine specimens stored dry (yellow), in ethanol (blue), and frozen (red). Outer bar graph showing the relative read abundance of each ZOTU sequence within the nine specimens stored dry (yellow), in ethanol (blue), and frozen (red). Axis for relative read abundance bar graph is reported as square root transformed to increase readability of low-abundant signals. Most frequently and abundant taxonomic groups are represented by silhouettes, including (a) Chondrichthyes, (b) Gadiformes, (c) Bathylagidae, (d) Nototheniidae, (e) Bathydraconidae, and (f) Channichthyidae.

```{code-block} R
```
