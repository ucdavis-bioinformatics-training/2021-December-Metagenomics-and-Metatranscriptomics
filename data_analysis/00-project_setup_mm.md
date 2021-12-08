
# The Dataset

The paper

Fuyong Li, et al. ["Comparative metagenomic and metatranscriptomic analyses reveal the breed effect on the rumen microbiome and its associations with feed efficiency in beef cattle."](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-019-0618-5#Sec2), Microbiome 7, 6 (2019)

The project on the NCBI Sequence Read Archive (SRA), [PRJNA448333](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA448333).

Relevant RNA sections of the paper

## Bovine models
Forty-eight steers were selected from a herd of 738 beef cattle that were born in 2014 and raised at the Roy Berg Kinsella Research Ranch, University of Alberta, according to their breeds and residual feed intake (RFI) ranking. These 48 steers belong to three breeds and two RFI groups (high RFI [H-RFI, inefficient] and low RFI [L-RFI, efficient]), including two purebreds (Angus [ANG]; H-RFI, n = 8; L-RFI, n = 8) and Charolais [CHAR]; H-RFI, n = 8; L-RFI, n = 8), and one crossbred (Kinsella composite hybrid [HYB]; H-RFI, n = 8; L-RFI, n = 8).

## Metagenomics and metatranscriptomics sequencing and data analysis
Total genomic DNA was isolated from rumen digesta using the repeated bead beating plus column (RBB + C) method. The quality and quantity of DNA was measured using a NanoDrop Spectrophotometer ND-1000 (Thermo Fisher Scientific Inc., Wilmington, DE, USA). Metagenome library was constructed using the TruSeq DNA PCR-Free Library Preparation Kit (Illumina, San Diego, CA, USA), and the quantity of each library was evaluated using a Qubit 2.0 fluorimeter (Invitrogen, Carlsbad, CA, USA). Sequencing of metagenome libraries was conducted at the McGill University and Génome Québec Innovation Centre (Montréal, QC, Canada) using Illumina HiSeq 2000 (100 bp paired-end sequencing of ~ 350 bp inserts).

Total RNA was extracted from rumen disgesta. The RNA yield was measured using a Qubit 2.0 fluorimeter (Invitrogen), and the RNA quality was measure using an Agilent 2200 TapeStation (Agilent Technologies, Santa Clara, CA, USA). Only samples with RNA integrity number (RIN) ≥ 7.0 were used to generate metatranscriptome libraries. In the current study, two types of metatranscriptome libraries were constructed: total RNA-based metatranscriptome libraries (T-metatranscriptome) and mRNA-enriched metatranscriptome libraries (M-metatranscriptome). For the M-metatranscriptome library construction, rRNA in each sample was depleted using the Ribo-Zero Gold rRNA Removal Kit (Epidemiology) (Illumina) according to the manufacturer’s instruction. Total RNA and enriched mRNA were used for T- and M-metatranscriptome library construction, respectively, using the TruSeq RNA Library Prep Kit v2 (Illumina). Sequencing of T- and M-metatranscriptome libraries was conducted at the McGill University and Génome Québec Innovation Centre (Montréal, QC, Canada) using Illumina HiSeq 2000 (100 bp paired-end sequencing of ~ 140 bp inserts) and 2500 (125 bp paired-end sequencing of ~ 140 bp inserts), respectively.

The quality control (QC) of each dataset was performed using Trimmomatic (version 0.35) to trim artificial sequences (adapters), cut low quality bases (quality scores < 20), and remove short reads (< 50 bp). The program SortMeRNA (version 1.9) was used to extract rDNA and rRNA reads from sequencing datasets. Non-rDNA/rRNA reads were then mapped to the bovine genome (UMD 3.1) using Tophat2 (version 2.0.9) to remove potential host DNA and RNA contaminations. Taxonomic profiles of the active rumen microbiota were generated using 16S rRNA extracted from T-metatranscriptomes. Briefly, post-QC bacterial and archaeal 16S rRNA reads were aligned to the V1-V3 region-enriched Greengenes database (version gg_13_8) and the V6-V8 region-enriched RIM-DB database, respectively. After that, mapped reads were taxonomically classified using the naive Bayesian approach in mothur.

To estimate rumen microbial functional profiles, non-rDNA sequences from all metagenomes (n = 48) were pooled, assembled, and annotated to create a functional reference database. Briefly, the pooled metagenomes were de novo assembled using Spherical program. Within Spherical, Velvet was set as the assembler with the kmer size of 31, Bowtie2 was set as the aligner, and 25% of total pooled sequences were subsampled as the input for each iteration of assembly with eight iterations in total. After the de novo assembly of pooled metagenomes reads, a total of 57,696,422 contigs with an average length of 144 bp (max 135,846 bp) and a N50 length of 140 bp were generated. Assembled contigs were then annotated using the blastx module in DIAMOND against the UniProt database, and only annotations with bitscore > 40 were kept for the downstream analysis. Overlapped annotations were filtered and converted to the GFF format using the MGKit package (https://bitbucket.org/setsuna80/mgkit). After discarding short contigs with length < 60 bp, 20,314,713 contigs (35.21%) were successfully annotated with an average length of 195 bp and a N50 length of 197 bp. To identify the functional categories of metagenomes, T-metatranscriptomes, and M-metatranscriptomes, non-rDNA/rRNA sequences were individually aligned to above annotated contigs using Bowtie2 and then were counted using HTSeq. Only reads mapped to contigs with eggNOG annotation information were further retrieved to calculate the abundances of genes and functional categories using MGKit.

## Statistical analysis
The comparison between efficient (L-RFI) and inefficient (H-RFI) animals were conducted using t test within each breed separately. Only microbial taxa with a relative abundance higher than 0.01% in at least 50% of individuals within each breed were considered as being observed and used for the analysis. Bacterial compositional profiles were summarized at phylum and genus levels, and archaeal communities were summarized at the species level. Relative abundances of microbial taxa were arcsine square root transformed, and then compared among breeds (using ANOVA) and between RFI groups within each breed (using t test).

Only functional categories and genes/transcripts with a minimum relative abundance of 0.01% in at least three samples within a dataset were considered as being detected. The abundance of each gene/transcript was then normalized into counts per million (cpm). Differential abundances of functional categories and genes (or transcripts) were compared among sequencing datasets, breeds, and RFI groups using DESeq2.

## Results
RFI values were not significantly different among the three beef cattle breeds (P = 0.73), but they were significant different between L- and H-RFI animals within each breed (P < 1.00e−5). 20,314,713 contigs (35.21%) from pooled metagenomes were successfully annotated based on the UniProt database. An average of 62.02 ± 0.56%, 33.04 ± 0.54%, and 32.19 ± 1.50% sequences from metagenomes, T-metatranscriptomes, and M-metatranscriptomes could be mapped back to these annotated contigs, respectively. The ratio of mapped metagenome reads to annotated genes (62.02%) is comparable with a recent rumen metagenomic study on dairy cattle (52.40%).

# Project Setup

Let's set up a project directory for the week, and talk a bit about project philosophy..

##  Creating a Project Directory

First, create a directory for you and the example project in the workshop share directory:

```bash
cd
mkdir -p /share/workshop/meta_workshop/$USER/meta_example
```

## Link raw fastq files

1. Next, go into that directory, create a raw data directory (we are going to call this 00-RawData) and cd into that directory. Lets then create symbolic links to the sample directories that contains the raw data.

    ```bash
    cd /share/workshop/meta_workshop/$USER/meta_example
    mkdir 00-RawData
    cd 00-RawData/
    ln -s /share/biocore/workshops/metagenomics_data/subset
    ```

    This directory now contains fastq files for all sample. The sequencing type is clearly indicated in the names of the fastq files.

1. We will create a directory to hold all of the scripts that will be used for this workshop, as well as a directory where we keep our reference sequence, genome, etc. The results of all our slurm script will output .out and .err files into the slurmout folder.

    ```bash
    cd ../
    mkdir -p scripts/slurmout
    mkdir References
    ```

1. Let's create a sample sheet for the project and store sample names in a file called samples.txt

    ```bash
    cd 00-RawData
    ls *_DNA_1.fastq |awk '{FS="_DNA"}{print $1}' - > ../scripts/allsamples.txt
    cat ../scripts/samples.txt
    ```

## Getting to know your data

1. Now, take a look at the raw data directory.

    ```bash
    ls /share/workshop/meta_workshop/$USER/meta_example/00-RawData
    ```

    You will see a list of the contents of each directory.

    ```bash
    ls *
    ```

    Lets get a better look at all the files in all of the directories.

    ```bash
    ls -lah */*
    ```

1. View the contents of the files using the 'less' command, when gzipped used 'zless' (which is just the 'less' command for gzipped files):

    ```bash
    cd 00-RawData
    zless ANG_301_DNA_1.fastq.gz
    ```

    Make sure you can identify which lines correspond to a read and which lines are the header, sequence, and quality values. Press 'q' to exit this screen.

1. Then, let's figure out the number of reads in this file. A simple way to do that is to count the number of lines and divide by 4 (because the record of each read uses 4 lines). In order to do this use cat to output the uncompressed file and pipe that to "wc" to count the number of lines:

    ```bash
    zcat ANG_301_DNA_1.fastq.gz | wc -l
    ```

    Divide this number by 4 and you have the number of reads in this file.

1. One more thing to try is to figure out the length of the reads without counting each nucleotide. First get the first 4 lines of the file (i.e. the first record):

    ```bash
    zcat ANG_301_DNA_1.fastq.gz  | head -2 | tail -1
    ```

    Note the header lines (1st and 3rd line) and sequence and quality lines (2nd and 4th) in each 4-line fastq block.

1. Then, copy and paste the DNA sequence line into the following command (replace [sequence] with the line):

    ```bash
    echo -n [sequence] | wc -c
    ```

    This will give you the length of the read. Also can do the bash one liner:

    ```bash
    echo -n $(zcat ANG_301_DNA_1.fastq.gz  | head -2 | tail -1) | wc -c
    ```

## Prepare our experiment folder for analysis

Now go back to your 'meta_example' directory and create one directory called '01-HTS_Preproc':

```bash
cd /share/workshop/meta_workshop/$USER/meta_example
mkdir References
mkdir 01-HTS_Preproc
```

The results of our preprocessing steps will be put into the 01-HTS_Preproc directory. The next step after that will go into a "02-..." directory, etc. You can collect scripts that perform each step, and notes and metadata relevant for each step, in the directory for that step. This way anyone looking to replicate your analysis has limited places to search for the commands you used. In addition, you may want to change the permissions on your original 00-RawData directory to "read only", so that you can never accidentally corrupt (or delete) your raw data. We won't worry about this here, because we've linked in sample folders.
Your directory should then look like the below:

```
$ ls
00-RawData  01-HTS_Preproc  References  scripts
```

### Questions you should now be able to answer.

1. How many reads are in the sample you checked?
2. How many basepairs is R1, how many is R2?
