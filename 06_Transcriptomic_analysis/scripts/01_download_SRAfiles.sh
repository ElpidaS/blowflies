#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o rna_download.$JOB_ID.log

# aim: look at differential gene expression (DGE) in blowfly embryos to identify sex determination genes (doublesex, transformer etc)

# original Script; Rob Braid
# Edit; Elpida Skarlou 
# date of Edit 18/5/2023

# In normal systems (without maternal control over sex determination), the maternal transcriptome (maternal transcripts deposited into eggs) should be the same for both sexes. However, this changes after the onset of ZGA (zygotic genome activation).
# DGE (differential gene expression) between the sexes in this data could be due to 1) different maternal transcriptomes (because we have maternal SD here) or 2) ZGA. Depends what time ZGA switches on.
# You could also download the other developmental stages to look at DGE across the stages if you want.

### PART 1; Download the reads from NCBI + dump the SRA files into fastq files

# in order to be able to use conda
source /ceph/users/eskarlou/miniconda3/etc/profile.d/conda.sh
conda activate sra-tools

set -e
SCRATCH=/scratch/$USER/$JOB_ID/download_rna
RAW_RNA=/data/ross/flies/raw/Chrysomya_rufifacies/embryo_RNAseq/01_raw_reads
FASTQ_RNA=/data/ross/flies/raw/Chrysomya_rufifacies/embryo_RNAseq/02_fastq_reads
TRANSCRIPTOME=/data/ross/flies/analyses/blowflies/06_Transcriptomic_analysis

# Getting SRA files from NCBI
# navigate to the SRA page, e.g. for the female embryos https://www.ncbi.nlm.nih.gov/sra/SRX2793047[accn]
# Click the SRR number under "Run"
# Click the 'data access' tab
# copy the AWS link (in this case https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5519687/SRR5519687)

### 1. download the reads to the cluster (do this within the directory you want to download them to):
# download all embryo replicates - looks like they have 3 for each sex 
# Instrument: Illumina HiSeq 2000

mkdir -p $SCRATCH 
cd $SCRATCH

### MALES ###

# embryos males (SRR5519706)
wget -O embryos_males_SRR5519706.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5519706/SRR5519706

# embryos males (SRR5519705)
wget -O embryos_males_SRR5519705.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5519705/SRR5519705

# embryos males (SRR5519698)
wget -O embryos_males_SRR5519698.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5519698/SRR5519698

### FEMALES ###

# embryos females (SRR5519687)
wget -O embryos_females_SRR5519687.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5519687/SRR5519687

# embryos females (SRR5519676)
wget -O embryos_females_SRR5519676.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5519676/SRR5519676

# embryos females (SRR5519675)
wget -O embryos_females_SRR5519675.sra https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5519675/SRR5519675


# use sra-tools to dump the SRA files into fastq files (https://anaconda.org/bioconda/sra-tools)
for file in $(ls embryos*)
do
	fastq-dump --gzip --skip-technical --readids --read-filter pass --dumpbase --clip ${file}
done
# this should output the forward and reverse files as something like SRR5519687.1.fq.gz, SRR5519687.2.fq.gz 
## The rest of the libraries were prepared as 100 bp single-end reads. -> those include the early embryo ones
## so I deleted  --split-3 which was the one that splits the reads into forward and reverse


# log file
mv $TRANSCRIPTOME/scripts/rna_download.$JOB_ID.log $TRANSCRIPTOME/logs

# raw RNA reads to the raw folder 
rsync -av *.sra $RAW_RNA

# for the trimmed reads
rsync -av *.fastq.gz $FASTQ_RNA

rm -rf $SCRATCH