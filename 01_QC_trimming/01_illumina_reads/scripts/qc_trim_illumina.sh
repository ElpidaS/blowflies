#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o fastp.$JOB_ID.log

'
Author: Elpida Skarlou
Date: 2023-03-22

Description: 
This script performs quality control on Illumina reads using FastQC and trims reads using Fastp. 
It takes as input gzipped paired-end Illumina reads located in $PIC_ILLUMINA
,and outputs FastQC reports for both pre- and post-trimmed reads, as well as trimmed gzipped reads. 
The outputs can be found in /data/ross/flies/analyses/blowflies/01_QC_trimming. 
'

# in order to be able to use fastp through conda
source /ceph/users/eskarlou/miniconda3/etc/profile.d/conda.sh
conda activate fastp

set -e
SCRATCH=/scratch/$USER/$JOB_ID/fastp_qc
PIC_ILLUMINA=/data/ross/sequencing/raw/blowflies/picard_lab_illumina

mkdir -p $SCRATCH
cd $SCRATCH

# sync files between $PIC_ILLUMINA and $SCRATCH 
rsync -av  $PIC_ILLUMINA/*.fastq.gz .

# fastqc pre-trim
fastqc -t 4 *.fastq.gz

# trim reads
for file in $(ls *.fastq.gz)
do
	# Extract the base name of the file (without the extension) and remove the "..." suffix
	# This base name will be used to construct the output file names for the trimmed and filtered reads
	base=$(basename $file "1_001.fastq.gz") 
  	fastp -i ${base}1_001.fastq.gz -I ${base}2_001.fastq.gz -o ${base}1_001.trimmed.fastq.gz -O ${base}2_001.trimmed.fastq.gz
done

# fastqc post-trim
fastqc -t 4 *.trimmed.fastq.gz

# syncing to final destinations

# log file
rsync -av fastp.$JOB_ID.log /data/ross/flies/analyses/blowflies/01_QC_trimming/01_illumina_reads/logs

# for the fastqc outputs (.html) / both pre and post trimming
rsync -av *.html /data/ross/flies/analyses/blowflies/01_QC_trimming/01_illumina_reads/outputs/01_QC

# for the trimmed reads
rsync -av *.trimmed.fastq.gz /data/ross/flies/analyses/blowflies/01_QC_trimming/01_illumina_reads/outputs/02_trimmed_illumina

rm *.gz
rm *.html
rm -rf /scratch/$USER/$JOB_ID
