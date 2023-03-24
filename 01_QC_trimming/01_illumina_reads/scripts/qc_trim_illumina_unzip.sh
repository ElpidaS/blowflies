#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o qc_fastp_unzip.$JOB_ID.log


# Author: Elpida Skarlou
# Date: 2024-03-22

# Description: 
#
# The same as qc_trim_illumina.sh, but I unzip the reads before using them
#
# This script performs quality control on Illumina reads using FastQC and trims reads using Fastp. 
# It takes as input unzipped paired-end Illumina reads located in $PIC_ILLUMINA
# ,and outputs FastQC reports for both pre- and post-trimmed reads, as well as trimmed unzipped reads. 
# The outputs can be found in /data/ross/flies/analyses/blowflies/01_QC_trimming. 


# in order to be able to use fastp through conda
source /ceph/users/eskarlou/miniconda3/etc/profile.d/conda.sh
conda activate fastp_qc

set -e
SCRATCH=/scratch/$USER/$JOB_ID/fastp_qc
PIC_ILLUMINA=/data/ross/sequencing/raw/blowflies/picard_lab_illumina
IL_TRIM=/data/ross/flies/analyses/blowflies/01_QC_trimming/01_illumina_reads/

mkdir -p $SCRATCH
cd $SCRATCH

# sync files between $PIC_ILLUMINA and $SCRATCH and unzip them
rsync -av $PIC_ILLUMINA/*.fastq.gz . && gunzip *.fastq.gz

# fastqc pre-trim
fastqc -t 4 *.fastq

# trim reads
for file in $(ls *1_001.fastq)
do
	# Extract the base name of the file (without the extension) and remove the "..." suffix
	# This base name will be used to construct the output file names for the trimmed and filtered reads
	base=$(basename $file "1_001.fastq") 
  	fastp -i ./${base}1_001.fastq -I ./${base}2_001.fastq -o ./${base}1_001.trimmed.fastq -O ./${base}2_001.trimmed.fastq
done

# fastqc post-trim
fastqc -t 4 *.trimmed.fastq

# zip all the fastq
gzip *.trimmed.fastq

# syncing to final destinations

# log file
mv $IL_TRIM/scripts/qc_fastp_unzip.$JOB_ID.log $IL_TRIM/logs

# for the fastqc outputs (.html) / both pre and post trimming
rsync -av ./*.html $IL_TRIM/outputs/01_QC

# for the trimmed reads
rsync -av ./*.trimmed.fastq.gz $IL_TRIM/outputs/02_trimmed_illumina

rm *.gz
rm *.html
rm -rf /scratch/$USER/$JOB_ID
