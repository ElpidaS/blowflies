#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o trimgalore.$JOB_ID.log


# Author: Elpida Skarlou
# Date: 16-5-2023

# Description: 


source /ceph/users/eskarlou/miniconda3/etc/profile.d/conda.sh
conda activate trim_galore

set -e
SCRATCH=/scratch/$USER/$JOB_ID/fastp_qc
PIC_ILLUMINA=/data/ross/sequencing/raw/blowflies/picard_lab_illumina
IL_TRIM=/data/ross/flies/analyses/blowflies/01_QC_trimming/01_illumina_reads/02_TrimGalore_trim

mkdir -p $SCRATCH
cd $SCRATCH

# sync files between $PIC_ILLUMINA and $SCRATCH 
rsync -av  $PIC_ILLUMINA/*.fastq.gz .

# fastqc pre-trim
fastqc -t 4 *.fastq.gz

# trim reads
for file in $(ls *1_001.fastq.gz)
do
	# Extract the base name of the file (without the extension) and remove the "..." suffix
	# This base name will be used to construct the output file names for the trimmed and filtered reads
	base=$(basename $file "1_001.fastq.gz") 
  	trim_galore --paired --quality 20 --length 50 --stringency 0  --illumina --fastqc --paired --output_dir ./ --gzip --output_file_1 /${base}1_001.trimmed.galore.fastq.gz --output_file_2 /${base}2_001.trimmed.galore.fastq.gz ./${base}1_001.fastq.gz ./${base}2_001.fastq.gz
done


# --paired: This option specifies that the reads are paired-end reads.
# --quality 20: This option specifies that bases with quality scores below 20 will be trimmed.
# --length 50: This option specifies that reads shorter than 50 bases will be trimmed.
# --stringency 0: This option specifies that TrimGalore will be less stringent in its trimming. This can be useful if you are working with low-quality reads.
# --illumina: This option specifies that the reads were sequenced on an Illumina platform.
# --fastqc: This option specifies that FastQC will be used to generate a report on the quality of the trimmed reads.
# --output_dir trimmed_reads: This option specifies the directory where the trimmed reads will be written.
# --gzip: This option specifies that the trimmed reads will be compressed using gzip.


# log file
mv $IL_TRIM/scripts/fastp.$JOB_ID.log $IL_TRIM/logs


# for the trimmed reads
rsync -av ./*.trimmed.galore.fastq.gz $IL_TRIM/outputs/


rm -rf /scratch/$USER/$JOB_ID