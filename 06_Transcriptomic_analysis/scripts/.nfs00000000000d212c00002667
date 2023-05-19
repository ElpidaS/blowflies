#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o rna_qc_trim.$JOB_ID.log

# aim: look at differential gene expression (DGE) in blowfly embryos to identify sex determination genes (doublesex, transformer etc)

# original Script; Rob Braid
# Edit; Elpida Skarlou 
# date of Edit 18/5/2023

# In normal systems (without maternal control over sex determination), the maternal transcriptome (maternal transcripts deposited into eggs) should be the same for both sexes. However, this changes after the onset of ZGA (zygotic genome activation).
# DGE (differential gene expression) between the sexes in this data could be due to 1) different maternal transcriptomes (because we have maternal SD here) or 2) ZGA. Depends what time ZGA switches on.
# You could also download the other developmental stages to look at DGE across the stages if you want.

### PART 2;  trim and QC the reads with fastqc + fastp
# the following is similar to the qc_trim_illumina.sh used to trimm Illumina reads

# in order to be able to use conda
source /ceph/users/eskarlou/miniconda3/etc/profile.d/conda.sh
conda activate fastp_qc

set -e
SCRATCH=/scratch/$USER/$JOB_ID/fastp_qc_RNA
RNA_GZ=/data/ross/flies/raw/Chrysomya_rufifacies/embryo_RNAseq/02_fastq_reads
RNA_TRIM=/data/ross/flies/raw/Chrysomya_rufifacies/embryo_RNAseq/03_trimmed_reads
QC_BEFORE=/data/ross/flies/analyses/blowflies/06_Transcriptomic_analysis/outputs/QC_before_trim
QC_AFTER=/data/ross/flies/analyses/blowflies/06_Transcriptomic_analysis/outputs/QC_after_trim
TRANSCRIPTOME=/data/ross/flies/analyses/blowflies/06_Transcriptomic_analysis

mkdir -p $SCRATCH
cd $SCRATCH

# sync files between $RNA_GZ and $SCRATCH 
rsync -av  $RNA_GZ/*.fastq.gz $SCRATCH

# fastqc pre-trim
fastqc -t 4 *.fastq.gz

# trim reads
for file in $(ls *.fastq.gz)
do
	# Extract the base name of the file (without the extension) and remove the "..." suffix
	# This base name will be used to construct the output file names for the trimmed and filtered reads
	base=$(basename $file ".fastq.gz") 
  	fastp -i ./${base}.fastq.gz -o ./${base}.trimmed.fastq.gz 
done

# fastqc post-trim
fastqc -t 4 *.trimmed.fastq.gz

# syncing to final destinations

# log file
mv $TRANSCRIPTOME/scripts/rna_qc_trim.$JOB_ID.log $TRANSCRIPTOME/logs

# after
rsync -av *trimmed_fastqc.html $QC_AFTER

# for the fastqc outputs (.html) 
# before trimming
rsync -av *pass_fastqc.html $QC_BEFORE

# for the trimmed reads
rsync -av *.trimmed.fastq.gz $RNA_TRIM

rm -rf $SCRATCH