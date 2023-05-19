#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o rna_assembly.$JOB_ID.log

# aim: look at differential gene expression (DGE) in blowfly embryos to identify sex determination genes (doublesex, transformer etc)

# original Script; Rob Braid
# Edit; Elpida Skarlou 
# date of Edit 18/5/2023

# In normal systems (without maternal control over sex determination), the maternal transcriptome (maternal transcripts deposited into eggs) should be the same for both sexes. However, this changes after the onset of ZGA (zygotic genome activation).
# DGE (differential gene expression) between the sexes in this data could be due to 1) different maternal transcriptomes (because we have maternal SD here) or 2) ZGA. Depends what time ZGA switches on.
# You could also download the other developmental stages to look at DGE across the stages if you want.

### PART 3;  Assemblying the trasncriptome 

# in order to be able to use conda
source /ceph/users/eskarlou/miniconda3/etc/profile.d/conda.sh
conda activate Trinity

set -e

SCRATCH=/scratch/$USER/$JOB_ID/assembly_RNA
RNA_TRIM=/data/ross/flies/raw/Chrysomya_rufifacies/embryo_RNAseq/03_trimmed_reads
OUT_ASS=/data/ross/flies/analyses/blowflies/06_Transcriptomic_analysis/outputs/rna_assemblies/
TRANSCRIPTOME=/data/ross/flies/analyses/blowflies/06_Transcriptomic_analysis

mkdir -p $SCRATCH
cd $SCRATCH

# sync all the trimmed rna reads
rsync -av  $RNA_TRIM/* $SCRATCH


### 3. assemble a transcriptome with Trinity (https://github.com/trinityrnaseq/trinityrnaseq)
for file in $(ls *.gz)
do
	base=$(basename $file "_pass.trimmed.fastq.gz")
	Trinity --seqType fq --single ${base}"_pass.trimmed.fastq.gz" --CPU 12 --max_memory 50G --output ${base}"_transcriptome.fasta"
done

mv $TRANSCRIPTOME/scripts/rna_assembly.$JOB_ID.log $TRANSCRIPTOME/logs

rsync -av *.fasta $OUT_ASS

rm -rf $SCRATCH