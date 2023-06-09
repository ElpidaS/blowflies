#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o abundance_transc.$JOB_ID.log

# aim: look at differential gene expression (DGE) in blowfly embryos to identify sex determination genes (doublesex, transformer etc)

# original Script; Rob Braid
# Edit; Elpida Skarlou 
# date of Edit 18/5/2023

# In normal systems (without maternal control over sex determination), the maternal transcriptome (maternal transcripts deposited into eggs) should be the same for both sexes. However, this changes after the onset of ZGA (zygotic genome activation).
# DGE (differential gene expression) between the sexes in this data could be due to 1) different maternal transcriptomes (because we have maternal SD here) or 2) ZGA. Depends what time ZGA switches on.
# You could also download the other developmental stages to look at DGE across the stages if you want.

# in order to be able to use conda
source /ceph/users/eskarlou/miniconda3/etc/profile.d/conda.sh
conda activate kallisto

set -e

SCRATCH=/scratch/$USER/$JOB_ID/abundance_transc
TRIMM=/data/ross/flies/raw/Chrysomya_rufifacies/embryo_RNAseq/03_trimmed_reads
TRINITY=/data/ross/flies/analyses/blowflies/06_Transcriptomic_analysis/outputs/trinity_rna_assemblies
CDS_DIR=/data/ross/flies/analyses/blowflies/06_Transcriptomic_analysis/outputs/transdecoder_filter
OUT=/data/ross/flies/analyses/blowflies/06_Transcriptomic_analysis/outputs/transcript_abundance

mkdir -p $SCRATCH
cd $SCRATCH

# sync all the trimmed rna reads
rsync -av  $TRIMM/* $SCRATCH

# sync the RNA assembly
rsync -av $TRINITY/Trinity.fasta $SCRATCH
rsync -av $CDS_DIR/transcripts.fa $SCRATCH

### PART 4; Quantify transcript abundance with Kallisto (https://github.com/pachterlab/kallisto)

# Whole fasta file

# build kallisto index
kallisto index -i trinity.kallisto.idx Trinity.fasta

# quantify
for file in $(ls *.fastq.gz)
do
	
  base=$(basename ${file} "_pass.trimmed.fastq.gz")
  kallisto quant --index=trinity.kallisto.idx --output-dir=${base}_out --threads=16 --single -l 100 -s 0.1 --plaintext ${base}_pass.trimmed.fastq.gz 
  mv ${base}_out/abundance.tsv ${base}_out/${base}.tsv 
done

#########3 Only for Cds #################

# we transfered the transcoder.cds to a fa file
# this .fa file will be used in kallisto to count the abundance of ONLY coding RNA

# build kallisto index
kallisto index -i transcripts.kallisto.idx transcripts.fa

# quantify
for file in $(ls *.fastq.gz)
do
	
  base=$(basename ${file} "_pass.trimmed.fastq.gz")
  kallisto quant --index=transcripts.kallisto.idx --output-dir=${base}_out --threads=16 --single -l 100 -s 0.1 --plaintext ${base}_pass.trimmed.fastq.gz 
  mv ${base}_out/abundance.tsv ${base}_out/${base}.transcript.tsv 
done

# for the -l and -s values you have to look at the fastq restults (after trimming)

# this will output counts in TPM (transcripts per million) in .tsv files - there will be one file per sample/replicate. 
# These can be read into R for plotting and DGE analyses.

mv ./*_out $OUT

rm -rf $SCRATCH