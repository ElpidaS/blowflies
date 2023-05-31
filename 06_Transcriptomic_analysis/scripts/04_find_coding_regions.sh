#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o coding_regions_find.$JOB_ID.log

# aim: look at differential gene expression (DGE) in blowfly embryos to identify sex determination genes (doublesex, transformer etc)

# original Script; Rob Braid
# Edit; Elpida Skarlou 
# date of Edit 22/5/2023

# In normal systems (without maternal control over sex determination), the maternal transcriptome (maternal transcripts deposited into eggs) should be the same for both sexes. However, this changes after the onset of ZGA (zygotic genome activation).
# DGE (differential gene expression) between the sexes in this data could be due to 1) different maternal transcriptomes (because we have maternal SD here) or 2) ZGA. Depends what time ZGA switches on.
# You could also download the other developmental stages to look at DGE across the stages if you want.

# in order to be able to use conda
source /ceph/users/eskarlou/miniconda3/etc/profile.d/conda.sh
# transdecoder (conda enviroment), icludes BLAST and transdecoder
conda activate transdecoder

set -e

SCRATCH=/scratch/$USER/$JOB_ID/coding_regions_find
TRINITY=/data/ross/flies/analyses/blowflies/06_Transcriptomic_analysis/outputs/trinity_rna_assemblies/
OUT=/data/ross/flies/analyses/blowflies/06_Transcriptomic_analysis/outputs/transdecoder_filter

# need to change the inputs according to my own

mkdir -p $SCRATCH
cd $SCRATCH

# sync all the trimmed rna reads
rsync -av  $TRINITY/Trinity.fasta $SCRATCH


### PART 3;  filter the transcriptome with transdecoder (find codding regions within transcriptome)

## with homology search:
TransDecoder.LongOrfs -m 50 -t Trinity.fasta

# you might need to make a blast database, in which case:
# makeblastdb -in ceph/users/rbaird/software/uniprot/uniprot_sprot.fasta -dbtype prot

blastp -query Trinity.fasta.transdecoder_dir/longest_orfs.pep -db /ceph/users/rbaird/software/uniprot/uniprot_sprot.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 > Trinity.fasta.blastp.outfmt6

TransDecoder.Predict -t Trinity.fasta --retain_blastp_hits Trinity.fasta.blastp.outfmt6

# move the output to the output directory
mv * $OUT

# delete the $SCRATCH file
rm -rf $SCRATCH
