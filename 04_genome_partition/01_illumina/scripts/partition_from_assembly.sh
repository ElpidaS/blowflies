#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o partition_from_assembly.$JOB_ID.log

# activate conda
source /ceph/users/eskarlou/miniconda3/etc/profile.d/conda.sh
conda activate seqtk

# PATHS
AF_assembly=/data/ross/flies/analyses/blowflies/02_genome_assembly/01_spades/outputs/AF_assembly/
TF_assembly=/data/ross/flies/analyses/blowflies/02_genome_assembly/01_spades/outputs/TF_assembly/
AF_part=/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/outputs/AF_out/partition/
TF_part=/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/outputs/TF_out/partition/

# AF
seqtk subseq $AF_assembly/contigs.fasta $AF_part/contig_keep_list_AF.txt > $AF_part/AF_contig.filtered.fasta

# TF
seqtk subseq $TF_assembly/contigs.fasta $TF_part/contig_keep_list_TF.txt > $TF_part/TF_contig.filtered.fasta