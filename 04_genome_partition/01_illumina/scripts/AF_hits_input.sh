#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o AF_hits.$JOB_ID.log

### Description 
# Main goal; 
# HITS FILES
# Produce oen of the 2 main files (coverage and hits) needed to run Blobtools

# More specifically; 
# Based on the script, it looks like the goal is 
# to prepare two types of files needed for running the software BlobTools: 
# a coverage file and a hits file. The coverage file contains information 
# about read coverage for each contig in the genome assembly, while the hits 
# file contains information about homology to other sequences in a reference database.

# Author; Elpida Skarlou
# Date; 01/05/2023

### Precode stuff

set -e # exit immediately on error

# PATHS used in this script
SCRATCH=/scratch/$USER/$JOB_ID/hits
AF_ASSEMBLY=/data/ross/flies/analyses/blowflies/02_genome_assembly/01_spades/outputs/AF_assembly/
TRIMMED_READS=/data/ross/flies/raw/Chrysomya_rufifacies/illumina/
GENOME_PART=/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/

# make $SCRATCH if it doesn't exist
mkdir -p $SCRATCH 

cd $SCRATCH

rsync -av $AF_ASSEMBLY/contigs.fasta $SCRATCH

### one (or more) hits file(s), e.g. example/blast.out (have to download)
	source /ceph/users/eskarlou/miniconda3/etc/profile.d/conda.sh
  conda activate BLAST

	blastn \
	 -query contigs.fasta \
	 -db nt \
	 -outfmt '6 qseqid staxids bitscore std' \
	 -max_target_seqs 1 \
	 -max_hsps 1 \
	 -evalue 1e-25 \
   -remote # perform the BLAST search in the NCBI servers
   
 # Post-code stuff

# remove all the input files 
rm -rf $SCRATCH/contigs.fasta


# log file
mv -av $GENOME_PART/scripts/AF_hits.$JOB_ID.log $GENOME_PART/logs

# rsync all the outputs 
rsync -av $SCRATCH $GENOME_PART/outputs/AF_out/