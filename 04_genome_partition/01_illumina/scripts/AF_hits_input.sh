#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o AF_hits.$JOB_ID.log
#$ -M elpida.skarlou@evobio.eu


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

# Set email address
EMAIL="elpidaskarlou@gmail.com"

# make $SCRATCH if it doesn't exist
mkdir -p $SCRATCH 

cd $SCRATCH

rsync -av $AF_ASSEMBLY/contigs.fasta $SCRATCH
rsync -av /ceph/software/databases/ncbi $SCRATCH

### one (or more) hits file(s), e.g. example/blast.out (have to download)
	source /ceph/users/eskarlou/miniconda3/etc/profile.d/conda.sh
  conda activate BLAST

cd ./ncbi

# remove the tar.gz.md5 files
rm -rf *.tar.gz.md5

# add execute persmission to all files
chmod a+x *

# combine all the databases
blastdb_aliastool -dblist "nt.00 nt.01 nt.02 nt.03 nt.04 nt.05 nt.06 nt.07 nt.08 nt.09 nt.10 nt.11 nt.12 nt.13 nt.14 nt.15 nt.16 nt.17 nt.18 nt.19 nt.20 nt.21 nt.22 nt.23 nt.24 nt.25 nt.26 nt.27 nt.28 nt.29 nt.30 nt.31 nt.32 nt.33 nt.34 nt.35 nt.36 nt.37 nt.38 nt.39 nt.40 nt.41 nt.42 nt.43 nt.44 nt.45 nt.46 nt.47 nt.48 nt.49 nt.50 nt.51 nt.52 nt.53 nt.54 nt.55 nt.56" -dbtype nucl -out nt_combined -title "combined database"

cd $SCRATCH

	blastn \
	 -query contigs.fasta \
   # this is the database allready available in the cluster
   -db ./ncbi/nt \
   -outfmt '6 qseqid staxids bitscore std' \
   -max_target_seqs 1 \
   -max_hsps 1 \
   -evalue 1e-25 \
   -num_threads 4
   
 # Post-code stuff

# remove all the input files 
rm -rf $SCRATCH/contigs.fasta


# log file
mv -av $GENOME_PART/scripts/AF_hits.$JOB_ID.log $GENOME_PART/logs

# rsync all the outputs 
rsync -av $SCRATCH $GENOME_PART/outputs/AF_out/



