#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o spades_AF.$JOB_ID.log

### DESCRIPTION 

# Date: 4.4.2023
# Author: Elpida Skarlou

# This script performs genome assembly for Illumina sequencing data using SPAdes software. It first transfers trimmed Illumina reads to a temporary
# directory, then performs one assembly for one sequencing library
# (AF). The output assembly is synced to a final destination and temporary 
# files are deleted. 
# The input files are; trimmed Illumina fastq files, and 
# the output files are; genome assembly files.

set -e # exit immediately on error

# PATHS used in this script
SCRATCH=/scratch/$USER/$JOB_ID/SPAdes 
TRIMMED=/data/ross/flies/raw/Chrysomya_rufifacies/illumina
SPADES=/data/ross/flies/analyses/blowflies/02_genome_assembly/01_spades

# get into the SPAdes conda enviroment
source /ceph/users/eskarlou/miniconda3/etc/profile.d/conda.sh
conda activate SPAdes

# make $SCRATCH if it doesn't exist
mkdir -p $SCRATCH 

cd $SCRATCH

# transfer $TRIMMED/*.trimmed.fastq.gz to the qurent dir
rsync -av $TRIMMED/*.trimmed.fastq.gz .


### run SPAdes with and MismatchCorrector ###

# spades vesrion 3.15.5
# we did not correct with BayesHammer because short reads are already pretty accurate 

#### for AF
AF="AF7_Chrysomya-rufifacies_S2_"

spades.py -1 ./${AF}R1_001.trimmed.fastq.gz -2 ./${AF}R2_001.trimmed.fastq.gz -o . --threads 16 --isolate --only-assembler  

# remove all the trimmed reads
rm *.trimmed.fastq.gz

# syncing to final destinations

# log file
rsync -av $SPADES/scripts/fastp.$JOB_ID.log $SPADES/logs

# the spades outputs 
rsync -av * $SPADES/outputs

rm -r *
rm -rf /scratch/$USER/$JOB_ID