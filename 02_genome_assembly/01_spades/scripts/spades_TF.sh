#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o spades_TF.$JOB_ID.log

### DESCRIPTION 

# Date: 18.4.2023
# Author: Elpida Skarlou

# Short summary #
# This shell script is designed to run the SPAdes genome assembler on Illumina sequencing data for the blowfly Chrysomya rufifacies, specifically two libraries called TF11 and TF19.

# Input for this script is a directory containing trimmed Illumina reads for Chrysomya rufifacies stored in the variable $TRIMMED. The script transfers the trimmed reads from the specified directory to a scratch directory ($SCRATCH) where the assembly will be performed.

# Output of this script is a genome assembly produced by SPAdes, which combines the reads from the TF11 and TF19 libraries, along with a log file and other intermediate files. The resulting assembly is transferred to the output directory specified in the variable $SPADES/outputs.

################################################################3333


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
rsync -av $TRIMMED/TF* .


### run SPAdes with and MismatchCorrector ###

# spades vesrion 3.15.5
# we did not correct with BayesHammer because short reads are already pretty accurate 

#### for TF - One assembly
TF11="TF11_Chrysomya-rufifacies_S3_"
TF19="TF19_Chrysomya-rufifacies_S4_"


#  create a single assembly using the reads from both libraries (TF11 and TF19)
spades.py --pe1-1 ./${TF19}R1_001.trimmed.fastq.gz --pe1-2 ./${TF19}R2_001.trimmed.fastq.gz --pe2-1 ./${TF11}R1_001.trimmed.fastq.gz --pe2-2 ./${TF11}R2_001.trimmed.fastq.gz -o . --threads 16 --only-assembler -m 400000


# remove all the trimmed reads
rm *.trimmed.fastq.gz

# syncing to final destinations

# log file
rsync -av $SPADES/scripts/spades_TF.$JOB_ID.log $SPADES/logs # be careful with the directories because if it crush there the files does not get transfered

# the spades outputs 
rsync -av * $SPADES/outputs

rm -r *
rm -rf /scratch/$USER/$JOB_ID