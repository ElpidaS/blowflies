#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o spades.$JOB_ID.log

### DESCRIPTION 

# Date: 23.03.2023
# Author: Elpida Skarlou

# The input of the script is trimmed Illumina paired-end reads in the TRIMMED directory. The script uses BayesHammer for error correction, and then it run
# SPAdes with the corrected reads and MismatchCorrector to perform genome assembly for both types of female blowflies. 
# The output of the script is two assembled genomes in two directories, TF_assembly (Male producing Females) and AF_assembly (Female producing Females). 
# The script also transfers the log file and output files to the SPADES directory and removes intermediate files from the scratch directory.

set -e # exit immediately on error

# PATHS used in this script
SCRATCH=/scratch/$USER/$JOB_ID/SPAdes 
TRIMMED=/data/ross/flies/analyses/blowflies/01_QC_trimming/01_illumina_reads/outputs/02_trimmed_illumina
SPADES=/data/ross/flies/analyses/blowflies/02_genome_assembly/01_spades

# get into the SPAdes conda enviroment
source /ceph/users/eskarlou/miniconda3/etc/profile.d/conda.sh
conda activate SPAdes

# make $SCRATCH if it doesn't exist
mkdir -p $SCRATCH 

cd $SCRATCH

# transfer $TRIMMED/*.trimmed.fastq.gz to the qurent dir
rsync -av $TRIMMED/*.trimmed.fastq.gz .


#### run BayesHammer for error correction ####

mkdir 01_corrected_reads

for file in $(ls *.trimmed.fastq.gz)
do
	base=$(basename $file "1_001.trimmed.fastq.gz") 
  	spades.py -1 ${base}1_001.trimmed.fastq.gz -2 ${base}2_001.trimmed.fastq.gz --only-error-correction -o ./01_corrected_reads/${base}1_001.trimmed.corrected.fastq.gz -O ./01_corrected_reads/${base}2_001.trimmed.corrected.fastq.gz
done

### run SPAdes with BayesHammer-corrected reads and MismatchCorrector ###

mkdir TF_assembly
mkdir AF_assembly

# for TF
TF11="TF11_Chrysomya-rufifacies_S3_"
TF19="TF19_Chrysomya-rufifacies_S4_"
spades.py --pe1-1 "${TF11}"R1_001.trimmed.corrected.fastq.gz --pe1-2 "${TF19}"R1_001.trimmed.corrected.fastq.gz --pe2-1 "${TF11}"R2_001.trimmed.corrected.fastq.gz --pe2-2 "${TF19}"R2_001.trimmed.corrected.fastq.gz -o TF_assembly --threads 16 --isolate --mismatch-correction

# for AF
AF="AF7_Chrysomya-rufifacies_S2_"
spades.py -1 "${AF}"R1_001.trimmed.corrected.fastq.gz -2 "${AF}"R2_001.trimmed.corrected.fastq.gz -o AF_assembly --threads 16 --isolate --mismatch-correction

# remove all the trimmed reads
rm *.trimmed.fastq.gz

# syncing to final destinations

# log file
rsync -av $SPADES/scripts/fastp.$JOB_ID.log $SPADES/logs

# for the fastqc outputs (.html) / both pre and post trimming
rsync -av * $SPADES/outputs

rm -r *
rm -rf /scratch/$USER/$JOB_ID