#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o kmer.$JOB_ID.log

################################ Description ####################################################

# This bash script performs k-mer analysis on Illumina sequencing reads of blowflies. 
#It uses KMC to count the k-mers and GenomeScope 2 to generate a report and informative plots.

#Inputs:

# Illumina sequencing (untrimmed) reads in .fastq.gz format located in 
# /data/ross/sequencing/raw/blowflies/picard_lab_illumina directory.

# Outputs:

# K-mer counts and histograms generated by KMC in the ./kmc directory.
# GenomeScope 2 report and informative plots in the ./genomescope directory.
# A log file in /data/ross/flies/analyses/blowflies/03_kmer/01_illumina_reads/logs directory.

# Author: Elpida Skarlou
# Date: 2023-03-27

#################################################################################################

set -e 

# PATHS

SCRATCH=/scratch/$USER/$JOB_ID/Kmer 
PIC_ILLUMINA=/data/ross/sequencing/raw/blowflies/picard_lab_illumina
KMER=/data/ross/flies/analyses/blowflies/03_kmer/01_illumina_reads/

# make $SCRATCH if it doesn't exist
mkdir -p $SCRATCH 

cd $SCRATCH

# transfer $TRIMMED/*.trimmed.fastq.gz to $PWD
rsync -av $PIC_ILLUMINA/*.fastq.gz .

# copy genomescope.R to SCRATCH (.)
cp /ceph/users/eskarlou/genomescope2.0/genomescope.R .
### KMC ### 
# count k-mers in FASTQ files

# initialize conda
eval "$(conda shell.bash hook)"

# activate your conda environment
conda activate KMC


for file in $(ls *1_001.fastq.gz)
do
  base=$(basename $file "1_001.fastq.gz")
  # combine fasta.gz
  cat ${base}1_001.fastq.gz ${base}2_001.fastq.gz > ${base}_files.fastq.gz
  # count  kmers
  kmc -k21 -t10 -m64 -ci1 -cs10000 -fq1 ${base}_files.fastq.gz ${base}_kmer_counts .  # use fq1 for fastq.gz input files
  kmc_tools transform ${base}_kmer_counts.kmc_pre_hist histogram ${base}_kmer_k21.histo -cx100000 # use counts.kmc_pre_hist file to produce .histo file


### Genome Scope 2 ### 
# produces a report and several informative plots describing the genome properties

# copy genomescope.R to SCRATCH (.)
cp /ceph/users/eskarlou/genomescope2.0/genomescope.R .

conda activate for_genomescope


# the command below follows the format bellow (suggested by the genimescope git hub)
# $ Rscript genomescope.R histogram_file k-mer_length read_length output_dir [kmer_max] [verbose]

# TF11
/ceph/users/eskarlou/miniconda3/envs/for_genomescope/bin/Rscript genomescope.R -i TF11_Chrysomya-rufifacies_S3_R_kmer_k21.histo -o ./ -k 21

# TF19
/ceph/users/eskarlou/miniconda3/envs/for_genomescope/bin/Rscript genomescope.R -i TF19_Chrysomya-rufifacies_S4_R_kmer_k21.histo -o ./ -k 21

# AF
/ceph/users/eskarlou/miniconda3/envs/for_genomescope/bin/Rscript genomescope.R -i AF7_Chrysomya-rufifacies_S2_R_kmer_k21.histo -o ./ -k 21

# syncing to final destinations #

# log file
mv $KMER/scripts/kmer.$JOB_ID.log $KMER/logs

rm -rf *.gz

rsync -av . $KMER/outputs/

rm -rf *
rm -rf /scratch/$USER/$JOB_ID