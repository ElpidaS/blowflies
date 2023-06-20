#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o sex_limited_g_extraction.$JOB_ID.log

# aim -> extract the sex limited genes (their names are found in the two txt files produced by the R script DGE_filtered) from the Trinity.fasta (transcriptome). 

# Author; Elpida Skarlou
# Date; 19/06/2023

set -e

SCRATCH=/scratch/$USER/$JOB_ID/sex_limited_g_extraction
TXT_FILES_IN=/data/ross/flies/analyses/blowflies/06_Transcriptomic_analysis/outputs/sex_limited_genes
TRANSCRIPTOME=/data/ross/flies/analyses/blowflies/06_Transcriptomic_analysis/outputs/transdecoder_filter
FASTA_OUT=/data/ross/flies/analyses/blowflies/06_Transcriptomic_analysis/outputs/sex_limited_genes
TRANSCRIPT=/data/ross/flies/analyses/blowflies/06_Transcriptomic_analysis

mkdir -p $SCRATCH
cd $SCRATCH

source /ceph/users/eskarlou/miniconda3/etc/profile.d/conda.sh



# sync input data
rsync -av  $TXT_FILES_IN/model01_98.txt $SCRATCH
rsync -av  $TRANSCRIPTOME/transcripts.fa $SCRATCH


########### model 01 ###########
assembly_file="transcripts.fa" # use that instead of Trinity.fasta, becasuse I used transcript.fa to find the genes and do the DGE analysis
header_file="model01_98.txt"

# Convert the header file to Unix line endings (LF only)
conda activate dos2unix
dos2unix "$header_file"

conda activate seqtk
# Read the header file line by line
while IFS=$'\t' read -r header _rest
do
    header=$(echo "$header" | tr -d '\r')
    echo "$header" > "${header}.txt"
    seqtk subseq "$assembly_file" "${header}.txt" > "${header}.fa"
    cat "${header}.txt"
    rm "${header}.txt"
done < "$header_file"

output_file="model01_98.fasta"

# Get a list of all FASTA files in the directory
fasta_files=$(find . -maxdepth 1 -type f -name "*.fa")

# Concatenate the files
for file in $fasta_files; do
    cat "$file" >> "$output_file"
done

rm *.fa

################ model 02 ################
rsync -av  $TRANSCRIPTOME/transcripts.fa $SCRATCH
rsync -av  $TXT_FILES_IN/model02_x2.txt $SCRATCH
header_file="model02_x2.txt"

conda activate dos2unix
dos2unix "$header_file"

conda activate seqtk
# Read the header file line by line
while IFS=$'\t' read -r header _rest
do
    header=$(echo "$header" | tr -d '\r')
    echo "$header" > "${header}.txt"
    seqtk subseq "$assembly_file" "${header}.txt" > "${header}.fa"
    cat "${header}.txt"
    rm "${header}.txt"
done < "$header_file"

output_file="model02_x2.fasta"

# Get a list of all FASTA files in the directory
fasta_files=$(find . -maxdepth 1 -type f -name "*.fa")

# Concatenate the files
for file in $fasta_files; do
    cat "$file" >> "$output_file"
done


rm *.fa

mv $TRANSCRIPT/scripts/sex_limited_g_extraction.$JOB_ID.log $TRANSCRIPT/logs

rsync -av model*.fasta $FASTA_OUT

rm -rf $SCRATCH
