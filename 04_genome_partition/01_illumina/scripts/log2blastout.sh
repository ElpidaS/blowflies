#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o logs2out.$JOB_ID.log

### Description
# For the BLAST search of both AF and TF there was no -out specified. So the output can be found in the .log files
# In order to create those blast.out files we are going to use grep

# AF
input_file_AF=/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/scripts/AF_hits.230274.log
output_hits_AF=/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/outputs/AF_out/hits_files
logs=/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/logs

grep "^NODE" $input_file_AF > $output_hits_AF/blastAF.out 

mv $input_file_AF $logs

# TF

input_file_TF=/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/scripts/TF_hits.230367.log
output_hits_TF=/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/outputs/TF_out/hits_files

grep "^NODE" $input_file_TF > $output_hits_TF/blastTF.out
  
mv $input_file_TF $logs

