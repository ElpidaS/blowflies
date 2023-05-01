#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o mv.$JOB_ID.log

rm -rf /data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/outputs/AF_out/* 

mv /scratch/eskarlou/204027/SPAdes/coverage_file /data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/outputs/AF_out/

rm -rf /scratch/eskarlou/204027