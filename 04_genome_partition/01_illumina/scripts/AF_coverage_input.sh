#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o AF_coverage.$JOB_ID.log

### Description 
# Main goal; 
# COVERAGE FILES
# Produce oen of the 2 main files (coverage and hits) needed to run Blobtools

# More specifically; 
# Based on the script, it looks like the goal is 
# to prepare two types of files needed for running the software BlobTools: 
# a coverage file and a hits file. The coverage file contains information 
# about read coverage for each contig in the genome assembly, while the hits 
# file contains information about homology to other sequences in a reference database.

# Author; Elpida Skarlou
# Date; 26/04/2023

### Precode stuff

set -e # exit immediately on error

# PATHS used in this script
SCRATCH=/scratch/$USER/$JOB_ID/Coverage
AF_ASSEMBLY=/data/ross/flies/analyses/blowflies/02_genome_assembly/01_spades/outputs/AF_assembly/
TRIMMED_READS=/data/ross/flies/raw/Chrysomya_rufifacies/illumina/
GENOME_PART=/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/

# make $SCRATCH if it doesn't exist
mkdir -p $SCRATCH 
mkdir $SCRATCH/coverage_file

# Some extra PATHS
COVERAGE=$SCRATCH/coverage_file


# rsync usefull input to working directory 
rsync -av $AF_ASSEMBLY/contigs.fasta $COVERAGE
rsync -av $TRIMMED_READS/AF* $COVERAGE

### TWO TYPES OF FILES NEEDED FOR THE BLOBTOOLS TO WORK ###

### one (or more) coverage file(s) e.g. example/mapping_1.bam
	
	# activate conda
	source /ceph/users/eskarlou/miniconda3/etc/profile.d/conda.sh
	conda activate bowtie2

	cd $COVERAGE

	## 1st Index the asssembly genome (contig.fasta) by using Bowtie2 
	# will output multiple assembly_index files

	bowtie2-build contigs.fasta AFcontigs_index

	## 2nd The second command aligns the reads in "reads.fastq" to the genome 
	# assembly index in "assembly_index", and outputs the resulting alignments 
	# in SAM format to "mappings.sam". 

	bowtie2 -x AFcontigs_index \
	-1 AF7_Chrysomya-rufifacies_S2_R1_001.trimmed.fastq.gz \
	-2 AF7_Chrysomya-rufifacies_S2_R2_001.trimmed.fastq.gz \
	-S AFmappings.sam

	# 3rd SAM coverted to BAM 
	conda activate samtools

	samtools view -b -o AFmappings.bam AFmappings.sam

	# 4th & 5th fourth and fifth commands sort 
	# and index the BAM file to create the coverage file "mapping.sorted.bam".

	samtools sort -o AFmapping.sorted.bam AFmappings.bam
	samtools index AFmapping.sorted.bam


# Post-code stuff

# remove all the input files 
rm -rf $COVERAGE/contigs.fasta
rm -rf $COVERAGE/*fastq.gz

# log file
mv -av $GENOME_PART/scripts/AF_coverage.$JOB_ID.log $GENOME_PART/logs

# rsync all the outputs 
rsync -av $COVERAGE $GENOME_PART/outputs/AF_out/

rm -r *
rm -rf /scratch/$USER/$JOB_ID