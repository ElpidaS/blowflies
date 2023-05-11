#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o TF_prepart.$JOB_ID.log



### Description 
# Execute Workflow A in https://blobtools.readme.io/docs/my-first-blobplot 
# Pre-Partition #
# Up Until step2 aka, the part that we are able to make a descission regarding 
# partitioning the data or not


#### 1 Create BlobDb #####

# Input; Hit files (TSV) + 
#        Assembly (FASTA) + 
#        Coverage files (BAM)

# Output; BlobDB (e.g. *.blobDB.json)


#### 2a Create BlobPlot + CovPlot #####

# Make partition descision #

# Input; BlobDB

# Output; BlobPlot + CovPlot


#### 2b Create Table #####

# Input; BlobDB

# Output; Table


# Author; Elpida Skarlou
# Date; 

### Precode stuff

set -e # exit immediately on error

# PATHS used in this script
SCRATCH=/scratch/$USER/$JOB_ID
GENOME_PART=/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/
TF_ASSEMBLY=/data/ross/flies/analyses/blowflies/02_genome_assembly/01_spades/outputs/TF_assembly/
TF_COV=/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/outputs/TF_out/coverage_file/
TF_HIT=/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/outputs/TF_out/hits_files/
TF_PART=/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/outputs/TF_out/pre_partition/
NCBI_TAX=/data/ross/flies/analyses/blowflies/05_NCBI_taxdump


# Needed steps to be able to use blobtools
# for the dependencies of blobtools 
source /ceph/users/eskarlou/miniconda3/etc/profile.d/conda.sh
conda activate blobtools_depend

# rsync ~/blobtools to $SCRATCH
rsync -av /ceph/users/eskarlou/blobtools/* $SCRATCH

##########################
#### 1 Create BlobDb #####
##########################

# rsyn input files #

# Assembly
rsync -av $TF_ASSEMBLY/contigs.fasta $SCRATCH
# coverage
rsync -av $TF_COV/TFmapping.sorted.bam.bai $SCRATCH
rsync -av $TF_COV/TFmapping.sorted.bam $SCRATCH # should include it so .bam.bai works 
# hits
rsync -av $TF_HIT/blastTF.out $SCRATCH
# --names
rsync -av $NCBI_TAX/names.dmp $SCRATCH
# --nodes
rsync -av $NCBI_TAX/nodes.dmp $SCRATCH

cd $SCRATCH

# using the files provided in the test_files/ folder 
./blobtools create -i contigs.fasta -b $SCRATCH/TFmapping.sorted.bam -t blastTF.out -o TF_blobplot --names names.dmp --nodes nodes.dmp

#######################################
#### 2a Create BlobPlot + CovPlot #####
#######################################

./blobtools plot -i TF_blobplot.blobDB.json -o ./

##########################
#### 2b Create Table #####
##########################

./blobtools view -i TF_blobplot.blobDB.json -o ./


### Post main code ###

# remove all the input files 
rm -rf $SCRATCH/contigs.fasta
rm -rf $SCRATCH/TFmapping.sorted.bam
rm -rf $SCRATCH/TFmapping.sorted.bam.bai
rm -rf $SCRATCH/blastTF.out
rm -rf $SCRATCH/blobtools

# move log file 
mv $GENOME_PART/scripts/TF/TF_prepart.$JOB_ID.log $GENOME_PART/logs

# rsync all the outputs 
rsync -av $SCRATCH/* $TF_PART

rm -rf *