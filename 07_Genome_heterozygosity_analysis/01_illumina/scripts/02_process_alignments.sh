#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o process_reads.$JOB_ID.log

# 2. Process alignments

# Use Rob's enviroment (up to step 4) as the packages downloading proccess is reall time consuming
# in order to be able to use conda
source /ceph/users/eskarlou/miniconda3/etc/profile.d/conda.sh
conda activate interspecific_alignments

# Here, we follow the GATK-4 best practices pipeline to pre-process the alignment files so they are ready for variant calling:
# https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows. 
# This involves use of picardtools.

set -e
SCRATCH=/scratch/$USER/$JOB_ID/process_alignments
GENOME_TF=/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/outputs/TF_out/partition/
GENOME_AF=/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/outputs/AF_out/partition/
STEP01=/data/ross/flies/analyses/blowflies/07_Genome_heterozygosity_analysis/01_illumina/outputs/01_step_align
HETERO=/data/ross/flies/analyses/blowflies/07_Genome_heterozygosity_analysis/01_illumina
OUT_PROC=/data/ross/flies/analyses/blowflies/07_Genome_heterozygosity_analysis/01_illumina/outputs

mkdir -p $SCRATCH
cd $SCRATCH

### SYNC GENOMES + sorted.bam ###

rsync -av  $GENOME_AF/AF_contig.filtered.fasta $SCRATCH
rsync -av  $GENOME_TF/TF_contig.filtered.fasta $SCRATCH

# should be two files 
rsync -av  $STEP01/*mapping.sorted.bam $SCRATCH 

#############################

# The samtools faidx command creates an index file for the reference genome in the .fai format. 
# This index file allows for fast random access to different regions of the genome, which is essential for many downstream analyses, 
# including variant calling.
samtools faidx AF_contig.filtered.fasta
samtools faidx TF_contig.filtered.fasta


# Define the sample names and corresponding input files
samples=("AF7" "TF11")
input_files=("AF7mapping.sorted.bam" "TF11mapping.sorted.bam")


# Iterate over each sample
for ((i=0; i<${#samples[@]}; i++))
do
    sample="${samples[$i]}"
    input_file="${input_files[$i]}"
    
    # Sort the bam file with picardtools (different to samtools sort)
    picard SortSam I="${input_file}" SORT_ORDER=coordinate O="${sample}.${sample}_contig_filtered.pic.sorted.bam" 
    
    # Mark and remove duplicates -- this is important to mitigate biases introduced by data generation steps such as PCR amplification.
    samtools index "${sample}.${sample}_contig_filtered.pic.sorted.bam" 
    
    java -Xmx20g -jar /ceph/users/eskarlou/software/picard/build/libs/picard.jar MarkDuplicates MAX_RECORDS_IN_RAM=400000 INPUT="${sample}.${sample}_contig_filtered.pic.sorted.bam" OUTPUT="${sample}.${sample}_contig_filtered.pic.rmdup.sorted.bam" M="${sample}.met.txt" REMOVE_DUPLICATES=true 
    
    # Assign unique read groups
    # Assign unique read groups - it is important your alignment has unique readgroups 
	# (RGID and RGSM should be different for each sample - this is not so important if you're just calling variants on one sample, 
	# but it does matter with multiple samples). One way to do this is to name RGID and RGSM as ${base}

    picard AddOrReplaceReadGroups INPUT="${sample}.${sample}_contig_filtered.pic.rmdup.sorted.bam" \
    OUTPUT="${sample}.${sample}_contig_filtered.pic.rmdup.rg.sorted.bam" \
    RGID="${sample}" \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=unit1 \
    RGSM="${sample}"
done


# After the data pre-processing, you can merge your bam files and call variants on the single merged file (assuming you have assigned distinct read group IDs to each sample). 
# Another option is to call variants on each sample separately via a loop - depends what you want to do with the outputs. 
# GATK-4 can take a few days to run on large alignment files, so I often submit variant calling jobs on different samples separately to make things faster.


# log file
mv $HETERO/scripts/process_reads.$JOB_ID.log $HETERO/logs

# rsync the output files 
rsync -av *_contig_filtered.pic.rmdup.rg.sorted.bam $OUT_PROC

# Going forward, you need the *.pic.rmdup.rg.sorted.bam file. Note: it can be useful to keep intermediate files in your scratch directory temporarily for debug purposes.