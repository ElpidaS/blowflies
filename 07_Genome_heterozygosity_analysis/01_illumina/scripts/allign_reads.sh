#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o allign_reads.$JOB_ID.log

### Analysis of heterozygosity across a genome

## Steps:
# 1. Align reads to the genome
# 2. Process alignments
# 3. Call, genotype and filter variants
# 4. Calculate heterozygosity
# 5. Plot

# original Script; Rob Braid
# Edit; Elpida Skarlou 
# date of Edit 18/5/2023

# Use Rob's enviroment (up to step 4) as the packages downloading proccess is reall time consuming
# in order to be able to use conda
source /ceph/users/eskarlou/miniconda3/etc/profile.d/conda.sh
conda activate interspecific_alignments


############################
########## STEP 1 ##########
############################


# 1. Align reads to the genome

# You can use any aligner really here - e.g. bowtie2 or bwa. I find bwa to have better mapping rates (but there's probably a trade-off with accuracy), 
# but bowtie2 is a bit faster. For the purpose of this analysis I would say it doesn't really matter. 
# If the files are large, maybe you would prefer bowtie2 so you're not waiting days for the job to run. Make sure to sort the bam file with samtools 
# afterwards.

set -e
SCRATCH=/scratch/$USER/$JOB_ID/allign_reads
GENOME_TF=/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/outputs/TF_out/partition/
GENOME_AF=/data/ross/flies/analyses/blowflies/04_genome_partition/01_illumina/outputs/AF_out/partition/
TRIMMED=/data/ross/flies/raw/Chrysomya_rufifacies/illumina/
OUT_ALLIGN=/data/ross/flies/analyses/blowflies/07_Genome_heterozygosity_analysis/01_illumina/outputs/01_step_align
HETERO=/data/ross/flies/analyses/blowflies/07_Genome_heterozygosity_analysis/01_illumina

mkdir -p $SCRATCH
cd $SCRATCH

### sync GENOMES ###
# We are using the genomes that have been "cleaned" from any contamination + include contigs > ~300 bp
# The TF Genome has been created using two libraries   

rsync -av  $GENOME_AF/AF_contig.filtered.fasta $SCRATCH
rsync -av  $GENOME_TF/TF_contig.filtered.fasta $SCRATCH

### ssync trimmes reads ###
rsync -av $TRIMMED/AF* $SCRATCH
# I am picking the TF11 instead of the TF19 for no particular reason I expect that there is going to be no difference :)
rsync -av $TRIMMED/TF11* $SCRATCH


##### START ALLINGING #####

# I am using the botwie2 because of time and the fact that it is speciallized for illumina data

# Array of sample prefixes and corresponding read files
sample_prefixes=("AF7" "TF11")
read_files=("AF7_Chrysomya-rufifacies_S2" "TF11_Chrysomya-rufifacies_S3")

# Create Bowtie2 index for the genome assembly (if not already done)
bowtie2-build AF_contig.filtered.fasta AF_genome_index
bowtie2-build TF_contig.filtered.fasta TF_genome_index

# Loop through the samples
for ((i=0; i<${#sample_prefixes[@]}; i++))
do
    # Align the reads to the indexed genome
    bowtie2 -x AF_genome_index -1 ${read_files[i]}_R1_001.trimmed.fastq.gz -2 ${read_files[i]}_R2_001.trimmed.fastq.gz -S ${sample_prefixes[i]}_aligned.sam

    # Convert SAM to BAM
    samtools view -b -o ${sample_prefixes[i]}mappings.bam ${sample_prefixes[i]}_aligned.sam

    # Sort the BAM file
    samtools sort -o ${sample_prefixes[i]}mapping.sorted.bam ${sample_prefixes[i]}mappings.bam
done

# log file
mv $HETERO/scripts/rna_download.$JOB_ID.log $HETERO/logs

# rsync the output files 
rsync -av * $OUT_ALLIGN

rm -rf $SCRATCH