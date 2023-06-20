#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -o cdstofa.$JOB_ID.log


CDS_dir=/data/ross/flies/analyses/blowflies/06_Transcriptomic_analysis/outputs/transdecoder_filter

awk '/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $CDS_dir/Trinity.fasta.transdecoder.cds >  $CDS_dir/transcripts.fa

