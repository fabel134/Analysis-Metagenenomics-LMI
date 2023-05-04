#!/bin/bash

#Version 1
#Script for to run an metagenomics analysis
#bash *.sh 
work_directory=$(pwd)

echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"
echo "||        Metagenome Binning             ||"
echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"

cd ${work_directory}/results/rawdata

for infile in *_R1.fastq.gz
do
	base=$(basename ${infile} _R1.fastq.gz)
	echo $base
	echo run_MaxBin.pl -thread 8 -contig ${work_directory}/results/assembly/${base}_assem/${base}_contigs.fasta \
	       	-reads ${work_directory}/results/trim/${base}_R1_val_1.fq.gz \
		-reads2 ${work_directory}/results/trim/${base}_R2_val_2.fq.gz \
		-out ${work_directory}/results/binning/${base}_binn


done
