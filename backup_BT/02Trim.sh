#!/bin/bash

#Version 1
#Script for to run an metagenomics analysis
#bash *.sh
#For to run this script is necessary activate the enviroment that content the trimgalore program: $ conda activate LIM-abel
work_directory=$(pwd)

echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"
echo "||        Trimming and Filtering         ||"
echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"

cd ${work_directory}/results/rawdata
for infile in *_R1.fastq.gz
do
base=$(basename ${infile} _R1.fastq.gz)
echo "Trimming for file: $base"

echo trim_galore --fastqc --clip_R1 15 --clip_R2 15 \
	--three_prime_clip_R1 5  --three_prime_clip_R2 5 \
	--paired --retain_unpaired ${work_directory}/results/rawdata/${base}_R1.fastq.gz ${work_directory}/results/rawdata/${base}_R2.fastq.gz \
	--output_dir ${work_directory}/results/trim/Trim_galore_${base}

echo mv ${work_directory}/results/trim/*${base}*/*val*.fq.gz ${work_directory}/results/trim
done

echo "Finished!"
