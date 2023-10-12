#!/bin/bash

#Version 1
#Script for to run an metagenomics analysis
#bash *.sh folder_raw_data (without diagonal)
work_directory=$(pwd)
echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"
echo "||      Setting programs and folders     ||"
echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"
raw=$1
[ -d results ] || mkdir results
mkdir -p results/rawdata
mkdir -p results/fastqc
mkdir -p results/trim/
mkdir -p results/assembly
mkdir -p results/binning
mkdir -p results/asstax
#Links simbolicos
cd ${work_directory}/results/rawdata
find ${work_directory}/${raw}/*/ -name "*fastq*" -exec ln -s {} . ';'
cd ${work_directory}
echo "Finished!"

echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"
echo "||        Assessing Read Quality         ||"
echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"

cd ${work_directory}/results/fastqc
echo fastqc ${work_directory}/results/rawdata/*.fastq* -o .
echo "Finished!"

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
	--paired --retain_unpaired ${base}_R1.fastq.gz ${base}_R2.fastq.gz \
	--output_dir ${work_directory}results/trim/Trim_galore_${base}

echo mv *.fq.gz ${work_directory}/results/trim
done

echo "Finished!"
