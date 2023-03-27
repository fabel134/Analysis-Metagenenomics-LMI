#!/bin/bash

#Version 1
#Script for to run an metagenomics analysis
#bash *.sh folder_raw_data (without diagonal)
echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"
echo "||      Setting programs and folders     ||"
echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"

work_directory=$(pwd)
raw=$1
[ -d results ] || mkdir results
mkdir -p results/rawdata
mkdir -p results/fastqc
mkdir -p results/fastqc/trim
mkdir -p results/trimming
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

trimmomatic PE ${infile} ${base}_R2.fastq.gz \
${work_directory}/results/trimming/${base}_R1.trim.fastq.gz ${work_directory}/results/trimming/${base}_R1un.trim.fastq.gz \
${work_directory}/results/trimming/${base}_R2.trim.fastq.gz ${work_directory}/results/trimming/${base}_R2un.trim.fastq.gz \
ILLUMINACLIP:${work_directory}/TruSeq3-PE.fa:1:30:1 LEADING:25 TRAILING:25 SLIDINGWINDOW:4:20 MINLEN:35

done
#FastqC again
#exit
cd ${work_directory}/results/fastqc/trim
fastqc ${work_directory}/results/trimming/*R1.trim.fastq* -o .
fastqc ${work_directory}/results/trimming/*R2.trim.fastq* -o .
echo "Finished!"
