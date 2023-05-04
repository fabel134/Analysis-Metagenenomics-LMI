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

