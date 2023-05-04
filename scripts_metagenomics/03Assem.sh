#!/bin/bash

#Version 1
#Script for to run an metagenomics analysis
#bash *.sh 
work_directory=$(pwd)

echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"
echo "||        Metagenome Assembly            ||"
echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"

cd ${work_directory}/results/trim
for infile in *_R1_val_1.fq.gz
do
base=$(basename ${infile} _R1_val_1.fq.gz)
echo "Assembly for file: $base"

echo mkdir ${work_directory}/results/assembly/${base}_assem

echo metaspades.py -1 ${base}_R1_val_1.fq.gz -2 ${base}_R2_val_2.fq.gz -o ${work_directory}/results/assembly/${base}_assem

cd ${work_directory}/results/assembly/${base}_assem
for name in *; 
	do 
		mv $name ${base}_$name; 
	done 

done

echo "Finished!"
