#!/bin/bash

#Version 1
#Script for to run an metagenomics analysis
#bash *.sh 
#remember that its necessary activate the environment 'metagenomics'
work_directory=$(pwd)

echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"
echo "||        Metagenome Assembly            ||"
echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"

cd ${work_directory}/results/trim
#for infile in LP23_R1_val_1.fq.gz LP21_R1_val_1.fq.gz
for infile in *_R1_val_1.fq.gz
do
base=$(basename ${infile} _R1_val_1.fq.gz)
echo "Assembly for file: $base"

mkdir ${work_directory}/results/assembly/${base}_assem

metaspades.py -1 ${work_directory}/results/trim/${base}_R1_val_1.fq.gz -2 ${work_directory}/results/trim/${base}_R2_val_2.fq.gz -o ${work_directory}/results/assembly/${base}_assem

cd ${work_directory}/results/assembly/${base}_assem
for name in *; 
	do 
	echo mv $name ${base}_$name; 
	mv $name ${base}_$name; 
	done 

done

echo "Finished!"
