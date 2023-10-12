#!/bin/bash

#Version 1
#Script for to run an metagenomics analysis
#bash *.sh 
work_directory=$(pwd)

echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"
echo "||        Taxonomic Assignment           ||"
echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"

cd ${work_directory}/results/trim

for infile in *_R1_val_1.fq.gz
do
	base=$(basename ${infile} _R1_val_1.fq.gz)
	echo $base
	echo "Taxonomic Assignment for reads"
	mkdir ${work_directory}/results/asstax/reads
	kraken2 --db /data/camda2023/krakenDB --threads 16 \
		--paired \
	       	${work_directory}/results/trim/${base}_R1_val_1.fq.gz \
		${work_directory}/results/trim/${base}_R2_val_2.fq.gz \
		--output ${work_directory}/results/asstax/reads \
		--report ${work_directory}/results/asstax/reads/${base}_reads.report

	cut -f2,3 ${work_directory}/results/asstax/reads/${base}_reads.kraken > ${work_directory}/results/asstax/reads/${base}_reads.krona.input
	
	#echo cut -f2,3 ${work_directory}/results/asstax/reads/${base}_reads.kraken \> ${work_directory}/results/asstax/reads/${base}_reads.krona.input

	echo "***ktImportTaxonomy***"
	ktImportTaxonomy ${work_directory}/results/asstax/reads/${base}_reads.krona.input -o ${work_directory}/results/asstax/reads/${base}_reads.krona.out.html

done
