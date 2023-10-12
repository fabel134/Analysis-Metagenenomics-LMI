#!/bin/bash

#Version 1
#Script for to run an metagenomics analysis
#bash *.sh 
work_directory=$(pwd)

echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"
echo "||        Taxonomic Assignment           ||"
echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"

cd ${work_directory}/results/trim

mkdir ${work_directory}/results/asstax/contigs

for infile in *_R1_val_1.fq.gz

do
	base=$(basename ${infile} _R1_val_1.fq.gz)
	echo $base

	echo "Taxonomic Assignment for contigs"
	kraken2 --db /data/camda2023/krakenDB --threads 18 \
		--output ${work_directory}/results/asstax/contigs/${base}_contigs.kraken \
	       --report ${work_directory}/results/asstax/contigs/${base}_contigs.report	\
	       ${work_directory}/results/assembly/${base}_contigs.fasta


	#echo cut -f2,3 ${work_directory}/results/asstax/contigs/${base}_contigs.kraken \> ${work_directory}/results/asstax/contigs/${base}_contigs.krona.input
	
	echo cut -f2,3 ${work_directory}/results/asstax/contigs/${base}_contigs.kraken > ${work_directory}/results/asstax/contigs/${base}_contigs.krona.input

	echo "***ktImportTaxonomy***"
	echo ktImportTaxonomy ${work_directory}/results/asstax/contigs/${base}_contigs.krona.input -o ${work_directory}/results/asstax/contigs/${base}_contigs.krona.out.html

done
