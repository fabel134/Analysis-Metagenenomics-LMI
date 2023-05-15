#!/bin/bash

#Version 1
#Script for to run an metagenomics analysis
#bash *.sh 
work_directory=$(pwd)

echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"
echo "||        Taxonomic Assignment           ||"
echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"

cd ${work_directory}/results/rawdata

for infile in *_R1.fastq.gz
do
	base=$(basename ${infile} _R1.fastq.gz)
	echo $base

	echo "Taxonomic Assignment for contigs"
	mkdir ${work_directory}/results/asstax/contigs
	echo kraken2 --db /datos/metamex/krakendb --threads 8 \
		--output ${work_directory}/results/asstax/contigs/${base}_contigs.kraken \
	       --report ${work_directory}/results/asstax/contigs/${base}_contigs.report	\
	       ${work_directory}/results/binning/${base}*.fasta


	echo cut -f2,3 ${work_directory}/results/asstax/contigs/${base}_contigs.kraken \> ${work_directory}/results/asstax/contigs/${base}_contigs.krona.input
	
	#echo cut -f2,3 ${work_directory}/results/asstax/contigs/${base}_contigs.kraken > ${work_directory}/results/asstax/contigs/${base}_contigs.krona.input

	echo "***ktImportTaxonomy***"
	echo ktImportTaxonomy ${work_directory}/results/asstax/contigs/${base}_contigs.krona.input -o ${work_directory}/results/asstax/contigs/${base}_contigs.krona.out.html

done
