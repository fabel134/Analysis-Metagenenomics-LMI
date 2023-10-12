#PBS -q ensam
#PBS -V
#PBS -N 05kraken
#PBS -e /LUSTRE/usuario/jlovaco/LIM/outs/05error.txt
#PBS -o /LUSTRE/usuario/jlovaco/LIM/outs/05out.txt
#PBS -l nodes=1:ppn=20,vmem=300gb,walltime=480:00:00

module load kraken/2.1.2

echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"
echo "||        Taxonomic Assignment           ||"
echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"

cd ${PBS_O_WORKDIR}/results/rawdata

for infile in *_R1.fastq.gz
do
	base=$(basename ${infile} _R1.fastq.gz)
	echo $base

	echo "Taxonomic Assignment for contigs"
	mkdir ${PBS_O_WORKDIR}/results/asstax/contigs
	kraken2 --db KRAKEN  --threads 12 \
		--output ${PBS_O_WORKDIR}/results/asstax/contigs/${base}_contigs.kraken \
	       --report ${PBS_O_WORKDIR}/results/asstax/contigs/${base}_contigs.report	\
	       ${PBS_O_WORKDIR}/results/binning/${base}*.fasta


	cut -f2,3 ${PBS_O_WORKDIR}/results/asstax/contigs/${base}_contigs.kraken > ${PBS_O_WORKDIR}/results/asstax/contigs/${base}_contigs.krona.input
	
	#echo cut -f2,3 ${PBS_O_WORKDIR}/results/asstax/contigs/${base}_contigs.kraken > ${PBS_O_WORKDIR}/results/asstax/contigs/${base}_contigs.krona.input

	echo "***ktImportTaxonomy***"
	echo ktImportTaxonomy ${PBS_O_WORKDIR}/results/asstax/contigs/${base}_contigs.krona.input -o ${PBS_O_WORKDIR}/results/asstax/contigs/${base}_contigs.krona.out.html

done
