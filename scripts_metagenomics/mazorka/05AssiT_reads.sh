#PBS -q ensam
#PBS -V
#PBS -N 05kraken-reads
#PBS -e /LUSTRE/usuario/jlovaco/LIM/outs/05READSerror.txt
#PBS -o /LUSTRE/usuario/jlovaco/LIM/outs/05READSout.txt
#PBS -l nodes=1:ppn=20,vmem=300gb,walltime=480:00:00

module load kraken/2.1.2


echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"
echo "||        Taxonomic Assignment           ||"
echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"

cd ${PBS_O_WORKDIR}/results/rawdata
mkdir ${PBS_O_WORKDIR}/results/asstax/reads

for infile in *_R1.fastq.gz
do
	base=$(basename ${infile} _R1.fastq.gz)
	
	echo "Taxonomic Assignment for reads:"$base
	kraken2 --db KRAKEN --threads 12 \
		--paired \
	       	${PBS_O_WORKDIR}/results/trim/${base}_R1_val_1.fq.gz \
		${PBS_O_WORKDIR}/results/trim/${base}_R2_val_2.fq.gz \
		--output ${PBS_O_WORKDIR}/results/asstax/reads/${base}_reads.kraken \
		--report ${PBS_O_WORKDIR}/results/asstax/reads/${base}_reads.report

	#echo cut -f2,3 ${PBS_O_WORKDIR}/results/asstax/reads/${base}_reads.kraken \> ${PBS_O_WORKDIR}/results/asstax/reads/${base}_reads.krona.input
	
	echo cut -f2,3 ${PBS_O_WORKDIR}/results/asstax/reads/${base}_reads.kraken > ${PBS_O_WORKDIR}/results/asstax/reads/${base}_reads.krona.input

	echo "***ktImportTaxonomy***"
	echo ktImportTaxonomy ${PBS_O_WORKDIR}/results/asstax/contigs/${base}_contigs.krona.input -o ${PBS_O_WORKDIR}/results/asstax/contigs/${base}_contigs.krona.out.html

done
