#PBS -q ensam
#PBS -V
#PBS -N 04Binning
#PBS -e /LUSTRE/usuario/jlovaco/LIM/outs/04error.txt
#PBS -o /LUSTRE/usuario/jlovaco/LIM/outs/04out.txt
#PBS -l nodes=1:ppn=8,vmem=150gb,walltime=480:00:00

##Version 1
##Script for to run an metagenomics analysis
##bash *.sh 

module load MaxBin/2.2.6

echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"
echo "||        Metagenome Binning             ||"
echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"

cd ${PBS_O_WORKDIR}/results/rawdata

for infile in *_R1.fastq.gz
do
	base=$(basename ${infile} _R1.fastq.gz)
	echo $base
	run_MaxBin.pl -thread 8 -contig ${PBS_O_WORKDIR}/results/assembly/${base}_assem/${base}_contigs.fasta \
	       	-reads ${PBS_O_WORKDIR}/results/trim/${base}_R1_val_1.fq.gz \
		-reads2 ${PBS_O_WORKDIR}/results/trim/${base}_R2_val_2.fq.gz \
		-out ${PBS_O_WORKDIR}/results/binning/${base}_binn


done
