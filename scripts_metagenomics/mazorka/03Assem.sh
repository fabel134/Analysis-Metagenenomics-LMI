#PBS -q ensam
#PBS -V
#PBS -N 03Assembly
#PBS -e /LUSTRE/usuario/jlovaco/LIM/outs/03error.txt
#PBS -o /LUSTRE/usuario/jlovaco/LIM/outs/03out.txt
#PBS -l nodes=1:ppn=20,vmem=200gb,walltime=480:00:00


module load SPAdes/3.15.2

echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"
echo "||        Metagenome Assembly            ||"
echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"

cd ${PBS_O_WORKDIR}/results/trim
#for infile in LP23_R1_val_1.fq.gz LP21_R1_val_1.fq.gz
for infile in *_R1_val_1.fq.gz
do
base=$(basename ${infile} _R1_val_1.fq.gz)
echo "Assembly for file: $base"

mkdir ${PBS_O_WORKDIR}/results/assembly/${base}_assem

metaspades.py -1 ${PBS_O_WORKDIR}/results/trim/${base}_R1_val_1.fq.gz -2 ${PBS_O_WORKDIR}/results/trim/${base}_R2_val_2.fq.gz -o ${PBS_O_WORKDIR}/results/assembly/${base}_assem

cd ${PBS_O_WORKDIR}/results/assembly/${base}_assem
for name in *; 
	do 
	echo mv $name ${base}_$name; 
	mv $name ${base}_$name; 
	done 

done

echo "Finished!"
