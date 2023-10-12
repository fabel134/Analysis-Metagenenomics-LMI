#PBS -N 02Trim
#PBS -q default
#PBS -l nodes=1:ppn=8,walltime=90:00:00,mem=8g,vmem=8g
#PBS -e /LUSTRE/usuario/jlovaco/LIM/outs/02error.txt
#PBS -o /LUSTRE/usuario/jlovaco/LIM/outs/02out.txt
#PBS -V

module load trim_galore/0.6.10
module load cutadapt/2.10

cd $PBS_O_WORKDIR

echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"
echo "||        Trimming and Filtering         ||"
echo "+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+:+"

#echo ${PBS_O_WORKDIR}/results
#echo $PBS_O_WORKDIR/results
#exit
cd ${PBS_O_WORKDIR}/results/rawdata
for infile in *_R1.fastq.gz
do
base=$(basename ${infile} _R1.fastq.gz)
echo "Trimming for file: $base"

trim_galore --fastqc --clip_R1 15 --clip_R2 15 \
	--three_prime_clip_R1 5  --three_prime_clip_R2 5 \
	--paired --retain_unpaired ${PBS_O_WORKDIR}/results/rawdata/${base}_R1.fastq.gz ${PBS_O_WORKDIR}/results/rawdata/${base}_R2.fastq.gz \
	--output_dir ${PBS_O_WORKDIR}/results/trim/Trim_galore_${base}

mv ${PBS_O_WORKDIR}/results/trim/*${base}*/*val*.fq.gz ${PBS_O_WORKDIR}/results/trim
done

echo "Finished!"
