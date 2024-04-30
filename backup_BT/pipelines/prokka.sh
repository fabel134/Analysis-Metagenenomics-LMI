#Run prokka
# Prokka  for  metagenomic's assamblies
inPath= "5_Assamblies" #this is the folder where the files will be analyzed
outPath= "7_Prokka" #the foler's name that will contain the outcomes files

#Define input files
inFiles="megahit_LP10_51K_contigs.fa"

#We start the cycle for input files
#for file in ${inFiles}; do
echo #"Start the cicle"
echo  "Original file"
#baseFile='basename $file'
echo "Real file $baseFile"
#Sintaxys to modify part of the content of a variable
#baseName=${baseFile/.prokka/}
echo "minm name": $baseName
#outFile="$baseName.fa"
echo $outFile

echo "The final name of the file $outFile"

#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP10_51K_contigs.fa --outdir 7_Prokka/megahit_LP10_51K  5_Assamblies/megahit_LP10_51K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP1_51K_contigs.fa --outdir 7_Prokka/megahit_LP1_51K  5_Assamblies/megahit_LP1_51K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP1_77K_contigs.fa --outdir 7_Prokka/megahit_LP1_77K  5_Assamblies/megahit_LP1_77K_contigs.fa

#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP3_51K_contigs.fa --outdir 7_Prokka/megahit_LP3_51K  5_Assamblies/megahit_LP3_51K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP3_55K.fa --outdir 7_Prokka/megahit_LP3_55K  5_Assamblies/megahit_LP5_55K.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP3_77K_contigs.fa --outdir 7_Prokka/megahit_LP3_77K  5_Assamblies/megahit_LP3_77K_contigs.fa

#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP4_51K_contigs.fa --outdir 7_Prokka/megahit_LP4_51K  5_Assamblies/megahit_LP4_51K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP4_77K_contigs.fa --outdir 7_Prokka/megahit_LP4_77K  5_Assamblies/megahit_LP4_77K_contigs.fa

#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP5_51K_contigs.fa --outdir 7_Prokka/megahit_LP5_51K  5_Assamblies/megahit_LP5_51K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP5_77K_contigs.fa --outdir 7_Prokka/megahit_LP5_77K  5_Assamblies/megahit_LP5_77K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP7_51K_contigs.fa --outdir 7_Prokka/megahit_LP7_51K  5_Assamblies/megahit_LP7_51K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP7_77K_contigs.fa --outdir 7_Prokka/megahit_LP7_77K  5_Assamblies/megahit_LP7_77K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP8_51K_contigs.fa --outdir 7_Prokka/megahit_LP8_51K  5_Assamblies/megahit_LP8_51K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP8_77K_contigs.fa --outdir 7_Prokka/megahit_LP8_77K  5_Assamblies/megahit_LP8_77K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP9_51K_contigs.fa --outdir 7_Prokka/megahit_LP9_51K  5_Assamblies/megahit_LP9_51K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP9_77K_contigs.fa --outdir 7_Prokka/megahit_LP9_77K  5_Assamblies/megahit_LP9_77K_contigs.fa

#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP10_51K_contigs.fa --outdir 7_Prokka/megahit_LP10_51K  5_Assamblies/megahit_LP10_51K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP10_77K_contigs.fa --outdir 7_Prokka/megahit_LP10_77K  5_Assamblies/megahit_LP10_77K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP11_51K_contigs.fa --outdir 7_Prokka/megahit_LP11_51K  5_Assamblies/megahit_LP11_51K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP11_77K_contigs.fa --outdir 7_Prokka/megahit_LP11_77K  5_Assamblies/megahit_LP11_77K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP12_51K_contigs.fa --outdir 7_Prokka/megahit_LP12_51K  5_Assamblies/megahit_LP12_51K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP12_77K_contigs.fa --outdir 7_Prokka/megahit_LP12_77K  5_Assamblies/megahit_LP12_77K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP13_51K_contigs.fa --outdir 7_Prokka/megahit_LP13_51K  5_Assamblies/megahit_LP13_51K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP13_77K_contigs.fa --outdir 7_Prokka/megahit_LP13_77K  5_Assamblies/megahit_LP13_77K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP13_777K_contigs.fa --outdir 7_Prokka/megahit_LP13_77K  5_Assamblies/megahit_LP13_777K_contigs.fa

#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP14_51K_contigs.fa --outdir 7_Prokka/megahit_LP14_51K  5_Assamblies/megahit_LP14_51K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP14_77K_contigs.fa --outdir 7_Prokka/megahit_LP14_77K  5_Assamblies/megahit_LP14_77K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP14_777K_contigs.fa --outdir 7_Prokka/megahit_LP14_777K  5_Assamblies/megahit_LP14_777K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP15_51K_contigs.fa --outdir 7_Prokka/megahit_LP15_51K  5_Assamblies/megahit_LP15_51K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP15_77K_contigs.fa --outdir 7_Prokka/megahit_LP15_77K  5_Assamblies/megahit_LP15_77K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP15_777K_contigs.fa --outdir 7_Prokka/megahit_LP15_777K  5_Assamblies/megahit_LP15_777K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP16_51K_contigs.fa --outdir 7_Prokka/megahit_LP16_51K  5_Assamblies/megahit_LP16_51K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP16_77K_contigs.fa --outdir 7_Prokka/megahit_LP16_77K  5_Assamblies/megahit_LP16_77K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP16_777K_contigs.fa --outdir 7_Prokka/megahit_LP16_777K  5_Assamblies/megahit_LP16_777K_contigs.fa

#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP17_51K_contigs.fa --outdir 7_Prokka/megahit_LP17_51K  5_Assamblies/megahit_LP17_51K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP17_77K_contigs.fa --outdir 7_Prokka/megahit_LP17_77K  5_Assamblies/megahit_LP17_77K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP17_777K_contigs.fa --outdir 7_Prokka/megahit_LP17_777K  5_Assamblies/megahit_LP17_777K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP18_51K_contigs.fa --outdir 7_Prokka/megahit_LP18_51K  5_Assamblies/megahit_LP18_51K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP18_77K_contigs.fa --outdir 7_Prokka/megahit_LP18_77K  5_Assamblies/megahit_LP18_77K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP18_777K_contigs.fa --outdir 7_Prokka/megahit_LP18_777K  5_Assamblies/megahit_LP18_777K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP19_51K_contigs.fa --outdir 7_Prokka/megahit_LP19_51K  5_Assamblies/megahit_LP19_51K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP19_77K_contigs.fa --outdir 7_Prokka/megahit_LP19_77K  5_Assamblies/megahit_LP19_77K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP19_777K_contigs.fa --outdir 7_Prokka/megahit_LP19_777K  5_Assamblies/megahit_LP19_777K_contigs.fa

#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP20_51K_contigs.fa --outdir 7_Prokka/megahit_LP20_51K  5_Assamblies/megahit_LP20_51K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP20_77K_contigs.fa --outdir 7_Prokka/megahit_LP20_77K  5_Assamblies/megahit_LP20_77K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP20_777K_contigs.fa --outdir 7_Prokka/megahit_LP20_777K  5_Assamblies/megahit_LP20_777K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP22_51K_contigs.fa --outdir 7_Prokka/megahit_LP22_51K  5_Assamblies/megahit_LP22_51K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP23_51K_contigs.fa --outdir 7_Prokka/megahit_LP23_51K  5_Assamblies/megahit_LP23_51K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP23_77K_contigs.fa --outdir 7_Prokka/megahit_LP23_77K  5_Assamblies/megahit_LP23_77K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP24_51K_contigs.fa --outdir 7_Prokka/megahit_LP24_51K  5_Assamblies/megahit_LP24_51K_contigs.fa
#prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP24_77K_contigs.fa --outdir 7_Prokka/megahit_LP24_77K  5_Assamblies/megahit_LP24_77K_contigs.fa
prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP25_51K_contigs.fa --outdir 7_Prokka/megahit_LP25_51K  5_Assamblies/megahit_LP25_51K_contigs.fa

prokka --metagenome --prefix prokka  --cpus 6 --locustag megahit_LP25_77K_contigs.fa --outdir 7_Prokka/megahit_LP25_77K  5_Assamblies/megahit_LP25_77K_contigs.fa
prokka --metagenome --prefix prokka  --cpus 6 --locustag spades_LP9_51K.fasta --outdir 7_Prokka/spades_LP9_51K 5_Assamblies/spades_LP9_51K.fasta
prokka --metagenome --prefix prokka  --cpus 6 --locustag spades_LP10_51K.fasta --outdir 7_Prokka/spades_LP10_51K 5_Assamblies/spades_LP10_51K.fasta
prokka --metagenome --prefix prokka  --cpus 6 --locustag spades_LP12_51K.fasta --outdir 7_Prokka/spades_LP12_51K 5_Assamblies/spades_LP12_51K.fasta
prokka --metagenome --prefix prokka  --cpus 6 --locustag spades_LP14_51K.fasta --outdir 7_Prokka/spades_LP14_51K 5_Assamblies/spades_LP14_51K.fasta
prokka --metagenome --prefix prokka  --cpus 6 --locustag spades_LP15_51K.fasta --outdir 7_Prokka/spades_LP15_51K 5_Assamblies/spades_LP15_51K.fasta
prokka --metagenome --prefix prokka  --cpus 6 --locustag spades_LP17_51K.fasta --outdir 7_Prokka/spades_LP17_51K 5_Assamblies/spades_LP17_51K.fasta
prokka --metagenome --prefix prokka  --cpus 6 --locustag spades_LP18_51K.fasta --outdir 7_Prokka/spades_LP18_51K 5_Assamblies/spades_LP18_51K.fasta
prokka --metagenome --prefix prokka  --cpus 6 --locustag spades_LP19_51K.fasta --outdir 7_Prokka/spades_LP19_51K 5_Assamblies/spades_LP19_51K.fasta
prokka --metagenome --prefix prokka  --cpus 6 --locustag spades_LP20_51K.fasta --outdir 7_Prokka/spades_LP20_51K 5_Assamblies/spades_LP20_51K.fasta
prokka --metagenome --prefix prokka  --cpus 6 --locustag spades_LP22_51K.fasta --outdir 7_Prokka/spades_LP22_51K 5_Assamblies/spades_LP22_51K.fasta
prokka --metagenome --prefix prokka  --cpus 6 --locustag spades_LP24_51K.fasta --outdir 7_Prokka/spades_LP24_51K 5_Assamblies/spades_LP24_51K.fasta
prokka --metagenome --prefix prokka  --cpus 6 --locustag spades_LP25_51K.fasta --outdir 7_Prokka/spades_LP25_51K 5_Assamblies/spades_LP25_51K.fasta




























































