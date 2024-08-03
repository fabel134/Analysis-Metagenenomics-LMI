#!/bin/bash
#Version 1.0 Modified from Flores-Nuñez et al., 2023
#Script must be run in the main_dir were data is storaged
#main_dir should contain a data folder
#TRIMMOMATIC, FLASH AND RSTUDIO ARE NEEDED

echo ###### ============================ ######
echo ###### Setting programs and folders ######
echo ###### ============================ ######

#Remove out files so they dont concatenate the same result
rm -r ./stats 

#Seting directories for each step in the analysis
mkdir -p ./trimmed
mkdir -p ./merged
mkdir -p ./process
mkdir -p ./result
mkdir -p ./stats

#set variables #avoid $o, $f, $n
#Trimmomatic
q=12 #minimum quality 
l=30 #minimun length

#Flash 
m=10 #minimum overlap
M=394 #maximum overlap

#Filtering #more can be added if needed
sr=18 #strip bases from right 
sl=39 #strip bases from left
mi=296 #minumum length
ma=394 #maximum length
ee=1 #expected error

#clustering
i=0.95 #threshold OTU for clustering

#otu table
b=0.80 #minimum boostrap value


#Create variables for programs 
#export trimmomatic=/home/fer/abel/analysis/db/Google-Drive/Trimmomatic-0.39/trimmomatic-0.39.jar
#export flash=/home/fer/abel/analysis/db/Google-Drive/FLASH-1.2.11/flash
#export vsearch=/home/fer/abel/analysis/db/Google-Drive/vsearch-2.18.0-macos-x86_64/bin/vsearch


echo ###### ============================ ######
echo ###### Q trimming and merging pairs ######
echo ###### ============================ ######


#TRIMM low quality bases from raw read
for f in $(find ./data/* | grep '_001.fastq.gz'); do
o=./trimmed/$(awk -F '/|_L001' '{print $3}' <<< "$f").fastq.gz
trimmomatic PE -basein $f -baseout $o -summary ./trimmed/out.txt LEADING:$q TRAILING:$q SLIDINGWINDOW:4:$q MINLEN:$l AVGQUAL:$q

echo $o >> ./stats/trimmed.out.txt
cat ./trimmed/out.txt >> ./stats/trimmed.out.txt
done

#MERGE paired reads
for f in $(find ./trimmed/* -type f | grep _1P.fastq.gz) ; do 
o=./merged/$(awk -F '/|_1P' '{print $3}' <<< "$f")
flash $f ${f/1P/2P} -o $o -m $m -M $M -z >> ./stats/merged.out.txt 
done

echo ###### ============================ ######
echo ###### Q filtering and cropping     ######
echo ###### ============================ ######


#RENAME; CROP to size and FILTER merged pairs
for f in $(find ./merged/* -type f | grep 'extendedFrags.fastq.gz') ; do 
n=$(awk -F '/|.extendedFrags' '{print $3}' <<< "$f")
o=./process/${n}_quality.fastq.gz 
vsearch --fastx_filter $f --relabel $n. --fastq_stripright $sr --fastq_stripleft $sl --fastq_maxlen $ma --fastq_minlen $mi --fastq_maxee $ee --fastqout $o --log ./process/out.txt 
cat ./process/out.txt >> ./stats/filter.out.txt
done

echo "###### ============================ ######"
echo "###### Dereplicating sequences      ######"
echo "###### ============================ ######"

#DEREPLICATE by sample
for f in $(find ./process/* | grep '_quality.fastq.gz') ; do 
o=./process/$(awk -F '/|_quality' '{print $3}' <<< "$f")_derep.fasta
vsearch --derep_fulllength $f --strand plus --output  $o --sizeout --fasta_width 0 --log ./process/out.txt 
cat ./process/out.txt >> ./stats/derep.sample.out.txt
done 
#DEREPLICATE across samples
cat ./process/*derep.fasta > ./result/all_derep.fasta
vsearch --derep_fulllength ./result/all_derep.fasta --strand plus --output ./result/all_unique.fasta --sizein --sizeout --fasta_width 0 --uc ./result/all-unique.uc --log ./stats/derep.all.out.txt 

echo "###### ============================ ######"
echo "###### Clustering OTUS              ######"
echo "###### ============================ ######"

#OTU CLUSTERING
vsearch --cluster_size ./result/all_unique.fasta --id $i --strand plus --sizein  --sizeout  --fasta_width 0 --centroids ./result/centroids.fasta --log ./result/out.txt
echo Clusters generated: $(grep -c "^>" ./result/centroids.fasta) >> ./stats/otu.out.txt
cat ./result/out.txt >> ./stats/otu.out.txt

#Remove SINGLETONS
vsearch --sortbysize ./result/centroids.fasta --sizein --sizeout  --fasta_width 0 --minsize 2 --output ./result/sorted.fasta --log ./result/out.txt
echo Clusters after removal of singletons: $(grep -c "^>" ./result/sorted.fasta) >> ./stats/otu.out.txt
cat ./result/out.txt >> ./stats/otu.out.txt

#Remove CHIMERAS de novo
vsearch --uchime_denovo ./result/sorted.fasta --sizein --sizeout --fasta_width 0 --qmask none --nonchimeras ./result/denovo.nonchimeras.fasta --log ./result/out.txt
echo Clusters after removal of chimeras de novo: $(grep -c "^>" ./result/denovo.nonchimeras.fasta) >> ./stats/otu.out.txt
cat ./result/out.txt >> ./stats/otu.out.txt

#RENAME OTUs
vsearch --fastx_filter ./result/denovo.nonchimeras.fasta --fasta_width 0 --relabel OTU --fastaout ./result/otus.fasta

echo "###### ============================ ######"
echo "###### Mapping reads and OTU table  ######"
echo "###### ============================ ######"

#9.1Construct concatenated SEMIRAW reads
for f in $(find ./merged/* -type f | grep extendedFrags.fastq.gz) ; do 
n=$(awk -F '/|.extendedFrags' '{print $3}' <<< "$f")
o=./process/${n}_renamed.fasta 
vsearch --fastx_filter $f --fastaout $o --relabel $n.
done

#9.2Concatenate
cat ./process/*renamed.fasta > ./result/all_semiraw.fasta
rm  ./process/*renamed.fasta

sed -i 's/-/_/g' ./result/all_semiraw.fasta #reformatting the header  

#Create OTU TABLE based on semiraw merged pairs
vsearch --usearch_global ./result/all_semiraw.fasta --db ./result/otus.fasta --id $i --strand plus --sizein --sizeout --fasta_width 0  --qmask none --dbmask none --otutabout ./result/otutab.txt --log ./stats/table.out.txt

#Create OTU TABLE based on dereplicated and filtered merged pairs
#$vsearch --usearch_global ./result/all_derep.fasta --db ./result/otus.fasta --id $i --strand plus --sizein --sizeout --fasta_width 0  --qmask none --dbmask none --otutabout ./result/otutab.txt --log ./stats/table.out.txt

echo "###### ============================ ######"
echo "###### Assign taxonomy to OTUs      ######"
echo "###### ============================ ######"


#Taxa classification. Make sure to use the appropiate database for 16S and ITS2
#DATABASE for ITS UNITe
vsearch --sintax ./result/otus.fasta --db /home/fer/abel/analysis/db/utax_reference_dataset_all_04.04.2024.fasta --tabbedout ./result/otus.sintax.uniteall --strand both --sintax_cutoff $b --log ./stats/taxauniteall.out.txt
#DATABASE 16s using SILVA 
vsearch --sintax ./result/otus.fasta --db /home/fer/abel/analysis/db/SILVA_138_16S.udb --tabbedout ./result/sotus.sintax --strand both --sintax_cutoff 0.8 --log ./stats/taxasilva.out
#DATABASE for 16S using UNITE
#$vsearch --sintax ./result/otus.fasta --db /home/fer/abel/analysis/db/UNITE_public_10.05.2021.fasta --tabbedout ./result/otus.unitepublic.sintax --strand both --sintax_cutoff 0.5 --log ./stats/taxaunitepublic.out.txt
#DATABASE for 16S using RDP
#$vsearch --sintax ./result/otus.fasta --db /home/fer/abel/analysis/db/rdp_its_v2.fasta --tabbedout ./result/otusrdp.sintax --strand both --sintax_cutoff $b --log ./stats/taxardp.out.txt


#12 Past OTU classification. 
#vsearch --sintax ./result/otus.fasta --db /home/fer/abel/analysis/db/otu_ITS.fasta --tabbedout ./result/potus.sintax --strand both --sintax_cutoff $b --log ./stats/ptaxa.out.txt


echo "###### ============================ ######"
echo "###### Stats                        ######"
echo "###### ============================ ######"

#We figure out that 16S sequencues needed to convert to reverse complemente.
#USING VSEARCH
#ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
#brew install wget
#wget https://github.com/torognes/vsearch/releases/download/v2.18.0/vsearch-2.18.0-macos-x86_64.tar.gz
#tar xzf vsearch-2.18.0-macos-x86_64.tar.gz

#export vsearch=/home/fer/abel/analysis/db/Google-Drive/vsearch-2.18.0-macos-x86_64/bin/vsearch
##vsearch --fastx_revcomp ./result/otusAMF.fasta --fastaout ./result/otusAMF.fasta
##vsearch --fastx_revcomp ./result/otus.daniel.fasta --fastaout ./result/otusdaniel.fasta
#$vsearch --fastx_revcomp ./result/otus.unidentified.fasta --fastaout ./result/otus.unidentified.fasta
#$vsearch --fastx_revcomp ./result/otus.unidentified2.fasta --fastaout ./result/otus.unidentified2.fasta


