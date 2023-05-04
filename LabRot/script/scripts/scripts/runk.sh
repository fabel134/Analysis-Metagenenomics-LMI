##Scripts para asignaciones taxonomicas con kraken usando contig (formato fasta)

##directory contig:lista de todos los archivos .fasta a analizar con kraken
##crea una carpeta llamada TAXONOMY_MAG
ls *.fasta | while read line
do
echo $line
#exit
name=$(echo $line | cut -d'.' -f1)
#echo $name
#exit
#echo kraken2 --db /home/betterlab/kraken2/database/db_kraken2 --threads 8 --output TAXONOMY_MAG/$name.kraken --report TAXONOMY_MAG/$name.report contig/$name.fasta

kraken2 --db /home/betterlab/kraken2/database/db_kraken2 --threads 8 --output TAXONOMY_MAG/$name.kraken --report TAXONOMY_MAG/$name.report $name.fasta
#exit
#exit
done
