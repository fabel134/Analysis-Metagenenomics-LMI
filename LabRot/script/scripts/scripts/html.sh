## Scripts para obtener los archivos html 
## origen: https://carpentries-incubator.github.io/metagenomics/

ls *.kraken | while read line
do
name=$(echo $line | cut -d'.' -f1 )
echo $name
cut -f2,3 $name.kraken > $name.krona.input
done

echo ***ktImportTaxonomy***

ls *.input | while read line
do
nameA=$(echo $line | cut -d'.' -f1)
echo $nameA
ktImportTaxonomy $nameA.krona.input -o $nameA.krona.out.html
done
