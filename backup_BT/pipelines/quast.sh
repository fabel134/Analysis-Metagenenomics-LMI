# This script is to run a qualylative analysis by quast
# Author: Zulia Fer
# Dataset: Metagenomes 

# For specific analyses fungus's sequences
mkdir 6_Quast_all_assamblies
quast.py -m 100 -t 4 5_Assamblies/Megahit/*.fa* -o 6_Quast_all_assamblies        
