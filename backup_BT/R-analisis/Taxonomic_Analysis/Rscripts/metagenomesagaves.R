#Exploring Taxonomy with R
library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("patchwork")
#Set directory
setwd("/home/abel/LIM/R-analisis/Taxonomic_Analysis/")

#import dataset
#agavesmetagenomes <- import_biom("ateq.biom")

#agavesmetagenomes <- import_biom("biom/LP_all_contigs_mazorka.biom") #contigs asignacion taxonomica mazorka
#agavesmetagenomes <- import_biom("biom/LP_all_reads_mazorka.biom") #reads asignacion taxonomica mazorka


#agavesmetagenomes <- import_biom("biom/LP_all_contigs_BT.biom") #contigs asignacion taxonomica BetterLab

#agavesmetagenomes <- import_biom("biom/LP_all_reads_ALN.biom") #reads asignacion taxonomica Alnitak
#agavesmetagenomes <- import_biom("biom/LP_all_bins_ALN.biom") #reads asignacion taxonomica Alnitak
agavesmetagenomes <- import_biom("biom/LP_all_contigs_ALN.biom") #reads asignacion taxonomica Alnitak

class(agavesmetagenomes)
#View Dataset
#View(agavesmetagenomes@tax_table@.Data)

#Put the names in the columns
agavesmetagenomes@tax_table@.Data <- substring(agavesmetagenomes@tax_table@.Data, 4)
colnames(agavesmetagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#Exploring the abudance table
#View(agavesmetagenomes@otu_table@.Data)

#colnames(agavesmetagenomes@otu_table@.Data) <- c("LP1", "LP2", "LP3", "LP4", "LP5", "LP6", "LP7", "LP8", "LP9", "LP10", 
#                                                 "LP11", "LP12", "LP13", "LP14", "LP15", "LP16", "LP17", "LP18", "LP19", "LP20",
#                                                 "LP21", "LP22", "LP23", "LP24", "LP25" )
#

otu.table <- agavesmetagenomes@otu_table@.Data
taxa = agavesmetagenomes@tax_table@.Data

#Checl how many reads by samples
sample_sums(agavesmetagenomes)

#Save data
write.table(otu.table, file = paste(getwd(), paste("tables/otu.table","txt",sep = "."),sep = "/"),
            row.names = T, quote = F, sep = "\t", col.names = T)

write.table(taxa,file=paste(getwd(),paste("tables/taxa.table","txt",sep = "."),sep = "/"),
            quote = F, col.names = T, row.names = T, sep = "\t")

#Load all necessary packages
library(stringr)
library(ggplot2)
library(vegan)
library(labdsv)
library(MASS)
library(Rmisc)
library(phyloseq)
library(tidyr)
library(dplyr)
library(VennDiagram)
library(ggvenn) 
library(plyr)
#All features to graph the data
tema=theme(#axis.text.x = element_text(color="black",size=12, angle=0,hjust=0.5,vjust=1.5),
  axis.text.x = element_text(color="black",size=12, angle=90,hjust=0.5,vjust=1.5),
  axis.text.y = element_text(color="black",size=12, vjust = 1.5),
  axis.title = element_text(color="black",size=12, face = "bold"),
  axis.title.x.bottom = element_blank(),
  panel.border =element_rect(color = "black", fill = NA),#element_blank(),
  strip.text.x = element_text(size=12, color="black",face="bold"),
  strip.text.y = element_text(size=12, color="black",face="bold"),
  strip.placement = "outside", strip.background = element_rect(fill = "white"), 
  panel.background = element_rect(fill = "white",colour = "white",size = 0.8, linetype = "solid"),
  panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
  legend.position = "right", legend.text = element_text(color = "black",size=12), legend.direction = "vertical",
  legend.title = element_text(color = "black",size=14, face = "bold"),
  legend.key.size = unit(0.3,"cm"))

#Set variables and data
ranks=c("Kingdom","Phylum","Class","Order","Family","Genus", "Species", "otu.id")
my_pal=colorRampPalette(c("#0C5C82", "#9CADBC", "#697987", "#7B4563","#B07696", "#6F502A", "grey80","grey20"))

#Load otu table, metadata table and taxa table
otu=read.delim(file=paste(getwd(),paste("tables/otu.table","txt",sep = "."),sep = "/"))
#met=read.delim(file=paste(getwd(),paste("metadata/metadata.table","txt",sep = "."),sep = "/"))
met=read.delim(file=paste(getwd(),paste("metadata/metadata.tableB","txt",sep = "."),sep = "/"))
rownames(met)=met$sample.ID
tax=read.delim(paste(getwd(),paste("tables/taxa.table","txt",sep = "."),sep = "/"))
rownames(tax)=tax$X
tax=tax[,-1]

## Para cuando los se tienen espacion en blanco en los datos
tax$Kingdom[tax$Kingdom == ""] <- "Unidentified"
tax$Phylum[tax$Phylum == ""] <- "Unidentified"
tax$Class[tax$Class == ""] <- "Unidentified"
tax$Order[tax$Order == ""] <- "Unidentified"
tax$Family[tax$Family == ""] <- "Unidentified"
tax$Genus[tax$Genus == ""] <- "Unidentified"
tax$Species[tax$Species == ""] <- "Unidentified"


## Para cuando los se tiene NA en los datos
tax$Kingdom[is.na(tax$Kingdom)] = "Unidentified"
tax$Phylum[is.na(tax$Phylum )] = "Unidentified"
tax$Class[is.na(tax$Class)] = "Unidentified"
tax$Order[is.na(tax$Order)] = "Unidentified"
tax$Family[is.na(tax$Family)] = "Unidentified"
tax$Genus[is.na(tax$Genus)] = "Unidentified" 
tax$Species[is.na(tax$species)] = "Unidentified"

 
#####BARPLOT OF RELATIVE ABUNDANCE BY SITE######
t=c("Kingdom","Phylum","Class","Order","Family","Genus", "Species")[4]

##Sumarise data
set.seed(23171341)
otu.s=t(rrarefy(t(otu), sample = min(colSums(otu))))


#Checklist
colSums(otu.s)
otu.s=as.data.frame(cbind(otu.s,tax))
otu.t=gather(data=otu.s, key = "sample", value = "abs", colnames(otu))

for (i in met$sample.ID){
  otu.t[otu.t$sample == i, c( "plant.specie", "treatment", "plant.compartment")] = met[i,c("plant.specie", "treatment", "plant.compartment")]
}


#Total reads by sample 
otu.t0 = otu.t %>% group_by(sample) %>% summarise(absab = sum(abs))

otu.t1 = otu.t %>% group_by(otu.t[,t],sample, plant.specie, plant.compartment, treatment)  %>% summarise(absab = sum(abs)) ; colnames(otu.t1)[1]=c("rank")
#otu.t1 = otu.t %>% group_by(otu.t[,t],sample)  %>% summarise(absab = sum(abs)) ; colnames(otu.t1)[1]=c("rank") #para agregar mas variables
otu.t1$relab=otu.t1$absab/otu.t0$absab*100

#Remove low abundant taxa
orden = otu.t1 %>% group_by(rank) %>% summarise(relab = mean(relab)) 
View(orden)
#You need to check the first 20 ranks to select the relab value. 
low = orden[orden$relab < 0.7, "rank"]
otu.t1[otu.t1$rank %in% low$rank, "rank" ] = "low abundant"
n=dim(orden)[1]-length(low$rank)+1 ; print(c(n,"taxa")) 

#Reorder taxa if needed
orden=orden[order(orden$relab,decreasing = T),]
niveles=orden[!orden$rank %in% low$rank, "rank"]$rank
otu.t1$rank=factor(otu.t1$rank, level = c(niveles[niveles!="Unclassified"],"Unclassified", "low abundant"))
otu.t1$sample=factor(otu.t1$sample) 

#df <- otu.t1 %>% filter(plant.specie == "Agave.tequilana")
#df <- otu.t1 %>% filter(plant.specie == "Agave.karwinskii")
#df <- otu.t1 %>% filter(plant.specie == "Agave.angustifolia")
#df <- otu.t1 %>% filter(plant.specie == "Agave.convallis")
df <- otu.t1 %>% filter(plant.specie == "Agave.potatorum")



#Save the table to percetage
#write.table(otu.t1, file = paste(paste(t,data,"taxas","txt",sep = "."),sep = "/"),
            #row.names = F, quote = F, sep = "\t", col.names = T)


#Graph for domain
#pdf(file = "figR/RA_BY_SAMPLE-Domain_betterlab_contigs.pdf", colormodel = "cmyk", width = 8.5, height = 11, compress = F)

a=ggplot(data=df, aes(y=relab, x=sample, fill=rank))+
  geom_bar(stat="identity",width = 0.95)+
  scale_fill_manual(values = my_pal(n))+
  labs(x = "Sample",y = "relative abundance (%)")+
  guides(fill=guide_legend(title = "Domain", ncol = 1))+tema

a
ggsave(filename = "figR/RA_BY_A.pot-Domain_alnitak_contigs.svg", plot = a , width = 20, height = 16, dpi = 300, units = "cm")
#print(a)
#dev.off()

