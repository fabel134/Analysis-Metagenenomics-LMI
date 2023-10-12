#Exploring Taxonomy with R
library("phyloseq")
library("ggplot2")
library("RColorBrewer")
library("patchwork")
library("plyr")
library("plotly") ##Graficas interactivas

#Set directory
setwd("/home/abel/LIM/R-analisis/Taxonomic_Analysis/")

#import dataset

#agavesmetagenomes <- import_biom("biom/LP_all_reads_ALN.biom") #reads asignacion taxonomica Alnitak
#agavesmetagenomes <- import_biom("biom/LP_all_bins_ALN.biom") #reads asignacion taxonomica Alnitak
#agavesmetagenomes <- import_biom("biom/LP_all_contigs_ALN.biom") #reads asignacion taxonomica Alnitak
#agavesmetagenomes <- import_biom("biom/SRR_all.biom") #reads asignacion taxonomica Alnitak datos VIC
agavesmetagenomes <- import_biom("biom/all_LIM_metagenomes.biom") #reads asignacion taxonomica Alnitak datos VIC and LP

class(agavesmetagenomes)

#View Dataset
#View(agavesmetagenomes@tax_table@.Data)

#Put the names in the columns
agavesmetagenomes@tax_table@.Data <- substring(agavesmetagenomes@tax_table@.Data, 4)
colnames(agavesmetagenomes@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

agavesmetagenomes <- subset_taxa(agavesmetagenomes, Family != "Hominidae") ##Delete Human contamintation

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
aes4 <- c("#024AAC","#1DACA8","#10B62F", "#E2E41C","#F48F06","#F2252C","#D140B7","grey20")

#Load otu table, metadata table and taxa table
otu=read.delim(file=paste(getwd(),paste("tables/otu.table","txt",sep = "."),sep = "/"))
#met=read.delim(file=paste(getwd(),paste("metadata/metadata.table","txt",sep = "."),sep = "/"))
#met=read.delim(file=paste(getwd(),paste("metadata/metadata.tableR","txt",sep = "."),sep = "/"))
#met=read.delim(file=paste(getwd(),paste("metadata/metadata_vic","txt",sep = "."),sep = "/"))
met=read.delim(file=paste(getwd(),paste("metadata/all_LIM_metagenomes","txt",sep = "."),sep = "/"))

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
otu.s=t(rrarefy(t(otu), sample = min(colSums(otu)))) ##Normalitzation de datos

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
low = orden[orden$relab < 0.9, "rank"]
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
#df <- otu.t1 %>% filter(plant.specie == "Agave.potatorum")



#Save the table to percetage
#write.table(otu.t1, file = paste(paste(t,data,"taxas","txt",sep = "."),sep = "/"),
            #row.names = F, quote = F, sep = "\t", col.names = T)


#Graph for domain
#pdf(file = "figR/RA_BY_SAMPLE-Domain_betterlab_contigs.pdf", colormodel = "cmyk", width = 8.5, height = 11, compress = F)

a=ggplot(data=otu.t1, aes(y=relab, x=sample, fill=rank))+
  geom_bar(stat="identity",width = 0.95)+
  scale_fill_manual(values = my_pal(n))+
  labs(x = "Sample",y = "relative abundance (%)")+
  guides(fill=guide_legend(title = "Domain", ncol = 1))+tema

a
#ggsave(filename = "figR/RA_BY_A.pot-Domain_alnitak_contigs.svg", plot = a , width = 20, height = 16, dpi = 300, units = "cm")
#print(a)
#dev.off()
################ Graficos interactivos######################################
# Convertir el gráfico a uno interactivo con plotly
p_interactivo <- ggplotly(a)
# Mostrar el gráfico interactivo
p_interactivo
######################################################

##############Alfa y Beta Diversity#############
## StackBar
percentages <- transform_sample_counts(agavesmetagenomes, function(x) x*100 / sum(x) ) # Normalizacion de datos
#percentages <- agavesmetagenomes
head(percentages@otu_table@.Data)

percentages_df <- psmelt(percentages)
## Para cuando los se tienen espacion en blanco en los datos
percentages_df$Kingdom[percentages_df$Kingdom == ""] <- "Unidentified"
percentages_df$Phylum[percentages_df$Phylum == ""] <- "Unidentified"
percentages_df$Class[percentages_df$Class == ""] <- "Unidentified"
percentages_df$Order[percentages_df$Order == ""] <- "Unidentified"
percentages_df$Family[percentages_df$Family == ""] <- "Unidentified"
percentages_df$Genus[percentages_df$Genus == ""] <- "Unidentified"
percentages_df$Species[percentages_df$Species == ""] <- "Unidentified"

## Para cuando los se tiene NA en los datos
percentages_df$Kingdom[is.na(percentages_df$Kingdom)] = "Unidentified"
percentages_df$Phylum[is.na(percentages_df$Phylum )] = "Unidentified"
percentages_df$Class[is.na(percentages_df$Class)] = "Unidentified"
percentages_df$Order[is.na(percentages_df$Order)] = "Unidentified"
percentages_df$Family[is.na(percentages_df$Family)] = "Unidentified"
percentages_df$Genus[is.na(percentages_df$Genus)] = "Unidentified" 
percentages_df$Species[is.na(percentages_df$species)] = "Unidentified"

# Agregar los metadatos
percentages@sam_data <- sample_data(met)

##obtener solo filosfera
# Crear un objeto para filtrar las muestras
filtro_phyllosphere <- sample_data(percentages)$plant.compartment == "phyllosphere"

# Aplicar el filtro a tu objeto phyloseq
percentages_phyllosphere <- subset_samples(percentages, filtro_phyllosphere)

#LIMININAR LAS MUESTRAS QUE LE METEN RUIDO A NUESTRA GRAFICA 
# Crear un objeto para filtrar las muestras que NO sean "LP10" ni "LP15"
#filtro_muestras <- !(sample_names(percentages_phyllosphere) %in% c("LP10_reads", "LP9_reads"))
#filtro_muestras <- !(sample_names(percentages_phyllosphere) %in% c("LP10_reads"))

# Aplicar el filtro a tu objeto phyloseq para eliminar las muestras deseadas
#percentages_phyllosphere_filtrado <- subset_samples(percentages_phyllosphere, filtro_muestras) #eliminar las no deseadas
percentages_phyllosphere_filtrado <- percentages_phyllosphere
#head(percentages_phyllosphere_filtrado@otu_table@.Data)

## diversidad beta
#meta_ord_phyll <- ordinate(physeq = percentages_phyllosphere, method = "NMDS", distance = "bray") 
meta_ord_phyll <- ordinate(physeq = percentages_phyllosphere_filtrado, method = "NMDS", distance = "bray") 

# Crear el gráfico de ordenación
phill <- plot_ordination(
  physeq = percentages_phyllosphere_filtrado,
  ordination = meta_ord_phyll,
  color = "plant.specie"
) +
  #geom_text(mapping = aes(label = colnames(percentages_phyllosphere_filtrado@otu_table@.Data)), size = 3, vjust = 1.5) +
  geom_point(aes(shape = treatment), size = 3) +  # Asignar las formas según la variable "treatment"
  tema
stress_plot <- stressplot(meta_ord_phyll)
phill

#ggsave(filename = "figR/NMDS_phyllophere.svg", plot = phill , width = 20, height = 16, dpi = 300, units = "cm")

##obtener solo rizosfera
# Crear un objeto para filtrar las muestras
filtro_rhizosphere <- sample_data(percentages)$plant.compartment == "rhizosphere"

# Aplicar el filtro a tu objeto phyloseq
percentages_rhizosphere <- subset_samples(percentages, filtro_rhizosphere)

#LIMININAR LAS MUESTRAS QUE LE METEN RUIDO A NUESTRA GRAFICA 
# Crear un objeto para filtrar las muestras que NO sean "LP10" ni "LP15"
#filtro_muestras <- !(sample_names(percentages_rhizosphere) %in% c("LP10_reads", "LP9_reads"))
#filtro_muestras <- !(sample_names(percentages_rhizosphere) %in% c("LP10_reads"))

# Aplicar el filtro a tu objeto phyloseq para eliminar las muestras deseadas
#percentages_rhizosphere_filtrado <- subset_samples(percentages_rhizosphere, filtro_muestras) #eliminar las no deseadas
percentages_rhizosphere_filtrado <- percentages_rhizosphere

## diversidad beta
meta_ord_rhizos <- ordinate(physeq = percentages_rhizosphere_filtrado, method = "NMDS", distance = "bray") 

# Crear el gráfico de ordenación
rhizos <- plot_ordination(
  physeq = percentages_rhizosphere_filtrado,
  ordination = meta_ord_rhizos,
  color = "plant.specie"
) +
  #geom_text(mapping = aes(label = colnames(percentages_rhizosphere_filtrado@otu_table@.Data)), size = 3, vjust = 1.5) +
  geom_point(aes(shape = treatment), size = 3) +  # Asignar las formas según la variable "treatment"
  tema
ggsave(filename = "figR/NMDS_rhizos.svg", plot = rhizos , width = 20, height = 16, dpi = 300, units = "cm")

rhizos
stress_plot <- stressplot(meta_ord_rhizos)

##obtener solo root.soil
# Crear un objeto para filtrar las muestras
filtro_root.zone.soil <- sample_data(percentages)$plant.compartment == "root.zone.soil"

# Aplicar el filtro a tu objeto phyloseq
percentages_root.zone.soil <- subset_samples(percentages, filtro_root.zone.soil)

#LIMININAR LAS MUESTRAS QUE LE METEN RUIDO A NUESTRA GRAFICA 
# Crear un objeto para filtrar las muestras que NO sean "LP10" ni "LP15"
#filtro_muestras <- !(sample_names(percentages_root.zone.soil) %in% c("LP10_reads", "LP9_reads"))
#filtro_muestras <- !(sample_names(percentages_root.zone.soil) %in% c("LP10_reads"))

# Aplicar el filtro a tu objeto phyloseq para eliminar las muestras deseadas
#percentages_root.zone.soil_filtrado <- subset_samples(percentages_root.zone.soil, filtro_muestras) #eliminar las no deseadas
percentages_root.zone.soil_filtrado <- percentages_root.zone.soil

## diversidad beta
meta_ord_soil <- ordinate(physeq = percentages_root.zone.soil_filtrado, method = "NMDS", distance = "bray") 

# Crear el gráfico de ordenación
soil <- plot_ordination(
  physeq = percentages_root.zone.soil_filtrado,
  ordination = meta_ord_soil,
  color = "plant.specie"
) +
  #geom_text(mapping = aes(label = colnames(percentages_root.zone.soil_filtrado@otu_table@.Data)), size = 3, vjust = 1.5) +
  geom_point(aes(shape = treatment), size = 3) +  # Asignar las formas según la variable "treatment"
  tema
#ggsave(filename = "figR/NMDS_soil.svg", plot = soil , width = 20, height = 16, dpi = 300, units = "cm")
soil
stress_plot <- stressplot(meta_ord_soil)

## NMDS para todos los datos 
percentages <- transform_sample_counts(agavesmetagenomes, function(x) x*100 / sum(x) ) # Normalizacion de datos
#percentages <- agavesmetagenomes
#head(percentages@otu_table@.Data)

#LIMININAR LAS MUESTRAS QUE LE METEN RUIDO A NUESTRA GRAFICA 
# Crear un objeto para filtrar las muestras que NO sean "LP10" ni "LP15"
filtro_muestras <- !(sample_names(percentages) %in% c("SRR5166355_reads"))

# Aplicar el filtro a tu objeto phyloseq para eliminar las muestras deseadas
percentages_all_filtrado <- subset_samples(percentages, filtro_muestras) #eliminar las no deseadas
#percentages_all_filtrado <- percentages
# Agregar los metadatos
percentages_all_filtrado@sam_data <- sample_data(met)
meta_ord_all <- ordinate(physeq = percentages_all_filtrado, method = "NMDS", distance = "bray") 

b = plot_ordination(physeq = percentages_all_filtrado, ordination = meta_ord_all, color = "plant.compartment") +
  #geom_text(mapping = aes(label = colnames(percentages_all_filtrado@otu_table@.Data)), size = 3, vjust = 1.5)+
  geom_point(aes(shape = treatment), size = 3) +
  tema
#ggsave(filename = "figR/NMDS_all.svg", plot = b , width = 20, height = 16, dpi = 300, units = "cm")
b
stress_plot <- stressplot(meta_ord_all)

##obtener solo los mocks
# Agregar los metadatos
percentages@sam_data <- sample_data(met)
# Crear un objeto para filtrar las muestras
filtro_mock <- sample_data(percentages)$treatment == "mock" | sample_data(percentages)$treatment == "vic"
# Aplicar el filtro a tu objeto phyloseq
percentages_mock <- subset_samples(percentages, filtro_mock)

#LIMININAR LAS MUESTRAS QUE LE METEN RUIDO A NUESTRA GRAFICA 
# Crear un objeto para filtrar las muestras que NO sean "LP10" ni "LP15"
#filtro_muestras <- !(sample_names(percentages_mock) %in% c("LP10_reads", "LP9_reads"))
#filtro_muestras <- !(sample_names(percentages_mock) %in% c("LP10_reads"))

# Aplicar el filtro a tu objeto phyloseq para eliminar las muestras deseadas
#percentages_root.zone.soil_filtrado <- subset_samples(percentages_root.zone.soil, filtro_muestras) #eliminar las no deseadas
percentages_mock_filtrado <- percentages_mock
## diversidad beta
meta_ord_mock <- ordinate(physeq = percentages_mock_filtrado, method = "NMDS", distance = "bray") 

# Crear el gráfico de ordenación
mock <- plot_ordination(
  physeq = percentages_mock_filtrado,
  ordination = meta_ord_mock,
  color = "plant.compartment"
) +
  #geom_text(mapping = aes(label = colnames(percentages_root.zone.soil_filtrado@otu_table@.Data)), size = 3, vjust = 1.5) +
  geom_point(aes(shape = treatment), size = 3) +  # Asignar las formas según la variable "treatment"
  tema
#ggsave(filename = "figR/NMDS_mock_all.svg", plot = mock , width = 20, height = 16, dpi = 300, units = "cm")
mock
stress_plot <- stressplot(meta_ord_mock)


##obtener solo la filosfera de los mocks
# Agregar los metadatos
percentages@sam_data <- sample_data(met)
# Crear un objeto para filtrar las muestras
filtro_mock_phyllosphere <- (sample_data(percentages)$treatment == "mock" | sample_data(percentages)$treatment == "vic") & sample_data(percentages)$plant.compartment == "phyllosphere"
percentages_mock_phyllosphere <- subset_samples(percentages, filtro_mock_phyllosphere)

#LIMININAR LAS MUESTRAS QUE LE METEN RUIDO A NUESTRA GRAFICA 
# Crear un objeto para filtrar las muestras que NO sean "LP10" ni "LP15"
#filtro_muestras <- !(sample_names(percentages_mock_phyllosphere) %in% c("LP10_reads", "LP9_reads"))
#filtro_muestras <- !(sample_names(percentages_mock_phyllosphere) %in% c("LP10_reads"))

# Aplicar el filtro a tu objeto phyloseq para eliminar las muestras deseadas
#percentages_mock_phyllosphere_filtrado <- subset_samples(percentages_root.zone.soil, filtro_muestras) #eliminar las no deseadas
percentages_mock_phyllosphere_filtrado <- percentages_mock_phyllosphere
## diversidad beta
meta_ord_mock_phyllosphere <- ordinate(physeq = percentages_mock_phyllosphere_filtrado, method = "NMDS", distance = "bray") 

# Crear el gráfico de ordenación
mock_phyllosphere <- plot_ordination(
  physeq = percentages_mock_phyllosphere_filtrado,
  ordination = meta_ord_mock_phyllosphere,
  color = "plant.specie"
) +
  #geom_text(mapping = aes(label = colnames(percentages_root.zone.soil_filtrado@otu_table@.Data)), size = 3, vjust = 1.5) +
  geom_point(aes(shape = treatment), size = 3) +  # Asignar las formas según la variable "treatment"
  tema
#ggsave(filename = "figR/NMDS_mock_phyllosphere.svg", plot = mock_phyllosphere , width = 20, height = 16, dpi = 300, units = "cm")
mock_phyllosphere
stress_plot <- stressplot(meta_ord_mock_phyllosphere)

##obtener solo la rizosfera de los mocks
# Agregar los metadatos
percentages@sam_data <- sample_data(met)
# Crear un objeto para filtrar las muestras
filtro_mock_rhizosphere <- (sample_data(percentages)$treatment == "mock" | sample_data(percentages)$treatment == "vic") & sample_data(percentages)$plant.compartment == "rhizosphere"
percentages_mock_rhizosphere <- subset_samples(percentages, filtro_mock_rhizosphere)

#LIMININAR LAS MUESTRAS QUE LE METEN RUIDO A NUESTRA GRAFICA 
# Crear un objeto para filtrar las muestras que NO sean "LP10" ni "LP15"
#filtro_muestras <- !(sample_names(percentages_mock_rhizosphere) %in% c("LP10_reads", "LP9_reads"))
#filtro_muestras <- !(sample_names(percentages_mock_rhizosphere) %in% c("LP10_reads"))

# Aplicar el filtro a tu objeto phyloseq para eliminar las muestras deseadas
#percentages_mock_rhizosphere_filtrado <- subset_samples(percentages_root.zone.soil, filtro_muestras) #eliminar las no deseadas
percentages_mock_rhizosphere_filtrado <- percentages_mock_rhizosphere
## diversidad beta
meta_ord_mock_rhizosphere <- ordinate(physeq = percentages_mock_rhizosphere_filtrado, method = "NMDS", distance = "bray") 

# Crear el gráfico de ordenación
mock_rhizosphere <- plot_ordination(
  physeq = percentages_mock_rhizosphere_filtrado,
  ordination = meta_ord_mock_rhizosphere,
  color = "plant.specie"
) +
  #geom_text(mapping = aes(label = colnames(percentages_mock_rhizosphere_filtrado@otu_table@.Data)), size = 3, vjust = 1.5) +
  geom_point(aes(shape = treatment), size = 3) +  # Asignar las formas según la variable "treatment"
  tema
#ggsave(filename = "figR/NMDS_mock_mock_rhizosphere.svg", plot = mock_rhizosphere , width = 20, height = 16, dpi = 300, units = "cm")
mock_rhizosphere
stress_plot <- stressplot(meta_ord_mock_rhizosphere)

##obtener solo root.zone.soil de los mocks
# Agregar los metadatos
percentages@sam_data <- sample_data(met)
# Crear un objeto para filtrar las muestras
filtro_mock_root.zone.soil <- (sample_data(percentages)$treatment == "mock" | sample_data(percentages)$treatment == "vic") & sample_data(percentages)$plant.compartment == "root.zone.soil"
percentages_mock_root.zone.soil <- subset_samples(percentages, filtro_mock_root.zone.soil)

#LIMININAR LAS MUESTRAS QUE LE METEN RUIDO A NUESTRA GRAFICA 
# Crear un objeto para filtrar las muestras que NO sean "LP10" ni "LP15"
#filtro_muestras <- !(sample_names(percentages_mock_rhizosphere) %in% c("LP10_reads", "LP9_reads"))
#filtro_muestras <- !(sample_names(percentages_mock_rhizosphere) %in% c("LP10_reads"))

# Aplicar el filtro a tu objeto phyloseq para eliminar las muestras deseadas
#percentages_mock_root.zone.soil_filtrado <- subset_samples(percentages_root.zone.soil, filtro_muestras) #eliminar las no deseadas
percentages_mock_root.zone.soil_filtrado <- percentages_mock_root.zone.soil
## diversidad beta
meta_ord_mock_root.zone.soil <- ordinate(physeq = percentages_mock_root.zone.soil_filtrado, method = "NMDS", distance = "bray") 

# Crear el gráfico de ordenación
mock_root.zone.soil <- plot_ordination(
  physeq = percentages_mock_root.zone.soil_filtrado,
  ordination = meta_ord_mock_root.zone.soil,
  color = "plant.specie"
) +
  #geom_text(mapping = aes(label = colnames(percentages_mock_rhizosphere_filtrado@otu_table@.Data)), size = 3, vjust = 1.5) +
  geom_point(aes(shape = treatment), size = 3) +  # Asignar las formas según la variable "treatment"
  tema
#ggsave(filename = "figR/NMDS_mock_root.zone.soil.svg", plot = mock_root.zone.soil , width = 20, height = 16, dpi = 300, units = "cm")
mock_root.zone.soil
stress_plot <- stressplot(meta_ord_mock_root.zone.soil)

shannon <- estimate_richness(agavesmetagenomes)

ggplot(shannon, aes(x = rownames(shannon) , y = Shannon)) +
  geom_bar(stat = "identity", fill = "green") +
  labs(x = "Position", y = "Depth", title = "Kellermania") +
  theme_minimal()

## veremos aqui solo el reino "Bacteria"
merge_Bacteria <-subset_taxa(agavesmetagenomes,Kingdom=="Bacteria")
merge_Bacteria_df <- psmelt(merge_Bacteria)
write.csv(merge_Bacteria_df, "tables/LIM_met_merge_Bacteria.csv")

## veremos aqui solo el reino "Eukaryota"
merge_Eukaryota <-subset_taxa(agavesmetagenomes,Kingdom=="Eukaryota")
merge_Eukaryota_df <- psmelt(merge_Eukaryota)
write.csv(merge_Eukaryota_df, "tables/LIM_met_merge_Eukaryota.csv")

## cuantos "Archaea"
sum(agavesmetagenomes@tax_table@.Data[,"Kingdom"] == "Archaea")
## veremos aqui solo el reino "Archaea"
merge_Archaea<-subset_taxa(agavesmetagenomes,Kingdom=="Archaea")
merge_Archaea_df <- psmelt(merge_Archaea)
write.csv(merge_Archaea_df, "tables/LIM_met_merge_Archaea.csv")

## cuantos "Viruses"
sum(agavesmetagenomes@tax_table@.Data[,"Kingdom"] == "Viruses")
## veremos aqui solo el reino "Viruses"
merge_Viruses<-subset_taxa(agavesmetagenomes,Kingdom=="Viruses")
merge_Viruses_df <- psmelt(merge_Viruses)
write.csv(merge_Viruses_df, "tables/LIM_met_merge_Viruses.csv")

index = estimate_richness(merge_Bacteria)
#index = estimate_richness(merge_Archaea)
#index = estimate_richness(merge_Viruses)
#index = estimate_richness(agavesmetagenomes)

index$Sample <- row.names(index) # cambiar las rowname por una columna y asignarle el nombre de Sample
rownames(index) <- NULL #Anular las Rowname
index
############################################

# Calcular la diversidad de Shannon para cada tratamiento


for (i in met$sample.ID){
  index[index$Sample == i, c( "plant.specie", "treatment", "plant.compartment")] = met[i,c("plant.specie", "treatment", "plant.compartment")]
}

#shannon_values <- index %>% filter(plant.compartment %in% c("phyllosphere", "rhizosphere"))

shannon_values <- index %>% filter(plant.compartment %in% c("phyllosphere"))

#shannon_values <- index %>% filter(plant.compartment %in% c("rhizosphere"))
View(shannon_values)

##Estadistica descriptiva por tratamientos (Tratamiento)

##Shannon index
tapply(shannon_values$Shannon,shannon_values$treatment,summary)
##sd
tapply(shannon_values$Shannon,shannon_values$treatment,sd)

###obtenemos datos por medición
library(Rmisc)

#####################Indice Shannon######################################
shannon=summarySE(shannon_values,measurevar=c("Shannon"),groupvars=c("treatment"))
shannon#n=6, ya sabemos que tenemos replicas

#Evaluamos la distribución y homocedasticidad en todos los tratamientos
#### Realizamos el test de normalidad
#Ho: no hay diferencia entre la distribucion normal y esta distribucion
#Ho se rechaza con p<0.05

tapply(shannon_values$Shannon,shannon_values$treatment,shapiro.test)
##p es menor a 0.05, se rechaza la Ho (no es distribucion normal)

# Realizar el test de Levene para verificar la igualdad de varianzas
library(car)
levene_test_result <- leveneTest(Shannon ~ as.factor(treatment), data = shannon_values)

print(levene_test_result)


###Usamos prueba no paramétrica

kruskal.test(Shannon ~ treatment, data = shannon_values)##muy Significativo

pairwise.wilcox.test(shannon_values$Shannon, shannon_values$treatment,
                     p.adjust.method = "BH")

### plots

library(ggplot2)
library(viridis)
#library(ggsignif)
##usamos reorder para ordenar ascendentemente
ggplot(shannon_values, aes(x=reorder(treatment,Shannon), y=Shannon,fill=reorder(treatment,Shannon))) + 
  geom_boxplot()+
  #geom_dotplot(binaxis='y', stackdir='center'
  #                           ,stackratio=0.7, dotsize=0.4,position=position_dodge())+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  #labs(title = substitute(bold(paste('Foliar area of six months ',italic('A. marmorata'),
  #                                   ' plants at greenhouse conditions.'))),
  #     subtitle = "Kruskal-Wallis chi-squared = 26.713, df = 5, p-value = 6.488e-05", x = "Treatment",
  #     y="area (cm²)",fill="Treatment")+
  theme(text = element_text(size = 11, face = "plain",family = "Gill_Sans_MT"),##Todo el texto
        strip.text.x = element_text(size = 10, face = "bold",hjust = 0.5), ##Texto de titulos de graficos
        plot.title = element_text(hjust = 0, size = 12),
        plot.title.position = "plot",
        plot.subtitle = element_text(hjust = 0, size = 11))+
  scale_fill_viridis(alpha = 1, begin = 0, end = 0.92, direction = -1,
                     discrete = T, option = "D")
