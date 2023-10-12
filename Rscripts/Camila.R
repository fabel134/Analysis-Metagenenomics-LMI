
library("phyloseq")
library("ggplot2")

setwd("/home/abel/LIM/R-analisis/Taxonomic_Analysis")

# lectura de datos en archivo .biom
# Camda2023 <- import_biom("/home/camila/GIT/ccm-bioinfomatica-lab/Hackaton_junio2023/PruebasHipotesis/camda23.biom")
# Assembly
# /home/haydee/CAMDA23/data/c23/taxonomy/camda23-assembly.biom
Camda2023 <- import_biom("biom/LP_all_reads_ALN.biom")

# ReadLevel
# /home/haydee/CAMDA23/data/c23/taxonomy/camda23-readlevel.biom
class(Camda2023)

# cambio de los nombres de los niveles taxonomicos
colnames(Camda2023@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# se elimina ruido del inicio de los nombres 
Camda2023@tax_table@.Data <- substring(Camda2023@tax_table@.Data,4)
# tabla de taxonomia
head(Camda2023@tax_table@.Data)
# tabla de otus
head(Camda2023@otu_table@.Data)
# colnames(Camda2023@otu_table@.Data) <- substring(colnames(Camda2023@otu_table@.Data),17)
# metadatos
# metadata_Camda2023 <- read.csv2("/home/camila/GIT/ccm-bioinfomatica-lab/Hackaton_junio2023/PruebasHipotesis/metadata_camda23.csv",header =  TRUE, row.names = 1, sep = ",")

metadata_Camda2023 <- read.csv2("metadata/metadata.tableR.txt",header =  TRUE, row.names = 1, sep = "\t")


#rownames(metadata_Camda2023) <- sample_names(metadata_Camda2023)
Camda2023@sam_data <- sample_data(metadata_Camda2023)
Camda2023 <- subset_taxa(Camda2023, Family != "Hominidae") ##Delete Human contamintation

Camda2023_df <- psmelt(Camda2023)
write.csv(Camda2023_df, "Camda2023_dataframe.csv")

## StackBar
percentages <- transform_sample_counts(Camda2023, function(x) x*100 / sum(x) )

head(percentages@otu_table@.Data)
percentages_df <- psmelt(percentages)

## Usaremos un comando llamado unique() para explorar cuántos reinos tenemos. 
unique(Camda2023@tax_table@.Data[,"Kingdom"])

## Subconjuntos separados por cada reino

## con esto podemos ver cuantos "Eukaryota" tenemos en "Kingdom"
sum(Camda2023@tax_table@.Data[,"Kingdom"] == "Eukaryota")
## veremos aqui solo el reino "Eukaryota"
merge_Eukaryota <-subset_taxa(Camda2023,Kingdom=="Eukaryota")
merge_Eukaryota_df <- psmelt(merge_Eukaryota)
write.csv(merge_Eukaryota_df, "Taxonomic_Analysis/tables/LIM_met_merge_Eukaryota.csv")

## cuantos "Bacteria"
sum(Camda2023@tax_table@.Data[,"Kingdom"] == "Bacteria")
## veremos aqui solo el reino "Bacteria"
merge_Bacteria <-subset_taxa(Camda2023,Kingdom=="Bacteria")
merge_Bacteria_df <- psmelt(merge_Bacteria)
write.csv(merge_Bacteria_df, "Taxonomic_Analysis/tables/LIM_met_merge_Bacteria.csv")

## cuantos "Archaea"
sum(Camda2023@tax_table@.Data[,"Kingdom"] == "Archaea")
## veremos aqui solo el reino "Archaea"
merge_Archaea<-subset_taxa(Camda2023,Kingdom=="Archaea")
merge_Archaea_df <- psmelt(merge_Archaea)
write.csv(merge_Archaea_df, "Taxonomic_Analysis/tables/LIM_met_merge_Archaea.csv")

## cuantos "Viruses"
sum(Camda2023@tax_table@.Data[,"Kingdom"] == "Viruses")
## veremos aqui solo el reino "Viruses"
merge_Viruses<-subset_taxa(Camda2023,Kingdom=="Viruses")
merge_Viruses_df <- psmelt(merge_Viruses)
write.csv(merge_Viruses_df, "Taxonomic_Analysis/tables/LIM_met_merge_Viruses.csv")

## funcion de aglomeracion por nivel taxonomico
glomToGraph<-function(phy,tax){
  ## creamos el subconjunto dependiendo del linaje taxonomico deseado
  glom <- tax_glom(phy, taxrank = tax)
  ## sacamos los porcentajes
  percentages <- transform_sample_counts(glom, function(x) x*100 / sum(x) )
  percentages_df <- psmelt(percentages)
  return(list(glom,percentages,percentages_df))
}
##########################################################
##########################################################
##########################################################

 percentages_df$Sample<-as.factor(percentages_df$Sample)
 percentages_df$plant.specie<-as.factor(percentages_df$plant.specie)
 # Ordenamos respecto a categoria 
 percentages_df<-percentages_df[order(percentages_df$plant.specie,percentages_df$Sample),]
 
 
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
 
 a=ggplot(data=percentages_df, aes_string(x='Sample', y='Abundance', fill='Order'))  +
   geom_bar(aes(), stat="identity", position="stack") +
   facet_wrap(~plant.specie, scales = "free")+
   #scale_x_discrete(limits = rev(levels(percentages_df$Category))) +
   labs(title = "Abundance", x='Sample', y='Abundance', color = 'plant.specie') +
   theme(legend.key.size = unit(0.2, "cm"),
         legend.key.width = unit(0.25,"cm"),
         legend.position = "bottom",
         legend.direction = "horizontal",
         legend.title=element_text(size=8, face = "bold"),
         legend.text=element_text(size=6),
         text = element_text(size=12),
         axis.text.x = element_text(angle=90, size=5, hjust=1, vjust=0.5))
a
#ggsave(filename = "Taxonomic_Analysis/figR/RA_BY_Phylum_betterlab_contigs.svg", plot = a , width = 20, height = 16, dpi = 300, units = "cm")
 
##########################################################
##########################################################
##########################################################

Abundance_barras <- function(phy,tax,attribute,abundance_percentage){
  ##llamar funcion de datos 
  Data <- glomToGraph(phy,tax)
  glom <- Data[[1]] #phyloseq
  percentages <- Data[[2]] #phyloseq
  percentages_df <- Data[[3]] # dataframe
  ## Graficamos para cada subconjunto las barras de abundancia
  plot_barras <- ggplot(data=percentages_df, aes_string(x='Sample', y='Abundance', fill=tax ,color=attribute)) + 
    scale_colour_manual(values=c('white','black')) +
    geom_bar(aes(), stat="identity", position="stack") +
    labs(title = "Abundance", x='Sample', y='Abundance', color = tax) +
    theme(legend.key.size = unit(0.2, "cm"),
          legend.key.width = unit(0.25,"cm"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title=element_text(size=8, face = "bold"),
          legend.text=element_text(size=6),
          text = element_text(size=12),
          axis.text.x = element_text(angle=90, size=5, hjust=1, vjust=0.5))
  percentages_df$tax<-percentages_df[,ncol(percentages_df)]
  percentages_df$tax[percentages_df$Abundance < abundance_percentage] <- "Others"
  percentages_df$tax <- as.factor(percentages_df$tax)
  plot_percentages <- ggplot(data=percentages_df, aes_string(x='Sample', y='Abundance', fill='tax' ,color=attribute))+
    scale_colour_manual(values=c('white','black')) +
    geom_bar(aes(), stat="identity", position="stack") +
    labs(title = "Abundance", x='Sample', y='Abundance', color = tax) +
    theme(legend.key.size = unit(0.3, "cm"),
          legend.key.width = unit(0.5,"cm"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title=element_text(size=10, face = "bold"),
          legend.text=element_text(size=8),
          text = element_text(size=12),
          axis.text.x = element_text(angle=90, size=5, hjust=1, vjust=0.5))
  return(list(plot_barras,plot_percentages))
}

## diversidad alfa
index = estimate_richness(Camda2023)
index
write.csv(index, "tables/indices_Alfa_reads.csv")

alpha_div_plot <- plot_richness(physeq = Camda2023, measures = c("Observed","Chao1","Shannon","simpson"))
alpha_div_plot
ggsave("Diversidad_Alfa_samples_reads.svg", plot = alpha_div_plot, path = "figR/" , width = 30, height = 15, dpi = 300, units = "cm")

alpha_div_plot_plant_compartment <- plot_richness(physeq = Camda2023, measures = c("Observed","Chao1","Shannon","simpson"),x = "plant.compartment", color = "plant.compartment") 
alpha_div_plot_plant_compartment
ggsave(filename = "figR/alpha_div_plot_plant_compartment_alnitak_reads.svg", plot = alpha_div_plot_plant_compartment , width = 20, height = 16, dpi = 300, units = "cm")

alpha_div_plot_site <- plot_richness(physeq = Camda2023, measures = c("Observed","Chao1","Shannon","simpson"),x = "Site", color = "Site") 
alpha_div_plot_site
ggsave(filename = "figR/alpha_div_plot_site_alnitak_reads.svg", plot = alpha_div_plot_plant_compartment , width = 20, height = 16, dpi = 300, units = "cm")

alpha_div_plot_specie <- plot_richness(physeq = Camda2023, measures = c("Observed","Chao1","Shannon","simpson"),x = "plant.specie", color = "plant.specie") 
alpha_div_plot_specie
ggsave(filename = "figR/alpha_div_plot_specie_alnitak_reads.svg", plot = alpha_div_plot_plant_compartment , width = 20, height = 16, dpi = 300, units = "cm")

# funcion automatizada de alfa diversidad
Alpha_diversity <- function(phy,tax,attribute){
  ## llamamos la funcion que crea los dataset
  Data <- glomToGraph(phy,tax)
  glom <- Data[[1]]
  percentages <- Data[[2]]
  percentages_df <- Data[[3]]
  ## Alfa diversidad
  plot_alpha <- plot_richness(physeq = glom, measures = c("Observed","Chao1","Shannon","simpson"),x = attribute, color = attribute) 
  return(plot_alpha)
}

## diversidad beta
meta_ord <- ordinate(physeq = percentages, method = "NMDS", distance = "bray") 

b = plot_ordination(physeq = percentages, ordination = meta_ord, color = "plant.compartment") +
  geom_text(mapping = aes(label = colnames(Camda2023@otu_table@.Data)), size = 3, vjust = 1.5)
ggsave(filename = "figR/diversidad_beta_plant.compartment_alnitak_reads.svg", plot = b , width = 20, height = 16, dpi = 300, units = "cm")

c = plot_ordination(physeq = percentages, ordination = meta_ord, color = "plant.specie") +
  geom_text(mapping = aes(label = colnames(Camda2023@otu_table@.Data)), size = 3, vjust = 1.5)
#ggsave(filename = "Taxonomic_Analysis/figR/diversidad_beta_plant.specie_betterlab_contigs.svg", plot = c , width = 20, height = 16, dpi = 300, units = "cm")

plot_ordination(physeq = percentages, ordination = meta_ord, color = "treatment") +
  geom_text(mapping = aes(label = colnames(Camda2023@otu_table@.Data)), size = 3, vjust = 1.5)
#ggsave("DiversidadBeta_Climate.png", plot = last_plot(), path = "/home/camila/GIT/ccm-bioinfomatica-lab/Hackaton_junio2023/PruebasHipotesis/Results_IMG"  , width = 30, height = 15, dpi = 300, units = "cm")


# funcion automatizada para diversidad beta
Beta_diversity <- function(phy,tax,attribute,distance){
  Data <- glomToGraph(phy,tax)
  glom <- Data[[1]]
  #CREAR UN GLOM AL 10%
  percentages <- Data[[2]]
  percentages_df <- Data[[3]]
  ## Beta diversidad 
  meta_ord <- ordinate(physeq = percentages, method = "NMDS", distance = distance) 
  plot_beta <- plot_ordination(physeq = percentages, ordination = meta_ord, color = attribute) +
    geom_text(mapping = aes(label = colnames(phy@otu_table@.Data)), size = 3, vjust = 1.5)
  return(plot_beta)
}

### Graficas por Familia y Género, para cada uno de los reinos

##-----------Eukaryota
###-----------Family
# Barras_Species <- Abundance_barras(merge_Eukaryota,'Family','ID_city',10.0)
# Barras_Species[1] 
# ggsave("Barras_Eukaryota_Family.png", plot = last_plot(), path = "/home/camila/GIT/ccm-bioinfomatica-lab/Hackaton_junio2023/PruebasHipotesis/Results_IMG"  , width = 30, height = 15, dpi = 300, units = "cm")
# Barras_Species[2]
# ggsave("Barras10_Eukaryota_Family.png", plot = last_plot(), path = "/home/camila/GIT/ccm-bioinfomatica-lab/Hackaton_junio2023/PruebasHipotesis/Results_IMG"  , width = 30, height = 15, dpi = 300, units = "cm")
Data<- glomToGraph(merge_Eukaryota , 'Family')
glom <- Data[[1]]
index_Eukaryota_Family <- estimate_richness(glom)
Alpha_diversity(merge_Eukaryota , 'Family' , 'Site')
ggsave("DiversidadAlfa_Eukaryota_Family.png", plot = last_plot(), path = "/home/camila/GIT/ccm-bioinfomatica-lab/Hackaton_junio2023/PruebasHipotesis/Results_IMG"  , width = 30, height = 15, dpi = 300, units = "cm")

Beta_diversity(merge_Eukaryota, 'Family' , 'ID_city', 'bray')
ggsave("DiversidadBeta_Eukaryota_Family.png", plot = last_plot(), path = "/home/camila/GIT/ccm-bioinfomatica-lab/Hackaton_junio2023/PruebasHipotesis/Results_IMG"  , width = 30, height = 15, dpi = 300, units = "cm")
###-----------Genero 
Barras_Species <- Abundance_barras(merge_Eukaryota,'Genus','ID_city',10.0)
Barras_Species[1] 
Barras_Species[2]
Alpha_diversity(merge_Eukaryota , 'Genus' , 'ID_city')
ggsave("DiversidadAlfa_Eukaryota_Genus.png", plot = last_plot(), path = "/home/camila/GIT/ccm-bioinfomatica-lab/Hackaton_junio2023/PruebasHipotesis/Results_IMG"  , width = 30, height = 15, dpi = 300, units = "cm")
Beta_diversity(merge_Eukaryota, 'Genus' , 'ID_city', 'bray')
ggsave("DiversidadBeta_Eukaryota_Genus.png", plot = last_plot(), path = "/home/camila/GIT/ccm-bioinfomatica-lab/Hackaton_junio2023/PruebasHipotesis/Results_IMG"  , width = 30, height = 15, dpi = 300, units = "cm")

##-----------Bacteria

###-----------Family
Barras_Species <- Abundance_barras(merge_Bacteria,'Family','ID_city',10.0)
Barras_Species[1] 
Barras_Species[2]
Data <- glomToGraph(merge_Bacteria , 'Family')
glom <- Data[[1]]
index_Bacteria_Family <- estimate_richness(glom)
Alpha_diversity(merge_Bacteria , 'Family' , 'ID_city')
ggsave("DiversidadAlfa_Bacteria_Family.png", plot = last_plot(), path = "/home/camila/GIT/ccm-bioinfomatica-lab/Hackaton_junio2023/PruebasHipotesis/Results_IMG"  , width = 30, height = 15, dpi = 300, units = "cm")
Beta_diversity(merge_Bacteria, 'Family' , 'ID_city', 'bray')
ggsave("DiversidadBeta_Bacteria_Family.png", plot = last_plot(), path = "/home/camila/GIT/ccm-bioinfomatica-lab/Hackaton_junio2023/PruebasHipotesis/Results_IMG"  , width = 30, height = 15, dpi = 300, units = "cm")
###-----------Genero 
Barras_Species <- Abundance_barras(merge_Bacteria,'Genus','ID_city',10.0)
Barras_Species[1] 
Barras_Species[2]
Alpha_diversity(merge_Bacteria , 'Genus' , 'ID_city')
ggsave("DiversidadAlfa_Bacteria_Genus.png", plot = last_plot(), path = "/home/camila/GIT/ccm-bioinfomatica-lab/Hackaton_junio2023/PruebasHipotesis/Results_IMG"  , width = 30, height = 15, dpi = 300, units = "cm")
Beta_diversity(merge_Bacteria, 'Genus' , 'ID_city', 'bray')
ggsave("DiversidadBeta_Bacteria_Genus.png", plot = last_plot(), path = "/home/camila/GIT/ccm-bioinfomatica-lab/Hackaton_junio2023/PruebasHipotesis/Results_IMG"  , width = 30, height = 15, dpi = 300, units = "cm")

##-----------Archaea
###-----------Family
Barras_Species <- Abundance_barras(merge_Archaea,'Family','ID_city',10.0)
Barras_Species[1] 
Barras_Species[2]
Alpha_diversity(merge_Archaea , 'Family' , 'ID_city')
ggsave("DiversidadAlfa_Archaea_Family.png", plot = last_plot(), path = "/home/camila/GIT/ccm-bioinfomatica-lab/Hackaton_junio2023/PruebasHipotesis/Results_IMG"  , width = 30, height = 15, dpi = 300, units = "cm")
Beta_diversity(merge_Archaea, 'Family' , 'ID_city', 'bray')
ggsave("DiversidadBeta_Archaea_Family.png", plot = last_plot(), path = "/home/camila/GIT/ccm-bioinfomatica-lab/Hackaton_junio2023/PruebasHipotesis/Results_IMG"  , width = 30, height = 15, dpi = 300, units = "cm")
###-----------Genero 
Barras_Species <- Abundance_barras(merge_Archaea,'Genus','ID_city',10.0)
Barras_Species[1] 
Barras_Species[2]
Alpha_diversity(merge_Archaea , 'Genus' , 'ID_city')
ggsave("DiversidadAlfa_Archaea_Genus.png", plot = last_plot(), path = "/home/camila/GIT/ccm-bioinfomatica-lab/Hackaton_junio2023/PruebasHipotesis/Results_IMG"  , width = 30, height = 15, dpi = 300, units = "cm")
Beta_diversity(merge_Archaea, 'Genus' , 'ID_city', 'bray')
ggsave("DiversidadBeta_Archaea_Genus.png", plot = last_plot(), path = "/home/camila/GIT/ccm-bioinfomatica-lab/Hackaton_junio2023/PruebasHipotesis/Results_IMG"  , width = 30, height = 15, dpi = 300, units = "cm")

##-----------Viruses
###-----------Family
Barras_Species <- Abundance_barras(merge_Viruses,'Family','ID_city',10.0)
Barras_Species[1] 
Barras_Species[2]
Alpha_diversity(merge_Viruses , 'Family' , 'ID_city')
ggsave("DiversidadAlfa_Viruses_Family.png", plot = last_plot(), path = "/home/camila/GIT/ccm-bioinfomatica-lab/Hackaton_junio2023/PruebasHipotesis/Results_IMG"  , width = 30, height = 15, dpi = 300, units = "cm")
Beta_diversity(merge_Viruses, 'Family' , 'ID_city', 'bray')
ggsave("DiversidadBeta_Viruses_Family.png", plot = last_plot(), path = "/home/camila/GIT/ccm-bioinfomatica-lab/Hackaton_junio2023/PruebasHipotesis/Results_IMG"  , width = 30, height = 15, dpi = 300, units = "cm")
###-----------Genero 
Barras_Species <- Abundance_barras(merge_Viruses,'Genus','ID_city',10.0)
Barras_Species[1] 
Barras_Species[2]
Alpha_diversity(merge_Viruses , 'Genus' , 'ID_city')
ggsave("DiversidadAlfa_Viruses_Genus.png", plot = last_plot(), path = "/home/camila/GIT/ccm-bioinfomatica-lab/Hackaton_junio2023/PruebasHipotesis/Results_IMG"  , width = 30, height = 15, dpi = 300, units = "cm")
Beta_diversity(merge_Viruses, 'Genus' , 'ID_city', 'bray')
ggsave("DiversidadBeta_Viruses_Genus.png", plot = last_plot(), path = "/home/camila/GIT/ccm-bioinfomatica-lab/Hackaton_junio2023/PruebasHipotesis/Results_IMG"  , width = 30, height = 15, dpi = 300, units = "cm")

