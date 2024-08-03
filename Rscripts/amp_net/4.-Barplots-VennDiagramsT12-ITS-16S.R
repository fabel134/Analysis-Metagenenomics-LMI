#Load all necessary libraries
library(stringr)
library(ggplot2)
library(vegan)
library(labdsv)
library(MASS)
library(Rmisc)
library(scales)
library(phyloseq)
library(tidyr)
library(dplyr)
library(plyr)
library(ggVennDiagram)
library(RColorBrewer)

#Choose all setting to plots
tema=theme(axis.text.x = element_text(color="black",size=12, angle=0,hjust=0.5,vjust=1.5, family = "sans" ),
           #axis.text.x = element_text(color="black",size=12, angle=90,hjust=0.5,vjust=1.5),
           axis.text.y = element_text(color="black",size=12, vjust = 1.5, family = "sans"),
           axis.title = element_text(color="black",size=12, face = "bold", family = "sans"),
           axis.title.x.bottom = element_blank(),
           panel.border =element_rect(color = "black", fill = NA),#element_blank(),
           strip.text.x = element_text(size=12, color="black",face="bold", family = "sans"),
           strip.text.y = element_text(size=12, color="black",face="bold", family = "sans"),
           strip.placement = "outside", strip.background = element_rect(fill = "white"), 
           panel.background = element_rect(fill = "white",colour = "white",size = 0.8, linetype = "solid"),
           panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
           legend.position = "right", legend.text = element_text(color = "black",size=12, family = "sans"), 
           legend.direction = "vertical", legend.title = element_text(color = "black",size=12, face = "bold", family = "sans"),
           legend.key.size = unit(0.4,"cm"))

# Set PARAMETERS
data=c("abel", "nestor", "all")[1]
amplicon=c("ITS2","16S")[1]
base=c("uniteall", "silva")[1]
my_pal=colorRampPalette(c("#024AAC","#1DACA8","#10B62F","#E2E41C","#F48F06","#F2252C","#D140B7", "grey80","dimgrey"))
pc=c("soil", "rhizosphere", "root.endosphere", "leaf.endosphere", "T0", "T6", "T12", "EV", "CAMP", "MOCK","PFC","PFCS","spores")
#Upload data R Base
if (amplicon == "ITS2"){
  setwd("~/abel/analysis/FITS_analy/result/") 
  tax=read.delim(paste(getwd(),paste("taxa.table",data,base,"txt",sep = "."),sep = "/"))
  #tax=read.delim(paste(getwd(),paste("taxa.table",data,"txt",sep = "."),sep = "/"))
  otu=read.delim(file=paste(getwd(),paste("otu.table",data,"txt",sep = "."),sep = "/"))
  colnames(otu)=gsub("X","",colnames(otu))
  met=read.delim(file=paste(getwd(),paste("metadata.table",data,"txt",sep = "."),sep = "/"))
  rownames(met)=met$ID_seq_ITS
  met=met[order(met$ID_seq_ITS),]
  met=met[order(met$ID_seq_ITS),]
  met=met[colnames(otu),]
  #Subset of data to process only at 12 months and soil
  met= subset(met, met$sampling != "T12")
  met= subset(met, met$sampling != "T0")
  #met=na.omit(met) ##delete rows withNA
  otu=otu[rownames(met)]
  ####Total reads by PLANT SPECIE AND PLANT COMPARTMENT
  t=c("domain","phylum","class","order","family","genus","otu.id")[5]
  ##Sumarise data
  set.seed(23171341)
  otu.s=t(rrarefy(t(otu), sample = min(colSums(otu))))
  #Checklist
  colSums(otu.s)
  #sum(rownames(otu.s)==row.names(tax))==dim(otu.s)[1]
  otu.s=as.data.frame(cbind(otu.s,tax))
  otu.t=gather(data=otu.s, key = "sample", value = "abs", colnames(otu))
  otu.t[is.na(otu.t[,t]),t]="Unclassified"
  for (i in met$ID_seq_ITS){
    otu.t[otu.t$sample == i, c("compartment","treatment","replica","site", "sampling", "ID_seq_ITS")] = met[i,c( "compartment","treatment","replica","site", "sampling", "ID_seq_ITS")]
  }
  
  #Total reads
  otu.t0 = otu.t %>% group_by(replica,treatment) %>% summarise(absab = sum(abs))
  otu.t1 = otu.t %>% group_by(otu.t[,t],replica,treatment) %>% summarise(absab = sum(abs))
  names(otu.t1)[1]=c("rank")
  sum(paste(otu.t1$replica,otu.t1$treatment,sep = ".")==paste(otu.t1$replica,otu.t1$treatment,sep = "."))==dim(otu.t1)[1]
  otu.t1$relab=otu.t1$absab/otu.t0$absab*100
  
  #Remove low abundant taxa
  orden = otu.t1 %>% group_by(rank) %>% summarise(relab = mean(relab)) 
  orden1 = otu.t1 %>% group_by(rank, treatment) %>% summarise(relab = mean(relab)) 
  #Check how many groups you can select for plots, modify relab
  #order
  low = orden[orden$relab < 0.96999746, "rank"]
  #phylum
  #low = orden[orden$relab < 0, "rank"]
  otu.t1[otu.t1$rank %in% low$rank, "rank" ] = "other"
  n=dim(orden)[1]-length(low$rank)+1 ; print(c(n,"taxa")) 
  #Reorder taxa if needed
  orden=orden[order(orden$relab,decreasing = T),]
  niveles=orden[!orden$rank %in% low$rank, "rank"]$rank
  otu.t1$rank=factor(otu.t1$rank, level = c(niveles[niveles!="Unclassified"],"Unclassified", "other"))
  otu.t1$replica=factor(otu.t1$replica, level = pc)
  otu.t1$treatment=factor(otu.t1$treatment, level = pc)
  
} else if ( amplicon == "16S"){
  setwd("~/abel/analysis/16s_analy/result/")
  tax=read.delim(paste(getwd(),paste("taxa.table",data,base,"txt",sep = "."),sep = "/"))
  otu=read.delim(file=paste(getwd(),paste("otu.table",data,"txt",sep = "."),sep = "/"))
  colnames(otu)=gsub("X","",colnames(otu))
  met=read.delim(file=paste(getwd(),paste("metadata.table",data,"txt",sep = "."),sep = "/"))
  rownames(met)=met$ID_seq_16s
  met=met[order(met$ID_seq_16s),]
  met=met[order(met$ID_seq_16s),]
  met=met[colnames(otu),]
  #Subset of data to process only at 12 months and soil
  #met= subset(met, met$treatment != "T0")
  #met= subset(met, met$treatment != "T6")
  otu=otu[rownames(met)]
  ####Total reads by PLANT SPECIE AND PLANT COMPARTMENT
  t=c("domain","phylum","class","order","family","genus","otu.id")[5]
  ##Sumarise data
  set.seed(23171341)
  otu.s=t(rrarefy(t(otu), sample = min(colSums(otu))))
  #Checklist
  colSums(otu.s)
  #sum(rownames(otu.s)==row.names(tax))==dim(otu.s)[1]
  otu.s=as.data.frame(cbind(otu.s,tax))
  otu.t=gather(data=otu.s, key = "sample", value = "abs", colnames(otu))
  otu.t[is.na(otu.t[,t]),t]="Unclassified"
  for (i in met$ID_seq_16s){
    otu.t[otu.t$sample == i, c("compartment","treatment","replica","site", "sampling", "ID_seq_ITS")] = met[i,c("compartment","treatment","replica","site", "sampling", "ID_seq_ITS")]
  }
  
  #Total reads
  otu.t0 = otu.t %>% group_by(plant.species,plant.compartment) %>% summarise(absab = sum(abs))
  otu.t1 = otu.t %>% group_by(otu.t[,t],plant.species,plant.compartment) %>% summarise(absab = sum(abs))
  names(otu.t1)[1]=c("rank")
  sum(paste(otu.t1$plant.species,otu.t1$plant.compartment,sep = ".")==paste(otu.t0$plant.species,otu.t0$plant.compartment,sep = "."))==dim(otu.t1)[1]
  otu.t1$relab=otu.t1$absab/otu.t0$absab*100
  #Remove low abundant taxa
  orden = otu.t1 %>% group_by(rank) %>% summarise(relab = mean(relab)) 
  orden1 = otu.t1 %>% group_by(plant.compartment, rank) %>% summarise(relab = mean(relab)) 
  #Check how many groups you can select for plots, modify relab
  #Order
  #low = orden[orden$relab < 1.17, "rank"]
  #Phylum or OTU
  low = orden[orden$relab < 0, "rank"]
  #Genus
  #low = orden[orden$relab < 5.256947e-01, "rank"]
  otu.t1[otu.t1$rank %in% low$rank, "rank" ] = "other"
  n=dim(orden)[1]-length(low$rank)+1 ; print(c(n,"taxa")) 
  
  #Reorder taxa if needed
  orden=orden[order(orden$relab,decreasing = T),]
  niveles=orden[!orden$rank %in% low$rank, "rank"]$rank
  otu.t1$rank=factor(otu.t1$rank, level = c(niveles[niveles!="Unclassified"],"Unclassified", "other"))
  otu.t1$plant.species=factor(otu.t1$plant.species, level = pc)
  otu.t1$plant.compartment=factor(otu.t1$plant.compartment, level = pc)

}

#Select the folder to save the figures
if (amplicon == "ITS2"){
  setwd("~/abel/analysis/FITS_analy/result/FIGURASS-ITS/")
} else if ( amplicon == "16S"){
  setwd("~/abel/analysis/16s_analy/result/FIGURASS-16S/")
}  

#FIGURES with labels by order
pdf(paste(amplicon, t, file = "RA BY PC-T12.pdf"), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
p=ggplot(data=otu.t1, aes(y=relab, x=replica, fill=rank))+
  geom_bar(stat="identity",width = 0.95)+
  scale_fill_manual(values = my_pal(n))+
  labs(x = "replica",y = "relative abundance (%)")+
  facet_grid(~treatment, scales = "free", space = "free_x")+
  guides(fill=guide_legend(title = "order", ncol = 1, keyheight = 0.7, default.unit="cm",
                           label.theme = element_text(size= 12, face = "plain", family = "sans")))+tema
p
dev.off()

#FIGURES without legend by order
pdf(paste(amplicon,t, file = "RA BY PC-T12-1.pdf"), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
p1=ggplot(data=otu.t1, aes(y=relab, x=compartment, fill=rank))+
  geom_bar(stat="identity",width = 0.95)+
  scale_fill_manual(values = my_pal(n))+
  labs(x = "compartment",y = "relative abundance (%)")+
  facet_grid(~treatment, scales = "free", space = "free_x")+tema+
  theme(legend.position = "none")
p1
dev.off()


#FIGURES with labels by plant species
pdf(paste(amplicon, t, file = "RA BY PS-T12.pdf"), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
p=ggplot(data=otu.t1, aes(y=relab, x=plant.compartment, fill=rank))+
  geom_bar(stat="identity",width = 0.95)+
  scale_fill_manual(values = my_pal(n))+
  labs(x = "Plant species",y = "relative abundance (%)")+
  facet_grid(~plant.species, scales = "free", space = "free_x")+
  guides(fill=guide_legend(title = "order", ncol = 1, keyheight = 0.7, default.unit="cm",
                           label.theme = element_text(size= 12, face = "plain", family = "sans")))+tema
p
dev.off()

#FIGURES without labels by plant species
pdf(paste(amplicon, t, file = "RA BY PS-T12-1.pdf"), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
p1=ggplot(data=otu.t1, aes(y=relab, x=plant.compartment, fill=rank))+
  geom_bar(stat="identity",width = 0.95)+
  scale_fill_manual(values = my_pal(n))+
  labs(x = "Plant species",y = "relative abundance (%)")+
  facet_grid(~plant.species, scales = "free", space = "free_x")+tema+
  theme(legend.position = "none")
p1
dev.off()



####SUBSET FOR Venn diagrams ONLY FUNGI OTUS IN SPORES#######
#At
ATspores= subset(otu.t1, otu.t1$plant.compartment == "spores", select = c("rank", "plant.compartment", "plant.species", "relab"))
ATspores= subset(ATspores, ATspores$plant.species =="Agave.tequilana", select = c("rank", "plant.compartment", "plant.species", "relab"))
ATspores=subset(ATspores, ATspores$relab >0.00)
ATspores= as.character(unique(ATspores$rank))
#As
ASspores= subset(otu.t1, otu.t1$plant.compartment == "spores", select = c("rank","plant.compartment", "plant.species","relab"))
ASspores= subset(ASspores, ASspores$plant.species =="Agave.salmiana", select = c("rank","plant.compartment", "plant.species","relab"))
ASspores=subset(ASspores, ASspores$relab >0.00)
ASspores= as.character(unique(ASspores$rank))
#Mg
Mgspores= subset(otu.t1, otu.t1$plant.compartment == "spores", select = c("rank","plant.compartment", "plant.species","relab"))
Mgspores= subset(Mgspores, Mgspores$plant.species =="Myrtillocactus.geometrizans", select = c("rank","plant.compartment", "plant.species", "relab"))
Mgspores=subset(Mgspores, Mgspores$relab >0.00)
Mgspores= as.character(unique(Mgspores$rank))
#
setlist <- list(AT=ATspores, AS=ASspores, MG=Mgspores)
#Venn diagrams of FUNGI OTUS in SPORES
pdf(paste(amplicon,t,file = "DV-OTUS IN SPORES.pdf"),colormodel = "cmyk", width = 11, height = 8.5, compress = F)
p=ggVennDiagram(setlist, category.names = c("A. tequilana","A. salmiana", "M. geometrizans"),
                label = "count", label_alpha = 0, edge_size = 1, label_size = 12, set_size = 9)+
  #scale_fill_distiller(palette = "Pastel1")+
  scale_fill_gradient2(low = "cornsilk3", mid = "grey", high = "dimgrey")+
  scale_color_manual(values = c("black", "black", "black"))+
  theme(legend.position = "none")
p
dev.off()


####SUBSET FOR Venn diagrams OF FUNGI OTUS IN ROOT ENDOSPHERE
#AT
ATRE= subset(otu.t1, otu.t1$plant.compartment == "root endosphere", select = c("rank","plant.compartment", "plant.species","relab"))
ATRE= subset(ATRE, ATRE$plant.species =="Agave.tequilana", select = c("rank","plant.compartment", "plant.species", "relab"))
ATRE=subset(ATRE, ATRE$relab >0.00)
ATRE= as.character(unique(ATRE$rank))
#AS
ASRE= subset(otu.t1, otu.t1$plant.compartment == "root endosphere", select = c("rank","plant.compartment", "plant.species","relab"))
ASRE= subset(ASRE, ASRE$plant.species =="Agave.salmiana", select = c("rank","plant.compartment", "plant.species", "relab"))
ASRE=subset(ASRE, ASRE$relab >0.00)
ASRE= as.character(unique(ASRE$rank))
#MG
MGRE= subset(otu.t1, otu.t1$plant.compartment == "root endosphere", select = c("rank","plant.compartment", "plant.species","relab"))
MGRE= subset(MGRE, MGRE$plant.species =="Myrtillocactus.geometrizans", select = c("rank","plant.compartment", "plant.species", "relab"))
MGRE=subset(MGRE, MGRE$relab >0.00)
MGRE= as.character(unique(MGRE$rank))

setlist <- list(AT=ATRE, AS=ASRE, MG=MGRE)

#Venn diagram of FUNGI OTUS in Root endosphere
pdf(paste(amplicon, t, file = "DV-OTUS IN ROOT ENDOSPHERE.pdf"), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
p=ggVennDiagram(setlist, category.names = c("A. tequilana","A. salmiana", "M. geometrizans"), 
                label = "count", label_alpha = 0, edge_size = 1, label_size = 12, set_size = 9)+
  #scale_fill_distiller(palette = "Pastel1")+
  scale_fill_gradient2(low = "cornsilk3", mid = "grey", high = "dimgrey")+
  scale_color_manual(values = c("black", "black", "black"))+
  theme(legend.position = "none")
p
dev.off()

####SUBSET fo Venn diagrams  BY PLANT COMPARTMENT##########
#SOIL
SOIL= subset(otu.t1, otu.t1$plant.compartment == "soil", select = c("rank", "plant.compartment", "plant.species", "relab"))
SOIL=subset(SOIL, SOIL$relab >0.00)
SOIL= as.character(unique(SOIL$rank))
#RE
RE= subset(otu.t1, otu.t1$plant.compartment == "root endosphere", select = c("rank", "plant.compartment", "plant.species", "relab"))
RE=subset(RE, RE$relab >0.00)
RE= as.character(unique(RE$rank))
#RHI
RHI= subset(otu.t1, otu.t1$plant.compartment == "rhizosphere", select = c("rank", "plant.compartment", "plant.species", "relab"))
RHI=subset(RHI, RHI$relab >0.00)
RHI= as.character(unique(RHI$rank))
#SPO
SPO= subset(otu.t1, otu.t1$plant.compartment == "spores", select = c("rank", "plant.compartment", "plant.species", "relab"))
SPO=subset(SPO, SPO$relab >0.00)
SPO= as.character(unique(SPO$rank))

setlist2 <- list(S=SOIL, RHI= RHI, RE=RE, SPO=SPO)

#Venn diagrams of FUNGI OTUS by plant compartment
pdf(paste(amplicon, t, file = "DV-OTUS IN PLANT COMPARTMENT.pdf"),  colormodel = "cmyk", width = 11, height = 8.5, compress = F)
p1=ggVennDiagram(setlist2, category.names = c("soil","rhizosphere", "root endosphere", "spores"), 
                 label = "count", label_alpha = 0, edge_size = 1, label_size = 12, set_size = 9)+
  #scale_fill_distiller(palette = "Pastel1")+
  scale_fill_gradient2(low = "cornsilk3", mid = "grey", high = "dimgrey")+
  scale_color_manual(values = c("black", "black", "black", "black"))+
  theme(legend.position = "none")
p1
dev.off()


##
####SUBSET for Venn diagrams by PLANT SPECIE########
#AT
AT= subset(otu.t1, otu.t1$plant.species == "Agave.tequilana", select = c("rank","plant.compartment", "plant.species","relab"))
AT=subset(AT, AT$relab >0.00)
AT= as.character(unique(AT$rank))
#AS
AS= subset(otu.t1, otu.t1$plant.species == "Agave.salmiana", select = c("rank","plant.compartment", "plant.species","relab"))
AS=subset(AS, AS$relab >0.00)
AS= as.character(unique(AS$rank))
#MG
MG= subset(otu.t1, otu.t1$plant.species == "Myrtillocactus.geometrizans", select = c("rank","plant.compartment", "plant.species","relab"))
MG=subset(MG, MG$relab >0.00)
MG= as.character(unique(MG$rank))

setlist3 <- list(AT, AS, MG)

#Venn diagrams of FUNGI OTUS in plant specie
pdf(paste(amplicon, t, file = "DV-OTUS IN PLANT SPECIE.pdf"),  colormodel = "cmyk", width = 11, height = 8.5, compress = F)
p2=ggVennDiagram(setlist3, category.names = c("A. tequilana","A. salmiana", "M. geometrizans"), 
                 label = "count", label_alpha = 0, edge_size = 1, label_size = 12, set_size = 9)+
  #scale_fill_distiller(palette = "Pastel1")+
  scale_fill_gradient2(low = "cornsilk3", mid = "grey", high = "dimgrey")+
  scale_color_manual(values = c("black", "black", "black"))+
  theme(legend.position = "none")
p2
dev.off()


#### GET UNIQUE OTUS from Spores
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/agave16s/")
g=VennDiagram::get.venn.partitions(setlist2)
g1=g$..values..$`8`
write.table(g1,
            file=paste(getwd(),paste("72OTUs-Spores",data,"txt",sep = "."),sep = "/"),
            quote = F, col.names = T, row.names = T, sep = "\t")
##424OTUs in Spores 
write.table(SPO,
            file=paste(getwd(),paste("OTUs-Spores",data,"txt",sep = "."),sep = "/"),
            quote = F, col.names = T, row.names = T, sep = "\t")

#spores=read.delim(file=paste(getwd(),paste("72OTUs-Spores",data,"txt",sep = "."),sep = "/"))
spores=read.delim(file=paste(getwd(),paste("OTUs-Spores",data,"txt",sep = "."),sep = "/"))

#spores=as.data.frame(SPO)
colnames(spores)= "OTUs-Spores"
rownames(spores)=spores$`OTUs-Spores`

#OBTAINF UNIQUE OTUS
library(Biostrings)
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/agave16s/16S")
#otueb=read.delim(file=paste(getwd(),paste("otu.table.EB",data,"txt",sep = "."),sep = "/"))
tax=read.delim(paste(getwd(),paste("taxa.table",data,base,"txt",sep = "."),sep = "/"))
#tax=read.delim(paste(getwd(),paste("taxa.BRE-MRE",data,"txt",sep = "."),sep = "/"))
#tax=tax[rownames(otueb),]
tax=tax[rownames(spores),]
tax1= subset(tax, tax$order == "o:Burkholderiales", select= c("domain","phylum","class","order","family","genus","otu.id"))
tax2= subset(tax, tax$order == "o:Mycoplasmatales", select= c("domain","phylum","class","order","family","genus","otu.id"))
tax3=rbind(tax1,tax2)
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/agave16s/")
write.table(tax3,
            file=paste(getwd(),paste("BRE-MRE",data,"txt",sep = "."),sep = "/"),
            quote = F, col.names = T, row.names = T, sep = "\t")
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/agave16s/16S")
SPO= readDNAStringSet("otus.daniel.fasta")
SPO=SPO[rownames(tax3),]
writeXStringSet(SPO, 
                file=paste(getwd(),paste("BRE-MRE","fasta",sep = "."),sep = "/"),
                format = "fasta")
SPO1=readDNAStringSet("OTUspores.fasta")
SPO2=readDNAStringSet("otus.CgMg.fasta")


write.table(met,
            file=paste(getwd(),paste("metadata-12",data,"txt",sep = "."),sep = "/"),
            quote = F, col.names = T, row.names = T, sep = "\t")



###SUBSET FOR OTUS in SPORES
EB= readDNAStringSet("otus.daniel.fasta")
EB=EB[rownames(tax),]
SPO=as.data.frame(SPO)
rownames(SPO)=SPO$SPO
tax=tax[rownames(SPO),]
otu=otu[rownames(SPO),]
writeXStringSet(EB, 
                file=paste(getwd(),paste("otus.AMF-SB","fasta",sep = "."),sep = "/"),
                format = "fasta")
write.table(otu,
            file=paste(getwd(),paste("otu.table.EB",data,"txt",sep = "."),sep = "/"),
            quote = F, col.names = T, row.names = T, sep = "\t")
