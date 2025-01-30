### Performs OTU enrichment analysis between treatments and plots ###
#Directories
setwd("~/abel/analysis/enrichment")
its=("~/abel/analysis/FITS_analy/result/") 
s16=("~/abel/analysis/16s_analy/result/")  

#Packages
library(stringr)
library(phyloseq)
library(ggplot2)
library(ggtern)
library(ggrepel)
library(vegan)
library(dplyr)
library(labdsv)
library(Rmisc)
library(tidyr)
library(FSA)
library(extrafont)


#Graphics
colores=c("black","#29854E","#74A557","#B8C466","#FFE081","#A02355","#C55DA1","#DB99EE","#234EA0","#008CCF","#00C6D9")
nivel=c("soil", "rhizosphere", "root.endosphere", "leaf.endosphere", "phyllosphere")
#etiquetas=c("Agave.tequilana", "Agave.marmorata", "Myrtillocactus.geometrizans")

tema=theme(axis.text.x = element_text(color="black",size=20, angle=0,hjust=0.5,vjust=1.5, family = "sans" ),
           #axis.text.x = element_text(color="black",size=12, angle=90,hjust=0.5,vjust=1.5),
           axis.text.y = element_text(color="black",size=20, vjust = 1.5, family = "sans"),
           axis.title = element_text(color="black",size=20, face = "bold", family = "sans"),
           axis.title.x.bottom = element_blank(),
           panel.border =element_rect(color = "black", fill = NA),#element_blank(),
           strip.text.x = element_text(size=20, color="black",face="bold", family = "sans"),
           strip.text.y = element_text(size=20, color="black",face="bold", family = "sans"),
           strip.placement = "outside", strip.background = element_rect(fill = "white"), 
           panel.background = element_rect(fill = "white",colour = "white",size = 0.8, linetype = "solid"),
           panel.grid.major.y = element_blank(),panel.grid.minor.x = element_blank(),
           legend.position = "right", legend.text = element_text(color = "black",size=14, family = "sans"), 
           legend.direction = "vertical", legend.title = element_text(color = "black",size=14, face = "bold", family = "sans"),
           legend.key.size = unit(0.4,"cm"))


# Set parameters
amplicon=c("ITS","16S")[1]
data=c("abel", "nestor", "all")[1]  #what sub sample would you need
base=c("uniteall","silva")[1]
t=c("domain","phylum","class","order","family","genus","otu.id")[4] #taxa level analysis
s=c("soil", "rhizosphere", "root.endosphere", "leaf.endosphere", "phyllosphere")[3] #Analyze only one treatment or all
#m=c("TO", "T6", "T12M","T12P") [2]
tr=c("MOCK","EV","PFC","PFCS","CAMP")[1]
test=c("kruskal")[1]

# Load metadata
#Upload data R Base
if (amplicon == "ITS"){
  setwd("~/abel/analysis/FITS_analy/result/") 
  tax=read.delim(paste(getwd(),paste("taxa.table",data,base,"txt",sep = "."),sep = "/"))
  otu=read.delim(file=paste(getwd(),paste("otu.table",data,"txt",sep = "."),sep = "/"))
  colnames(otu)=gsub("X","",colnames(otu))
  met=read.delim(file=paste(getwd(),paste("metadata.table",data,"txt",sep = "."),sep = "/"))
  #met=read.delim(file=paste(getwd(),paste("EV_metadata.table",data,"txt",sep = "."),sep = "/"))
  #met=read.delim(file=paste(getwd(),paste("PFC_metadata.table",data,"txt",sep = "."),sep = "/"))
  #met=read.delim(file=paste(getwd(),paste("PFCS_metadata.table",data,"txt",sep = "."),sep = "/"))
  rownames(met)=met$ID_seq_ITS
  met=met[order(met$ID_seq_ITS),]
  met=met[order(met$ID_seq_ITS),]
  met=met[colnames(otu),]
  # Removing unwanted samples
  met= subset(met, met$sampling != "T0")
  #met= subset(met, met$sampling != "T12M")
  #met= subset(met, met$sampling != "T12P")
  #met= subset(met, met$sampling != "T6")
  #met= subset(met, met$compartment != "rhizosphere")
  met= subset(met, met$compartment != "root.endosphere")
  met= subset(met, met$compartment != "leaf.endosphere")
  met= subset(met, met$compartment != "phyllosphere")
  #met= subset(met, met$compartment != "soil")
  #met= subset(met, met$treatment != "MOCK")
  met= subset(met, met$treatment != "EV")
  met= subset(met, met$treatment != "PFC")
  met= subset(met, met$treatment != "PFCS")
  met= subset(met, met$treatment != "CAMP")
  otu=otu[rownames(met)]
  # Sub sample reads 
  set.seed(23171341)
  otu.s=t(rrarefy(t(otu), sample = min(colSums(otu)))) 
  #Checklist
  colSums(otu.s)
  #otu.r=t(t(otu.s)/colSums(otu.s)*100)
  otu.r=cbind(otu.s,tax)
  otu.t=gather(data=otu.r, key = "sample", value = "abs", colnames(otu.s))
  for (i in met$ID_seq_ITS){
    otu.t[otu.t$sample == i, c("compartment","treatment","sampling","site", "ID_seq_ITS")] = met[i,c( "compartment","treatment","sampling","site", "ID_seq_ITS")]
    }
  
  otu.t[is.na(otu.t[,t]),t]="Unclassified"
  #Calculates relative abundance for each taxa (OTU) in each treatment
  otu.t0 = otu.t %>% group_by(treatment,compartment) %>% dplyr::summarise(absab = sum(abs))
  #Subset to eliminate unclassified
  otu.t= subset(otu.t, otu.t$class != "Unclassified")
  #otu.t= subset(otu.t, otu.t$family != "Unclassified")
  #contnue
  otu.t1 = otu.t %>% group_by(otu.t[,t],sample,sampling,compartment,treatment) %>% dplyr::summarise(absab = sum(abs)) #T12M YT12P
  names(otu.t1)[1]=c("rank")
  sum(paste(otu.t1$sampling,otu.t1$treatment,otu.t1$compartment,sep = ".")==paste(otu.t1$sampling,otu.t1$treatment,otu.t1$compartment,sep = "."))==dim(otu.t1)[1] #T12M YT12P
  #sum(paste(otu.t1$sample,otu.t1$treatment,otu.t1$compartment,sep = ".")==paste(otu.t1$sample,otu.t1$treatment,otu.t1$compartment,sep = "."))==dim(otu.t1)[1]
  otu.t1$relab=otu.t1$absab/otu.t0$absab*100
  
  #Subset the data 
  #otu.t1=otu.t1[otu.t1$sample %in% met[met$compartment %in% s, "ID_seq_ITS"],]
  otu.t1 <- subset(otu.t1, otu.t1$compartment=='soil' | otu.t1$compartment=='rhizosphere' & otu.t1$sampling=='T6') #to analyze soil samples
  #otu.t1 <- subset(otu.t1, otu.t1$treatment=='MOCK') #to analyze mock samples
  taxon=sort(unique(otu.t1$rank))
  

} else if ( amplicon == "16S"){
  setwd("~/abel/analysis/16s_analy/result/")
  tax=read.delim(paste(getwd(),paste("taxa.table",data,base,"txt",sep = "."),sep = "/"))
  otu=read.delim(file=paste(getwd(),paste("otu.table",data,"txt",sep = "."),sep = "/"))
  colnames(otu)=gsub("X","",colnames(otu))
  met=read.delim(file=paste(getwd(),paste("metadata.table",data,"txt",sep = "."),sep = "/"))
  #met=read.delim(file=paste(getwd(),paste("EV_metadata.table",data,"txt",sep = "."),sep = "/"))
  #met=read.delim(file=paste(getwd(),paste("PFC_metadata.table",data,"txt",sep = "."),sep = "/"))
  #met=read.delim(file=paste(getwd(),paste("PFCS_metadata.table",data,"txt",sep = "."),sep = "/"))
  #met=read.delim(file=paste(getwd(),paste("CAMP_metadata.table",data,"txt",sep = "."),sep = "/"))
  rownames(met)=met$ID_seq_16s
  met=met[order(met$ID_seq_16s),]
  met=met[order(met$ID_seq_16s),]
  met=met[colnames(otu),]
  # Removing unwanted samples
  #met= subset(met, met$sampling != "T0")
  #met= subset(met, met$sampling != "T12M")
  #met= subset(met, met$sampling != "T12P")
  #met= subset(met, met$sampling != "T6")
  #met= subset(met, met$compartment != "rhizosphere")
  #met= subset(met, met$compartment != "root.endosphere")
  #met= subset(met, met$compartment != "leaf.endosphere")
  #met= subset(met, met$compartment != "phyllosphere")
  #met= subset(met, met$compartment != "soil")
  #met= subset(met, met$treatment != "MOCK")
  #met= subset(met, met$treatment != "EV")
  #met= subset(met, met$treatment != "PFC")
  #met= subset(met, met$treatment != "PFCS")
  #met= subset(met, met$treatment != "CAMP")
  otu=otu[rownames(met)]
  #ONLY SPORES OTUS
 #otueb=read.delim(file=paste(getwd(),paste("otu.table.EB",data,"txt",sep = "."),sep = "/"))
  # Sub sample reads 
  set.seed(23171341)
  otu.s=t(rrarefy(t(otu), sample = min(colSums(otu)))) 
  #Checklist
  colSums(otu.s)
  #otu.r=t(t(otu.s)/colSums(otu.s)*100)
  otu.r=cbind(otu.s,tax)
  otu.t=gather(data=otu.r, key = "sample", value = "abs", colnames(otu.s))
  for (i in met$ID_seq_16s){
    otu.t[otu.t$sample == i, c("compartment","treatment","sampling","site", "ID_seq_16s")] = met[i,c( "compartment","treatment","sampling","site", "ID_seq_16s")]
  }
  otu.t[is.na(otu.t[,t]),t]="Unclassified"
  #Calculates relative abundance for each taxa (OTU) in each treatment
  otu.t0 = otu.t %>% group_by(sample,sampling,compartment,treatment) %>% dplyr::summarise(absab = sum(abs))
  #otu.t0 = otu.t %>% group_by(sample,compartment, treatment) %>% dplyr::summarise(absab = sum(abs))
    #Subset to eliminate unclassified
  #otu.t= subset(otu.t, otu.t$class != "Unclassified")
  otu.t= subset(otu.t, otu.t$order != "Unclassified")
  #contnue
  otu.t1 = otu.t %>% group_by(otu.t[,t],sample,sampling,compartment,treatment) %>% dplyr::summarise(absab = sum(abs)) #T12M YT12P
  #otu.t1 = otu.t %>% group_by(otu.t[,t],sample,compartment,treatment) %>% dplyr::summarise(absab = sum(abs))
  names(otu.t1)[1]=c("rank")
  sum(paste(otu.t1$sampling,otu.t1$compartment,otu.t1$treatment,sep = ".")==paste(otu.t1$sampling,otu.t1$compartment,otu.t1$treatment,sep = "."))==dim(otu.t1)[1] #T12M YT12P
  #sum(paste(otu.t1$compartment,otu.t1$treatment,sep = ".")==paste(otu.t1$compartment,otu.t1$treatment,sep = "."))==dim(otu.t1)[1]
  otu.t1$relab=otu.t1$absab/otu.t0$absab*100
  ####################
  #Subset the data 
  #otu.t1=otu.t1[otu.t1$sample %in% met[met$compartment %in% s, "ID_seq_16s"],]
  otu.t1 <- subset(otu.t1, otu.t1$compartment=='soil' | otu.t1$compartment=='rhizosphere' & otu.t1$sampling=='T6') #to analyze soil samples
  #otu.t1 <- subset(otu.t1, otu.t1$treatment=='EV') #to analyze mock samples
  #################
  taxon=sort(unique(otu.t1$rank))
}

#Test
taxon=as.character(sort(unique(otu.t1$rank)))

if(test=="kruskal"){
kruskal=data.frame(row.names = taxon, rank=taxon, p.value=rep(1, length(taxon)[1]))

  for (id in taxon){
    sub1=otu.t1[otu.t1$rank==id,]
    #kruskal[id,"p.value"]=kruskal.test(relab ~ treatment, data = sub1)$p.value
    kruskal[id,"p.value"]=kruskal.test(relab ~ sampling, data = sub1)$p.value
    kruskal$p.value=format(kruskal$p.value, scientific = F)
  }

  if (sum(kruskal$p.value<=0.05, na.rm = T)>0){
    
    dunn=data.frame(row.names = taxon)
    
    for (h in taxon){
      sub1=otu.t1[otu.t1$rank==h,]
    
      if (sum(sub1$relab)>0){
        #test1=dunnTest(relab ~ treatment, data = sub1, method = "bh")
        test1=dunnTest(relab ~ sampling, data = sub1, method = "bh")
        p.value=test1$res$P.adj
        #dunn[h,1:choose(length(unique(sub1$treatment)),2)]=p.value
        dunn[h,1:choose(length(unique(sub1$sampling)),2)]=p.value
      } else if (sum(sub1$relab)==0){
        #p.value=rep(1,choose(length(unique(sub1$treatment)),2))
        p.value=rep(1,choose(length(unique(sub1$sampling)),2))
        #dunn[h,1:choose(length(unique(sub1$treatment)),2)]=p.value
        dunn[h,1:choose(length(unique(sub1$sampling)),2)]=p.value
      }
  }
}
      
  colnames(dunn)=test1$res$Comparison
  dunn[is.na(dunn)]=1
  pv=dunn
  #pv=select(pv,"EV - MOCK","EV - PFC","MOCK - PFC","EV - PFCS","MOCK - PFCS","PFC - PFCS")
  #pv=select(pv,"CAMP - EV","CAMP - MOCK","EV - MOCK","CAMP - PFC","EV - PFC","MOCK - PFC",
   #         "CAMP - PFCS","EV - PFCS","MOCK - PFCS","PFC - PFCS") #treatment filosfera
  #pv=select(pv,"CAMP - MOCK","EV - MOCK","MOCK - PFC","MOCK - PFCS") #treatment filosfera
#pv=select(pv,"EV - MOCK","MOCK - PFC","MOCK - PFCS") #TREATMENTS
  #pv=select(pv, "T0 - T6","T0 - T12M","T0 - T12P") #sampling
  pv=select(pv, "T12M - T12P","T12M - T6") #sampling soil
  #pv=dunn[,grep("Agave.",colnames(dunn))]
  #colnames(pv)=gsub("T0 - ","",colnames(pv)) #sampling T0-T12
  colnames(pv)=gsub("T12M - ","",colnames(pv)) #sampling soil
#colnames(pv)=gsub(" - MOCK","",colnames(pv))
#colnames(pv)=gsub("MOCK - ","",colnames(pv))
  
  for (h in colnames(pv)){
    print(h)
    print(sum(pv[,h]<=0.05))
  }
}

#Calculate log fold for plotting 
#otu.d=summarySE(otu.t1, measurevar = "relab", groupvars = c("rank","treatment")) #mean abundance of each OTU
otu.d=summarySE(otu.t1, measurevar = "relab", groupvars = c("rank","sampling")) #mean abundance of each OTU
en=data.frame(row.names=taxon)
for (i in colnames(pv)){
  #logfold=log2(otu.d[otu.d$treatment==i,"relab"]/otu.d[otu.d$treatment=="MOCK","relab"])
  #logfold=log2(otu.d[otu.d$sampling==i,"relab"]/otu.d[otu.d$sampling=="T0","relab"])
  logfold=log2(otu.d[otu.d$sampling==i,"relab"]/otu.d[otu.d$sampling=="T12M","relab"])
  
  en=cbind(en, logfold)
}
colnames(en)=colnames(pv)
en=en[rownames(pv), colnames(pv)]
en[is.infinite(as.matrix(en))]=0
en[is.nan(as.matrix(en))]=0

#Select any colors 
colores=c("#A74476","#F56AB0","#804795","#DF87FF",
          "#524CA1","#9A94EF","#2C79AE","#6EC4FF","#369DA4","#6EF6FF",
          "#29A876","#43FFB5","#55A054","#A5F77A","#A5A53F","#F1F18A",
          "#AA723D","#FFBB7A","#A12F2F","#FF4848", "#A74476","#F56AB0","#804795","#DF87FF",
          "#524CA1","#9A94EF","#2C79AE","#6EC4FF","#369DA4","#6EF6FF",
          "#29A876","#43FFB5","#55A054","#A5F77A","#A5A53F","#F1F18A",
          "#AA723D","#FFBB7A","#A12F2F","#FF4848", "#A74476","#F56AB0","#804795","#DF87FF",
          "#524CA1","#9A94EF","#2C79AE","#6EC4FF","#369DA4","#6EF6FF",
          "#29A876","#43FFB5","#55A054","#A5F77A","#A5A53F","#F1F18A",
          "#AA723D","#FFBB7A","#A12F2F","#FF4848", "#A74476","#F56AB0","#804795","#DF87FF",
          "#524CA1","#9A94EF","#2C79AE","#6EC4FF","#369DA4","#6EF6FF",
          "#29A876","#43FFB5","#55A054","#A5F77A","#A5A53F","#F1F18A",
          "#AA723D","#FFBB7A","#A12F2F","#FF4848", "#F56AB0")

#Select the folder to save the figures
if (amplicon == "ITS"){
  setwd("~/abel/analysis/enrichment/FITS/")
} else if ( amplicon == "16S"){
  setwd("~/abel/analysis/enrichment/16s/")
}

# Generates a plot per each comparison 

for (w in colnames(pv)){
  
  # threshold for enrichment
  #w="EV"
  w="T12P"
  #w="T6"
  thf=1
  #thf=2.5
  thp=0.05
  lab="order"

  #ab=otu.d[otu.d$plant.species == "Agave.salmiana",]
  #ab=otu.d[otu.d$treatment == w,]
  ab=otu.d[otu.d$sampling == w,]
  plot.data=data.frame(row.names = rownames(pv), pvalue=pv[,w], logf=en[,w], relab=ab$relab)

  
  for (i in rownames(plot.data)){
    plot.data[i, c("domain","phylum","class","order","family","genus","otu.id")]=unique(tax[tax[,t] %in% i, c("domain","phylum","class","order","family","genus","otu.id")])}
  
  
sig=plot.data[plot.data$pvalue > 0.05, "order"]
  
 
  #plot.data$family="differential"
  #plot.data[plot.data$order %in% sig,"order"] ="no differencial"
  plot.data[plot.data$order %in% sig,"order"] =""
  #plot.data$order=factor(plot.data$order, levels = c(niveles,"Unclassified","Low.abundant"))
  niveles=unique(plot.data[plot.data$order %in% plot.data$order, "order"])
  col=data.frame(color=c(colores[1:length(niveles)],"grey50"),taxa=c(niveles,"no differencial"))

  plot.data$label=plot.data[,lab]
  #plot.data[plot.data$pvalue>thp | abs(plot.data$logf)<thf,"label"]="_"
  
}

plot.data$treatment <- "MOCK"
plot.data$sampling <- w
write.csv(plot.data, file = "soil plot.data.csv", row.names = FALSE)

library(dplyr)
tabla1 <- read.csv("T6-MOCK plot.data.csv")  # Reemplaza con la ruta de tu primer archivo CSV
tabla2 <- read.csv("T6-PFCS plot.data.csv")  # Reemplaza con la ruta de tu segundo archivo CSV
tabla3 <- read.csv("T12P-MOCK plot.data.csv")  # Reemplaza con la ruta de tu primer archivo CSV
tabla4 <- read.csv("T12P-PFC plot.data.csv")  # Reemplaza con la ruta de tu segundo archivo CSV

# Juntar las dos tablas usando rbind (asumiendo que las columnas tienen el mismo nombre y orden)
plot.data.1 <- rbind(tabla1, tabla2,tabla3, tabla4)
sig=plot.data.1[plot.data.1$pvalue > 0.05, "order"]
plot.data.1[plot.data.1$order %in% sig,"order"] =""
#plot.data$order=factor(plot.data$order, levels = c(niveles,"Unclassified","Low.abundant"))
niveles=unique(plot.data.1[plot.data.1$order %in% plot.data.1$order, "order"])
col=data.frame(color=c(colores[1:length(niveles)],"grey50"),taxa=c(niveles,"no differencial"))

pdf(paste(amplicon,s,tr,t,w, file = "V4_enrichment.pdf"),  colormodel = "srgb", width = 15, height = 16, compress = F)
p=ggplot(data=plot.data.1, aes(x=logf, y=-log10(pvalue), color=order))+
  facet_wrap(~sampling~treatment)+
    geom_point(size=4)+
    scale_color_manual(values=c(col[col$taxa %in% unique(plot.data.1$order),"color"]))+
    geom_text_repel(aes(label=label),size=5, color="black", max.overlaps = 102, force_pull = 2, force = 2)+
    geom_hline(yintercept = -log10(thp), linetype=2)+tema+
    #geom_vline(xintercept = c(-thf,thf), linetype=3)+
    ggtitle(paste("T0","vs",w))+
    #ggtitle(paste(w,"vs", "MOCK"))+
    scale_x_continuous(limits = c(-9,6), breaks = c(seq(-9,6,3)))
p
dev.off()

pdf(paste(amplicon,s,tr,t,w, file = "enrichment-V4.pdf"),  colormodel = "srgb", width = 15, height = 16, compress = F)
p=ggplot(data=plot.data.1, aes(x=logf, y=-log10(pvalue), color=order))+
  facet_wrap(~sampling~treatment)+
  geom_point(size=4)+
  scale_color_manual(values=c(col[col$taxa %in% unique(plot.data.1$order),"color"]))+
  geom_text_repel(aes(label=label),size=5, color="black", max.overlaps = 42, force_pull = 2, force = 2)+
  geom_hline(yintercept = -log10(thp), linetype=3)+tema+
  #geom_vline(xintercept = c(-thf,thf), linetype=3)+
  ggtitle(paste("T0","vs",w))+
  scale_x_continuous(limits = c(-9,6), breaks = c(seq(-9,6,3)))+
  theme(legend.position = "none")
p
dev.off()

################SOIL      ############################
pdf(paste(amplicon,"phyl-mock",tr,t,w, file = "V4_enrichment.pdf"),  colormodel = "srgb", width = 15, height = 8.5, compress = F)
p=ggplot(data=plot.data, aes(x=logf, y=-log10(pvalue), color=order))+
  #facet_wrap(~sampling~treatment)+
  geom_point(size=4)+
  scale_color_manual(values=c(col[col$taxa %in% unique(plot.data$order),"color"]))+
  geom_text_repel(aes(label=label),size=5, color="black", max.overlaps = 102, force_pull = 2, force = 2)+
  geom_hline(yintercept = -log10(thp), linetype=2)+tema+
  #geom_vline(xintercept = c(-thf,thf), linetype=3)+
  ggtitle(paste("T12M","vs",w))+
  #ggtitle(paste(w,"vs", "MOCK"))+
  scale_x_continuous(limits = c(-9,6), breaks = c(seq(-9,6,3)))
p
dev.off()

pdf(paste(amplicon,"soil",tr,t,w, file = "enrichment-V4.pdf"),  colormodel = "srgb", width = 15, height = 8.5, compress = F)
p=ggplot(data=plot.data, aes(x=logf, y=-log10(pvalue), color=order))+
  #facet_wrap(~sampling~treatment)+
  geom_point(size=4)+
  scale_color_manual(values=c(col[col$taxa %in% unique(plot.data$order),"color"]))+
  geom_text_repel(aes(label=label),size=5, color="black", max.overlaps = 42, force_pull = 2, force = 2)+
  geom_hline(yintercept = -log10(thp), linetype=3)+tema+
  #geom_vline(xintercept = c(-thf,thf), linetype=3)+
  ggtitle(paste("T12M","vs",w))+
  scale_x_continuous(limits = c(-9,6), breaks = c(seq(-9,6,3)))+
  theme(legend.position = "none")
p
dev.off()
