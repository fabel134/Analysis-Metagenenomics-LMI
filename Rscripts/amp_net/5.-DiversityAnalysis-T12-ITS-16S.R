#Load all necessary libraries
library(stringr)
library(ggplot2)
library(vegan)
library(labdsv)
library(MASS)
library(Rmisc)
library(phyloseq)
library(tidyr)
library(gridExtra)
library(dplyr)
library(FSA)

# Choose all setting to plots
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
pc=c("soil","rhizosphere", "root.endosphere", "leaf.endosphere","phyllosphere", "T0", "T6", "T12M","T12P", 
     "MOCK", "EV","PFC","PFCS","CAMP", "uno","dos","tres","malla.sombra","mono", "poli")

#Upload data R Base
if (amplicon == "ITS2"){
  setwd("~/abel/analysis/FITS_analy/result/") 
  tax=read.delim(paste(getwd(),paste("taxa.table",data,base,"txt",sep = "."),sep = "/"))
  otu=read.delim(file=paste(getwd(),paste("otu.table",data,"txt",sep = "."),sep = "/"))
  colnames(otu)=gsub("X","",colnames(otu))
  met=read.delim(file=paste(getwd(),paste("metadata.table",data,"txt",sep = "."),sep = "/"))
  rownames(met)=met$ID_seq_ITS
  met=met[order(met$ID_seq_ITS),]
  met=met[order(met$ID_seq_ITS),]
  met=met[colnames(otu),]
  #Subset of data to process only at 12 months and soil
  #met= subset(met, met$sampling != "T0")
  #met= subset(met, met$sampling != "T12")
  #met= subset(met, met$sampling != "T6")
  #met= subset(met, met$compartment != "rhizosphere")
  #met= subset(met, met$compartment != "root.endosphere")
  #met= subset(met, met$compartment != "leaf.endosphere")
  #met= subset(met, met$compartment != "phyllosphere")
  met= subset(met, met$compartment != "soil")
  #met= subset(met, met$treatment != "MOCK")
  #met= subset(met, met$treatment != "EV")
  #met= subset(met, met$treatment != "PFC")
  #met= subset(met, met$treatment != "PFCS")
  #met= subset(met, met$treatment != "CAMP")
  otu=otu[rownames(met)]
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
  #met= subset(met, met$treatment != "T4")
  #met= subset(met, met$sampling != "T6")
  #met= subset(met, met$sampling != "T0")
  #met= subset(met, met$sampling != "T12")
  #met= subset(met, met$compartment != "rhizosphere")
  #met= subset(met, met$compartment != "root.endosphere")
  #met= subset(met, met$compartment != "leaf.endosphere")
  met= subset(met, met$compartment != "soil")
  #met= subset(met, met$compartment != "phyllosphere")
  #met= subset(met, met$treatment != "MOCK")
  #met= subset(met, met$treatment != "EV")
  #met= subset(met, met$treatment != "PFC")
  #met= subset(met, met$treatment != "PFCS")
  #met= subset(met, met$treatment != "CAMP")
  
  otu=otu[rownames(met)]
}

com=c("soil", "phyllosphere", "rhizosphere", 
      "leaf.endosphere", "root.endosphere","all")[6]

#com=c("EV", "PFC", "PFCS","CAMP", "MOCK")[4]
#Rarefy otu table
set.seed(23171341)
#all fungi
ric=(rarefy(t(otu), sample = min(colSums(otu))))
#all fungi
otu.s=t(rrarefy(t(otu), sample = min(colSums(otu))))
#Checklist
colSums(otu.s)
shan=diversity(t(otu.s),index="shannon")
#Generate plot data
alpha=data.frame(richness=ric, shannon=shan)
alpha=alpha[met$ID_seq_ITS,]
#alpha=alpha[met$ID_seq_16s,]
alpha=cbind(alpha, met[,c( "compartment","sampling", "ID_seq_16s", "treatment","site")])
alpha$sampling=factor(alpha$sampling, level = pc)
alpha$compartment=factor(alpha$compartment, level = pc)
alpha$treatment=factor(alpha$treatment, level = pc)
alpha$site=factor(alpha$site, level = pc)

#KW by sampling
kruskal.test(richness~sampling, data = alpha) ; dunnTest(richness~sampling, data = alpha, method = "bh")
kruskal.test(shannon~sampling, data = alpha) ; dunnTest(shannon~sampling, data = alpha, method = "bh")

#KW by compartment
kruskal.test(richness~compartment, data = alpha) ; dunnTest(richness~compartment, data = alpha, method = "bh")
kruskal.test(shannon~compartment, data = alpha) ; dunnTest(shannon~compartment, data = alpha, method = "bh")

#KW by treatment
kruskal.test(richness~treatment, data = alpha) ; dunnTest(richness~treatment, data = alpha, method = "bh")
kruskal.test(shannon~treatment, data = alpha) ; dunnTest(shannon~treatment, data = alpha, method = "bh")

#KW by site
#kruskal.test(richness~site, data = alpha) ; dunnTest(richness~site, data = alpha, method = "bh")
#kruskal.test(shannon~site, data = alpha) ; dunnTest(shannon~site, data = alpha, method = "bh")

#ANOVA by sampling
test=aov(richness~sampling, data = alpha) ; print(test) ; TukeyHSD(test)
test=aov(shannon~sampling, data = alpha) ; print(test) ; TukeyHSD(test)

#ANOVA by compartment
test=aov(richness~compartment, data = alpha) ; print(test) ; TukeyHSD(test)
test=aov(shannon~compartment, data = alpha) ; print(test) ; TukeyHSD(test)

#ANOVA by site
#test=aov(richness~site, data = alpha) ; print(test) ; TukeyHSD(test)
#test=aov(shannon~site, data = alpha) ; print(test) ; TukeyHSD(test)

#ANOVA by treatment
#test=aov(richness~treatment, data = alpha) ; print(test) ; TukeyHSD(test)
#test=aov(shannon~treatment, data = alpha) ; print(test) ; TukeyHSD(test)

#Select the folder to save the figures
if (amplicon == "ITS2"){
  setwd("~/abel/analysis/FITS_analy/result/FIGURAS-ITS/")
} else if ( amplicon == "16S"){
  setwd("~/abel/analysis/16s_analy/result/FIGURAS-16S/")
} 
#dataframe by build mean of richhnes and shannon boxplot
mean.r=summarySE(data = alpha, measurevar = "richness", groupvars = c("compartment","sampling","ID_seq_16s"))
#mean.r=summarySE(data = alpha, measurevar = "richness", groupvars = c("compartment","sampling"))
#mean.r$treatment=factor(mean.r$treatment, level = pc)
mean.r$compartment=factor(mean.r$compartment, level = pc)
mean.s=summarySE(data = alpha, measurevar = "shannon", groupvars = c("compartment","sampling","ID_seq_16s"))
#mean.s=summarySE(data = alpha, measurevar = "shannon", groupvars = c("compartment","sampling"))
#mean.s$treatment=factor(mean.s$treatment, level = pc)
mean.s$compartment=factor(mean.s$compartment, level = pc)

#Analisis estadisticos diversidad  
#RICHENESS BY 
pdf(paste(amplicon,"MOCK-",com,file = "sampling-richness.pdf"), colormodel = "cmyk", width = 15, height =8.5, compress = F)
p=ggplot(data=mean.r, aes(x=sampling, y=richness, fill=sampling))+
  #geom_bar(stat="identity", colour="black",size=1)+
  #geom_boxplot(outlier.size = 0, aes(fill=treatment))+
  geom_boxplot(outlier.size = 0)+ ##delete outliers
  #geom_errorbar(aes(ymin=richness-sd, ymax=richness+sd),width=0.2,position=position_dodge(1.5))+
  geom_jitter(aes(colour = sampling),shape=16, position=position_jitter(0.2))+
  #geom_dotplot(data=alpha,binaxis='y', stackdir='center', stackratio=2, dotsize=0.8, aes(fill=treatment))+
  scale_fill_manual(values = my_pal(9)[3:8])+
  scale_color_manual(values = c("T0"="gray80","T6"="gray60","T12M"="gray40","T12P"="gray30"))+
  facet_wrap(~ compartment, nrow = 1)+tema
p
dev.off()

#Shannon BY site and soil
pdf(paste(amplicon,"MOCK-",com,file= "sampling.shannon.pdf"), colormodel = "cmyk", width = 15, height =8.5, compress = F)
p1=ggplot(data=mean.s, aes(x=sampling, y=shannon, fill = sampling))+
  #geom_bar(stat="identity", colour="black",size=1)+
  #geom_boxplot(outlier.size = 0, aes(fill=treatment))+
  geom_boxplot(outlier.size = 0)+ ##delete outliers
  #geom_errorbar(aes(ymin=richness-sd, ymax=richness+sd),width=0.2,position=position_dodge(1.5))+
  geom_jitter(aes(colour = sampling),shape=16, position=position_jitter(0.2))+
  #geom_dotplot(data=alpha,binaxis='y', stackdir='center', stackratio=2, dotsize=0.8, aes(fill=treatment))+
  scale_fill_manual(values = my_pal(9)[3:8])+
  scale_color_manual(values = c("T0"="gray80","T6"="gray60","T12M"="gray40","T12P"="gray30"))+
  facet_wrap(~compartment, nrow = 1)+tema
p1
dev.off()

#dataframe by build mean of richhnes and shannon barplot
mean.r=summarySE(data = alpha, measurevar = "richness", groupvars = c("compartment","sampling"))
mean.r$compartment=factor(mean.r$compartment, level = pc)
mean.s=summarySE(data = alpha, measurevar = "shannon", groupvars = c("compartment","sampling"))
mean.s$compartment=factor(mean.s$compartment, level = pc)

#RICHENESS BY site and treantment
pdf(paste(amplicon, "Bar-MOCK-",com,file = ".richeness.pdf"), colormodel = "cmyk", width = 15, height =8.5, compress = F)
p=ggplot(data=mean.r, aes(x=sampling, y=richness, fill = sampling ))+
  geom_bar(stat="identity", colour="black",size=1)+
  facet_wrap(~compartment, nrow = 1)+
  geom_errorbar(aes(ymin=richness-sd, ymax=richness+sd),width=0.2,position=position_dodge(1.5))+
  scale_fill_manual(values = my_pal(9)[4:8])+tema
p
dev.off()

#Shannon BY site and treantment
pdf(paste(amplicon, "Bar-MOCK-",com,file = ".shannon.pdf"), colormodel = "cmyk", width = 15, height =8.5, compress = F)
p1=ggplot(data=mean.s, aes(x=sampling, y=shannon, fill = sampling))+
  geom_bar(stat="identity", colour="black",size=1)+
  #geom_bar(stat="identity", fill=c("gray100","gray80","gray70","gray60","gray40","gray30"), colour="black",size=1)+
  facet_wrap(~compartment, nrow = 1)+
  geom_errorbar(aes(ymin=shannon-sd, ymax=shannon+sd),width=0.2,position=position_dodge(1.5))+
    #geom_dotplot(data=alpha,binaxis='y', stackdir='center', stackratio=2, dotsize=0.8, aes(fill=treatment))+
  scale_fill_manual(values = my_pal(9)[4:8])+tema
p1
dev.off()

######################
#Richness BY PLANT.treatment
pdf(paste(amplicon, file = "site.richnessplot.pdf"), colormodel = "cmyk", width = 8.5, height = 8.5, compress = F)
p=ggplot(data=mean.r, aes(x=site, y=richness))+
  geom_bar(stat="identity",fill=c("gray100","gray80"),colour="black",size=1)+
  #geom_bar(stat="identity",fill=c("gray100","gray80"),colour="black",size=1)+
  geom_errorbar(aes(ymin=shannon-sd, ymax=shannon+sd),width=0.2,position=position_dodge(1.5))+
  #geom_dotplot(data=alpha,binaxis='y', stackdir='center', stackratio=2, dotsize=0.8, aes(fill=treatment))+
  scale_fill_manual(values = my_pal(9)[1:4])+tema
p
dev.off()

#SHANNON BY PLANT.treatment
pdf(paste(amplicon, file = "site.shannonplot.pdf"), colormodel = "cmyk", width = 8.5, height = 8.5, compress = F)
p1=ggplot(data=mean.s, aes(x=site, y=shannon))+
  geom_bar(stat="identity",fill=c("gray100","gray80"),colour="black",size=1)+
  #geom_bar(stat="identity",fill=c("gray100","gray80"),colour="black",size=1)+
  geom_errorbar(aes(ymin=shannon-sd, ymax=shannon+sd),width=0.2,position=position_dodge(1.5))+
  #geom_dotplot(data=alpha,binaxis='y', stackdir='center', stackratio=2, dotsize=0.8, aes(fill=treatment))+
  scale_fill_manual(values = my_pal(9)[1:4])+tema
p1
dev.off()
 ##PLOTS_diversidad

####NMDS#####

#Beta diversity
otu.n=as.data.frame(t(t(otu.s)/colSums(otu.s)*100)) #relative abundance
otu.n=log10(otu.n+1)
otu.n=sqrt(otu.n)

#Select the folder to save the figures
if (amplicon == "ITS2"){
  setwd("~/abel/analysis/FITS_analy/result/FIGURAS-ITS/")
} else if ( amplicon == "16S"){
  setwd("~/abel/analysis/16s_analy/result/FIGURAS-16S/")
} 
########## generate table for TDA analysis############
#tax <- apply(tax, 1, function(row) paste(row, collapse = "_")) #abel
#tax <- as.data.frame(tax) #abel
#otu.n=as.data.frame(cbind(otu.n,tax))
#write.table(otu.n,
#            file=paste(getwd(),paste("otu.table",data,"_tmap","txt",sep = "."),sep = "/"),
#            quote = F, col.names = T, row.names = T, sep = "\t") #abel
###################################################################

#calculate dimensions
scaling=vegdist(t(otu.n), method = "bray", binary = T, na.rm = FALSE) #calculate distance
#scaling2=isoMDS(scaling) ; scaling2$stress #create NMDS 
#scaling2=metaMDS(otu.n, distance = "bray") ; scaling2$stress
scaling2=monoMDS(scaling); scaling2$stress
scaling3=data.frame(scaling2$points) #select cordenates
scaling3=cbind(scaling3,alpha)

#plot
#plant.treatment
pdf(paste(amplicon,com,file = "NMDS.comparment.pdf"), colormodel = "srgb", width = 11, height = 8.5, compress = F)
n=ggplot(data=scaling3, aes(x=MDS1, y=MDS2, colour=compartment, shape=sampling))+
  geom_hline(yintercept = 0, linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_point(size=9)+
  #scale_colour_manual(values = my_pal(9)[2:5])+
  scale_colour_manual(values = c("#000000","#EB73B3","#21B087","#2CB5C0","#F64971"))+
  scale_alpha_continuous(range=c(0.4,1))+
  scale_shape_manual(values=c(15:20))+
  labs(x = "NMDS1",y = "NMDS2")+tema
n
dev.off()

#sampling
pdf(paste(amplicon,file = "NMDS.treatment-1.pdf"), colormodel = "srgb", width = 11, height = 8.5, compress = F)
m=ggplot(data=scaling3, aes(x=MDS1, y=MDS2, colour=compartment, shape=sampling))+
  geom_hline(yintercept = 0, linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_point(size=9)+
  #scale_colour_manual(values = my_pal(9)[2:5])+
  scale_colour_manual(values = c("#000000","#EB73B3","#21B087","#2CB5C0","#F64971"))+
  scale_alpha_continuous(range=c(0.4,1))+
  scale_shape_manual(values=c(15:18))+
  labs(x = "NMDS1",y = "NMDS2")+tema+
  theme(legend.position = "none")
print(m)
dev.off()

#site
pdf(paste(amplicon,com,file = "NMDS.sampling.pdf"), colormodel = "srgb", width = 11, height = 8.5, compress = F)
a=ggplot(data=scaling3, aes(x=MDS1, y=MDS2, colour=sampling, shape=compartment))+
  geom_hline(yintercept = 0, linetype=2)+
  geom_vline(xintercept = 0,linetype=2)+
  geom_point(size=9)+
  #scale_colour_manual(values = my_pal(9)[2:5])+
  scale_colour_manual(values = c("#000000","#EB73B3","#21B087","#2CB5C0","#F64971"))+
  scale_alpha_continuous(range=c(0.4,1))+
  scale_shape_manual(values=c(15:20))+
  labs(x = "NMDS1",y = "NMDS2")+tema
print(a)
dev.off()

mod = with(met, betadisper(scaling, sampling))
#mod = with(met, betadisper(scaling, site))
alpha$distance_to_centroid=mod$distances
mean.c=summarySE(data = alpha, measurevar = "distance_to_centroid", groupvars = "sampling")
#mean.c=summarySE(data = alpha, measurevar = "distance_to_centroid", groupvars = "site")
mean.c$sampling=factor(mean.c$sampling, level = pc)
#mean.c$site=factor(mean.c$site, level = pc)

######
pdf(paste(amplicon, com,file = "treatment.distances.pdf"),colormodel = "srgb", width = 8.5, height = 8.5, compress = F)
d=ggplot(data=mean.c, aes(x=sampling, y=distance_to_centroid, fill=sampling))+
  geom_bar(stat="identity",colour="black",size=1, fill=c("#000000","#EB73B3","#21B087","#2CB5C0"))+
  geom_errorbar(size=1,aes(ymin=distance_to_centroid-sd, ymax=distance_to_centroid+sd),width=0.2,position=position_dodge(1.5))+
  scale_fill_manual(values = my_pal(12)[c(1,3,5,7)])+
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2))+tema
d
dev.off()

TukeyHSD(mod)
#KW by distance_to_centroid
kruskal.test(distance_to_centroid~treatment, data = alpha) ; dunnTest(distance_to_centroid~treatment, data = alpha, method = "bh")
kruskal.test(distance_to_centroid~site, data = alpha) ; dunnTest(distance_to_centroid~treatment, data = alpha, method = "bh")
kruskal.test(distance_to_centroid~sampling, data = alpha) ; dunnTest(distance_to_centroid~sampling, data = alpha, method = "bh")

plot(mod)

###PERMANOVA ANALYSIS
#Beta diversity
otu.n=as.data.frame(t(t(otu.s)/colSums(otu.s)*100)) #relative abundance
otu.n=log10(otu.n+1)
otu.n=sqrt(otu.n)
#
p=t(otu.n)
set.seed(173612)
scaling=vegdist(p, method = "bray", binary = T) 
adonis2(scaling~treatment, data = met, permutations = 1000)
adonis2(scaling~sampling, data = met, permutations = 1000)
adonis2(scaling~compartment, data = met, permutations = 1000)
adonis2(scaling~sampling*compartment*treatment, data = met, permutations = 1000)

