#Combine script
{# Definir paquetes
  # Repositorio CRAN
  cran_packages <- c("bookdown", "knitr", "devtools", "igraph", "qgraph", "RColorBrewer", "tidyverse", "network", "sna", "grid", "gridExtra")
  # Repositorio Bioconductor
  bioc_packages <- c("phyloseq", "microbiome", "genefilter")
  # Repositorio GitHub
  git_source <- c("zdk123/SpiecEasi", "hallucigenia-sparsa/seqtime", "briatte/ggnet") # fuente/nombre del paquete
  git_packages <- c("SpiecEasi", "seqtime", "ggnet", "pulsar") # nombre del paquete
  
  # Instalar paquetes CRAN
  .inst <- cran_packages %in% installed.packages()
  if(any(!.inst)) {
    install.packages(cran_packages[!.inst])
  }
  # Intalar paquetes BioConductor
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  .inst <- bioc_packages %in% installed.packages()
  if(any(!.inst)) {
    BiocManager::install(bioc_packages[!.inst])
  }
  # Instalar paquetes GitHub
  .inst <- git_source %in% installed.packages()
  if(any(!.inst)) {
    devtools::install_github(git_source[!.inst])
  }
  
  # Cargar paquetes
  sapply(c(cran_packages, bioc_packages, git_packages), require, character.only = TRUE)
  #Installing packages
# install.packages("pulsar")

  ## this version works too
  #install.packages("huge")
  library(devtools)
  install_github("kingaa/pomp")
  install_github("zdk123/SpiecEasi")

#It is important to install all repositories needded for SpiecEasi such as gfrotran 12.1 verion
}

#Load necessarry libraries
{
library(devtools)
library(phyloseq)
library(microbiome)
library(genefilter)
library(huge)
library(pulsar)
library(MASS)
library(seqtime)
library(remotes)
library(SpiecEasi)
library(igraph)
library(qgraph)
library(ggnet)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(network)
library(sna)
library(ggplot2)
library(ggnetwork)
library(tidyr)
library(Hmisc)
library(rgexf)
library(GGally)
library(Rmisc)
library(FSA)
library(intergraph)
}

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
t=c("Agave.tequilana", "Agave.marmorata", "Myrtillocactus.geometrizans","spores")[2]
colores=c("#29854E","#74A557","#B8C466","#FFE081",
          "#A02355","#C55DA1","#DB99EE",
          "#234EA0","#008CCF","#00C6D9")
colnet=c("goldenrod3","#8c96c6", "grey80")

{
  ## Generate OTU table
  setwd("~/abel/analysis/network")
  data.dir=("~/abel/analysis/network")
  main.dir=getwd()
  data=c("abel","nestor", "all")[1]  #what sub sample would you need
  s16=("~/abel/analysis/16s_analy/result/")
  its=("~/abel/analysis/FITS_analy/result/")  
  met=read.delim(file=paste(its,paste("metadata.table",data,"txt",sep = "."),sep = "/"))
  #met= subset(met, met$compartment != "rhizosphere")
  met= subset(met, met$compartment != "leaf.endosphere")
  met= subset(met, met$compartment != "root.endosphere")
  met= subset(met, met$compartment != "soil")
  met= subset(met, met$compartment != "phyllosphere")
  f="AMF"
  
  #ITS
  fotu=read.delim(file=paste(its,paste("otu.table",data,"net.txt",sep = "."),sep = "/"))
  colnames(fotu)=gsub("X","",colnames(fotu))
  rownames(met)=met$met$ID_seq_ITS
  ftax=read.delim(paste(its,paste("taxa.table",data,"uniteall","txt",sep = "."),sep = "/"))
  ##SUBSET GLOMEROMYCOTA
  #ftax= subset(ftax, ftax$phylum == "p:Glomeromycota", select= c("domain","phylum","class","order","family","genus","otu.id"))
  #fotu=fotu[rownames(ftax),]
  
  #16S OTUS in Spores
  potu=read.delim(file=paste(s16,paste("otu.table",data,"txt",sep = "."),sep = "/"))
  #potu=read.delim(file=paste(s16,paste("otu.table.EB1",data,"txt",sep = "."),sep = "/"))
  #All bacteria
  #potu=read.delim(file=paste(s16,paste("otu.table",data,"txt",sep = "."),sep = "/"))
  colnames(potu)=gsub("X","",colnames(potu))
  rownames(met)=met$ID_seq_16s
  ptax=read.delim(paste(s16,paste("taxa.table",data,"silva","txt",sep = "."),sep = "/"))
  #SUBSET ONLY BACTERIA OTUS IN SPORES
  #ptax=subset(ptax[rownames(potu),])
  
  #rearrage taxonomy
  rownames(ftax)=gsub("OTU","FOTU",rownames(ftax))
  rownames(ptax)=gsub("OTU","POTU",rownames(ptax))
  ftax$otu.id=rownames(ftax)
  ptax$otu.id=rownames(ptax)
  
  tax=rbind(ptax,ftax)
  #tax=rbind(ptax)
  tax$domain[is.na(tax$domain)] = "d:Unidentified"
  tax$phylum[is.na(tax$phylum)] = "p:Unidentified"
  tax$class[is.na(tax$class)] = "c:Unidentified"
  tax$order[is.na(tax$order)] = "o:Unidentified"
  tax$family[is.na(tax$family)] = "f:Unidentified"
  tax$genus[is.na(tax$genus)] = "g:Unidentified" 
  
  #meake sure to have the same samples
  rownames(fotu)=gsub("OTU","FOTU",rownames(fotu))
  rownames(potu)=gsub("OTU","POTU",rownames(potu))
  
  #Unir fotu y potu 
  potu=potu[,colnames(potu) %in% colnames(fotu)]
  fotu=fotu[,colnames(fotu) %in% colnames(potu)]
  otu=rbind(potu,fotu)
  #otu=rbind(potu)
  
  #generate metadata
  met.2=met[colnames(potu),]
  #limit= 1 ; thr=0.80 #limit = minimum read count for abundant OTUs
  limit=3
  ms=3 #minimum sample frequency for abundant OTUs 
  pc=c("soil", "rhizosphere", "root.endosphere", "leaf.endosphere", "phyllosphere")[2]
  tr=c("T0", "T6", "T12M","T12P")[2]
  ps=c("EV", "MOCK", "PFC", "PFCS","CAMP")[1:4]
  
  for ( t in pc){
  samples=rownames(met.2[met.2$treatment %in% ps & met.2$compartment %in% t & met.2$sampling %in% tr,])
  }
  #for ( t in pc){
  #samples=rownames(met.2[met.2$treatment %in% ps & met.2$compartment %in% t,])
  #}
  otus=otu[,samples]
  otus=otus[rowSums(otus>=limit)>=ms,]
  #otus=otus[rowSums(otus)>=length(samples)*thr,]
  rowSums(otus)
  #Generate matrices
  otu=as.matrix(otus)
  taxon=as.matrix(tax)
  
  #Save Data
  #exit_name = "FUNGI"
  exit_name = "Bacteria"
  save(otu, file = paste(main.dir,paste(t,exit_name,"OTU_table",sep = "_"),sep = "/"))
  save(met.2, file = paste(main.dir,paste(t,exit_name,"metadata_table",sep = "_"),sep = "/"))
  save(taxon, file = paste(main.dir,paste(t,exit_name,"taxonomy_table",sep = "_"),sep = "/"))
  
  Bacteria=phyloseq(otu_table(otu, taxa_are_rows = T),tax_table(taxon),sample_data(met))
  save(Bacteria, file = paste(getwd(),paste(t,exit_name,"phyloseq",sep="_"),sep = "/"))

}

#Load Phyloseq objet
#load("~/abel/analysis/network/rhizosphere_Bacteria_phyloseq")
#load("~/abel/analysis/network/root.endosphere_Bacteria_phyloseq")
#load("~/abel/analysis/network/leaf.endosphere_Bacteria_phyloseq")
#load("~/abel/analysis/network/soil_Bacteria_phyloseq")
load("~/abel/analysis/network/phyllosphere_Bacteria_phyloseq")

Bacteria

tax_tbl <- as.data.frame(Bacteria@tax_table@.Data)

#Using spiec.easi

#using default parameter for spores (rizosfer:1e-4,)
se_mb= spiec.easi(Bacteria, method='mb',lambda.min.ratio=1e-4, nlambda=100, #1e-2
                  pulsar.params=list(rep.num=50, ncores=6))
 
getStability(se_mb)

saveRDS(se_mb,"T6_phyllosphere.ed_mb.RDS")
#saveRDS(se_mb,"T6_rhizosphere.ed_mb.RDS")


#3.1 SPIEC-EASI
setwd("~/abel/analysis/network")
# Read network, select the plant 
#se_mb <- readRDS(file = "ateqAMF_mb.RDS")
#se_mb <- readRDS(file = "asalAMF_mb.RDS")
#se_mb <- readRDS(file = "mgeoAMF_mb.RDS")
#se_mb <- readRDS(file = "allspo_mb.RDS")
#se_mb <- readRDS(file = "allrhizo_mb.RDS")
#se_mb <- readRDS(file = "all_root.ed_mb.RDS")
#se_mb <- readRDS(file = "all_leaf.ed_mb.RDS")
#se_mb <- readRDS(file = "all_phyllosphere.ed_mb.RDS")
#We have now fit a network, but since we have only a rough, discrete sampling of networks along the lambda path, we should check how far we are from the target stability threshold (0.05).  
getStability(se_mb)

{
# Add OTU names to rows and columns
# Create igraph objects
se_net <- adj2igraph(getRefit(se_mb),
                     rmEmptyNodes = TRUE, diag = FALSE, 
                     vertex.attr = list(name = taxa_names(Bacteria)))

  #save igraph
save(se_net, file = paste(paste(t,"network","igraph",sep = "."),sep = "/"))
#Plot the network, you can choose the color by domain, phylum, order etc etc
pdf(file=paste(paste(t,"network.pdf")), colormodel = "cmyk", width = 8.5, height = 11, compress = F)
a=plot_network(se_net, Bacteria, type = "taxa", color = "phylum", shape = "domain", label = NULL, title = t)
print(a)
dev.off()

#Using ggnet2
# Extract matrix of adjacency
net_class <- as_adjacency_matrix(se_net, type = "both")
# Generate network object
net_class <- network(as.matrix(net_class), 
                     vertex.attrnames = taxa_names(Bacteria), 
                     matrix.type = "adjacency", directed = F)
#using ggnet2
ggnet2(net_class)

# Copy se_net to the object
#Metrics
dens=c();diam=c();modu=c();tran=c();aspa=c();node=c();edge=c();node.names=c();pos.edges=c();
neg.edges=c();to.edges=c();bf=c();bfp=c();bb=c();bbp=c();ff=c();ffp=c()
#calculate metrics
net <- se_net
nodes <- V(net) # Nodes
edges <- E(net) # Edges
node.names <- V(net)$name
num.nodes <- vcount(net) 
num.edges <- ecount(net)
node.names <- V(net)$name
dens=c(dens,edge_density(net, loops = FALSE)) 
diam=c(diam,diameter(net, directed=F, weights=NA))
modu=c(modu,modularity(net, membership(cluster_fast_greedy(net,weights = abs(E(net)$weight))), weights = abs(E(net)$weight)))
tran=c(tran,transitivity(net, type = "undirected", vids = NULL, weights = abs(E(net)$weight), isolates = c("NaN")))
aspa=c(aspa,mean_distance(net, directed = TRUE, unconnected = TRUE))
node=c(node,length(V(net)))
edge=c(edge,length(E(net)))

# Extract the regression coefficients of se_mb
betaMat <- as.matrix(symBeta(getOptBeta(se_mb)))
# Calculate the number of positive and negative edges in the network
positive <- length(betaMat[betaMat>0])/2 
negative <- length(betaMat[betaMat<0])/2 
total <- length(betaMat[betaMat!=0])/2
#Determine bacterial-fungal edges
ver=as_data_frame(net, what = c("edges"))
one=ver[grep("POTU_",ver$from),]
two=dim(one[grep("FOTU_",one$to),])[1]
three=ver[grep("FOTU_",ver$from),]
four=dim(three[grep("POTU_",three$to),])[1]
bf=c(bf,c(two+four))
bfp=c(bfp,c((two+four)/dim(ver)[1]*100))

#Determine bacterial-bacteria edges
one=ver[grep("POTU_",ver$from),]
two=dim(one[grep("POTU_",one$to),])[1]
bb=c(bb,c(two))
bbp=c(bbp,c((two)/dim(ver)[1]*100))

#Determine fungal-fungal edges
one=ver[grep("FOTU_",ver$from),]
two=dim(one[grep("FOTU_",one$to),])[1]
ff=c(ff,c(two))
ffp=c(ffp,c((two)/dim(ver)[1]*100))

metrics=data.frame(network=t,density=dens,diameter=diam,modularity=modu,clustering=tran,av.short.phat=aspa,
                   vertex=node,edges=edge, p.e.edges=bfp,p.p.edges=bbp,e.e.edges=ffp,pos.edges=positive,
                   neg.edges=negative,total.edges=total)

#Save metric table 
write.table(metrics, file = paste(paste(t,"plantsmetrics","txt",sep = "."),sep = "/"),
            row.names = F, quote = F, sep = "\t", col.names = T)

# The first step is to extract the signs of the regression coefficients from the regression coefficient matrix
tax_ids <- taxa_names(Bacteria)
edges <- E(net) # edges
edge_colors <- c()
for(e_index in 1:length(edges)){
  adj_nodes <- ends(net,edges[e_index])
  xindex <- which(tax_ids==adj_nodes[1])
  yindex <- which(tax_ids==adj_nodes[2])
  beta <- betaMat[xindex,yindex]
  if(beta>0){
    edge_colors=append(edge_colors,"darkgreen") # positive
  }else if(beta<0){
    edge_colors=append(edge_colors,"red3") # negative
  }
}
E(net)$color <- edge_colors

# Extract matrix of adjacency
net_class <- as_adjacency_matrix(net, type = "both")
# Generate network object
net_class <- network(as.matrix(net_class), 
                     vertex.attrnames = taxa_names(Bacteria), 
                     matrix.type = "adjacency", directed = F)

ggnet2(net_class, edge.color = E(net)$color)

}

#Extract taxonomy table
tax_tbl <- as.data.frame(Bacteria@tax_table@.Data)
tax_tbl$domain=gsub("d:unidentified", "d:Unidentified", tax_tbl$domain)
tax_tbl$domain=gsub("d:Unidentified", "d:Unidentified", tax_tbl$domain)
#Replace the name of each node with its Phylum, you can set each level taxa
V(net)$name <- as.character(getTaxonomy(V(net)$name, tax_tbl, level = "phylum", useRownames = TRUE))
#V(net)$name <- as.character(getTaxonomy(V(net)$name, tax_tbl, level = "domain", useRownames = TRUE))
#V(net)$name <- as.character(getTaxonomy(V(net)$name, tax_tbl, level = "order", useRownames = TRUE))

#Save Phylum list per node in 'nodenames
nodenames <- V(net)$name
# Extract matrix of adjacency
net_class <- as_adjacency_matrix(net, type = "both")
# Generate network object
net_class <- network(as.matrix(net_class), 
                     vertex.attrnames = taxa_names(Bacteria), 
                     matrix.type = "adjacency", directed = F)

#Check unique phylum
unique(tax_tbl$phylum)
#unique(tax_tbl$order)


#color set for Phylum
colorsphylum=c("p:Proteobacteria"="greenyellow", "p:Bacteroidota"="grey80","p:Firmicutes"="gold", 
          "p:Unidentified"= "grey20","p:Actinobacteriota"= "pink","p:Entotheonellaeota"="grey80",
          "p:Chloroflexi"= "grey80","p:Gemmatimonadota"="grey80","p:Crenarchaeota"="grey80", 
          "p:Acidobacteriota"= "grey80","p:Cyanobacteria" = "grey80","p:Myxococcota"="grey80",
          "p:Verrucomicrobiota"="grey80", "p:Planctomycetota"="grey80","p:Deinococcota"= "grey",
          "p:Nitrospirota"="grey80","p:Armatimonadota"="orange", "p:Bdellovibrionota"="grey80", "p:Methylomirabilota" = "grey80",
          "p:Sumerlaeota"="grey80", "p:Abditibacteriota"="grey80","p:Deinococcota"= "grey80",
          "p:Patescibacteria"="grey80", "p:NB1-j"="grey80", "p:Thermoplasmatota"="grey80", "p:Fibrobacterota"="grey80",
          "s:uncultured_bacterium"="grey80","f:uncultured"="grey80","p:Elusimicrobiota"="grey80", "p:Latescibacterota"="grey80", 
          "p:Desulfobacterota"="grey80", "p:Dependentiae"="grey80", "p:RCP2--54"="grey80", "p:WS2"="grey80", "p:RCP2-54"="grey80",
          "g:uncultured"="grey80","p:Ascomycota"="red","p:Basidiomycota"="blue","p:Mucoromycota"="grey80",
          "p:Chytridiomycota"="grey80","p:Glomeromycota"="grey80","p:Aphelidiomycota"="grey80","p:Mortierellomycota"="grey80","p:Spirochaetota"="grey80" )

#color set for order subset AMF
#colorsorder=c("o:Pseudonocardiales"="grey80", "o:Micrococcales"="grey80","o:Cytophagales"="grey80", "o:Pasteurellales"="grey80",
#              "o:Mycoplasmatales"="gold","o:Rhizobiales"="grey80","o:Pseudomonadales"="grey80", "o:Burkholderiales"="greenyellow",
 #             "o:Unidentified"="grey80","o:Bacillales"="grey80", "o:Micromonosporales"="grey80","o:Corynebacteriales"="grey80",         
  #            "o:Rhodocyclales"="grey80", "o:Rhodobacterales"="grey80", "o:Neisseriales"="grey80","o:Propionibacteriales"="grey80",
   #           "o:Caulobacterales"="grey80","o:Streptomycetales"="grey80", "o:S0134_terrestrial_group" ="grey80","o:Xanthomonadales"="grey80", 
  #          "o:Selenomonadales"="grey80", "o:Clostridiales"="grey80", "o:Myxococcales"="grey80","o:AT425-EubC11_terrestrial_group"="grey80",
   #          "o:Nitriliruptorales"="grey80","o:Fusobacteriales" ="grey80","o:Flavobacteriales"="grey80", "o:Frankiales"="grey80",  
  #          "o:Lactobacillales"="grey80","o:Bacteroidales"="grey80","o:Enterobacteriales" ="grey80","o:Streptosporangiales"="grey80",
  #            "o:Acidimicrobiales"="grey80","o:Oligoflexales"="grey80","o:Rhodospirillales"="grey80","o:Anaerolineales"="grey80",            
   #           "o:Actinomycetales"="grey80", "o:Sphingomonadales"="grey80","o:Campylobacterales"="grey80", "o:Aeromonadales"="grey80","o:Chromatiales"="grey80",
    #          "o:Archaeosporales"="violetred","o:Glomerales"="sienna1","o:Paraglomerales"="royalblue","o:Diversisporales"="orchid")

colorsorder=c("o:Cytophagales"="grey80","o:Lachnospirales"="grey80","o:Xanthomonadales"="grey80","o:Sphingomonadales"="#48c9b0","o:Chitinophagales"="grey80","o:Unidentified"="grey80","o:Microtrichales"="grey80","o:Rhodobacterales"="grey80","o:Burkholderiales"="grey80","o:Entotheonellales"="grey80","o:S085"="grey80","o:Longimicrobiales"="grey80","o:Nitrososphaerales"="grey80","o:Bryobacterales"="grey80","o:Rhizobiales"="purple","o:Gemmatimonadales"="grey80","o:Paenibacillales"="grey80","o:Spirochaetales"="grey80","o:Thermoanaerobaculales"="grey80",
"o:Solirubrobacterales"="grey80","o:Gaiellales"="grey80","o:Azospirillales"="grey80","o:Pseudomonadales"="grey80","o:KD4-96"="grey80","o:Gammaproteobacteria_Incertae_Sedis"="grey80","o:Acetobacterales"="grey80","o:uncultured"="grey80","o:SAR202_clade"="grey80","o:Nannocystales"="grey80","o:Rubrobacterales"="grey80","o:Haliangiales"="grey80","o:Caulobacterales"="grey80","o:Chthoniobacterales"="grey80","o:Elsterales"="grey80","o:Pirellulales"="grey80","o:Corynebacteriales"="grey80","o:Isosphaerales"="grey80","o:IMCC26256"="grey80","o:Tistrellales"="grey80",
"o:Vicinamibacterales"="grey80","o:Propionibacteriales"="brown","o:0319-7L14"="grey80","o:Pyrinomonadales"="grey80","o:Polyangiales"="grey80","o:S0134_terrestrial_group"="grey80","o:Pedosphaerales"="grey80","o:Alicyclobacillales"="grey80","o:Tepidisphaerales"="grey80","o:Nitrospirales"="grey80","o:Gemmatales"="grey80","o:bacteriap25"="grey80","o:Thermomicrobiales"="grey80","o:Cytophagales"="grey80","o:Lachnospirales"="grey80","o:Xanthomonadales"="grey80","o:Chitinophagales"="grey80","o:Unidentified"="grey80",
"o:Microtrichales"="grey80","o:Rhodobacterales"="grey80","o:Burkholderiales"="grey80","o:Entotheonellales"="grey80","o:S085"="grey80","o:Longimicrobiales"="grey80","o:Nitrososphaerales"="grey80","o:Bryobacterales"="grey80","o:Rhizobiales"="grey80","o:Gemmatimonadales"="grey80","o:Paenibacillales"="grey80","o:Spirochaetales"="grey80","o:Thermoanaerobaculales"="grey80","o:Solirubrobacterales"="grey80","o:Gaiellales"="grey80","o:Azospirillales"="grey80","o:Pseudomonadales"="grey80","o:KD4-96"="grey80","o:Gammaproteobacteria_Incertae_Sedis"="grey80",
"o:Acetobacterales"="grey80","o:uncultured"="grey80","o:SAR202_clade"="grey80","o:Nannocystales"="grey80","o:Rubrobacterales"="grey80","o:Haliangiales"="grey80","o:Caulobacterales"="grey80","o:Chthoniobacterales"="grey80","o:Elsterales"="grey80","o:Pirellulales"="grey80","o:Corynebacteriales"="grey80","o:Isosphaerales"="grey80","o:IMCC26256"="grey80","o:Tistrellales"="grey80","o:Vicinamibacterales"="grey80","o:Propionibacteriales"="grey80","o:0319-7L14"="grey80","o:Pyrinomonadales"="grey80","o:Polyangiales"="grey80","o:S0134_terrestrial_group"="grey80",
"o:Pedosphaerales"="grey80","o:Alicyclobacillales"="grey80","o:Tepidisphaerales"="grey80","o:Nitrospirales"="grey80","o:Gemmatales"="grey80","o:bacteriap25"="grey80","o:Thermomicrobiales"="grey80","o:MB-A2-108"="grey80","o:JG30-KF-CM66"="grey80","o:Streptosporangiales"="grey80","o:Alteromonadales"="grey80","o:Micrococcales"="pink","o:Kallotenuales"="grey80","o:Micromonosporales"="grey80","o:Bacteroidales"="grey80","o:Peptostreptococcales-Tissierellales"="grey80","o:Pseudonocardiales"="grey80","o:Planctomycetales"="grey80","o:Sphingobacteriales"="red",
"o:Bacillales"="yellow","o:Fibrobacterales"="grey80","o:Gitt-GS-136"="grey80","o:TK10"="grey80","o:Opitutales"="grey80","o:Streptomycetales"="grey80","o:Brevibacillales"="grey80","o:Flavobacteriales"="grey80","o:Verrucomicrobiales"="grey80","o:Chloroflexales"="grey80","o:Veillonellales-Selenomonadales"="grey80","o:Euzebyales"="grey80","o:Steroidobacterales"="grey80","o:C0119"="grey80","o:PLTA13"="grey80","o:BD2-11_terrestrial_group"="grey80","o:Thermoactinomycetales"="grey80","o:Subgroup_7"="grey80","o:Subgroup_17"="grey80","o:Rokubacteriales"="grey80",
"o:Clostridiales"="grey80","o:Frankiales"="grey80","o:R7C24"="grey80","o:Solibacterales"="grey80","o:mle1-27"="grey80","o:Enterobacterales"="grey80","o:Saccharimonadales"="grey80","s:uncultured_bacterium"="grey80","o:Obscuribacterales"="grey80","o:Thalassobaculales"="grey80","o:Salinisphaerales"="grey80","o:Glycomycetales"="grey80","o:Erysipelotrichales"="grey80","o:Blastocatellales"="grey80","o:Kineosporiales"="grey80","o:Myxococcales"="grey80","o:Eurotiales"="grey80","o:Agaricales"="grey80","o:Verrucariales"="grey80","o:Sordariales"="grey80",
"o:Cantharellales"="grey80","o:Pleosporales"="grey80","o:Cephalothecales"="grey80","o:Hypocreales"="grey80","o:Mucorales"="grey80","o:Orbiliales"="grey80","o:Cystobasidiomycetes_ord_Incertae_sedis"="grey80","o:Xylariales"="grey80","o:Chaetothyriales"="grey80","o:Glomerellales"="grey80","o:Onygenales"="grey80","o:Wallemiales"="grey80","o:Amphisphaeriales"="grey80","o:Pezizales"="grey80","o:Mycosphaerellales"="grey80","o:Branch06"="grey80","o:Microascales"="grey80","o:Botryosphaeriales"="grey80","o:Geminibasidiales"="grey80","o:Spizellomycetales"="grey80",
"o:Glomerales"="grey80","o:Leotiales"="grey80","o:Dothideomycetes_ord_Incertae_sedis"="grey80","o:Helotiales"="grey80","o:Trichosporonales"="grey80","o:Filobasidiales"="grey80","o:Phaeomoniellales"="grey80","o:Mortierellales"="grey80","o:Tremellales"="grey80","o:MB-A2-108"="grey80","o:JG30-KF-CM66"="grey80","o:Streptosporangiales"="grey80","o:Alteromonadales"="grey80","o:Micrococcales"="grey80","o:Kallotenuales"="grey80","o:Micromonosporales"="grey80","o:Bacteroidales"="grey80","o:Peptostreptococcales-Tissierellales"="grey80","o:Pseudonocardiales"="grey80",
"o:Planctomycetales"="grey80","o:Sphingobacteriales"="grey80","o:Bacillales"="grey80","o:Fibrobacterales"="grey80","o:Gitt-GS-136"="grey80","o:TK10"="grey80","o:Opitutales"="grey80","o:Streptomycetales"="grey80","o:Brevibacillales"="grey80","o:Flavobacteriales"="grey80","o:Verrucomicrobiales"="grey80","o:Chloroflexales"="grey80","o:Veillonellales-Selenomonadales"="grey80","o:Euzebyales"="grey80","o:Steroidobacterales"="grey80","o:C0119"="grey80","o:PLTA13"="grey80","o:BD2-11_terrestrial_group"="grey80","o:Thermoactinomycetales"="grey80","o:Subgroup_7"="grey80",
"o:Subgroup_17"="grey80","o:Rokubacteriales"="grey80","o:Clostridiales"="grey80","o:Frankiales"="grey80","o:R7C24"="grey80","o:Solibacterales"="grey80","o:mle1-27"="grey80","o:Enterobacterales"="grey80","o:Saccharimonadales"="grey80","s:uncultured_bacterium"="grey80","o:Obscuribacterales"="grey80","o:Thalassobaculales"="grey80","o:Salinisphaerales"="grey80","o:Glycomycetales"="grey80","o:Erysipelotrichales"="grey80","o:Blastocatellales"="grey80","o:Kineosporiales"="grey80","o:Myxococcales"="grey80","o:Eurotiales"="grey80","o:Agaricales"="grey80",
"o:Verrucariales"="grey80","o:Sordariales"="grey80","o:Cantharellales"="grey80","o:Pleosporales"="grey80","o:Cephalothecales"="grey80","o:Hypocreales"="grey80","o:Mucorales"="grey80","o:Orbiliales"="grey80","o:Cystobasidiomycetes_ord_Incertae_sedis"="grey80","o:Xylariales"="grey80","o:Chaetothyriales"="grey80","o:Glomerellales"="grey80","o:Onygenales"="grey80","o:Wallemiales"="grey80","o:Amphisphaeriales"="grey80","o:Pezizales"="grey80","o:Mycosphaerellales"="grey80","o:Branch06"="grey80","o:Microascales"="grey80","o:Botryosphaeriales"="grey80",
"o:Geminibasidiales"="grey80","o:Spizellomycetales"="grey80","o:Glomerales"="grey80","o:Leotiales"="grey80","o:Dothideomycetes_ord_Incertae_sedis"="grey80","o:Helotiales"="grey80","o:Trichosporonales"="grey80","o:Filobasidiales"="grey80","o:Phaeomoniellales"="grey80","o:Mortierellales"="grey80","o:Tremellales"="grey80","p:Bacteroidota"="grey80","p:Firmicutes"="grey80","p:Proteobacteria"="grey80","p:Unidentified"="grey80","p:Actinobacteriota"="grey80","p:Entotheonellaeota"="grey80","p:Chloroflexi"="grey80","p:Gemmatimonadota"="grey80",
"p:Crenarchaeota"="grey80","p:Acidobacteriota"="grey80","p:Spirochaetota"="grey80","p:Myxococcota"="grey80","p:Verrucomicrobiota"="grey80","p:Planctomycetota"="grey80","p:Nitrospirota"="grey80","p:Fibrobacterota"="grey80","p:Methylomirabilota"="grey80","p:Armatimonadota"="grey80","p:Patescibacteria"="grey80","p:Cyanobacteria"="grey80","g:uncultured"="grey80","p:Ascomycota"="grey80","p:Basidiomycota"="grey80","p:Mucoromycota"="grey80","p:Chytridiomycota"="grey80","p:Glomeromycota"="grey80","p:Mortierellomycota"="grey80","p:Bdellovibrionota"="grey80",
"p:Rs-K70_termite_group"="grey80", "p:Deinococcota"="grey80", "p:Sumerlaeota"="grey80", "o:Oligoflexales"="grey80", "o:Ardenticatenales"="grey80", "o:Deinococcales"="grey80", "o:Armatimonadales"="grey80", "o:Bacteriovoracales"="grey80", "o:Cellvibrionales"="grey80", "o:Blfdi19"="grey80", "o:Rhodothermales"="grey80", "o:Chloroplast"="grey80", "o:Dothideales"="grey80", "o:Agaricostilbales"="grey80", "o:Russulales"="grey80", "o:Hysteriales"="grey80","o:Staphylococcales"="grey80", "o:SBR1031"="grey80", "o:Nitrosococcales"="grey80", "o:11-24"="grey80", 
"o:Oxyphotobacteria_Incertae_Sedis"="grey80", "o:Cyanobacteriales"="grey80", "o:Sumerlaeales"="grey80", "o:Caldilineales"="grey80", "o:Reyranellales"="grey80", "o:Aeromonadales"="grey80", "o:Sebacinales"="grey80", "o:Lichenostigmatales"="grey80", "o:Ascomycota_ord_Incertae_sedis"="grey80", "o:Cladosporiales"="grey80", "o:Rhizophlyctidales"="grey80")



#color set for domain/kingdom
#unique(tax_tbl$domain)
#colorsdomain=c("k:Bacteria"= "royalblue", "d:Fungi"= "gold3","d:Unidentified"="grey20")

#Network for order
pdf(file=paste(paste(t,"order-1.pdf")), colormodel = "srgb", width = 11, height = 8.5, compress = F)
a=ggnet2(net_class, color = nodenames, palette = colorsorder, edge.color = E(net)$color, 
       edge.size = 0.5, title = t)
print(a)
dev.off()

#Network for phlyum
pdf(file=paste(paste(t,"phylum.pdf")), colormodel = "srgb", width = 11, height = 8.5, compress = F)
a1=ggnet2(net_class, size= 6, color = nodenames, palette = colorsphylum, edge.color = E(net)$color, edge.size = 0.5, tittle = t)
print(a1)
dev.off()
#network for domain
#pdf(file=paste(paste(t,"domain.pdf")), colormodel = "srgb", width = 8, height = 8, compress = F)
#a2=ggnet2(net_class, size= 3, color = nodenames, palette = colorsdomain, edge.color = E(net)$color, edge.size = 0.5, title = t)
#print(a2)
#dev.off()

#5 network statistics
#Create a copy of se_net
net <- se_net
# Assign coordinates to the network design
net$layout <- array(1:40, dim = c(20, 2))
# Assign the layout as a fixed attribute of the network
net$layout <- layout.fruchterman.reingold(net)

# #Select the most connected OTUS
deg <- igraph::degree(net, mode = "all")
#degree distribution
deg.dist <- degree_distribution(net, mode = "all", cumulative = F)
#Plot of degree distribution
plot(deg.dist, xlab = "Nodes degree", ylab = "Probability")
lines(deg.dist)
#Select the most central OTUS
clos <- igraph::closeness(net, mode = "all")
betw <- igraph::betweenness(net, v = V(net))

pdf(file=paste(paste(t,"centralityplot.pdf")), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a3=centralityPlot(net, include = c("Betweenness", "Closeness", "Degree")) + 
  theme(axis.text.y = element_blank())
print(a3)
dev.off()

# Para calcular el valor de ANND para todos los nodos en la red, usamos el argumento vids = V(net)
net.knn <- knn(net, vids = V(net))

# Esta función intenta detectaar sub-redes densamente conectadas, usando 'random walks'
# 'random walks' se refiere a "recorrer" la red de forma aleatoria
#wt <- fastgreedy.community(net)
wt <- cluster_fast_greedy(net)

# Consultar mimebros de cada nodo
pdf(file=paste(paste(t,"dendogram.pdf")), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a4=igraph::plot_dendrogram(wt)
print(a4)
dev.off()
# Plot
pdf(file=paste(paste(t,"Clusters.pdf")), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a5=plot(wt, net,  main=paste(paste(t,"clusters")))
print(a5)
dev.off()

#KEYSTONES
#degree
deg <- igraph::degree(net)
#Sort nodes by degree from highest to lowest
deg_sort <- sort(deg, decreasing = TRUE)
# Transform 'deg_sort' to 'data frame' neccesary to graph the histogram.
deg_df <- as.data.frame(deg_sort)


# Histogram
ggplot(deg_df, aes(x = deg_sort)) + 
  geom_histogram(binwidth = 1) + 
  scale_x_continuous(breaks=c(10,20,30,40,50)) + 
  geom_vline(aes(xintercept=mean(deg_sort)),
             color="black", linetype="dashed", size=1) + 
  theme_minimal() + 
  labs(x = "Degree", y = "Node count")

#Node centrality
#betweenness
bn <- igraph::betweenness(net)
# Sort nodes according to betweenness
bn_sort <- sort(bn, decreasing = TRUE)
bn_df <- as.data.frame(bn_sort)
# Graficar histograma
ggplot(bn_df, aes(x = bn_sort)) + 
  geom_histogram(binwidth = 30) +  
  geom_vline(aes(xintercept=mean(bn_sort)), 
             color="black", linetype="dashed", size=1) + 
  theme_minimal() + 
  labs(x = "Betweenness", y = "Node count")

# Usamos la función 'grid.arrange' del paquete 'gridExtra' para visualizar ambos histogramas juntos

pdf(file=paste(paste(t,"Network Nodes Centrality:Betweenness.pdf")), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
pdeg <- ggplot(deg_df, aes(x = deg_sort)) + 
  geom_histogram(binwidth = 1) + 
  scale_x_continuous(breaks=c(10,20,30,40,50)) + 
  geom_vline(aes(xintercept=mean(deg_sort)),
             color="black", linetype="dashed", size=1) + 
  theme_minimal() + 
  labs(x = "Degree", y = "Node count", title = "Network Nodes Degree")
pbn <- ggplot(bn_df, aes(x = bn_sort)) + 
  geom_histogram(binwidth = 30) +  
  geom_vline(aes(xintercept=mean(bn_sort)), 
             color="black", linetype="dashed", size=1) + 
  theme_minimal() + 
  labs(x = "Betweenness", y = "Node count", title = "Network Nodes Centrality: Betweenness")
a6=grid.arrange(pdeg, pbn, ncol = 2)
print(a6)
dev.off()

#Edit degree table. You must check if you are using differentes variables. 
deg_df$TaxID <- row.names(deg_df)
#Filter taxax with degree >6 (ateq), 7 (asal), 7(spores). 6 (mgeo)
deg_high_df <- dplyr::filter(deg_df, deg_sort > 30)
head(deg_high_df)
#Edit betweenness table
bn_df$TaxID <- row.names(bn_df)
# Filtramos las taxas con betweenness > 700 (ateq), >900 (asal), >600 (spores), 500(mgeo)
bn_high_df <- dplyr::filter(bn_df, bn_sort > 2000)
head(bn_high_df)
# #Select the hub OTUS and las tablas de degree y betweenness filtradas
keystone <- merge(deg_high_df, bn_high_df, all.x = FALSE)
print(keystone)

# Usamos ggplot2
ggplot(keystone, aes(x = bn_sort, y = deg_sort, label = TaxID)) + 
  scale_y_continuous(limits = c(1,30), breaks = c(10,20,30)) + 
  scale_x_continuous(limits = c(8000, 3000), breaks = c(2000, 4000, 6000,8000)) + 
  geom_point(alpha = 4, color = "#829FD9", size = 6) + 
  geom_text(size = 4) + 
  theme_minimal() + 
  labs(x = "Betweenness", y = "Degree", 
       title = "Nodes With Highest Degree And Betweenness: Keystone Nodes")

# Graficar red chagen, color of domain or order
ggnet2(se_net, mode = net$layout, edge.size = 0.2,
       color = nodenames, palette = colorsphylum, 
       node.alpha = 0.7, 
       size = 4,
       edge.color = edge_colors, 
       label = keystone$TaxID, label.size = 3.5)
       #legend.position = "none")

#Alternative plot network for Figure 4
#Save previous data
grafica=data.frame(degree=deg <- igraph::degree(net), betweenesscentrality= bn <- igraph::betweenness(net))
grafica=cbind(grafica,tax_tbl[rownames(grafica),])
grafica$hub=0 ; grafica$central=0 ; grafica$grado=0
grafica[keystone$TaxID,"hub"]=1 ;grafica[bn_high_df$TaxID,"central"]=1 ;grafica[deg_high_df$TaxID,"grado"]=1
#Table of nodes
write.table(grafica, file = paste(paste(t,"spores_nodes","txt",sep = "."),sep = "/"),
            row.names = F, quote = F, sep = "\t" ,col.names = T)

###alternative form to graphic hubs.
  file = paste(paste(t,"network","igraph",sep = "."),sep = "/")
  print(t)
  load(file = paste(paste(t,"network","igraph",sep = "."),sep = "/"))
  nodos = read.delim(paste(paste(t,"spores_nodes","txt",sep = "."),sep = "/"))
  nodos$compartment=t
  vertex_attr(se_net, name="kingdom", index = V(se_net)) <- tax_tbl[V(se_net)$name,4]
  vertex_attr(se_net, name="hub", index = V(se_net)) <- ifelse(names(V(se_net)) %in% nodos[nodos$hub==1,"otu.id"],"hub","no_hub")
  colnet=c("goldenrod3","#8c96c6", "grey80")
  edge_attr(se_net, name="peso",index=E(se_net)) <- ifelse(E(se_net)$weight>0, "darkgreen", "red3")

  #table of hub
  hub=nodos[nodos$hub==1,]
  write.table(hub, file = paste(paste(t,"hub","txt",sep = "."),sep = "/"),
              row.names = F, quote = F, sep = "\t", col.names = T)
  
#graph domain
pdf(file=paste(paste(t,"Network-Hubs--Bacteria-L.pdf")), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a7=ggnet2(se_net, size= 4, node.color = "kingdom",  edge.color = edge_colors, node.size = "hub",edge.size = 0.2,
         mode="fruchtermanreingold", label = keystone$TaxID, label.size = 3, legend.position = "none")+
    geom_point(aes(fill=color,size=size),shape=21)+
    scale_fill_manual(values =  colorsorder)+
    scale_size_manual(values = c(5,3))+
    ggtitle(t)
print(a7)
dev.off()

pdf(file=paste(paste(t,"Network-Hubs-Bacteria-NL.pdf")), colormodel = "cmyk", width = 110, height = 80.5, compress = F)
a8=ggnet2(se_net, node.color = "kingdom",  edge.color = edge_colors, node.size = "hub", edge.size = 0.2,
          mode="fruchtermanreingold", label = F, label.size = 3)+
  geom_point(aes(fill=color,size=size),shape=21)+
  scale_fill_manual(values =  colorsorder)+
  scale_size_manual(values = c(5,3))+
  ggtitle(t)
print(a8)
dev.off()
    
# graph red##
vertex_attr(se_net, name="domain", index = V(se_net)) <- tax_tbl[V(se_net)$name,1]
vertex_attr(se_net, name="phylum", index = V(se_net)) <- tax_tbl[V(se_net)$name,4]
se_net1 <- vertex_attr(se_net, name="hub", index = V(se_net)) <- ifelse(names(V(se_net)) %in% nodos[nodos$hub==1,"otu.id"],"hub","no_hub")

#Check if you choose domain-phylum, orden. etc etc 
pdf(file=paste(paste(t,"Network-Hubs-O.pdf")), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a81=ggnet2(se_net, mode = net$layout, edge.size = 0.2, node.size = "hub", 
          node.color = "phylum", palette = colorsorder,edge.color = edge_colors, 
          label = keystone$TaxID, label.size = 2.5, legend.position = "none")+ 
          scale_size_manual(values = c(5, 3))+
          scale_fill_manual(values =  colorsorder)
print(a81)
dev.off()

pdf(file=paste(paste(t,"Network-Hubs-O-NL.pdf")), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a82=ggnet2(se_net, size = 6, mode = net$layout, edge.size = 0.2, node.size = "hub",
       node.color = "phylum", palette = colorsorder, legend.position = "none",
       edge.color = edge_colors, label.size = 2.5,  label = T)+
      scale_size_manual(values = c(5, 3))
print(a82)
dev.off()

##other mode
pdf(file=paste(paste(t,"Network-Hubs-O-NL1.pdf")), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a82=ggnet2(se_net, size = 5, mode = "fruchtermanreingold", edge.size = 0.2, node.size = "hub",
           node.color = "phylum", palette = colorsorder, legend.position = "none",
           edge.color = edge_colors, label.size = 2.5,  label = F, alpha = 0.95, 
           node.shape = "domain",layout.exp = 0.1,edge.alpha = 0.4)+
  scale_size_manual(values = c(5, 3))
print(a82)
dev.off()

pdf(file=paste(paste(t,"Network-Hubs-O-NL-1.1.pdf")), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a82=ggnet2(se_net, size = 7, mode = "fruchtermanreingold", edge.size = 0.2, node.size = "hub",
           node.color = "phylum", palette = colorsorder, legend.position = "none",
           edge.color = edge_colors, label.size = 1.5,  label = T, alpha = 0.95, 
           node.shape = "domain")+
  scale_size_manual(values = c(5, 3))
print(a82)
dev.off()

##EXPORT DATA

#write.csv(a82, "spores.cvs")
#write_graph(a82, "spores.cvs", format = c("gml"))

###Subnetoworks ONLY PHYLUM DATA
#Extract nodes of each of 10 gruops.
{
m1 <- V(net)[wt$membership == 1]
m2 <- V(net)[wt$membership == 2]
m3 <- V(net)[wt$membership == 3]
m4 <- V(net)[wt$membership == 4]
m5 <- V(net)[wt$membership == 5]
m6 <- V(net)[wt$membership == 6]
m7 <- V(net)[wt$membership == 7]


m1_subnet <- induced_subgraph(net, m1)
m2_subnet <- induced_subgraph(net, m2)
m3_subnet <- induced_subgraph(net, m3)
m4_subnet <- induced_subgraph(net, m4)
m5_subnet <- induced_subgraph(net, m5)
m6_subnet <- induced_subgraph(net, m6)
m7_subnet <- induced_subgraph(net, m7)


###SUBNETM1
# Asignar coordenadas al diseño de la red
m1_subnet$layout <- array(2:50, dim = c(20, 2))
# Asignar el diseño como un atributo fijo de la red
m1_subnet$layout <- layout.fruchterman.reingold(m1_subnet)

# Extraer matriz de adyacencia
m1_subnetclass <- as_adjacency_matrix(m1_subnet, type = "both")
# Generar objeto de clase 'network'
m1_subnetclass <- network(as.matrix(m1_subnetclass), 
                          matrix.type = "adjacency", directed = F)
#Plot
ggnet2(m1_subnetclass, mode = m1_subnet$layout)
#Extract taxa names
m1_subnetnodenames <- as.character(getTaxonomy(V(m1_subnet)$name, tax_tbl, level = "phylum", useRownames = TRUE))
# Identificar y asignar color a los edges según si representan una relación positiva o negativa
edges1 <- E(m1_subnet)
edge_colors <- c()
for(e_index in 1:length(edges1)){
  adj_nodes <- ends(m1_subnet,edges1[e_index])
  xindex <- which(tax_ids==adj_nodes[1])
  yindex <- which(tax_ids==adj_nodes[2])
  beta <- betaMat[xindex,yindex]
  if(beta>0){
    edge_colors=append(edge_colors,"darkgreen") # positive
  }else if(beta<0){
    edge_colors=append(edge_colors,"red3") # negative
  }
}
E(m1_subnet)$color <- edge_colors
vertex_attr(m1_subnet, name="phylum", index = V(m1_subnet)) <- tax_tbl[V(m1_subnet)$name,2]
vertex_attr(m1_subnet, name="hub", index = V(m1_subnet)) <- ifelse(names(V(m1_subnet)) %in% nodos[nodos$hub==1,"otu.id"],"hub","no_hub")

#Subnet
pdf(file=paste(paste(t,"sub-Network-M1.pdf")), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
b1=ggnet2(m1_subnet, node.color = "phylum", mode = m1_subnet$layout, node.size = "hub",
       color = selegoG_nodenames, palette = colorsphylum, edge.size = 0.4,
       edge.color = E(m1_subnet)$color, label = T, label.size = 4)+
       scale_size_manual(values = c(6, 3))
print(b1)
dev.off()

###SUBNETM2
# Asignar coordenadas al diseño de la red
m2_subnet$layout <- array(1:40, dim = c(20, 2))
# Asignar el diseño como un atributo fijo de la red
m2_subnet$layout <- layout.fruchterman.reingold(m2_subnet)

# Extraer matriz de adyacencia
m2_subnetclass <- as_adjacency_matrix(m2_subnet, type = "both")
# Generar objeto de clase 'network'
m2_subnetclass <- network(as.matrix(m2_subnetclass), 
                          matrix.type = "adjacency", directed = F)
#Plot
ggnet2(m2_subnetclass, mode = m2_subnet$layout)
#Extract taxa names
m2_subnetnodenames <- as.character(getTaxonomy(V(m2_subnet)$name, tax_tbl, level = "phylum", useRownames = TRUE))
# Identificar y asignar color a los edges según si representan una relación positiva o negativa
edges1 <- E(m2_subnet)
edge_colors <- c()
for(e_index in 1:length(edges1)){
  adj_nodes <- ends(m2_subnet,edges1[e_index])
  xindex <- which(tax_ids==adj_nodes[1])
  yindex <- which(tax_ids==adj_nodes[2])
  beta <- betaMat[xindex,yindex]
  if(beta>0){
    edge_colors=append(edge_colors,"darkgreen") # positive
  }else if(beta<0){
    edge_colors=append(edge_colors,"red3") # negative
  }
}
E(m2_subnet)$color <- edge_colors
vertex_attr(m2_subnet, name="phylum", index = V(m2_subnet)) <- tax_tbl[V(m2_subnet)$name,2]
vertex_attr(m2_subnet, name="hub", index = V(m2_subnet)) <- ifelse(names(V(m2_subnet)) %in% nodos[nodos$hub==1,"otu.id"],"hub","no_hub")

#Subnet
pdf(file=paste(paste(t,"sub-Network-M2.pdf")), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
b2=ggnet2(m2_subnet, node.color = "phylum", mode = m2_subnet$layout, node.size = "hub",
       color = selegoG_nodenames, palette = colorsphylum, edge.size = 0.4,
       edge.color = E(m2_subnet)$color, label = T, label.size = 4)+
      scale_size_manual(values = c(6, 3))
print(b2)
dev.off()

###SUBNETM3
# Asignar coordenadas al diseño de la red
m3_subnet$layout <- array(1:40, dim = c(20, 2))
# Asignar el diseño como un atributo fijo de la red
m3_subnet$layout <- layout.fruchterman.reingold(m3_subnet)

# Extraer matriz de adyacencia
m3_subnetclass <- as_adjacency_matrix(m3_subnet, type = "both")
# Generar objeto de clase 'network'
m3_subnetclass <- network(as.matrix(m3_subnetclass), 
                          matrix.type = "adjacency", directed = F)
#Plot
ggnet2(m3_subnetclass, mode = m3_subnet$layout)
#Extract taxa names
m3_subnetnodenames <- as.character(getTaxonomy(V(m3_subnet)$name, tax_tbl, level = "phylum", useRownames = TRUE))
# Identificar y asignar color a los edges según si representan una relación positiva o negativa
edges1 <- E(m3_subnet)
edge_colors <- c()
for(e_index in 1:length(edges1)){
  adj_nodes <- ends(m3_subnet,edges1[e_index])
  xindex <- which(tax_ids==adj_nodes[1])
  yindex <- which(tax_ids==adj_nodes[2])
  beta <- betaMat[xindex,yindex]
  if(beta>0){
    edge_colors=append(edge_colors,"darkgreen") # positive
  }else if(beta<0){
    edge_colors=append(edge_colors,"red3") # negative
  }
}
E(m3_subnet)$color <- edge_colors
vertex_attr(m3_subnet, name="phylum", index = V(m3_subnet)) <- tax_tbl[V(m3_subnet)$name,2]
vertex_attr(m3_subnet, name="hub", index = V(m3_subnet)) <- ifelse(names(V(m3_subnet)) %in% nodos[nodos$hub==1,"otu.id"],"hub","no_hub")

#SubnetM3
pdf(file=paste(paste(t,"sub-Network-M3.pdf")), colormodel = "cmyk", width = 14, height = 12, compress = F)
b3=ggnet2(m3_subnet, node.color = "phylum", mode = m3_subnet$layout, node.size = "hub",
          color = selegoG_nodenames, palette = colorsphylum, edge.size = 0.4,
          edge.color = E(m3_subnet)$color, label = T, label.size = 4)+
  scale_size_manual(values = c(6, 3))
print(b3)
dev.off()


###SUBNETM4
# Asignar coordenadas al diseño de la red
m4_subnet$layout <- array(1:40, dim = c(20, 2))
# Asignar el diseño como un atributo fijo de la red
m4_subnet$layout <- layout.fruchterman.reingold(m4_subnet)

# Extraer matriz de adyacencia
m4_subnetclass <- as_adjacency_matrix(m4_subnet, type = "both")
# Generar objeto de clase 'network'
m4_subnetclass <- network(as.matrix(m4_subnetclass), 
                          matrix.type = "adjacency", directed = F)
#Plot
ggnet2(m4_subnetclass, mode = m4_subnet$layout)
#Extract taxa names
m4_subnetnodenames <- as.character(getTaxonomy(V(m4_subnet)$name, tax_tbl, level = "phylum", useRownames = TRUE))
# Identificar y asignar color a los edges según si representan una relación positiva o negativa
edges1 <- E(m4_subnet)
edge_colors <- c()
for(e_index in 1:length(edges1)){
  adj_nodes <- ends(m4_subnet,edges1[e_index])
  xindex <- which(tax_ids==adj_nodes[1])
  yindex <- which(tax_ids==adj_nodes[2])
  beta <- betaMat[xindex,yindex]
  if(beta>0){
    edge_colors=append(edge_colors,"darkgreen") # positive
  }else if(beta<0){
    edge_colors=append(edge_colors,"red3") # negative
  }
}
E(m4_subnet)$color <- edge_colors
vertex_attr(m4_subnet, name="phylum", index = V(m4_subnet)) <- tax_tbl[V(m4_subnet)$name,2]
vertex_attr(m4_subnet, name="hub", index = V(m4_subnet)) <- ifelse(names(V(m4_subnet)) %in% nodos[nodos$hub==1,"otu.id"],"hub","no_hub")

#SubnetM4
pdf(file=paste(paste(t,"sub-Network-M4.pdf")), colormodel = "cmyk", width = 14, height = 12, compress = F)
b4=ggnet2(m4_subnet, node.color = "phylum", mode = m4_subnet$layout, node.size = "hub",
          color = selegoG_nodenames, palette = colorsphylum, edge.size = 0.4,
          edge.color = E(m4_subnet)$color, label = T, label.size = 4)+
          scale_size_manual(values = c(6, 3))
print(b4)
dev.off()


###SUBNETM5
# Asignar coordenadas al diseño de la red
m5_subnet$layout <- array(1:40, dim = c(20, 2))
# Asignar el diseño como un atributo fijo de la red
m5_subnet$layout <- layout.fruchterman.reingold(m5_subnet)

# Extraer matriz de adyacencia
m5_subnetclass <- as_adjacency_matrix(m5_subnet, type = "both")
# Generar objeto de clase 'network'
m5_subnetclass <- network(as.matrix(m5_subnetclass), 
                          matrix.type = "adjacency", directed = F)
#Plot
ggnet2(m5_subnetclass, mode = m5_subnet$layout)
#Extract taxa names
m5_subnetnodenames <- as.character(getTaxonomy(V(m5_subnet)$name, tax_tbl, level = "family", useRownames = TRUE))
# Identificar y asignar color a los edges según si representan una relación positiva o negativa
edges1 <- E(m5_subnet)
edge_colors <- c()
for(e_index in 1:length(edges1)){
  adj_nodes <- ends(m5_subnet,edges1[e_index])
  xindex <- which(tax_ids==adj_nodes[1])
  yindex <- which(tax_ids==adj_nodes[2])
  beta <- betaMat[xindex,yindex]
  if(beta>0){
    edge_colors=append(edge_colors,"darkgreen") # positive
  }else if(beta<0){
    edge_colors=append(edge_colors,"red3") # negative
  }
}
E(m5_subnet)$color <- edge_colors
vertex_attr(m5_subnet, name="phylum", index = V(m5_subnet)) <- tax_tbl[V(m5_subnet)$name,2]
vertex_attr(m5_subnet, name="hub", index = V(m5_subnet)) <- ifelse(names(V(m5_subnet)) %in% nodos[nodos$hub==1,"otu.id"],"hub","no_hub")

#SubnetM5
pdf(file=paste(paste(t,"sub-Network-M5.pdf")), colormodel = "cmyk", width = 14, height = 12, compress = F)
b5=ggnet2(m5_subnet, node.color = "phylum", mode = m5_subnet$layout, node.size = "hub",
          color = selegoG_nodenames, palette = colorsphylum, edge.size = 0.4,
          edge.color = E(m5_subnet)$color, label = T, label.size = 4)+
    scale_size_manual(values = c(6, 3))
print(b5)
dev.off()


###SUBNETM6
# Asignar coordenadas al diseño de la red
m6_subnet$layout <- array(1:40, dim = c(20, 2))
# Asignar el diseño como un atributo fijo de la red
m6_subnet$layout <- layout.fruchterman.reingold(m6_subnet)

# Extraer matriz de adyacencia
m6_subnetclass <- as_adjacency_matrix(m6_subnet, type = "both")
# Generar objeto de clase 'network'
m6_subnetclass <- network(as.matrix(m6_subnetclass), 
                          matrix.type = "adjacency", directed = F)
#Plot
ggnet2(m6_subnetclass, mode = m6_subnet$layout)
#Extract taxa names
m6_subnetnodenames <- as.character(getTaxonomy(V(m6_subnet)$name, tax_tbl, level = "phylum", useRownames = TRUE))
# Identificar y asignar color a los edges según si representan una relación positiva o negativa
edges1 <- E(m6_subnet)
edge_colors <- c()
for(e_index in 1:length(edges1)){
  adj_nodes <- ends(m6_subnet,edges1[e_index])
  xindex <- which(tax_ids==adj_nodes[1])
  yindex <- which(tax_ids==adj_nodes[2])
  beta <- betaMat[xindex,yindex]
  if(beta>0){
    edge_colors=append(edge_colors,"darkgreen") # positive
  }else if(beta<0){
    edge_colors=append(edge_colors,"red3") # negative
  }
}
E(m6_subnet)$color <- edge_colors
vertex_attr(m6_subnet, name="phylum", index = V(m6_subnet)) <- tax_tbl[V(m6_subnet)$name,2]
vertex_attr(m6_subnet, name="hub", index = V(m6_subnet)) <- ifelse(names(V(m6_subnet)) %in% nodos[nodos$hub==1,"otu.id"],"hub","no_hub")

#SubnetM6
pdf(file=paste(paste(t,"sub-Network-M6.pdf")), colormodel = "cmyk", width = 14, height = 12, compress = F)
b6=ggnet2(m6_subnet, node.color = "phylum", mode = m6_subnet$layout, node.size = "hub",
          color = selegoG_nodenames, palette = colorsphylum, edge.size = 0.4,
          edge.color = E(m6_subnet)$color, label = T, label.size = 4)+
        scale_size_manual(values = c(6, 3))
print(b6)
dev.off()


###SUBNETM7
# Asignar coordenadas al diseño de la red
m7_subnet$layout <- array(1:40, dim = c(20, 2))
# Asignar el diseño como un atributo fijo de la red
m7_subnet$layout <- layout.fruchterman.reingold(m7_subnet)

# Extraer matriz de adyacencia
m7_subnetclass <- as_adjacency_matrix(m7_subnet, type = "both")
# Generar objeto de clase 'network'
m7_subnetclass <- network(as.matrix(m7_subnetclass), 
                          matrix.type = "adjacency", directed = F)
#Plot
ggnet2(m7_subnetclass, mode = m7_subnet$layout)
#Extract taxa names
m6_subnetnodenames <- as.character(getTaxonomy(V(m7_subnet)$name, tax_tbl, level = "phylum", useRownames = TRUE))
# Identificar y asignar color a los edges según si representan una relación positiva o negativa
edges1 <- E(m7_subnet)
edge_colors <- c()
for(e_index in 1:length(edges1)){
  adj_nodes <- ends(m7_subnet,edges1[e_index])
  xindex <- which(tax_ids==adj_nodes[1])
  yindex <- which(tax_ids==adj_nodes[2])
  beta <- betaMat[xindex,yindex]
  if(beta>0){
    edge_colors=append(edge_colors,"darkgreen") # positive
  }else if(beta<0){
    edge_colors=append(edge_colors,"red3") # negative
  }
}
E(m7_subnet)$color <- edge_colors
vertex_attr(m7_subnet, name="phylum", index = V(m7_subnet)) <- tax_tbl[V(m7_subnet)$name,2]
vertex_attr(m7_subnet, name="hub", index = V(m7_subnet)) <- ifelse(names(V(m7_subnet)) %in% nodos[nodos$hub==1,"otu.id"],"hub","no_hub")

#SubnetM6
pdf(file=paste(paste(t,"sub-Network-M6.pdf")), colormodel = "cmyk", width = 14, height = 12, compress = F)
b7=ggnet2(m7_subnet, node.color = "phylum", mode = m7_subnet$layout, node.size = "hub",
          color = selegoG_nodenames, palette = colorsphylum, edge.size = 0.4,
          edge.color = E(m7_subnet)$color, label = T, label.size = 4)+
  scale_size_manual(values = c(6, 3))
print(b7)
dev.off()

}

# make a list of the names of the nodes of interest #Allspores
nodes_of_interest <- c("POTU_706")
# make a list of the names of the nodes of interest #Allspores
#nodes_of_interest <- c("FOTU_1134", "FOTU_1278")

# select the nodes having these names
selnodes <- V(se_net)[wt$membership == 
                          wt$membership[which(wt$names %in% nodes_of_interest)]]


# get their network neighborhood 
selegoV <- ego(se_net, order=1, nodes = selnodes, mode = "all", mindist = 0)
#selegoV <- make_ego_graph(se_net, nodes = selnodes, mode = "all")


# turn the returned list of igraph.vs objects into a graph
selegoG_subnet <- induced_subgraph(se_net,unlist(selegoV))
#selegoG_subnet <- induced_subgraph(se_net,vids = selnodes)


# plot the subgraph
plot(selegoG_subnet,vertex.label=V(selegoG_subnet)$name)

# Asignar coordenadas al diseño de la red
selegoG_subnet$layout <- array(1:30, dim = c(50, 3))
# Asignar el diseño como un atributo fijo de la red
selegoG_subnet$layout <- layout.fruchterman.reingold(selegoG_subnet)
  
#Extrac matrix od adjacency
selegoGclass<- as_adjacency_matrix(selegoG_subnet, type = "both")
  # Generar objeto de clase 'network'
selegoGclass <- network(as.matrix(selegoGclass), 
                            matrix.type = "adjacency", directed = F)


ggnet2(selegoGclass, mode= selegoG_subnet$layout)

#Extract taxa names
selegoG_nodenames <- as.character(getTaxonomy(V(selegoG_subnet)$name, tax_tbl, level = "order", useRownames = TRUE))
#selegoG_nodenames <- as.character(getTaxonomy(V(selegoG_subnet)$name, tax_tbl, level = "domain", useRownames = TRUE))

ggnet2(selegoGclass, mode = selegoG_subnet$layout, color = selegoG_nodenames)

unique(selegoG_nodenames)

# Identificar y asignar color a los edges según si representan una relación positiva o negativa
edges1 <- E(selegoG_subnet)
edge_colors <- c()
for(e_index in 1:length(edges1)){
  adj_nodes <- ends(selegoG_subnet,edges1[e_index])
  xindex <- which(tax_ids==adj_nodes[1])
  yindex <- which(tax_ids==adj_nodes[2])
  beta <- betaMat[xindex,yindex]
  if(beta>0){
    edge_colors=append(edge_colors,"darkgreen") # positive
  }else if(beta<0){
    edge_colors=append(edge_colors,"red3") # negative
  }
}
E(selegoG_subnet)$color <- edge_colors
vertex_attr(selegoG_subnet, name="phylum", index = V(selegoG_subnet)) <- tax_tbl[V(selegoG_subnet)$name,4]
vertex_attr(selegoG_subnet, name="domain", index = V(selegoG_subnet)) <- tax_tbl[V(selegoG_subnet)$name,1]
vertex_attr(selegoG_subnet, name="hub", index = V(selegoG_subnet)) <- ifelse(names(V(selegoG_subnet)) %in% nodos[nodos$hub==1,"otu.id"],"hub","no_hub")

#Final graphic
#pdf(file=paste(paste(t,"sub-Network-Hubs-ARCHA")), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
#a9=ggnet2(selegoG_subnet, node.color = "phylum", mode = "fruchtermanreingold", node.size = "hub",
       #color = selegoG_nodenames, palette = colorsorder, edge.size = 0.2,
       #edge.color = E(selegoG_subnet)$color, label = keystone$TaxID, label.size = 3)+
    #scale_size_manual(values = c(6, 3))
#print(a9)
#dev.off()

pdf(file=paste(paste(t,"sub-Network-Hubs-ARCHA.pdf")), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a9=ggnet2(selegoG_subnet, size= 8, node.color = "phylum", mode = selegoG_subnet$layout, node.size = "hub",
          color = selegoG_nodenames, palette = colorsorder, edge.size = 0.2, alpha = 0.95,
          node.alpha = 1.2, legend.position = "none", node.shape = "domain",
          edge.color = E(selegoG_subnet)$color, label = keystone$TaxID, label.size = 3)+
  scale_size_manual(values = c(6, 3))
print(a9)
dev.off()

pdf(file=paste(paste(t,"sub-Network-Hubs-Nolabel-ARCHA.pdf")), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a10=ggnet2(selegoG_subnet, size = 8, node.color = "phylum", mode = "fruchtermanreingold", node.size = "hub",
          color = selegoG_nodenames, palette = colorsorder, edge.size = 0.2,alpha = 0.95,
          node.alpha = 1.2, legend.position = "none",
          edge.color = E(selegoG_subnet)$color, node.shape = "domain", label.size = 3,
          label = T)+
  scale_size_manual(values = c(6, 3))
print(a10)
dev.off()

pdf(file=paste(paste(t,"sub-Network-Hubs-Label.pdf")), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a11=ggnet2(selegoG_subnet, node.color = "domain", mode = selegoG_subnet$layout, #node.size = "hub",
       color = selegoG_nodenames, palette = colorsdomain, edge.size = 0.5,
       edge.color = E(selegoG_subnet)$color, label = T, label.size = 3)+
  scale_size_manual(values = c(6, 3))
print(a11)
dev.off()

pdf(file=paste(paste(t,"sub-Network-Hubs-NoLabel1.pdf")), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a11=ggnet2(selegoG_subnet, node.color = "domain", mode = selegoG_subnet$layout, #node.size = "hub",
           color = selegoG_nodenames, palette = colorsdomain, edge.size = 0.5,
           edge.color = E(selegoG_subnet)$color, label = F, label.size = 3)+
  scale_size_manual(values = c(6, 3))
print(a11)
dev.off()




plot_network(selegoG_subnet, Bacteria, type = "taxa", color = "order", shape = "domain", 
             label = "selegoG_nodenames", title = t, point_size = 7)

