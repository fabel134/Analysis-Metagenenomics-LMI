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
  install.packages("pulsar")

  ## this version works too
  install.packages("huge")
  library(devtools)
  install_github("kingaa/pomp")
  install_github("zdk123/SpiecEasi")
  
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
library(tidyr)
library(Hmisc)
library(rgexf)
library(GGally)
library(Rmisc)
library(FSA)
library(intergraph)
}


#Choose all setting to plots
{tema=theme(axis.text.x = element_text(color="black",size=12, angle=0,hjust=0.5,vjust=1.5, family = "sans" ),
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
t=c("Agave.tequilana", "Agave.potatorum", "Agave.angustifolia","Agave.karwinskii")[1:4]
colores=c("#29854E","#74A557","#B8C466","#FFE081",
          "#A02355","#C55DA1","#DB99EE",
          "#234EA0","#008CCF","#00C6D9")
colnet=c("goldenrod3","#8c96c6", "grey80") #Color for domain hubs 
}

{
  ## Generate OTU table
  setwd("/home/abel/LIM/R-analisis/networks")
  data.dir=("/home/abel/LIM/R-analisis/networks")
  main.dir=getwd()
  #data=c("syncoms","daniel", "all")[2]  #what sub sample would you need
  #s16=("~/LIM/R-analisis")
  #its=("~/Library/Mobile Documents/com~apple~CloudDocs/Documentos/agaveITS/ITS2_1/")  
  #met=read.delim(file=paste(getwd(),paste("metadata.tableC.txt",sep = "."),sep = "/"))
  #met= subset(met, met$treatment != "mock")
  #met= subset(met, met$treatment != "ethyl.isovalerate")
  #f="AMF"
  
  #ITS
  #fotu=read.delim(file=paste(its,paste("otu.table",data,"txt",sep = "."),sep = "/"))
  #colnames(fotu)=gsub("X","",colnames(fotu))
  #rownames(met)=met$AIM.sample.id
  #ftax=read.delim(paste(its,paste("taxa.table",data,"uniteall","txt",sep = "."),sep = "/"))
  ##SUBSET GLOMEROMYCOTA
  #ftax= subset(ftax, ftax$phylum == "p:Glomeromycota", select= c("domain","phylum","class","order","family","genus","otu.id"))
  #fotu=fotu[rownames(ftax),]
  
  #16S OTUS in Spores
  potu=read.delim(file=paste(paste("otu.table","txt",sep = "."),sep = "/"))
  
  #All bacteria
  #potu=read.delim(file=paste(s16,paste("otu.table",data,"txt",sep = "."),sep = "/"))
  #colnames(potu)=gsub("X","",colnames(potu))
  #rownames(met)=met$Sample.id
  ptax=read.delim(paste(paste("taxa.table","txt",sep = "."),sep = "/"))
  
  #LOAD METADATA
  met=read.delim(file=paste(getwd(),paste("metadata.tableC.txt",sep = "."),sep = "/"))
  rownames(met)=met$sample.ID ### hacer que el LP1_contigs sea rowname
  
  #rearrage taxonomy
  #rownames(ftax)=gsub("OTU","FOTU",rownames(ftax))
  #rownames(ptax)=gsub("OTU","POTU",rownames(ptax))
  #ftax$otu.id=rownames(ftax)
  #ptax$otu.id=rownames(ptax)
  
  #tax=rbind(ptax,ftax) # unir dos dataframe 16s/ITS
  ## Para cuando los se tiene NA en los datos
  #ptax$Kingdom[is.na(ptax$Kingdom)] = "Unidentified"
  #ptax$Phylum[is.na(ptax$Phylum )] = "Unidentified"
  #ptax$Class[is.na(ptax$Class)] = "Unidentified"
  #ptax$Order[is.na(ptax$Order)] = "Unidentified"
  #ptax$Family[is.na(ptax$Family)] = "Unidentified"
  #ptax$Genus[is.na(ptax$Genus)] = "Unidentified" 
  #ptax$Species[is.na(ptax$species)] = "Unidentified"
  
  ## Para cuando los se tienen espacion en blanco en los datos
  
  ptax$Kingdom[ptax$Kingdom == ""] <- "K:Unidentified"
  ptax$Phylum[ptax$Phylum == ""] <- "P:Unidentified"
  ptax$Class[ptax$Class == ""] <- "C:Unidentified"
  ptax$Order[ptax$Order == ""] <- "O:Unidentified"
  ptax$Family[ptax$Family == ""] <- "F:Unidentified"
  ptax$Genus[ptax$Genus == ""] <- "G:Unidentified"
  ptax$Species[ptax$Species == ""] <- "S:Unidentified"
  
  ##Agregar a las rowname la palabra OTU
  #rownames(potu) <- sub("^", "POTU_", rownames(potu))
  
  
  #meake sure to have the same samples
  #rownames(fotu)=gsub("OTU","FOTU",rownames(fotu))
  #rownames(potu)=gsub("OTU","POTU",rownames(potu))
  
  #potu=potu[,colnames(potu) %in% colnames(fotu)]
  #fotu=fotu[,colnames(fotu) %in% colnames(potu)]
  #otu=rbind(potu,fotu)
  
  #generate metadata
  met.2=met[colnames(potu),]
  #limit= 1 ; thr=0.70 #limit = minimum read count for abundant OTUs
  limit=1
  ms=2 #minimum sample frequency for abundant OTUs 
  pc=c("root.zone.soil", "rhizosphere", "phyllosphere")[3]

  ps=c("Agave.tequilana", "Agave.potatorum", "Agave.angustifolia","Agave.karwinskii",
       "Agave.convallis") [1:5]
  
  #for ( t in ps){
  #samples=rownames(met.2[met.2$plant.compartment %in% pc &  met.2$plant.specie %in% t,])
  #}
  
  
  for ( t in pc){
  samples=rownames(met.2[met.2$plant.specie %in% ps & met.2$plant.compartment %in% t,])
  }
  
  otus=potu[,samples]
  otus=otus[rowSums(otus>=limit)>=ms,] # de las sumatorioas sea mayor al limite (>83)
 
  rowSums(otus)
  ##Generate matrices
  otu=as.matrix(otus)
  taxon=as.matrix(ptax)
  
  #Save Data
  exit_name = "All.phill"
  save(otu, file = paste(main.dir,paste(t,exit_name,"OTU_table",sep = "_"),sep = "/"))
  save(met.2, file = paste(main.dir,paste(t,exit_name,"metadata_table",sep = "_"),sep = "/"))
  save(taxon, file = paste(main.dir,paste(t,exit_name,"taxonomy_table",sep = "_"),sep = "/"))
  
  All.phill=phyloseq(otu_table(otu, taxa_are_rows = T),tax_table(taxon),sample_data(met))
  
  save(All.phill, file = paste(getwd(),paste(t,exit_name,"phyloseq",sep="_"),sep = "/"))

}

#Load Phyloseq objet
load("phyllosphere_All.phill_phyloseq")

All.phill

#Using spiec.easi
#Using default parameters for A. tequilana, M. geometrizans
#pargs <- list(seed=10010)

#se_mb= spiec.easi(A.teq, method='mb',lambda.min.ratio=1e-1, nlambda=100, 
#                    pulsar.params=pargs, ncores=6)

#using defaut parameter for A. salmiana
se_mb= spiec.easi(All.phill, method='mb',lambda.min.ratio=1e-2, nlambda=5, 
                  pulsar.params=list(rep.num=20, ncores=6))

#se_mb= spiec.easi(A.teq, method='mb',lambda.min.ratio=1e-2, nlambda=20, 
#                  pulsar.params=list(rep.num=50, ncores=6))

getStability(se_mb)

#Save the objet'se_mb'
#saveRDS(se_mb, "ateqAMF_mb.RDS")
#saveRDS(se_mb, "asalAMF_mb.RDS")
#saveRDS(se_mb, "mgeoAMF_mb.RDS")
saveRDS(se_mb, "All.phill_mb.RDS")


#3.1 SPIEC-EASI
setwd("/home/abel/LIM/R-analisis/networks")
# Read network, select the plant 
#se_mb <- readRDS(file = "ateqAMF_mb.RDS")
#se_mb <- readRDS(file = "asalAMF_mb.RDS")
#se_mb <- readRDS(file = "mgeoAMF_mb.RDS")
se_mba <- readRDS(file = "All.phill_mb.RDS")

#We have now fit a network, but since we have only a rough, discrete sampling of networks along the lambda path, we should check how far we are from the target stability threshold (0.05).  
getStability(se_mba)

{
# Add OTU names to rows and columns
# Create igraph objects
se_net <- adj2igraph(getRefit(se_mba), 
                     rmEmptyNodes = TRUE, diag = FALSE, 
                     vertex.attr = list(name = taxa_names(All.phill))) # Usamos el ID de las taxas para nombrar los vértices o nodos de la red
#save igraph
save(se_net, file = paste(paste(t,"network","igraph",sep = "."),sep = "/"))
#Plot the network, you can choose the color by domain, phylum, order etc etc
pdf(file=paste(paste(t,"network.pdf")), colormodel = "cmyk", width = 8.5, height = 11, compress = F)
a=plot_network(se_net, All.phill, type = "Kingdom", color = NULL, shape = NULL, label = NULL, title = t) #print(a)
print(a)

#dev.off()

#Using ggnet2
# Extract matrix of adjacency
net_class <- as_adjacency_matrix(se_net, type = "both")
# Generate network object
net_class <- network(as.matrix(net_class), 
                     vertex.attrnames = taxa_names(All.phill), 
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
betaMat <- as.matrix(symBeta(getOptBeta(se_mba)))
# Calculate the number of positive and negative edges in the network
positive <- length(betaMat[betaMat>0])/2 
negative <- length(betaMat[betaMat<0])/2 
total <- length(betaMat[betaMat!=0])/2
#Determine bacterial-fungal edges
#ver=as_data_frame(net, what = c("edges"))
#one=ver[grep("POTU_",ver$from),]
#two=dim(one[grep("FOTU_",one$to),])[1]
#three=ver[grep("FOTU_",ver$from),]
#four=dim(three[grep("POTU_",three$to),])[1]
#bf=c(bf,c(two+four))
#bfp=c(bfp,c((two+four)/dim(ver)[1]*100))

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
tax_ids <- taxa_names(All.phill) #File phyloseq
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
                     vertex.attrnames = taxa_names(All.phill), 
                     matrix.type = "adjacency", directed = F)

ggnet2(net_class, edge.color = E(net)$color)

}

#Extract taxonomy table
tax_tbl <- as.data.frame(All.phill@tax_table@.Data)

#Replace the name of each node with its Phylum, you can set each level taxa

V(net)$name <- as.character(getTaxonomy(V(net)$name, tax_tbl, level = "phylum", useRownames = TRUE))


#Save Phylum list per node in 'nodenames
nodenames <- V(net)$name
# Extract matrix of adjacency
net_class <- as_adjacency_matrix(net, type = "both")
# Generate network object
net_class <- network(as.matrix(net_class), 
                     vertex.attrnames = taxa_names(AMF), 
                     matrix.type = "adjacency", directed = F)

#Check unique phylum
unique(tax_tbl$phylum)
#color set for Phylum
colorsphylum=c("Proteobacteria"="#6972D8", "Actinobacteria"="#bf360c","Tenericutes"="#007FD8", 
          "Bacteroidetes"= "#00837A","Acidobacteria"= "#F9F871","Gemmatimonadetes"="#FFCC5E",
          "Fusobacteria"= "#009688","Chloroflexi"="#FF9800","Verrucomicrobia"="#795548", 
          "Planctomycetes"= "#91003f","Saccharibacteria" = "#4a148c","Deinococcus-Thermus"="#845A12",
          "Gemmatimonadetes"="#03A9F4", "Hydrogenedentes"="#9C27B0","Glomeromycota"= "#10441A",
          "P:Unidentified" = "grey20", "Firmicutes"= "#00BCD4", "Chordata"="#2196F3", "Cyanobacteria" = "#FFEB3B", 
          "Bacteroidota"= "#607D8B","Planctomycetota"= "#F44336", "Spirochaetes"= "#c7e9b4", "Thermodesulfobacteria" = "#4CAF50",
          "Thermodesulfobacteria" = "#6EFF33" , "Uroviricota" = "#9599CD")

#color set for domain/kingdom
unique(tax_tbl$domain)
colorsdomain=c("k:Bacteria"= "royalblue", "d:Fungi"= "gold3","d:Unidentified"="grey20")

#Network for phlyum
pdf(file=paste(paste(t,"phylum.pdf")), colormodel = "srgb", width = 11, height = 8.5, compress = F)
a1=ggnet2(net_class, size= 6, color = nodenames, palette = colorsphylum, edge.color = E(net)$color, edge.size = 0.5, tittle = t)
print(a1)
dev.off()
#network for domain
pdf(file=paste(paste(t,"domain.pdf")), colormodel = "srgb", width = 8, height = 8, compress = F)
a2=ggnet2(net_class, size= 3, color = nodenames, palette = colorsdomain, edge.color = E(net)$color, edge.size = 0.5, title = t)
print(a2)
dev.off()

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
wt <- fastgreedy.community(net)
# Consultar membresía de cada nodo
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
  scale_x_continuous(breaks=c(0,2,4,6,8,10)) + 
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
  scale_x_continuous(breaks=c(0,2,4,6,8,10)) + 
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
deg_high_df <- dplyr::filter(deg_df, deg_sort > 7)
head(deg_high_df)
#Edit betweenness table
bn_df$TaxID <- row.names(bn_df)
# Filtramos las taxas con betweenness > 700 (ateq), >900 (asal), >600 (spores), 500(mgeo)
bn_high_df <- dplyr::filter(bn_df, bn_sort > 550)
head(bn_high_df)
# #Select the hub OTUS and las tablas de degree y betweenness filtradas
keystone <- merge(deg_high_df, bn_high_df, all.x = FALSE)
print(keystone)

# Usamos ggplot2
ggplot(keystone, aes(x = bn_sort, y = deg_sort, label = TaxID)) + 
  scale_y_continuous(limits = c(4,10), breaks = c(4,6,8,10)) + 
  scale_x_continuous(limits = c(800, 3700), breaks = c(800, 1600, 2400, 3200, 3700)) + 
  geom_point(alpha = 4, color = "#829FD9", size = 6) + 
  geom_text(size = 4) + 
  theme_minimal() + 
  labs(x = "Betweenness", y = "Degree", 
       title = "Nodes With Highest Degree And Betweenness: Keystone Nodes")

# Graficar red
ggnet2(se_net, mode = net$layout, edge.size = 0.4,
       color = nodenames, palette = colorsdomain, 
       node.alpha = 0.9, 
       edge.color = edge_colors, 
       label = keystone$TaxID, label.size = 6)

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
  nodos$plant=t
  vertex_attr(se_net, name="kingdom", index = V(se_net)) <- tax_tbl[V(se_net)$name,1]
  vertex_attr(se_net, name="hub", index = V(se_net)) <- ifelse(names(V(se_net)) %in% nodos[nodos$hub==1,"otu.id"],"hub","no_hub")
  colnet=c("goldenrod3","#8c96c6", "grey80")
  edge_attr(se_net, name="peso",index=E(se_net)) <- ifelse(E(se_net)$weight>0, "darkgreen", "red3")

  #table of hub
  hub=nodos[nodos$hub==1,]
  write.table(hub, file = paste(paste(t,"hub","txt",sep = "."),sep = "/"),
              row.names = F, quote = F, sep = "\t", col.names = T)
#graph domain
pdf(file=paste(paste(t,"Network-Hubs-D.pdf")), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a7=ggnet2(se_net, size= 4, node.color = "kingdom",  edge.color = edge_colors, node.size = "hub",edge.size = 0.4,
         mode="fruchtermanreingold", label = keystone$TaxID, label.size = 4)+
    geom_point(aes(fill=color,size=size),shape=21)+
    scale_fill_manual(values =  colnet)+
    scale_size_manual(values = c(5,3))+
    ggtitle(t)
print(a7)
dev.off()

ggnet2(se_net, mode = net$layout, edge.size = 0.4, node.size = "hub", 
       node.color = "kingdom", palette = colorsdomain ,edge.color = edge_colors, 
       label = keystone$TaxID, label.size = 4, node.shape = "kingdom")+ 
  scale_size_manual(values = c(5, 3))+
  scale_fill_manual(values =  colnet)


pdf(file=paste(paste(t,"Network-Hubs-D-NL.pdf")), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a8=ggnet2(se_net, node.color = "kingdom",  edge.color = edge_colors, node.size = "hub", edge.size = 0.4,
          mode="fruchtermanreingold", label = F)+
  geom_point(aes(fill=color,size=size),shape=21)+
  scale_fill_manual(values =  colnet)+
  scale_size_manual(values = c(5,3))+
  ggtitle(t)
print(a8)
dev.off()
    
# graph red##
#vertex_attr(se_net, name="phylum", index = V(se_net)) <- tax_tbl[V(se_net)$name,2]
vertex_attr(se_net, name="phylum", index = V(se_net)) <- tax_tbl[V(se_net)$name,1]
vertex_attr(se_net, name="hub", index = V(se_net)) <- ifelse(names(V(se_net)) %in% nodos[nodos$hub==1,"otu.id"],"hub","no_hub")
#Check if you choose domain-phylum, orden. etc etc 
pdf(file=paste(paste(t,"Network-Hubs-P.pdf")), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a81=ggnet2(se_net, mode = net$layout, edge.size = 0.4, node.size = "hub", 
          node.color = "kingdom", palette = colorsdomain,edge.color = edge_colors, 
          label = keystone$TaxID, label.size = 4)+ 
          scale_size_manual(values = c(5, 3))+
          scale_fill_manual(values =  colnet)
print(a81)
dev.off()

pdf(file=paste(paste(t,"Network-Hubs-P-NL.pdf")), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a82=ggnet2(se_net, size = 12, mode = net$layout, edge.size = 0.5, node.size = "hub",
       node.color = "kingdom", palette = colorsdomain,
       edge.color = edge_colors,  label = F)+
      scale_size_manual(values = c(5, 3))
print(a82)
dev.off()


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
# make a list of the names of the nodes of interest #Ateq
#nodes_of_interest <- c("FOTU_87")
# make a list of the names of the nodes of interest #Asal
#nodes_of_interest <- c("FOTU_140")
# make a list of the names of the nodes of interest #Mgeo
#nodes_of_interest <- c("FOTU_490")
# make a list of the names of the nodes of interest #Allspores OTU140
nodes_of_interest <- c("FOTU_140", "POTU_6345")
# make a list of the names of the nodes of interest #Allspores
#nodes_of_interest <- c("POTU_1263", "POTU_3160")


# select the nodes having these names
selnodes <- V(se_net)[wt$membership == 
                          wt$membership[which(wt$names %in% nodes_of_interest)]]


# get their network neighborhood 
selegoV <- ego(se_net, order=1, nodes = selnodes, mode = "all", mindist = 0)


# turn the returned list of igraph.vs objects into a graph
#selegoG_subnet <- induced_subgraph(se_net,unlist(selegoV))
selegoG_subnet <- induced_subgraph(se_net,vids = selnodes)


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
#selegoG_nodenames <- as.character(getTaxonomy(V(selegoG_subnet)$name, tax_tbl, level = "phylum", useRownames = TRUE))
selegoG_nodenames <- as.character(getTaxonomy(V(selegoG_subnet)$name, tax_tbl, level = "domain", useRownames = TRUE))

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
#vertex_attr(selegoG_subnet, name="phylum", index = V(selegoG_subnet)) <- tax_tbl[V(selegoG_subnet)$name,2]
vertex_attr(selegoG_subnet, name="domain", index = V(selegoG_subnet)) <- tax_tbl[V(selegoG_subnet)$name,1]
vertex_attr(selegoG_subnet, name="hub", index = V(selegoG_subnet)) <- ifelse(names(V(selegoG_subnet)) %in% nodos[nodos$hub==1,"otu.id"],"hub","no_hub")

#Final graphic
pdf(file=paste(paste(t,"sub-Network-Hubs.pdf")), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a9=ggnet2(selegoG_subnet, node.color = "domain", mode = selegoG_subnet$layout, node.size = "hub",
       color = selegoG_nodenames, palette = colorsdomain, edge.size = 0.4,
       edge.color = E(selegoG_subnet)$color, label = keystone$TaxID, label.size = 3)+
    scale_size_manual(values = c(6, 3))
print(a9)
dev.off()

pdf(file=paste(paste(t,"sub-Network-Hubs-Nolabel.pdf")), colormodel = "cmyk", width = 11, height = 8.5, compress = F)
a10=ggnet2(selegoG_subnet, node.color = "domain", mode = selegoG_subnet$layout, #node.size = "hub",
          color = selegoG_nodenames, palette = colorsdomain, 
          node.alpha = 1.2, 
          edge.color = E(selegoG_subnet)$color, 
          label = F)+
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



plot_network(selegoG_subnet, AMF, type = "taxa", color = "domain", shape = "domain", 
             label = "selegoG_nodenames", title = t, point_size = 7)

