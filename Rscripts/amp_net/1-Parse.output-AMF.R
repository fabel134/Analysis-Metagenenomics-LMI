setwd(getwd())
#setwd("~/Google-Drive")

library(stringr)
library(splitstackshape)

#Stats by sample
data=readLines("~/abel/analysis/16s_analy/stats/trimmed.out.txt")
out=data.frame(aim.id=gsub(".fastq.gz","",gsub("~/abel/analysis/16s_analy/trimmed/","", gsub("\\./trimmed/", "",grep("fastq.gz",data, value = T)))),
               input.pairs=str_extract(grep("Input Read Pairs",data, value = T),"\\d+"),
               surviving.pairs=str_extract(grep("Both Surviving Reads",data,value = T),"\\d+"),
               surviving.pct=str_extract(grep("Both Surviving Read Percent",data, value = T),"\\d+"))
out <- unique(out[,])

data=readLines("~/abel/analysis/16s_analy/stats/merged.out.txt")
sub=data.frame(sample=str_extract(grep("extendedFrags",data, value = T), "EC24-034-(\\d{4}-16S-block_[A-Z]\\d{2}_IP\\d+_S\\d+)|(\\d+-[A-Z0-9-]+_[A-Z0-9]+_[A-Z0-9]+_S\\d+)"),
#sub=data.frame(sample=str_extract(grep("extendedFrags",data, value = T), "EC24-034-(\\d{4}-FITS_[A-Z]\\d{2}_IP\\d+_S\\d+)|(\\d+-[A-Z0-9-]+_[A-Z0-9]+_[A-Z0-9]+_S\\d+)"),
               merged.pairs=str_extract(grep("Combined pairs:", data, value = T), "\\d+"),
               merged.pct=str_extract(grep("Percent combined:", data, value = T), "\\d+"))
out[out$aim.id %in% sub$sample, c("merged.pairs","merged.pct")]=sub[,c("merged.pairs","merged.pct")]

data=readLines("~/abel/analysis/16s_analy/stats/filter.out.txt")
kept.pairs=str_extract(grep("truncated", data, value = T), "^\\d+")
sub=data.frame(sample=str_extract(grep("EC24", data, value = TRUE), "EC24-034-(\\d{4}-16S-block_[A-Z]\\d{2}_IP\\d+_S\\d+)|(\\d+-[A-Z0-9-]+_[A-Z0-9]+_[A-Z0-9]+_S\\d+)"),
#sub=data.frame(sample=str_extract(grep("EC24", data, value = TRUE), "EC24-034-(\\d{4}-FITS_[A-Z]\\d{2}_IP\\d+_S\\d+)|(\\d+-[A-Z0-9-]+_[A-Z0-9]+_[A-Z0-9]+_S\\d+)"),
               kept.pairs=kept.pairs,
               kept.pct=round(as.numeric(kept.pairs)/as.numeric(sub$merged.pairs)*100,0))
out[out$aim.id %in% sub$sample, c("kept.pairs","kept.pct")]=sub[,c("kept.pairs","kept.pct")]


data=readLines("~/abel/analysis/16s_analy/stats/derep.sample.out.txt")
unique=str_extract(grep("unique", data, value = T),"^\\d+")

sub=data.frame(sample=str_extract(grep("EC24",data, value = T),"EC24-034-(\\d{4}-16S-block_[A-Z]\\d{2}_IP\\d+_S\\d+)|(\\d+-[A-Z0-9-]+_[A-Z0-9]+_[A-Z0-9]+_S\\d+)"),
               unique.pairs=unique,
               unique.pct=round(as.numeric(unique)/as.numeric(sub$kept.pairs)*100,0))
out[out$aim.id %in% sub$sample, c("unique.pairs","unique.pct")]=sub[,c("unique.pairs","unique.pct")]

samples=length(unique)

#rownames(out)=out[,1] ; out=out[,-1]
#out=t(out)

write.table(out,file = "~/abel/analysis/16s_analy/stats/processing.stats.txt", row.names = F, col.names = T, quote = F,
            sep = "\t")

#stats of otus

data=readLines("~/abel/analysis/16s_analy/stats/derep.all.out.txt")
total=str_extract(str_extract(grep("seqs",data, value = T),"\\d+ seqs"),"\\d+")
unique=str_extract(grep("unique sequences",data, value = T),"^\\d+")
unique.pct=round(as.numeric(unique)/as.numeric(total)*100,0)

data=readLines("~/abel/analysis/16s_analy/stats/otu.out.txt")
clusters=str_extract(grep("Clusters",data, value = T),"\\d+")

clusters.i=clusters[1]
clusters.s=clusters[2]
clusters.n=clusters[3]
clusters.b=clusters[4]

data=readLines("~/abel/analysis/16s_analy/stats/table.out.txt")
matching=str_extract(grep("Matching total",data, value = T), "\\d+ \\([^()]+\\)")
to.match=str_extract(matching,"^\\d+")
matching.pct=str_extract(str_extract(matching,"\\([^()]+\\)"),"\\d+\\.\\d+")

data=readLines("~/abel/analysis/16s_analy/stats/taxauniteall.out.txt")
classi=str_extract(grep("Classified ",data, value = T), "\\d+ sequences \\([^()]+\\)")
classi.pct=str_extract(str_extract(classi,"\\([^()]+\\)"),"\\d+\\.\\d+")

data=readLines("~/abel/analysis/16s_analy/stats/taxasilva.out")
pclassi=str_extract(grep("Classified ",data, value = T), "\\d+ sequences \\([^()]+\\)")
pclassi.pct=str_extract(str_extract(pclassi,"\\([^()]+\\)"),"\\d+\\.\\d+")

out.con=data.frame(
  samples=samples,
  unique.per.sample=total,
  unique.total=unique,
  unique.pct=unique.pct,
  otus.total=clusters.i,
  otus.no.singles=clusters.s,
  otus.no.chimeras=clusters.n,
  otus.final=clusters.b,
  reads.to.match=to.match,
  matching.pct=matching.pct,
  classified.otu.pct=classi.pct#,
  #classified.potu.pct=pclassi.pct
)

out.con=t(out.con)

write.table(out.con,file = "~/abel/analysis/16s_analy/stats/clustering.stats.txt", row.names = T, col.names = F, quote = F,
            sep = "\t")

#Taxonomy

#Sintax derived otus
#Formating taxa OTUs already bootstrap filtered by --sintax

sintax=read.delim(file=paste("~/abel/analysis/16s_analy/result","otus.sintax.uniteall",sep = "/"), header = F, stringsAsFactors = F)
sintax=sintax[,-c(2,3,5)]
sintax[,1]=as.numeric(gsub("OTU","",sintax[,1]))
sintax=sintax[order(sintax[,1]),]
sintax[,1]=paste("otu",sintax[,1],sep = "_")

otus=splitstackshape::cSplit(indt = sintax, splitCols = "V4", sep = ",", type.convert = F)
colnames(otus)=c("otu.id","domain","phylum","class","order","family","genus") #you can substract the d: p: c: o: f: g: using gsub but I thing it is usefull for counting and grepping 

taxa=data.frame(
  
#how many different taxa: 
unique.domain=length(unique(otus$domain))-1,
unique.phylum=length(unique(otus$phylum))-1,
unique.class=length(unique(otus$class))-1,
unique.order=length(unique(otus$order))-1,
unique.family=length(unique(otus$family))-1,
unique.genus=length(unique(otus$genus))-1,

#How many 

domain.unclassified=round(sum(is.na(otus$domain))/dim(otus)[1]*100,2),
domain.eukaryota=round(sum(otus$domain=="d:Eukaryota", na.rm = T)/dim(otus)[1]*100,2),
domain.bacteria=round(sum(otus$domain=="d:Bacteria", na.rm = T)/dim(otus)[1]*100,2),
domain.archaea=round(sum(otus$domain=="d:Archaea", na.rm = T)/dim(otus)[1]*100,2),

#Percentage of unclassified OTUS

phylum.unclassified=round(sum(is.na(otus[,phylum])) /dim(otus)[1]*100,2) ,
class.unclassified=round(sum(is.na(otus[,class])) /dim(otus)[1]*100,2),
order.unclassified=round(sum(is.na(otus[,order])) /dim(otus)[1]*100,2) ,
family.unclassified=round(sum(is.na(otus[,family])) /dim(otus)[1]*100,2) ,
genus.unclassified=round(sum(is.na(otus[,genus])) /dim(otus)[1]*100,2) 

)

taxa=t(taxa)
write.table(taxa,file = "~/abel/analysis/16s_analy/classification.stats.txt", row.names = T, col.names = F, quote = F,
            sep = "\t")
##################################################################################
##############################ITS#################################################
##################################################################################
setwd(getwd())

library(stringr)
library(splitstackshape)

#Stats by sample
data=readLines("~/abel/analysis/FITS_analy/stats/trimmed.out.txt")
out=data.frame(aim.id=gsub(".fastq.gz","",gsub("~/abel/analysis/FITS_analy/trimmed/","", gsub("\\./trimmed/", "",grep("fastq.gz",data, value = T)))),
               input.pairs=str_extract(grep("Input Read Pairs",data, value = T),"\\d+"),
               surviving.pairs=str_extract(grep("Both Surviving Reads",data,value = T),"\\d+"),
               surviving.pct=str_extract(grep("Both Surviving Read Percent",data, value = T),"\\d+"))
out <- unique(out[,])

data=readLines("~/abel/analysis/FITS_analy/stats/merged.out.txt")
sub=data.frame(sample=str_extract(grep("extendedFrags",data, value = T), "EC24-034-(\\d{4}-FITS_[A-Z]\\d{2}_IP\\d+_S\\d+)|(\\d+-[A-Z0-9-]+_[A-Z0-9]+_[A-Z0-9]+_S\\d+)"),
               #sub=data.frame(sample=str_extract(grep("extendedFrags",data, value = T), "EC24-034-(\\d{4}-FITS_[A-Z]\\d{2}_IP\\d+_S\\d+)|(\\d+-[A-Z0-9-]+_[A-Z0-9]+_[A-Z0-9]+_S\\d+)"),
               merged.pairs=str_extract(grep("Combined pairs:", data, value = T), "\\d+"),
               merged.pct=str_extract(grep("Percent combined:", data, value = T), "\\d+"))
out[out$aim.id %in% sub$sample, c("merged.pairs","merged.pct")]=sub[,c("merged.pairs","merged.pct")]

data=readLines("~/abel/analysis/FITS_analy/stats/filter.out.txt")
kept.pairs=str_extract(grep("truncated", data, value = T), "^\\d+")
sub=data.frame(sample=str_extract(grep("EC24", data, value = TRUE), "EC24-034-(\\d{4}-FITS_[A-Z]\\d{2}_IP\\d+_S\\d+)|(\\d+-[A-Z0-9-]+_[A-Z0-9]+_[A-Z0-9]+_S\\d+)"),
               #sub=data.frame(sample=str_extract(grep("EC24", data, value = TRUE), "EC24-034-(\\d{4}-FITS_[A-Z]\\d{2}_IP\\d+_S\\d+)|(\\d+-[A-Z0-9-]+_[A-Z0-9]+_[A-Z0-9]+_S\\d+)"),
               kept.pairs=kept.pairs,
               kept.pct=round(as.numeric(kept.pairs)/as.numeric(sub$merged.pairs)*100,0))
out[out$aim.id %in% sub$sample, c("kept.pairs","kept.pct")]=sub[,c("kept.pairs","kept.pct")]

data=readLines("~/abel/analysis/FITS_analy/stats/derep.sample.out.txt")
unique=str_extract(grep("unique", data, value = T),"^\\d+")

sub=data.frame(sample=str_extract(grep("EC24",data, value = T),"EC24-034-(\\d{4}-FITS_[A-Z]\\d{2}_IP\\d+_S\\d+)|(\\d+-[A-Z0-9-]+_[A-Z0-9]+_[A-Z0-9]+_S\\d+)"),
               unique.pairs=unique,
               unique.pct=round(as.numeric(unique)/as.numeric(sub$kept.pairs)*100,0))
out[out$aim.id %in% sub$sample, c("unique.pairs","unique.pct")]=sub[,c("unique.pairs","unique.pct")]

samples=length(unique)

#rownames(out)=out[,1] ; out=out[,-1]
#out=t(out)

write.table(out,file = "~/abel/analysis/FITS_analy/stats/processing.stats.txt", row.names = F, col.names = T, quote = F,
            sep = "\t")

#stats of otus

data=readLines("~/abel/analysis/FITS_analy/stats/derep.all.out.txt")
total=str_extract(str_extract(grep("seqs",data, value = T),"\\d+ seqs"),"\\d+")
unique=str_extract(grep("unique sequences",data, value = T),"^\\d+")
unique.pct=round(as.numeric(unique)/as.numeric(total)*100,0)

data=readLines("~/abel/analysis/FITS_analy/stats/otu.out.txt")
clusters=str_extract(grep("Clusters",data, value = T),"\\d+")

clusters.i=clusters[1]
clusters.s=clusters[2]
clusters.n=clusters[3]
clusters.b=clusters[4]

data=readLines("~/abel/analysis/FITS_analy/stats/table.out.txt")
matching=str_extract(grep("Matching total",data, value = T), "\\d+ \\([^()]+\\)")
to.match=str_extract(matching,"^\\d+")
matching.pct=str_extract(str_extract(matching,"\\([^()]+\\)"),"\\d+\\.\\d+")

data=readLines("~/abel/analysis/FITS_analy/stats/taxauniteall.out.txt")
classi=str_extract(grep("Classified ",data, value = T), "\\d+ sequences \\([^()]+\\)")
classi.pct=str_extract(str_extract(classi,"\\([^()]+\\)"),"\\d+\\.\\d+")

data=readLines("~/abel/analysis/FITS_analy/stats/taxasilva.out")
pclassi=str_extract(grep("Classified ",data, value = T), "\\d+ sequences \\([^()]+\\)")
pclassi.pct=str_extract(str_extract(pclassi,"\\([^()]+\\)"),"\\d+\\.\\d+")

out.con=data.frame(
  samples=samples,
  unique.per.sample=total,
  unique.total=unique,
  unique.pct=unique.pct,
  otus.total=clusters.i,
  otus.no.singles=clusters.s,
  otus.no.chimeras=clusters.n,
  otus.final=clusters.b,
  reads.to.match=to.match,
  matching.pct=matching.pct,
  classified.otu.pct=classi.pct,
  classified.potu.pct=pclassi.pct
)

out.con=t(out.con)

write.table(out.con,file = "~/abel/analysis/FITS_analy/stats/clustering.stats.txt", row.names = T, col.names = F, quote = F,
            sep = "\t")

#Taxonomy

#Sintax derived otus
#Formating taxa OTUs already bootstrap filtered by --sintax

sintax=read.delim(file=paste("~/abel/analysis/FITS_analy/result","otus.sintax.uniteall",sep = "/"), header = F, stringsAsFactors = F)
sintax=sintax[,-c(2,3,5)]
sintax[,1]=as.numeric(gsub("OTU","",sintax[,1]))
sintax=sintax[order(sintax[,1]),]
sintax[,1]=paste("otu",sintax[,1],sep = "_")

otus=splitstackshape::cSplit(indt = sintax, splitCols = "V4", sep = ",", type.convert = F)
colnames(otus)=c("otu.id","domain","phylum","class","order","family","genus", "species") #you can substract the d: p: c: o: f: g: using gsub but I thing it is usefull for counting and grepping 

taxa=data.frame(
  
  #how many different taxa: 
  unique.domain=length(unique(otus$domain))-1,
  unique.phylum=length(unique(otus$phylum))-1,
  unique.class=length(unique(otus$class))-1,
  unique.order=length(unique(otus$order))-1,
  unique.family=length(unique(otus$family))-1,
  unique.genus=length(unique(otus$genus))-1,
  
  #How many 
  
  domain.unclassified=round(sum(is.na(otus$domain))/dim(otus)[1]*100,2),
  domain.eukaryota=round(sum(otus$domain=="d:Eukaryota", na.rm = T)/dim(otus)[1]*100,2),
  domain.bacteria=round(sum(otus$domain=="d:Bacteria", na.rm = T)/dim(otus)[1]*100,2),
  domain.archaea=round(sum(otus$domain=="d:Archaea", na.rm = T)/dim(otus)[1]*100,2),
  
  #Percentage of unclassified OTUS
  
  phylum.unclassified=round(sum(is.na(otus[,phylum])) /dim(otus)[1]*100,2) ,
  class.unclassified=round(sum(is.na(otus[,class])) /dim(otus)[1]*100,2),
  order.unclassified=round(sum(is.na(otus[,order])) /dim(otus)[1]*100,2) ,
  family.unclassified=round(sum(is.na(otus[,family])) /dim(otus)[1]*100,2) ,
  genus.unclassified=round(sum(is.na(otus[,genus])) /dim(otus)[1]*100,2) 
  
)

taxa=t(taxa)
write.table(taxa,file = "~/abel/analysis/FITS_analy/classification.stats.txt", row.names = T, col.names = F, quote = F,
            sep = "\t")
