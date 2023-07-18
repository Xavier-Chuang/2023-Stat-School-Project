rm(list = ls())
library(stringi)
library(seqinr)
library(stringr)
library(ape)
library(philentropy)

setwd("G:/COVID-19 Genomic data anaylsis")


## read reference MN908947 ==========
ref <- read.fasta("data/MN908947.fas",as.string = T)
nchar(ref) # 29903

## read aligned fasta data
fas <- read.fasta("data/n1932_29903.fas",as.string = T)
length(fas) # 1932 strains
table(sapply(1:length(fas), function(x) nchar(fas[[x]])))


## extract 1st site
do.call(c,stri_sub_all(fas, from = 1  , length = 1))

## snv frequency
snv_freq <- sapply(1:nchar(ref), function(x) sum(do.call(c,stri_sub_all(fas, from = x  , length = 1))!= stri_sub_all(ref, from = x  , length = 1)))

plot(1:nchar(ref),snv_freq/length(fas),type = "l", xlab="Position",ylab="Frequency", bty = 'n')

head(snv_freq)
table(do.call(c,stri_sub_all(fas, from = 1  , length = 1)))

summary(snv_freq/length(fas))
sum(snv_freq/length(fas)!=0) #2629

## gap frequency
gap_freq <- sapply(1:nchar(ref), function(x) sum(do.call(c,stri_sub_all(fas, from = x  , length = 1))=="-"))

plot(1:nchar(ref),gap_freq/length(fas),type = "l", xlab="Position",ylab="Frequency", bty = 'n')

## extract snv
idx <- which(snv_freq/length(fas)>=0.01)
idx <- setdiff(idx,c(1:266,29675:29903))

length(idx) # 2135

## extract SNP
data <- do.call(rbind,stri_sub_all(fas, from = idx  , length = 1))
dim(data) # 1932 2135
row.names(data) <- names(fas)
colnames(data) <- paste0("pos",idx)

## compared to reference
data_binary <- t(t(data)!=do.call(c,stri_sub_all(ref, from = idx  , length = 1)))*1
data_ambig <- 1-t(t(data)=="a")-t(t(data)=="t")-t(t(data)=="c")-t(t(data)=="g")
data_binary_rmambig <- data_binary - data_ambig
row.names(data_binary_rmambig) <- row.names(data)
colnames(data_binary_rmambig) <- paste0("pos",idx)

## Pairwise Distance ==========

# names of implemented distance/similarity functions
getDistMethods()


## pearson correlation 
dist <- distance(data_binary_rmambig, method = "pearson")

## jaccard distance
# dist <- distance(data_binary_rmambig, method = "jaccard")


## Hierarchical Clustering ==========
hc <- hclust(as.dist(dist), method = "average")


## Visualization ==========

phylo <- as.phylo(hc)
phylo$tip.label <- row.names(data_binary_rmambig)[hc$order]


data_binary_rmambig <- data_binary_rmambig[match(phylo$tip.label,row.names(data_binary_rmambig)),]

dir.create("result")
pdf("result/plot.pdf", width = 11.7 , height = 8.3)
layout(matrix(c(1,2),ncol=2,byrow = T),widths = c(0.15,0.85))
par(oma=c(0,0,1,0),mar=c(0,0,0,0))


## tree =====

plot(phylo,
     use.edge.length = T,
     show.tip.label = F,
     cex = 0.3,align.tip.label= F,edge.width=0.001,edge.color = "black")

plot(0,0,type="n",xlim=c(0,29903),ylim=c(1,length(fas)), xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '',xaxs = "i")

for(i in 1:dim(data_binary_rmambig)[1]){
  for(j in which(data_binary_rmambig[i,]%in%1)){
    nt <- colnames(data_binary_rmambig)[j]
    nt <- as.numeric(gsub("pos","",nt))
    rect(nt-10,i-0.5,nt+10,(i+0.5),col="black",border = "transparent")
  }
}

dev.off()

