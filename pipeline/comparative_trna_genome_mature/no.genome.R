d<-read.table("yes.mod.cca.genomic",sep="\t")
head(d)
exp<-rowMeans(d[,3:20])
len<-unlist(lapply(as.character(d$V2),function(x){
  length(unlist(strsplit(x,split="")))
}))

data <- data.frame(average=exp,size=len)
library(ggplot2)
ggplot(data,aes(x=len,y=exp))+
  geom_point()
tot<-read.table("fasta/define.cluster.Dec/res/ann.tab",header=T,sep="\t")
tot.counts<-colSums(tot[grepl("tRNA",tot$ann),3:30])
sum(exp)/mean(tot.counts)
summary(exp)