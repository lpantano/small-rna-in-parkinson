library(data.table)
library(DESeq2)
library(qvalue)
library(ggplot2)
library("RColorBrewer")
library("gplots")

get_accuracy<-function(hr,c,t){
  clus<-cutree(hr,2)
  g1<-names(clus[clus==1])
  ng1.cc<-sum(sub("[0-9]","",(g1))== c )
  ng1.ct<-sum(sub("[0-9]","",(g1))== t )
  ng1<-(max(ng1.cc,ng1.ct))/length(g1)
  g2<-names(clus[clus==2])
  ng2.cc<-sum(sub("[0-9]","",(g2))== c )
  ng2.ct<-sum(sub("[0-9]","",(g2))== t )
  ng2<-(max(ng2.cc,ng2.ct))/length(g2)
  a<-min(ng1,ng2)
  return(a)
}


table<-read.table("../fasta/define.cluster.Dec/res/ann.tab",sep="\t",header=T,row.names=1)
ann<-table[,2]
names(ann)<-row.names(table)
table<-table[,2:(ncol(table)-1)]

idx<-read.table("../fasta/define.cluster.Dec/res/clusters.clean2.tab",header=T,row.names=1)
table<-table[row.names(idx),]

con<-"cc"
treat<-"ct"
ini<-15
end<-28
name<-"clus"

#########################################################
design<-data.frame(condition=c(rep(con,7),rep(treat,7)))
row.names(design)<-names(table)[ini:end]

mirna<-table[ini:end]
mirna<-mirna[apply(mirna>10,1,sum)>=5,]
mirna<-mirna[!is.na(mirna[,1]),]

#write.table(row.names(mirna),"tables.res/ids.ccvsct.filter.txt",row.names=F,quote=F,col.names=F)

dds <- DESeqDataSetFromMatrix(countData = mirna,
                              colData = design,
                              design = ~ condition)
dds <- estimateSizeFactors( dds )
save(dds,file=paste(treat,"vs",con,".dss",sep="\t"))

dds<-DESeq(dds)

res<-results(dds,independentFiltering=FALSE)
res.all<-mcols(dds,use.names=TRUE)
res.dt<-as.data.frame(res)
res.dt<-res.dt[!is.na(res.dt$padj),]
pdf(paste("heatmap.res/",unique(design[,1])[2],"vs",unique(design[,1])[1],sep="",
          ".pvalue.clu.pdf"))
hist(res.dt$pvalue,main="pvalue distribution")
dev.off()

head(res.dt[order(res.dt$pvalue),],13)
q<-qvalue(res.dt$pvalue[!is.na(res.dt$pvalue)],
          fdr.level=0.3)
#summary(q)
rld <- rlogTransformation(dds, blind=TRUE)
save(rld,file=paste(treat,"vs",con,".rld",sep="\t"))

idx.lim<-max(which(sort(res.dt$pvalue)<0.06))
hmcol<- colorRampPalette(brewer.pal(9,"YlGnBu"))(256)
max.acc<-0
for (lim in seq(10,idx.lim,5)){
  select<-row.names((res.dt[order(res.dt$pvalue),])[1:lim,])
  cor_t <- 1 - cor((assay(rld)[select,]))
  hr<-hclust(as.dist(cor_t),method="ward")
  cor_tc <- 1 - cor(t((assay(rld)[select,])))
  hc<-hclust(as.dist(cor_tc),method="ward")
  par(mar=c(3,5,5,2),cex=.5)
  pdf(paste("heatmap.res/",unique(design[,1])[2],"vs",unique(design[,1])[1],sep="",
           ".top",lim,".",name,".pdf"))
  heatmap.2(assay(rld)[select,], col = hmcol,
            Colv=as.dendrogram(hr), Rowv=as.dendrogram(hc),
            scale="row",
            dendrogram="col", trace="none",cexRow=.7,
            ColSideColors=c(rep("gray",7),rep("yellow",7)),
            ,margins=c(7,7),keysize=1.5,main=paste("top: ",lim))
  dev.off()
  acc<-get_accuracy(hr,con,treat)
  if (acc>=max.acc){
    max.heatmap<-lim
    max.acc<-acc
  }
}

write.table(cbind(res.dt,assay(rld)[row.names(res.dt),],ann[row.names(res.dt)]),paste("tables.res/DEseq.",name,".",unique(design[,1])[2],"vs",unique(design[,1])[1],sep="",".txt"),quote=F,sep="\t")
realvalue<-max.acc

#permutation of micros but same DE: 22 87%
acc<-vector()
for (perm in 1:400){
  select<-row.names((res.dt[sample(1:nrow(res.dt),max.heatmap),])[1:max.heatmap,])  
  cor_t <- 1 - cor((assay(rld)[select,]))
  hr<-hclust(as.dist(cor_t),method="ward")
  acc<-c(get_accuracy(hr,con,treat),acc)
}
pval=sum(acc>=realvalue)/401
d<-density(acc,bw=0.1)
pdf(paste("heatmap.res/",unique(design[,1])[2],"vs",unique(design[,1])[1],sep="",
          ".signif.random.",name,".pdf"))
plot(d, main="density random score",xlab="clustering scores",xlim=c(0,1.1))
polygon(d, col="steelblue", border="blue")
abline(v=realvalue,col="red",lwd=2)
text(1,3,labels=paste("p-value:",round(pval,digits=2)),cex=0.7)
dev.off()


rand1<-combinations(7,4)
rand2<-combinations(7,3)
acc<-vector()
idx<-0
statistics<-vector("list")
for (row in sample(1:35,20)){
  for (row2 in sample(1:35,20)){
    idx<-idx+1
    design.r<-design
    design.r$lib<-1    
    design.r$condition<-treat
    design.r$condition[rand1[row,]]<-con
    design.r$condition[rand2[row2,]+7]<-con
    row.names(design.r[design.r$condition==con,])<-paste(con,1:7,sep="")
    row.names(design.r[design.r$condition==treat,])<-paste(treat,1:7,sep="")
    
    dds <- DESeqDataSetFromMatrix(countData = mirna,
                                  colData = design.r,
                                  design = ~ condition)
    dds<-DESeq(dds,quiet=T)
    res<-results(dds,independentFiltering=FALSE)
    res.dt<-as.data.frame(res)
    res.dt<-res.dt[!is.na(res.dt$padj),]
    statistics[[idx]]<-(mcols(dds,use.names=TRUE))[,17]
    rld <- rlogTransformation(dds, blind=TRUE)
    select<-row.names((res.dt[order(res.dt$pvalue),])[1:max.heatmap,])  
    cor_t <- 1 - cor((assay(rld)[select,]))
    hr<-hclust(as.dist(cor_t),method="ward")
    acc<-c(get_accuracy(hr,con,treat),acc)
    
  }
}
save(statistics,file=paste("tables.res/",unique(design[,1])[2],"vs",unique(design[,1])[1],sep="",
                           name,".DE.statistics.Robj"))

pval=sum(acc>=realvalue)/401
d<-density(acc,bw=0.1)
pdf(paste("heatmap.res/",unique(design[,1])[2],"vs",unique(design[,1])[1],sep="",
          ".signif.random.",name,".DE.pdf"))
plot(d, main="density random score",xlab="clustering scores",xlim=c(0,1.1))
polygon(d, col="steelblue", border="blue")
abline(v=realvalue,col="red",lwd=2)
text(1,3,labels=paste("p-value:",round(pval,digits=2)),cex=0.7)
dev.off()
