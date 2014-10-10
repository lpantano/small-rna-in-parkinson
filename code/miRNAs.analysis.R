library(data.table)
library(DESeq2)
library(qvalue)
library(ggplot2)
library("RColorBrewer")
library("gplots")
library(gtools)

heatmap.dir<-"heatmaps.res"
tables.dir<-"tables.res"

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

config<-read.table("config.2",sep="\t")
table<-data.frame(mir="hsa")
#for (a in  list.files(".",pattern="mirna.out.f")){
for (a in 1:nrow(config)) {
  file<-paste(config[a,1],".clean.fastq.ad.fa.mirna.out.f",sep="")
  d<-read.table(file,sep="\t",header=F,skip=1)
  d.f<-d[d$V14==1,c(2,3,4)]
  #d.f<-d.f[d.f$V3>=5,]  
  table.dt<-data.table(d)
  tabmir<-as.data.frame(table.dt[,list(freq=sum(V3)),by=V4])
  d.f.m<-merge(d.f[,c(3,1,2)],tabmir,by=1)
  #d.f.m$r<-d.f.m$V3/d.f.m$freq
  #d.f.m<-d.f.m[d.f.m$r>=0.1,]
  tabmir.f<-as.data.frame(data.table(d.f.m)[,list(freq=sum(V3)),by=V4])
  names(tabmir.f)<-c("mirna",as.character(config[a,2]))
  table<-merge(table,tabmir.f,by=1,all=TRUE)
}
table<-table[2:nrow(table),]
row.names(table)<-table[,1]
table<-table[,2:ncol(table)]
table[is.na(table)]<-0

con<-"pc"
treat<-"pt"
ini<-1
end<-14
#########################################################
design<-data.frame(condition=c(rep(con,7),rep(treat,7)))
row.names(design)<-names(table)[ini:end]
mirna<-table[,ini:end]
mirna<-mirna[apply(mirna>10,1,sum)>=5,]
mirna<-mirna[!is.na(mirna[,1]),]
dds <- DESeqDataSetFromMatrix(countData = mirna,
                              colData = design,
                              design = ~ condition)
dds <- estimateSizeFactors( dds )

dds<-DESeq(dds)
#dds<-nbinomWaldTest(dds,cooksCutoff=FALSE)
#dds<-replaceOutliersWithTrimmedMean(dds)

res<-results(dds,independentFiltering=FALSE)
#res<-results(dds)
res.all<-mcols(dds,use.names=TRUE)
#names(res.all)
res.dt<-as.data.frame(res)
res.dt<-res.dt[!is.na(res.dt$padj),]
head(res.dt[order(res.dt$pvalue),],21)
q<-qvalue(res.dt$pvalue[!is.na(res.dt$pvalue)],
          fdr.level=0.3)
#summary(q)
pdf(paste(heatmap.dir,"/",unique(design[,1])[2],"vs",unique(design[,1])[1],sep="",
          ".pvalue.mir.pdf"))
hist(res.dt$pvalue,main="pvalue distribution")
dev.off()
#plotDispEsts(dds)
#plotMA(dds)
rld <- rlogTransformation(dds, blind=TRUE)
save(rld,file=paste(treat,"vs",con,".rld",sep="\t"))
idx.lim<-max(which(sort(res.dt$pvalue)<0.06))

#head(res.dt)
#head(assay(rld))
hmcol<- colorRampPalette(brewer.pal(9,"YlGnBu"))(256)
max.acc<-0
for (lim in seq(10,idx.lim,2)){
  select<-row.names((res.dt[order(res.dt$pvalue),])[1:lim,])
  cor_t <- 1 - cor((assay(rld)[select,]))
  hr<-hclust(as.dist(cor_t),method="ward")
  cor_tc <- 1 - cor(t((assay(rld)[select,])))
  hc<-hclust(as.dist(cor_tc),method="ward")
  par(mar=c(3,5,5,2),cex=.5)
  pdf(paste(heatmap.dir,"/",unique(design[,1])[2],"vs",unique(design[,1])[1],sep="",
            ".top",lim,".mirna.pdf"))
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

write.table(cbind(res.dt,assay(rld)[row.names(res.dt),]),paste(tables.dir,"/DEseq.miraligner.mirna.",unique(design[,1])[2],"vs",unique(design[,1])[1],sep="",".txt"),quote=F,sep="\t")
realvalue<-max.acc
  
#permutation of micros but same DE: 22 87%
acc<-vector()
for (perm in 1:400){
  select<-row.names((res.dt[sample(1:nrow(res.dt),max.heatmap),]))  
  cor_t <- 1 - cor((assay(rld)[select,]))
  hr<-hclust(as.dist(cor_t),method="ward")
  acc<-c(get_accuracy(hr,con,treat),acc)
}
pval=sum(acc>=realvalue)/401
d<-density(acc,bw=0.1)
pdf(paste(heatmap.dir,"/",unique(design[,1])[2],"vs",unique(design[,1])[1],sep="",
          ".signif.random.mirna.pdf"))
 plot(d, main=paste("density random score with: ",max.heatmap),xlab="clustering scores",xlim=c(0,1.1))
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
   #dds<-nbinomWaldTest(dds)
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
save(statistics,file=paste(tables.dir,"/",unique(design[,1])[2],"vs",unique(design[,1])[1],sep="",
                           "mirna.DE.statistics.Robj"))

pval=sum(acc>=realvalue)/401
d<-density(acc,bw=0.1)
pdf(paste(heatmap.dir,"/",unique(design[,1])[2],"vs",unique(design[,1])[1],sep="",
          ".signif.random.DE.pdf"))
plot(d, main="density random score",xlab="clustering scores",xlim=c(0,1.1))
polygon(d, col="steelblue", border="blue")
abline(v=realvalue,col="red",lwd=2)
text(1,3,labels=paste("p-value:",round(pval,digits=2)),cex=0.7)
dev.off()


###########isomirs#################
###########isomirs#################
###########isomirs#################

config<-read.table("config.2",sep="\t")
table<-data.frame(mir="hsa")
ids<-read.table("ids.isomirs2")
ids<-ids[order(ids$V1),]
ids<-ids[!duplicated(ids$V1),]
row.names(ids)<-ids$V1
ids<-ids[,2:3]
#for (a in 1:nrow(config)) {
for (a in 1:28) {
  print(a)
  file<-paste(config[a,1],".clean.fastq.ad.fa.mirna.out.f",sep="")
  d<-read.table(file,sep="\t",header=F,skip=1)
  d.f<-d[d$V14==1,]
  d.f<-d.f[d.f$V3>=3,]  
  d.f$id<-paste(d.f$V4,d.f$V7,d.f$V8,d.f$V9,d.f$V10,sep=".")
  tabmir<-d.f[,c("V1","V3")]
  #table.dt<-data.table(d)
  #tabmir<-as.data.frame(table.dt[,list(freq=sum(V3)),by=V4])
  #d.f.m<-merge(d.f[,c(3,1,2)],tabmir,by=1)
  #d.f.m$r<-d.f.m$V3/d.f.m$freq
  #d.f.m<-d.f.m[d.f.m$r>=0.1,]
  #tabmir.f<-as.data.frame(data.table(d.f.m)[,list(freq=sum(V3)),by=V4])
  names(tabmir)<-c("isomir",as.character(config[a,2]))
  table<-merge(table,tabmir,by=1,all=TRUE)
}
table<-table[2:nrow(table),]
row.names(table)<-table[,1]
table<-table[,2:ncol(table)]
table[is.na(table)]<-0

con<-"cc"
treat<-"ct"
ini<-15
end<-28

#########################################################
design<-data.frame(condition=c(rep(con,7),rep(treat,7)))
row.names(design)<-names(table)[ini:end]
mirna<-table[,ini:end]
mirna<-mirna[apply(mirna>10,1,sum)>=5,]
mirna<-mirna[!is.na(mirna[,1]),]

dds <- DESeqDataSetFromMatrix(countData = mirna,
                              colData = design,
                              design = ~ condition)
dds <- estimateSizeFactors( dds )

dds<-DESeq(dds)
#dds<-nbinomWaldTest(dds,cooksCutoff=FALSE)
#res<-results(dds,independentFiltering=FALSE,cooksCutoff=FALSE)
res<-results(dds,independentFiltering=FALSE)
res.all<-mcols(dds,use.names=TRUE)
#names(res.all)
res.dt<-as.data.frame(res)
res.dt<-res.dt[!is.na(res.dt$padj),]
res.dt$gene<-paste(ids[row.names(res.dt),1],ids[row.names(res.dt),2])
head(res.dt[order(res.dt$pvalue),],21)
q<-qvalue(res.dt$pvalue[!is.na(res.dt$pvalue)],
          fdr.level=1, pi0.method="smoother",robust=TRUE)
#summary(q)
pdf(paste(heatmap.dir,"/",unique(design[,1])[2],"vs",unique(design[,1])[1],sep="",
          ".pvalue.isomir.pdf"))
hist(res.dt$pvalue,main="pvalue distribution")
dev.off()
#plotDispEsts(dds)
#plotMA(dds)
rld <- rlogTransformation(dds, blind=TRUE)
#vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
idx.lim<-max(which(sort(res.dt$pvalue)<=0.05))

#head(res.dt)
#head(assay(rld))
hmcol<- colorRampPalette(brewer.pal(9,"YlGnBu"))(256)
max.acc<-0
for (lim in seq(10,idx.lim,5)){
  select<-row.names((res.dt[order(res.dt$pvalue),])[1:lim,])
  cor_t <- 1 - cor((assay(rld)[select,]))
  hr<-hclust(as.dist(cor_t),method="ward")
  cor_tc <- 1 - cor(t((assay(rld)[select,])))
  hc<-hclust(as.dist(cor_tc),method="ward")
  par(mar=c(3,5,5,2),cex=.5)
  pdf(paste(heatmap.dir,"/",unique(design[,1])[2],"vs",unique(design[,1])[1],sep="",
            ".top",lim,".isomir.pdf"))
  heatmap.2(assay(rld)[select,], col = hmcol,
            Colv=as.dendrogram(hr),Rowv=as.dendrogram(hc), scale="row",labRow="",
            dendrogram="col", trace="none",cexRow=.7,
            ColSideColors=c(rep("gray",7),rep("yellow",7)),
            ,margins=c(7,7),keysize=1.5,main=paste("top: ",lim))
  dev.off()
  ids.sort<-(paste(ids[select,1],ids[select,2]))[hc$order]
  write.table(ids.sort,paste(tables.dir,"/ids.sort.miraligner.isomir.",lim,".",unique(design[,1])[2],"vs",unique(design[,1])[1],sep="",".txt"),quote=F,sep="\t")
  
  acc<-get_accuracy(hr,con,treat)
  if (acc>=max.acc){
    max.heatmap<-lim
    max.acc<-acc
  }
}
realvalue<-max.acc

write.table(cbind(ids[row.names(res.dt),],res.dt,assay(rld)[row.names(res.dt),]),paste(tables.dir,"/DEseq.miraligner.isomir.",unique(design[,1])[2],"vs",unique(design[,1])[1],sep="",".txt"),quote=F,sep="\t")

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
pdf(paste(heatmap.dir,"/",unique(design[,1])[2],"vs",unique(design[,1])[1],sep="",
          ".signif.random.isomir.pdf"))
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
    #dds<-nbinomWaldTest(dds,cooksCutoff=FALSE)
    res<-results(dds,independentFiltering=FALSE)
    res.dt<-as.data.frame(res)
    res.dt<-res.dt[!is.na(res.dt$padj),]    
    res.dt<-as.data.frame(res)
    statistics[[idx]]<-(mcols(dds,use.names=TRUE))[,17]
    rld <- rlogTransformation(dds, blind=TRUE)
    select<-row.names((res.dt[order(res.dt$pvalue),])[1:max.heatmap,])  
    cor_t <- 1 - cor((assay(rld)[select,]))
    hr<-hclust(as.dist(cor_t),method="ward")
    acc<-c(get_accuracy(hr,con,treat),acc)
    
  }
}
save(statistics,file=paste(tables.dir,"/",unique(design[,1])[2],"vs",unique(design[,1])[1],sep="",
                           "isomir.DE.statistics.Robj"))

pval=sum(acc>=realvalue)/(length(acc)+1)
d<-density(acc,bw=0.1)
pdf(paste(heatmap.dir,"/",unique(design[,1])[2],"vs",unique(design[,1])[1],sep="",
          ".isomir.signif.random.DE.pdf"))
plot(d, main="density random score",xlab="clustering scores")
polygon(d, col="steelblue", border="blue")
abline(v=realvalue,col="red",lwd=2)
text(1,3,labels=paste("p-value:",round(pval,digits=2)),cex=0.7)
dev.off()



####list mirna with eu list
####list mirna with eu list
####list mirna with eu list
####list mirna with eu list
####list mirna with eu list

for (i in select){
  table.counts<-melt(assay(rld)[i,])
  table.counts$group<-design$condition
  p<-ggplot(table.counts, aes(x=factor(group),y=value,fill=factor(group) ) )+
    geom_boxplot(outlier.size = 0)   + geom_jitter(position=position_jitter(width=0.2),aes(factor(group), value,colour=group)) +
    scale_fill_manual("",values=c("tomato2","SlateBlue3"))+
    theme_bw(base_size = 14, base_family = "") +
    theme(axis.text.x=element_text(angle = 60, hjust = 1))+
    labs(list(title=i,y="log2(norm_counts)",x=""))
  png(paste(sep="",i,".cc.vs.ct.png"))
  print(p)
  dev.off()
}

de<-read.table("../list.pd.vs.c.eulalia")
#select<-setdiff(intersect(as.character(de[,1]),row.names(assay(rld))),c("hsa-miR-7-5p","hsa-miR-146a-5p","hsa-miR-21-5p","hsa-miR-577","hsa-miR-584-5p","hsa-miR-181b-5p","hsa-miR-181c-5p","hsa-miR-181a-5p"))
select<-intersect(as.character(de[,1]),row.names(assay(rld)))
cor_t <- 1 - cor((assay(rld)[select,]))
hr<-hclust(as.dist(cor_t),method="ward")
heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv=F,Colv=as.dendrogram(hr), scale="row",
          dendrogram="col", trace="none")

for (i in select){
  table.counts<-melt(assay(rld)[i,])
  table.counts$group<-design$condition
  p<-ggplot(table.counts, aes(x=factor(group),y=value,fill=factor(group) ) )+
    geom_boxplot(outlier.size = 0)   + geom_jitter(position=position_jitter(width=0.2),aes(factor(group), value,colour=group)) +
    scale_fill_manual("",values=c("tomato2","SlateBlue3"))+
    theme_bw(base_size = 14, base_family = "") +
    theme(axis.text.x=element_text(angle = 60, hjust = 1))+
    labs(list(title=i,y="log2(norm_counts)",x=""))
  png(paste(sep="","/home/lpantano/crickhome/parkinson/micros/seqbuster.eu/",i,".cc.vs.ct.png"))
  print(p)
  dev.off()
}

###########
head(table.n[,1:7],2)
head(table[,1:7],2)
table.n<-t(t(table)/(colSums(table)))*1000000

idx<-row.names(assay(rld))
par(mfrow=c(4,7),mar=c(1,1,2,2))
for (i in 1:28){
 plot(assay(rld)[,i],log2(table.n[idx,i]+0.1),pch=20)
}


plot(log2(mirna.g10.cf10[idx,1]),log2(mirna.g10[idx,1]),pch=20)
