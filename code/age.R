######################
#######normalize with deseq2######
library(data.table)

do_analysis<-function(rld,de,top,name){
  
  index.de<-intersect(row.names((de[order(de$pvalue),])[1:top,]),row.names(assay(rld)))
  ma.de<-assay(rld[index.de,])
  cor_t <- 1 - cor(ma.de)
  hc<-hclust(as.dist(cor_t),method="ward")
  cor_tr <- 1 - cor(t(ma.de))
  hr<-hclust(as.dist(cor_tr),method="ward")
  
  #pdf(paste("../figures/",name,sep="",
  #        ".aging.isomir.pdf"))
  par(mar=c(5,5,5,5))
  heatmap.2(ma.de,
            Colv=as.dendrogram(hc),
            Rowv=as.dendrogram(hr),
            scale="row",labRow="",
            dendrogram="col", trace="none",cexRow=.7,
            ColSideColors=cols,
            keysize=1.5)
  #dev.off()
  
  realvalue<-get_accuracy_age( hc)
  
  acc<-vector()
  hcl[["real"]]<-hc
  for (perm in 1:400){
    select<-sample(setdiff(row.names(ma),index.de),top)
    ma.de<-assay(rld[select,])
    cor_t <- 1 - cor(ma.de)
    hcl[[perm]]<-hclust(as.dist(cor_t),method="ward")
    
  }
  
  acc<-sapply(hcl,get_accuracy_age)
  #dis<-cl_dissimilarity(cl_ensemble(list=hcl))
  #plot(hclust(dis,method="ward"))
  
  pval=sum(acc[names(acc)!="real"]>=acc["real"])/401
  d<-density(acc,bw=0.1)
  #pdf(paste("../figures/age.",name,sep="",
  #          ".signif.random.isomir.pdf"))
  plot(d, main="density random score",xlab="clustering scores")
  polygon(d, col="steelblue", border="blue")
  abline(v=realvalue,col="red",lwd=2)
  text(1,1,labels=paste("p-value:",round(pval,digits=2)),cex=0.7)
  #dev.off()
  
}

get_accuracy_age<-function(m){
  lab<-cutreeDynamicTree(m,minModuleSize=4)
  c1<-length(unique(lab[1:3]))
  c2<-length(unique(lab[4:7]))
  c3<-length(unique(lab[8:12]))
  v<-(3/(c1+c2+c3))
  return(v)
}

age<-c("y0d2","y0d4","y0d34","y0d204","y8d2","y13d360","y25d152","y53d112","y66d0","y80d0","y88d0","y98d0")
table<-data.frame(mir="hsa")
for (a in age[1:length(age)]){
  d<-read.table(paste("~/crickhome/isomirs/mirbase20/",a,".hsa.fa.ad.new.mirna",sep=""),sep="\t",skip=1)

  d.f<-d[d$V14==1,c(1,3,4)]
  #d.f<-d.f[d.f$V3>=5,]  
  table.dt<-data.table(d)
  tabmir<-as.data.frame(table.dt[,list(freq=sum(V3)),by=V4])
  d.f.m<-merge(d.f[,c(3,1,2)],tabmir,by=1)
  d.f.m$r<-d.f.m$V3/d.f.m$freq
  d.f.m<-d.f.m[d.f.m$r>=0.1,]
  
  tabmir.f<-as.data.frame(data.table(d.f.m)[,list(freq=sum(V3)),by=V4])
  
  names(tabmir.f)<-c("mirna",a)
  table<-merge(table,tabmir.f,by=1,all=TRUE)
}

table.a<-table[2:nrow(table),]
row.names(table.a)<-table.a[,1]
table.a<-table.a[,2:ncol(table.a)]
table.a[is.na(table.a)]<-0

ma<-table.a[apply(table.a<5,1,sum)<=6,]


design2<-data.frame(condition=c(rep(1,length(age)/2),rep(2,length(age)/2)))
row.names(design2)<-age

dds <- DESeqDataSetFromMatrix(countData = data.frame(ma),
                              colData = design2,
                              design = ~ condition)

#dds <- estimateSizeFactors( dds )
#sizeFactors(dds)
#ma.n<-t(t(ma)*sizeFactors(dds))
rld <- rlogTransformation(dds,blind=TRUE)
#vsd <- varianceStabilizingTransformation(dds,blind=TRUE)

cols<-c(rep("slategray3",3),rep("SlateGray4",4),rep("Grey",5))
colfunc <-colorRampPalette(c("green", "blue"))
cols<-colfunc(12)
#list in de
de<-read.table("tables.res/DEseq.miraligner.mirna.ptvspc.txt",header=T)

index.de<-intersect(row.names((de[order(de$pvalue),])[1:10,]),row.names(ma))

ma.de<-assay(rld[index.de,])
cor_t <- 1 - cor(ma.de)
hc<-hclust(as.dist(cor_t),method="ward")
cor_tr <- 1 - cor(t(ma.de))
hr<-hclust(as.dist(cor_tr),method="ward")
pdf(paste("heatmaps.res/ptvspc",sep="",
          ".aging.mirna.pdf"))
heatmap.2(ma.de, 
          Colv=as.dendrogram(hc),
          Rolv=as.dendrogram(hr),
          scale="row",
          dendrogram="col", trace="none",cexRow=.7,
          ColSideColors=cols,margins=c(10,10),
          keysize=1.5)
dev.off()

###########isomirs#################
###########isomirs#################
###########isomirs#################
table<-data.frame(mir="hsa")
for (a in age[1:length(age)]){
  d<-read.table(paste("~/crickhome/isomirs/mirbase20/",a,".hsa.fa.ad.new.mirna",sep=""),sep="\t",header=T)
  
  d.f<-d[d$ambiguity==1,]
  d.f<-d.f[d.f$freq>=3,]  
  d.f$id<-paste(d.f$mir,d.f$mism,d.f$add,d.f$t5,d.f$t3,sep=".")
  tabmir<-d.f[,c("seq","freq")]
  names(tabmir)<-c("isomir",a)
  table<-merge(table,tabmir,by=1,all=TRUE)
}

table<-table[2:nrow(table),]
row.names(table)<-table[,1]
table.a<-table[,2:ncol(table)]
table.a[is.na(table.a)]<-0

ma<-table.a[apply(table.a<10,1,sum)<=3,]

design2<-data.frame(condition=age)
row.names(design2)<-age

dds <- DESeqDataSetFromMatrix(countData = data.frame(ma),
                              colData = design2,
                              design = ~ condition)
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds,fitType="local")
plotDispEsts(dds)
rld <- rlogTransformation(dds,blind=FALSE)

cols<-c(rep("slategray3",3),rep("SlateGray4",4),rep("Grey",5))
colfunc <-colorRampPalette(c("green", "blue"))
cols<-colfunc(12)
#list in de
de<-read.table("../micros/mirlaigner.sept10.2013/tables.res/DEseq.miraligner.isomir.ptvspc.txt",header=T,sep="\t")
top<-85
do_analysis(rld,de,top,"ptvspc")
de<-read.table("../micros/mirlaigner.sept10.2013/tables.res/DEseq.miraligner.isomir.ctvscc.txt",header=T,sep="\t")
top<-233
do_analysis(rld,de,top,"ctvscc")
