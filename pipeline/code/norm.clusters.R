library(data.table)
library(DESeq2)
library(qvalue)
library(ggplot2)
library("RColorBrewer")
library("gplots")
library(gtools)


table<-read.table("../../fasta/prepare_18_35_10/res/counts.tsv",sep="\t",header=T,row.names=1)
ann<-table[,1]
names(ann)<-row.names(table)
table<-table[,2:(ncol(table))]
table <- table[,mixedsort(names(table))]

table<-read.table("../../fasta/define.cluster.Dec/res/ann.tab",sep="\t",header=T,row.names=1)
ann<-table[,1]
names(ann)<-row.names(table)
table<-table[,2:(ncol(table)-1)]
table <- table[,mixedsort(names(table))]


#########################################################
get_norm_values <- function(table,con,treat,ini,end){
    design<-data.frame(condition=c(rep(con,7),rep(treat,7)))
    row.names(design)<-names(table)[ini:end]
    
    mirna<-table[,ini:end]
    mirna<-mirna[rowSums(mirna>10)>=5,]
    
    dds <- DESeqDataSetFromMatrix(countData = mirna,
                                  colData = design,
                                  design = ~ condition)
    dds <- estimateSizeFactors( dds )
    #summary(q)
    rlogTransformation(dds, blind=FALSE)
}


rld1 <- get_norm_values(table,"cc","ct",1,14)
trna_ids <- names(ann)[grepl("tRNA",ann)]
idx <- intersect(row.names(assay(rld1)),trna_ids)
trna_counts_clininc <- assay(rld1)[idx,]


rld1 <- get_norm_values(table,"pc","pt",15,28)
trna_ids <- names(ann)[grepl("tRNA",ann)]
idx <- intersect(row.names(assay(rld1)),trna_ids)
trna_counts_preclininc <- assay(rld1)[idx,]


library(devtools)
devtools::load_all("~/repos/isomiRs")
counts<-setClass("DataSeq",
                 slots=c(counts="matrix",
                         normcounts="matrix",
                         design="data.frame"
                 ))


do_pls<-function(counts)
{
    obj<-counts()
    obj@normcounts<-as.matrix(counts)
    obj@design<-data.frame(g=gsub("[0-9]+","",colnames(obj@normcounts)), b=1)
    pls <- isoPLSDA(obj,"g",nperm = 400)    
    pls
}

cc_trna <- do_pls(trna_counts_clininc)
pc_tnra <- do_pls(trna_counts_preclininc)

