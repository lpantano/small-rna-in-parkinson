#!/bin/bash

set -e

MIRA="/home/lpantano/soft/srnabench/sRNAbench.jar"

OPATH=miranalyzer_output_trimmed
CONFIG=$1
DB=~/soft/srnabench-db/sRNAbenchDB
BPATH=/home/lpantano/soft/bowtie-1.0.1
TLIBS=hsa_tRNA

mkdir -p $OPATH
cd $OPATH

while read line; do

 NAME=`echo $line| cut -d " " -f 2`   
 INPUTFASTA=`echo $line| cut -d " " -f 1`   
 
 echo "$NAME $INPUTFASTA"

 INPUT=$NAME.rc

 mkdir -p $NAME
 cd $NAME
 #create rc count
 INPUT=input.rc 
 sed 's/_x/\t/' ../../fasta/trimmed/$INPUTFASTA.trimmed.fa | awk '{if ($0~/>/){counts=$2}else{print $0"\t"counts}}' > $INPUT

 #run MIRA
 java -jar $MIRA dbPath=$DB microRNA=hsa input=$INPUT  output=out maxReadLength=40 isoMiR=true libs=hg19-tRNAs.fa,hsa_Rfam.fa

 cd ..

done < $CONFIG
