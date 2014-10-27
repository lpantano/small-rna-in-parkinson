#!/bin/bash

set -e


OPATH=flaimapper_out
CONFIG=$1
BPATH=/home/lpantano/soft/bowtie-1.0.1
GENOME="/home/lpantano/projects/index/ody_copy/hg19"
MYPATH=$(which ./myscript)

mkdir -p $OPATH
cd $OPATH
LIT_F=""


while read line; do    

 NAME=`echo $line| cut -d " " -f 2`   
 INPUTFASTA=`echo $line| cut -d " " -f 1`   
 
 echo "$NAME $INPUTFASTA"

 INPUT=$NAME.uncollapse.fa

 mkdir -p $NAME
 cd $NAME
 if [ ! -e $NAME.bam ] ;  then
  #create rc count
  python $MYPATH/uncollapse.py ../../raw/$INPUTFASTA > $INPUT
  #run bowtie
  $BPATH/bowtie -a --best --strata -f -S -m 5000  $GENOME $INPUT| samtools view -Shb /dev/stdin | samtools sort -o /dev/stdin temp > $NAME.bam
 fi
 cd ..

 LIST_F=`echo $LIST_F $NAME/$NAME.bam`

done < $CONFIG


mkdir -p output
flaimapper -f 1 -o output/counts.tab -m ../../data/ncrnadb09_hg19.gtf -r $GENOME.fa  $LIST_F

