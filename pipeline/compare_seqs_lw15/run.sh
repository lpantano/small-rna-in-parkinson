#!/bin/bash

#analyze sequences lower than 15 nucleotides

set -e

FOLDER=$1
RUN=$2
INDEX=../miranalyzer/Human/miRanalyzerDB/bowtie/translibs/hsa_tRNA

cd $FOLDER

rm -f  config

for dir in `dir` ; do
 if [ "$dir" != "prepare" ] ; then
  echo "$dir"
  if [ "$RUN" == "get15" ] ; then
   FILE=`ls $dir/*.ad`
   awk '{l=length($1);if (l <18 && l > 14){i=i+1;print ">seq_"i"_x"$2"\n"$1}}' $FILE > $dir/$dir.lw15.fa
  fi

  if [ "$RUN" == "prepare" ] ; then
   printf "$dir/$dir.lw15.fa\t$dir\n" >>config
  fi
 fi
done

if [ "$RUN" == "prepare" ] ; then
 echo "seqcluster prepare"
 PATH=~/soft/bcbiome/anaconda/bin:$PATH
 mkdir -p prepare
 seqcluster prepare -c config -o prepare
fi

if [ "$RUN" == "bowtie" ] ; then
 echo "bowtie"
 bowtie -f -a --best --strata -S -m 500 $INDEX prepare/seqs.fa > prepare/seqs.sam
fi

if [ "$RUN" == "parse" ] ; then
 echo "parse"
 python ../pipeline/compare_seqs_lw15/counts_mapped.py prepare/seqs.sam prepare/seqs.ma prepare/stats.tsv
fi

cd ..

