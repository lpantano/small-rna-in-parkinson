DB=../../data
seqcluster cluster -a seqs.sort.bam -m seqs.ma -o res -b $DB/tRNA.bed,$DB/miRNA.bed,$DB/pirna2.bed,$DB/rmsk.bed,$DB/snoRNA.bed,$DB/wgRNA.bed,$DB/refGene.bed,$DB/simpleR.bed
