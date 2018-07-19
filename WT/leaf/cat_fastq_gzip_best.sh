#!/bin/bash

## ************** TEST BEFORE USE ON IMPORTANT FASTQ FILES ******************
# Example usage via condor submission system on hydrogen node7:
# csmit -m 10G -c 1 "bash cat_fastq_gzip_best.sh SRR4204558 SRR4204559 SRR4204560 SRR4204561 WT_RNASeq_leaf_Rep1_pooled"

run1=$1
run2=$2
run3=$3
run4=$4
name=$5

if [ ! -f "$name.fastq.gz" ]; then 
    cat $run1".fastq" $run2".fastq" $run3".fastq" $run4".fastq" \
    | gzip -c -k --best > $name.fastq.gz;
else 
    echo "skipping $name"
fi

