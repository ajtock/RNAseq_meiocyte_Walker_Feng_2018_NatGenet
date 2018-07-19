#!/bin/bash

## ************** TEST BEFORE USE ON IMPORTANT FASTQ FILES ******************
# Example usage via condor submission system on hydrogen node7:
# csmit -m 10G -c 1 "bash fastq_gzip_best.sh SRR4204566 WT_RNASeq_leaf_Rep3_SRR4204566"

run1=$1
name=$2

if [ ! -f "$name.fastq.gz" ]; then 
    gzip -c -k --best $run1".fastq" > $name.fastq.gz;
else 
    echo "skipping $name"
fi

