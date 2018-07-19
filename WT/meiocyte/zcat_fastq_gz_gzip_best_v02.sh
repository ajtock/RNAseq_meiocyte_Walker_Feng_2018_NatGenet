#!/bin/bash

## ************** TEST BEFORE USE ON IMPORTANT FASTQ FILES ******************
# Example usage via condor submission system on hydrogen node7:
# csmit -m 10G -c 1 "bash zcat_fastq_gz_gzip_best_v02.sh SRR4204536 SRR4204538 SRR4204539 SRR4204540 WT_RNAseq_meiocyte_Rep3_pooled"

run1=$1
run2=$2
run3=$3
run4=$4
name=$5

if [ ! -f "$name.fastq.gz" ]; then 
    zcat $run1".fastq.gz" $run2".fastq.gz" $run3".fastq.gz" $run4".fastq.gz" \
    | gzip -c -k --best > $name.fastq.gz;
else 
    echo "skipping $name"
fi

