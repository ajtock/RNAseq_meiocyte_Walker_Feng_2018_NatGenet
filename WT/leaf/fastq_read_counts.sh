#!/bin/bash

# Example usage via condor submission system on hydrogen node7
# csmit -m 1G -c 3 "bash fastq_read_counts.sh WT_RNAseq_leaf_Rep1_pooled WT_RNAseq_leaf_Rep2_pooled WT_RNAseq_leaf_Rep3_SRR4204566" 

# Output files with the suffix "fastq.stats" contain read counts
for i in "$@"
do
( gunzip -k ${i}.fastq.gz
  cat ${i}.fastq | echo $((`wc -l`/4)) > ${i}.fastq.stats ) &
done
wait
