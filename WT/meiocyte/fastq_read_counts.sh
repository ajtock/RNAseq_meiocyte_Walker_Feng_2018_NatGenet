#!/bin/bash

# Example usage via condor submission system on hydrogen node7
# csmit -m 1G -c 3 "bash fastq_read_counts.sh WT_RNAseq_meiocyte_Rep1_SRR4204534 WT_RNAseq_meiocyte_Rep2_SRR4204535 WT_RNAseq_meiocyte_Rep3_pooled" 

# Output files with the suffix "fastq.stats" contain read counts
for i in "$@"
do
( gunzip -k ${i}.fastq.gz
  cat ${i}.fastq | echo $((`wc -l`/4)) > ${i}.fastq.stats ) &
done
wait
