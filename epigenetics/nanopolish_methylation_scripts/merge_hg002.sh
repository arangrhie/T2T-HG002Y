#!/bin/bash

# PWH 10-19-2022
# Purpose of this script is to merge aligned nanoNOMe data

# Run in conda environment with winnowmap (v2.03), meryl (v1.3), and samtools (v1.9) installed

# Universal variables
out=/pym/Data/paul/hg002_y/meth_calls
out_bams=${out}/bams

# chm13 v2.0 merge
## variables
prefix_20="CHM13_v2.0"

## merge
echo `ls ${out_bams}/*${prefix_20}*.bam`
samtools merge -@48 ${out_bams}/sorted_merged_${prefix_20}_HG002_nanonome.bam ${out_bams}/*${prefix_20}*.bam
samtools index -@48 ${out_bams}/sorted_merged_${prefix_20}_HG002_nanonome.bam


# chm13 v2.7 merge
## variables
prefix_27="CHM13_v2.7"

## merge
echo `ls ${out_bams}/*${prefix_27}*.bam`
samtools merge -@48 ${out_bams}/sorted_merged_${prefix_27}_HG002_nanonome.bam ${out_bams}/*${prefix_27}*.bam
samtools index -@48 ${out_bams}/sorted_merged_${prefix_27}_HG002_nanonome.bam

exit
