#!/bin/bash

# PWH 10-19-2022
# The purpose of this script is to filter the combined BAMs for primary alignments

# Run in conda environment with winnowmap (v2.03), meryl (v1.3), and samtools (v1.9) installed


# variables
outdir=/pym/Data/paul/hg002_y/meth_calls/bams
bam20=sorted_merged_CHM13_v2.0_HG002_nanonome.bam
bam27=sorted_merged_CHM13_v2.7_HG002_nanonome.bam

# CHM13 v2.0
## keep only primary alignments 
samtools view -@48 -h -b -F 256 -F 2048 ${outdir}/${bam20} > ${outdir}/filtered_${bam20}
samtools index -@48 ${outdir}/filtered_${bam20}

# CHM13 v2.7
## keep only primary alignments 
samtools view -@48 -h -b -F 256 -F 2048 ${outdir}/${bam27} > ${outdir}/filtered_${bam27}
samtools index -@48 ${outdir}/filtered_${bam27}
