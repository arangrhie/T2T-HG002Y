#!/bin/bash 

# PWH 10-19-2022
# The purpose of this script is to call methylation for a specific fastq file (given on the command line) for the CHM13 v2.7 genome. See "nanonome_cpggpc_call_meth.sh" for implementation of this script

# command line input
fastq=$1

# variables
basename=cpggpc_nanonome_CHM13v2.7_HG002_nanonome
ref27=/pym/Data/paul/hg002_y/meth_calls/references/chm13v1.1_hg002XYv2.7.fasta
bam27=/pym/Data/paul/hg002_y/meth_calls/bams/filtered_sorted_merged_CHM13_v2.7_HG002_nanonome.bam
fastq_dir=/uru/Data/old_atium/Data/Nanopore/projects/HG002_nanonome/fastq
repo=/home/phook2/software/nanopolish_nanonome
outdir=/pym/Data/paul/hg002_y/meth_calls/calls/nanonome

# print some things to know we are good
echo $fastq
echo $ref27
echo $bam27
${repo}/nanopolish --version

# calling meth
${repo}/nanopolish call-methylation --progress -b $bam27 -r ${fastq_dir}/$fastq -g $ref27 -q cpggpc -t 48 > ${outdir}/${basename}_methylation_${fastq%.*}.tsv

exit
