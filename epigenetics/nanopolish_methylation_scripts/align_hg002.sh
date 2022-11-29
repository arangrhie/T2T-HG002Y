#!/bin/bash

# PWH 10-19-2022
# The purpose of this script is to prep the reference genomes for alignment and then align the nanoNOMe data

# Run in conda environment with winnowmap (v2.03), meryl (v1.3), and samtools (v1.9) installed

# Universal variables
fastq=/uru/Data/old_atium/Data/Nanopore/projects/HG002_nanonome/fastq
files=$(ls ${fastq}/*HG002*.fastq)
out=/pym/Data/paul/hg002_y/meth_calls
out_bams=${out}/bams
out_references=${out}/references

echo $files

# Get genomes
wget -nc -P $out_references  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
wget -nc -P $out_references  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/chm13v1.1_hg002XYv2.7.fasta

# Run chm13 v2.0 preparation and alignment
## variables
prefix_20="CHM13_v2.0"
ref_20=${out_references}/chm13v2.0.fa.gz

## making the kmers necessary for alignment
meryl count threads=48 k=15 output merylDB_20 $ref_20
meryl print greater-than distinct=0.9998 merylDB_20 > ${out_references}/${prefix_20}_repetitive_k15.txt

## Aligning, sorting, indexing
for i in $files; do
	base_20=$(basename "$i" .fastq)
	echo $base_20
	winnowmap -t 48 -W ${out_references}/${prefix_20}_repetitive_k15.txt -ax map-ont $ref_20 $i | samtools view -@48 -Sb | samtools sort -@48 -o ${out_bams}/sorted_${prefix_20}_${base_20}.bam &> ~/${prefix_20}_${base_20}.log
	samtools index -@48 ${out_bams}/sorted_${prefix_20}_${base_20}.bam
done

# Run chm13 v2.7 preparation and alignment
## variables
prefix_27="CHM13_v2.7"
ref_27=${out}/references/chm13v1.1_hg002XYv2.7.fasta

## making the kmers necesarry for alignment
meryl count threads=48 k=15 output merylDB_27 $ref_27
meryl print greater-than distinct=0.9998 merylDB_27 > ${out_references}/${prefix_27}_repetitive_k15.txt

## Aligning, sorting, indexing
for i in $files; do
        base_27=$(basename "$i" .fastq)
        echo $base_27
        winnowmap -t 48 -W ${out_references}/${prefix_27}_repetitive_k15.txt -ax map-ont $ref_27 $i | samtools view -@48 -Sb | samtools sort -@48 -o ${out_bams}/sorted_${prefix_27}_${base_27}.bam &> ~/${prefix_20}_${base_20}.log
        samtools index -@48 ${out_bams}/sorted_${prefix_27}_${base_27}.bam
done

exit
