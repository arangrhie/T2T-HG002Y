#!/bin/bash

set -e
set -x

module add conda-modules-py37
conda activate winnowmap
module load bedtools-2.26.0 seqtk-1.0

ampliconName=( "blue" "gray" "green" "red" "teal" "yellow" )

#map each amplicon type to the T2T Y chromosome (chrY_hg002_v2.7.fasta)
for color in "${ampliconName[@]}"
do
	winnowmap -W repetitive_k15.txt --sv-off -ax map-pb chrY_hg002_v2.7.fasta ${color}.hg38.default.fa | samtools view -h -F 2048 >${color}.default.sam

	samtools view -bhS ${color}.default.sam >${color}.default.bam
	samtools sort ${color}.default.bam >${color}.default.sorted.bam
	samtools index ${color}.default.sorted.bam
	bedtools bamtobed -i ${color}.default.sorted.bam >${color}.coordinates.HG002.txt #these are the coordinates to be used/potentially manually curated
done