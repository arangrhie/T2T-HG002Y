#!/bin/bash

# PWH 10-19-2022
# The purpose of this script is to process nanopolish output into a form that can be analyzed

# variables
indir=/pym/Data/paul/hg002_y/meth_calls/calls/nanonome_2.0
basename=cpggpc_nanonome_CHM13v2.0_HG002
outdir=/pym/Data/paul/hg002_y/meth_calls/calls/nanonome_2.0/processed

bam=/pym/Data/paul/hg002_y/meth_calls/bams/filtered_sorted_merged_CHM13_v2.0_HG002_nanonome.bam

repo=/home/phook2/software/nanopore-methylation-utilities

ref20=/pym/Data/paul/hg002_y/meth_calls/references/chm13v2.0.fa

nanopolish=/home/phook2/software/nanopolish_nanonome

# Processing all data
## Combining nanopolish meth probability tsvs from all the fastqs
awk 'FNR==1 && NR!=1 {next} {print}' ${indir}/cpggpc_nanonome_CHM13v2.0_HG002_nanonome_methylation*.tsv > ${outdir}/${basename}_methylation.tsv
echo "Meth combined: `date`"

## Creating a list of reads with primary alignment and length >20kb. This comes from grabbing the read ids from the filterd merged BAMs
samtools view -@48 $bam | awk 'length($10) > 20000' | cut -f1 > ${outdir}/chm13v2.0_20kb_readIDs.txt
echo "20kb reads determined: `date`"

## Grabbing the header from the combined tsv files and removing the columns that display the 4-state probabilities. These columns do not play well with downstream processing and are not needed. There are better ways to do this, but I grabbed the header because in the processing steps below, the header will be lost.
head -n 1 ${outdir}/${basename}_methylation.tsv | awk -v FS='\t' -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}'> ${outdir}/cpg_nanonome_CHM13v2.0_header.txt
echo "got header: `date`"

## Filter for methylation probabilities from reads with primary alignments and >20kb in length
awk 'NR==FNR {a[$1]++; next}a[$5]' ${outdir}/chm13v2.0_20kb_readIDs.txt ${outdir}/${basename}_methylation.tsv > ${outdir}/20kb_${basename}_methylation.tsv
echo "Reads filtered: `date`"

# Processing GpC data
## Filter for GpC. I not only filter for just CG calls here, I also only keep the relevant columns, removing the 4-state probability columns which are not necessary for downstream analysis
cat ${outdir}/20kb_${basename}_methylation.tsv | awk -v FS='\t' -v OFS='\t' '$11=="GC" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}'> ${outdir}/gpc_filtered_header_20kb_${basename}_methylation.tsv
echo "filtered for gpc: `date`"

## Adding header back in
cat ${outdir}/cpg_nanonome_CHM13v2.0_header.txt ${outdir}/gpc_filtered_header_20kb_${basename}_methylation.tsv > ${outdir}/header_gpc_filtered_header_20kb_${basename}_methylation.tsv

## Process to tmp file. This converts tsvs to bedgraphs and sets a LLR threshold for methylation calling
python3 ${repo}/mtsv2bedGraph.py -q gpc -c 1.0 --nome -i ${outdir}/header_gpc_filtered_header_20kb_${basename}_methylation.tsv -g $ref20 > ${outdir}/gpc_hg002_chm13v2.0_meth.tmp
echo "gpc processed tmp file: `date`"

## Sort and bgzip the temp file
sort ${outdir}/gpc_hg002_chm13v2.0_meth.tmp -k1,1 -k2,2n | bgzip > ${outdir}/gpc_hg002_chm13v2.0_meth.bed.gz
echo "cpg sorted and bgzipped: `date`"

## Index the bedgraph
tabix -p bed ${outdir}/gpc_hg002_chm13v2.0_meth.bed.gz
echo "cpg indexed: `date`"

## Convert to an aggregated frequency file. This file can go forward and be used in downstream analyses.
python3 ${repo}/parseMethylbed.py frequency -i ${outdir}/gpc_hg002_chm13v2.0_meth.bed.gz -m gpc > ${outdir}/gpc_hg002_chm13v2.0_meth.freq
echo "cpg frequency conversion finished: `date`"

# Processing CpG data
## Filter for CpG. I not only filter for just CG calls here, I also only keep the relevant columns, removing the 4-state probability columns which are not necessary for downstream analysis
cat ${outdir}/20kb_${basename}_methylation.tsv | awk -v FS='\t' -v OFS='\t' '$11=="CG" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}'> ${outdir}/cpg_filtered_header_20kb_${basename}_methylation.tsv
echo "filtered for cpg: `date`"

## Adding header back in
cat ${outdir}/cpg_nanonome_CHM13v2.0_header.txt ${outdir}/cpg_filtered_header_20kb_${basename}_methylation.tsv > ${outdir}/header_cpg_filtered_header_20kb_${basename}_methylation.tsv

## Process to tmp file. This converts tsvs to bedgraphs and sets a LLR threshold for methylation calling
python3 ${repo}/mtsv2bedGraph.py -q cpg -c 1.5 --nome -i ${outdir}/header_cpg_filtered_header_20kb_${basename}_methylation.tsv -g $ref20 > ${outdir}/cpg_hg002_chm13v2.0_meth.tmp
echo "cpg processed tmp file: `date`"

## Sort and bgzip the temp file
sort ${outdir}/cpg_hg002_chm13v2.0_meth.tmp -k1,1 -k2,2n | bgzip > ${outdir}/cpg_hg002_chm13v2.0_meth.bed.gz
echo "cpg sorted and bgzipped: `date`"

## Index the bedgraph
tabix -p bed ${outdir}/cpg_hg002_chm13v2.0_meth.bed.gz
echo "cpg indexed: `date`"

## Convert to an aggregated frequency file. This file can go forward and be used in downstream analyses.
python3 ${repo}/parseMethylbed.py frequency -i ${outdir}/cpg_hg002_chm13v2.0_meth.bed.gz > ${outdir}/cpg_hg002_chm13v2.0_meth.freq
echo "cpg frequency conversion finished: `date`"

exit
