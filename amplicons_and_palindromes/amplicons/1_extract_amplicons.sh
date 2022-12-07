#!/bin/bash

set -e
set -x

#create a meryl database
meryl count k=15 output merylDB chrY_hg002_v2.7.fasta
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt

#use hg38 Y chromosome and the coordinates from Teitz et al. 2018 to extract amplicons in the fasta format
bedtools getfasta -fi chrY.fa -bed hg38.amplicons.colors.bed -name -s -fullHeader >hg38.colors.fasta
csplit hg38.colors.fasta /\>.*/ {*} 
for a in x*; do echo $a; mv $a $(head -1 $a).fa; done;

#default orientation because bedtools was run with -s parameter
cat hg38.BLUE.1::chrY:21925349-22093051-.fa hg38.BLUE.2::chrY:22493356-22661025+.fa hg38.BLUE.3::chrY:23535162-23702730-.fa hg38.BLUE.4::chrY:25967344-26135454+.fa >blue.hg38.default.fa
cat hg38.GRAY.1::chrY:23359751-23474022-.fa hg38.GRAY.2::chrY:26196591-26311317+.fa >gray.hg38.default.fa
cat hg38.GREEN.1::chrY:22746680-23061370-.fa hg38.GREEN.2::chrY:24380446-24695116-.fa hg38.GREEN.3::chrY:24974970-25289694+.fa >green.hg38.default.fa
cat hg38.RED.1::chrY:23061754-23208205-.fa hg38.RED.2::chrY:23210394-23359020+.fa hg38.RED.3::chrY:24695499-24822589-.fa hg38.RED.4::chrY:24824725-24974626+.fa >red.hg38.default.fa
cat hg38.TEAL.1::chrY:22093051-22208739-.fa hg38.TEAL.2::chrY:22377680-22493356+.fa >teal.hg38.default.fa
cat hg38.YELLOW.1::chrY:23702730-24380446-.fa hg38.YELLOW.2::chrY:25289694-25967344+.fa >yellow.hg38.default.fa