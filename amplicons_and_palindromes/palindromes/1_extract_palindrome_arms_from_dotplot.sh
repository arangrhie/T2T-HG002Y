module add conda-modules-py37
module add bedtools-2.26.0
module load emboss-6.5.7

#Start with manually extracted coordinates of arms from gepard dot-lot in the file P4-P5.hg002_v2.7.initial.bed

bedtools slop -i P4-P5.hg002_v2.7.bed -g chrY.2.7.genome.txt -b 10 >P4-P5.hg002_v2.7.extended.bed
bedtools getfasta -fi chrY_hg002_v2.7.fasta -bed P4-P5.hg002_v2.7.extended.bed -name >P4-P5.hg002_v2.7.fasta

