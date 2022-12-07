module add conda-modules-py37
module add bedtools-2.26.0
module load emboss-6.5.7

bedtools getfasta -fi chrY_hg002_v2.7.fasta -bed P4-P5.hg002_v2.7.fixed.bed -name >P4-P5.hg002_v2.7.fasta
csplit P4-P5.hg002_v2.7.fasta /\>P.*/ {*} 
for a in x*; do echo $a; mv $a $(head -1 $a).txt; done; 