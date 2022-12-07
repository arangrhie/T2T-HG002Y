module add conda-modules-py37
module add bedtools-2.26.0
conda activate winnowmap

#prepare meryl database
meryl count k=15 output merylDB chrY_hg002_v2.7.fasta
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt

#create overlapping windows of the size 5k in the bed format
bedtools makewindows -g chrY.2.7.genome.txt -w 5000 -s 1000 -i srcwinnum >chrY.2.7.windows.5000.bed

#generate windows of the size 5k using a hardmasked version of the Y chromosome
bedtools getfasta -fi chrY.fasta -bed chrY.2.7.windows.5000.bed >chrY.2.7.hardmasked.windows.5000.fasta
winnowmap --MD -W repetitive_k15.txt -ax map-ont chrY_hg002_v2.7.fasta chrY.2.7.hardmasked.windows.5000.fasta > chrY.2.7.hardmasked.windows.sam

samtools view -bh -F 4 chrY.2.7.hardmasked.windows.sam >chrY.2.7.hardmasked.windows.bam
samtools sort chrY.2.7.hardmasked.windows.bam -o chrY.2.7.hardmasked.windows.sorted.bam
samtools index chrY.2.7.hardmasked.windows.sorted.bam
samtools view  chrY.2.7.hardmasked.windows.sorted.bam | grep -v "^@" | cut -f1-6 >chrY.2.7.hardmasked.windows.sorted.bam.f1-6.txt

#now calculate the identities between the windows and the reference
module load python36-modules-gcc
python parse_bam.py chrY.2.7.hardmasked.windows.sorted.bam



