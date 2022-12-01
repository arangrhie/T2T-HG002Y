#!/bin/sh
#SBATCH --cpus-per-task=16
#SBATCH --mem=24g
#SBATCH --time=1:00:00
#SBATCH -a 1-32
#SBATCH -J map64chm13

# --cpus-per-task=16 --mem=24g --time=1:00:00 -a 1-32 -J map64 run.sh

part=$( expr $SLURM_ARRAY_TASK_ID - 1 )
part=$( printf %02d $part )

echo ""
echo "INPUTS for part $part"
echo ""

ls -l release213/bacteria/bacteria.${part}??.genomic.fna.gz

echo ""
echo "START"
echo ""

mkdir chm13v2

../meryl/build/bin/position-lookup \
  -s genomes/chm13v2.chrY.fasta \
  -m genomes/chm13v2.chrY.k64.meryl \
  -hpq chm13v2/hpq.$part \
  -mpb chm13v2/mpb.$part \
  -qpb chm13v2/qpb.$part \
  release213/bacteria/bacteria.${part}??.genomic.fna.gz
> chm13v2/err.$part 2>&1

