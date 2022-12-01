#!/bin/sh
#SBATCH --cpus-per-task=16
#SBATCH --mem=24g
#SBATCH --time=1:00:00
#SBATCH -a 1-32
#SBATCH -J map64grch38

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

mkdir grch38

../meryl/build/bin/position-lookup \
  -s genomes/grch38.chrY.fasta \
  -m genomes/grch38.chrY.k64.meryl \
  -hpq grch38/hpq.$part \
  -mpb grch38/mpb.$part \
  -qpb grch38/qpb.$part \
  release213/bacteria/bacteria.${part}??.genomic.fna.gz
> grch38/err.$part 2>&1

