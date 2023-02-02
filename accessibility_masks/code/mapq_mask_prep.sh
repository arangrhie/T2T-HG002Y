#!/bin/bash

#SBATCH --job-name=mapq_prep
#SBATCH -N 1
#SBATCH --cpus-per-task 48
#SBATCH --time=48:00:00
#SBATCH --array=1-24

#====================#
# Set up environment #
#====================#

source ~/.bashrc
source activate acc_mask


#===============#
# Set filepaths #
#===============#

projectDir=/scratch4/mschatz1/mschatz/T2T-chrY_mask

cramDir=$projectDir/data/T2T-CHM13v2_crams
fullRef=$projectDir/data/references/chm13v2.0.fasta
XXRef=$projectDir/data/references/chm13v2.0.XX.fasta
XYRef=$projectDir/data/references/chm13v2.0.XY.fasta
sampleInfo=$projectDir/data/samples.info


#=================#
# Get chrom masks #
#=================#

cd $cramDir

if [[ ${SLURM_ARRAY_TASK_ID} -eq 23 ]]
then
  chrom="chrX"
elif [[ ${SLURM_ARRAY_TASK_ID} -eq 24 ]]
then
  chrom="chrY"
else
  chrom="chr${SLURM_ARRAY_TASK_ID}"
fi

# If not chrY...
if [ $chrom != "chrY" ]; then
  samtools merge \
    -R $chrom \
    -h HG00448.cram \
    -O CRAM \
    --reference $fullRef \
    -@ 4 \
    ${chrom}_merged.cram \
    HG*.cram

  samtools index -@ 16 ${chrom}_merged.cram
  exit 0
fi


# If chrY...

# Extract chrY alignments from each sample
for sample in `cat $sampleInfo | tr -s ' ' '\t' | cut -f 1`
do
  sampleSex=`grep -m1 $sample $sampleInfo | tr -s ' ' '\t' | cut -f 2`
  if [[ $sampleSex == "M" ]]
  then
    ref=$XYRef
  elif [[ $sampleSex == "F" ]]
  then
    ref=$XXRef
  fi
  samtools view -@ 8 -O BAM -T $ref -o ${sample}.chrY.bam ${sample}.cram chrY
  samtools index -@ 8 ${sample}.chrY.bam
done

# Merge Y chromosome
samtools merge \
  -R chrY \
  -h HG00448.cram \
  -O CRAM \
  --reference $fullRef \
  -@ 4 \
  chrY_merged.cram \
  HG*.chrY.bam

samtools index -@ 16 chrY_merged.cram

rm HG*.chrY.bam
