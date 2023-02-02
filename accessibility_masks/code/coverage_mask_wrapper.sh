#!/bin/bash

#SBATCH --job-name=cov_mask
#SBATCH -N 1
#SBATCH --cpus-per-task 48
#SBATCH --time=48:00:00
#SBATCH --array=1-24

#====================#
# Set up environment #
#====================#

source ~/.bashrc
source activate pysam


#===============#
# Set filepaths #
#===============#

projectDir=/scratch4/mschatz1/mschatz/T2T-chrY_mask
outDir=${projectDir}/masks/chrom_coverage

if [ ! -d $outDir ]; then
  mkdir -p $outDir
fi


#===============================#
# Generate median coverage mask #
#===============================#

cd $projectDir

if [[ ${SLURM_ARRAY_TASK_ID} -eq 23 ]]
then
  chrom="chrX"
elif [[ ${SLURM_ARRAY_TASK_ID} -eq 24 ]]
then
  chrom="chrY"
else
  chrom="chr${SLURM_ARRAY_TASK_ID}"
fi

if [[ $chrom == "chrX" ]]
then
  python code/coverage_mask.py ${chrom}_nonPAR  $projectDir > $outDir/coverage_mask_${chrom}_unsorted.txt
  python code/coverage_mask.py ${chrom}_PAR $projectDir >> $outDir/coverage_mask_${chrom}_unsorted.txt
  bedtools sort -i $outDir/coverage_mask_${chrom}_unsorted.txt > $outDir/coverage_mask_${chrom}.txt
  rm $outDir/coverage_mask_${chrom}_nonPAR.txt $outDir/coverage_mask_${chrom}_PAR.txt
elif [[ $chrom == "chrY" ]]
then
  python code/coverage_mask.py ${chrom}_nonPAR $projectDir > $outDir/coverage_mask_${chrom}.txt
else
  python code/coverage_mask.py ${chrom} $projectDir > $outDir/coverage_mask_${chrom}.txt
fi