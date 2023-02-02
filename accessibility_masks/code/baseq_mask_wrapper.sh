#!/bin/bash

#SBATCH --job-name=baseq_mask
#SBATCH -N 1
#SBATCH --cpus-per-task 48
#SBATCH --time=48:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=dtaylo95@jhu.edu
#SBATCH --array=1-24%24

#====================#
# Set up environment #
#====================#

source ~/.bashrc
source activate pysam


#===============#
# Set filepaths #
#===============#

projectDir=/scratch4/mschatz1/mschatz/T2T-chrY_mask
outDir=${projectDir}/masks/baseq


#=================#
# Get chrom masks #
#=================#

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

python code/baseq_mask.py ${chrom} ${projectDir} > ${outDir}/baseq_mask_${chrom}.txt