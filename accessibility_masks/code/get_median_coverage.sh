#!/bin/bash

#SBATCH --job-name=median_cov
#SBATCH -N 1
#SBATCH --cpus-per-task 48
#SBATCH --time=2:00:00
#SBATCH --array=1-29

#====================#
# Set up environment #
#====================#

source ~/.bashrc
source activate acc_mask


#===============#
# Set filepaths #
#===============#

projectDir=/scratch4/mschatz1/mschatz/T2T-chrY_mask

assembly_fasta=$projectDir/data/references/chm13v2.0.XY.fasta
autosomes_bed=$projectDir/data/references/autosomes.bed
sample_list=$projectDir/data/sample.info

SAMPLE=`sed "${SLURM_ARRAY_TASK_ID}q;d" $sample_list | awk '{print $1}'`
SAMPLE_CRAM=$projectDir/data/T2T-CHM13v2_crams/${SAMPLE}.cram

outDir=${projectDir}/masks/median_coverage

if [ ! -d $outDir ]; then
	mkdir -p $outDir
fi


#===============================#
# Generate median coverage mask #
#===============================#

cd $projectDir

SAMPLE=`sed "${SLURM_ARRAY_TASK_ID}q;d" data/sample_list.txt | awk '{print $1}'`

mosdepth \
	-n \
	-t 48 \
	--use-median \
	-f $assembly_fasta \
	-b $autosomes_bed \
	$SAMPLE \
	$SAMPLE_CRAM

cat ${SAMPLE}.mosdepth.region.dist.txt | awk '{if ($1 == "total" && $3 > 0.5) print}' | head -1 | awk -v "sample=${SAMPLE}" '{print sample"\t"$2}' > ${outDir}/${SAMPLE}.median

rm ${SAMPLE}*.txt ${SAMPLE}*.csi ${SAMPLE}*.gz