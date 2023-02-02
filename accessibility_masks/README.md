# Generating the short-read SNV accessibility mask

The approach below generates three genome masks, according to the three criteria described in the following article: https://www.illumina.com/science/genomics-research/articles/identifying-genomic-regions-with-high-quality-single-nucleotide-.html

These masks define regions of the genome where short-read SNV discovery and genotyping is relatively reliable. The three masks are then intersected to create one combined mask. Note that the results uploaded here are for the `T2T-CHM13v2.0` assembly, though the code could be generalized to any assembly version.

## 0. Getting set up
If you wish to run this pipeline on your own system, there are some steps you'll need to take to make sure everything is set up properly:
### Loading the required packages
There are four packages you'll need in order to run these scripts: `samtools`, `mosdepth`, `bedtools`, and `pysam`. Luckily, these can all be easily installed with conda. You can load them yourself, or alternatively, you can use the YAML file we've provided in this repo:
```
conda env create -f acc_mask.yml
```

### SLURM
Many of the scripts in this repo use SLURM array jobs for parallelization. The SLURM parameters are in the scripts themselves, and can be modified based on your system's capabilities.

### The file tree
The scripts in this repo assume files are organized in a specific manner. You can simply copy the structure of this repo to your local system. Scripts are located in the `code` directory, and any necessary data files are located in the `data` directory. The `masks` directory will be created and filled as you run the pipeline.

We've provided as many small data files as we can within this repo, but you will need to download the larger data files (i.e. the reference genome and sample crams) yourself. The file tree below describes where to put them.
```
projectDir/
┣━ code/
┗━ data/
   ┣━ references/
   ┃  ┣━ chm13v2.0.XX.fasta
   ┃  ┣━ chm13v2.0.XX.fasta.fai
   ┃  ┣━ chm13v2.0.XY.fasta
   ┃  ┗━ chm13v2.0.XY.fasta.fai
   ┗━ T2T-CHM13v2_crams/
      ┣━ HG*.cram
      ┗━ HG*.cram.crai
```
Finally, within the bash scripts, you'll need to change the `projectDir` variable to reflect the location of the working directory on your own system. The other file paths are defined relative to this path, so you shouldn't have to change them.

That said, if you are running this pipeline on a different reference with different files, OR if you rename any of the files, you will have to change the file names with the scripts as appropriate.

## 1. Normalized depth mask

### Compute median autosomal coverage
```
sbatch code/get_median_coverage.sh
```
When complete:
```
cat masks/median_coverage/*.median > masks/median_coverage/median_cov.txt
```
Technically, the above value may differ slightly from the median, as `mosdepth` does not report the median itself, but the cumulative distribution of coverages. We extract the first reported coverage for which the reported cumulative proportion of total bases exceeds 0.5. 

### Compare observed coverage to median autosomal coverage

We next generate masks for each chromosome, based on median autosomal coverage across samples.

For autosomes, we define "high quality" bases as those where the mean normalized coverage across samples falls within 25% of the median autosomal coverage.

For sex chromosomes, we adjust the median expectation based on the dosage of the chromosome for that sample:
* In XX samples, the expected coverage on the X chromosome is the median autosomal coverage
* In XY samples, the median autosomal coverage is divided by 2 for the non-PAR portions of both chrX and chrY
* Because the Y-PAR is masked in XY samples during alignment, the X-PAR is expected to have similar coverage as the autosomes, and we use the median autosomal coverage for the X-PAR in these samples
* Coverage for female samples is not considered when producing the coverage mask for the Y chromosome.

```
sbatch coverage_mask_wrapper.sh
```

When complete:
```
# concatenate and bedtools merge individual chromosome masks
for i in {1..22} X Y;
do
    cat masks/chrom_coverage/coverage_mask_chr${i}.txt | bedtools merge >> masks/coverage_mask.bed;
done

bgzip masks/coverage_mask.bed
```

## 2. Mapping quality mask
For both this mask and the base quality mask, we need our alignments partitioned by chromosome. For each chromosome, we extract the corresponding alignments from each sample CRAM, and merge these into a chromosome CRAM.
```
sbatch mapq_mask_prep.sh
```
When this is complete, we can then use these to generate the mapping quality mask:
```
sbatch mapq_mask_wrapper.sh
```
When complete:
```
# concatenate and bedtools merge individual chromosome masks
for i in {1..22} X Y;
do
    cat masks/mapq/mapq_mask_chr${i}.txt | bedtools merge >> masks/mapq_mask.bed;
done

bgzip masks/mapq_mask.bed
```

## 3. Base quality mask
Once the `mapq_mask_prep.sh` script above has finished and we have the chromosome-wise alignments, we can start generating this mask. There is no need to wait for `mapq_mask_wrapper.sh` to finish. 
```
sbatch baseq_mask_wrapper.sh
```
When complete:
```
# concatenate and bedtools merge individual chromosome masks
for i in {1..22} X Y;
do
    cat masks/baseq/baseq_mask_chr${i}.txt | bedtools merge >> masks/baseq_mask.bed;
done

bgzip masks/baseq_mask.bed
```

## 4. Combined mask
Once the above 3 steps are finished, and we have the three separate mask files in `masks`: `coverage_mask.bed.gz`, `mapq_mask.bed.gz` and `baseq_mask.bed.gz`, run the following to merge these masks into the final mask:
```
bedtools multiinter \
  -header \
  -i masks/baseq_mask.bed masks/coverage_mask.bed masks/mapq_mask.bed \
  -names baseq coverage mapq |\
  awk '{if ($4 == 3) print $1"\t"$2"\t"$3}' > masks/combined_mask.bed

bgzip masks/combined_mask.bed
```
