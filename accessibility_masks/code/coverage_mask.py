import sys
import pysam
import pysamstats
import itertools
import statistics

#=============#
# Set globals #
#=============#

chrom = sys.argv[1]
min_norm_cov = 0.75
max_norm_cov = 1.25


#===============#
# Set filepaths #
#===============#

projectDir = sys.argv[2]
XX_ref = f"{projectDir}/data/references/chm13v2.0.XX.fasta"
XY_ref = f"{projectDir}/data/references/chm13v2.0.XY.fasta"
medianCovFile = f"{projectDir}/masks/median_coverage/median_cov.txt"
sampleInfoFile = f"{projectDir}/data/samples.info"
cramDir = f"{projectDir}/data/T2T-CHM13v2_crams"


#===============================#
# Generate median coverage mask #
#===============================#

# load dictionary of median autosomal coverage per sample
med_cov = {}
cov_file = open(medianCovFile, 'r')
Lines = cov_file.readlines()
for line in Lines:
	line = line.strip().split()
	med_cov[line[0]] = int(line[1])

samples = list(med_cov.keys())

# adjust expectations if sex chromosome
sample_sex = {}
sex_file = open(sampleInfoFile, 'r')
Lines = sex_file.readlines()
for line in Lines:
	line = line.strip().split()
	sample_sex[line[0]] = line[1]


if (chrom == "chrX_nonPAR"):
	for sample in samples:
		if (sample_sex[sample] == "F"):
			continue
		elif (sample_sex[sample] == "M"):
			med_cov[sample] = med_cov[sample] / 2

if (chrom == "chrX_PAR"):
	for sample in list(sample_sex.keys()):
		if (sample_sex[sample] == "F"):
			continue
		elif (sample_sex[sample] == "M"):
			continue

if (chrom == "chrY_nonPAR"):
	for sample in list(sample_sex.keys()):
		if (sample_sex[sample] == "F"):
			samples.remove(sample);
		elif (sample_sex[sample] == "M"):
			med_cov[sample] = med_cov[sample] / 2

# If you are using a different reference, you will have to adjus the "start"
# and "stop" values in the below code based on the locations of the chrX
# and chrY PARs (see "data/references/CHM13v2.0.PAR.bed")
pys_list = []
for idx, sample in enumerate(samples):
	if sample_sex[sample] == 'F':
		ref = XX_ref
	elif sample_sex[sample] == 'M':
		ref = XY_ref
	alignFile = pysam.AlignmentFile(f"{cramDir}/{sample}.cram", "rb", reference_filename = ref)
	if chrom == "chrX_nonPAR":
		sampleStats = pysamstats.stat_coverage(alignFile, chrom = "chrX", start = 2394410, end = 153925834, truncate = True, pad = True)
	elif chrom == 'chrX_PAR':
		sampleStats1 = pysamstats.stat_coverage(alignFile, chrom = "chrX", start = 0, end = 2394410, truncate = True, pad = True)
		sampleStats2 = pysamstats.stat_coverage(alignFile, chrom = "chrX", start = 153925834, end = 154259566, truncate = True, pad = True)
		sampleStats = itertools.chain(sampleStats1, sampleStats2)
	elif chrom == 'chrY_nonPAR':
		sampleStats = pysamstats.stat_coverage(alignFile, chrom = "chrY", start = 2458320, end = 62122809, truncate = True, pad = True)
	else:
		sampleStats = pysamstats.stat_coverage(alignFile, chrom = chrom, truncate = True, pad = True)
	pys_list.append(sampleStats)

for rec in zip(*pys_list):
	norm_cov = []
	for idx, pys in enumerate(rec):
		sample = samples[idx]
		norm_cov_sample = pys['reads_all'] / med_cov[sample]
		norm_cov.append(norm_cov_sample)
	mean_norm_cov = statistics.mean(norm_cov)
	if (mean_norm_cov > 0.75 and mean_norm_cov < 1.25):
		print(pys['chrom'], pys['pos'], pys['pos'] + 1, mean_norm_cov, sep = "\t")