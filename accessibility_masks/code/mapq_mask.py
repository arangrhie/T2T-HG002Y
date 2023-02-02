import sys
import pysam
import statistics

#=============#
# Set globals #
#=============#

chrom = sys.argv[1]


#===============#
# Set filepaths #
#===============#

projectDir = sys.argv[2]
ref = f"{projectDir}/data/references/chm13v2.0.fasta"
cramDir = f"{projectDir}/data/T2T-CHM13v2_crams"


#====================#
# Generate mapq mask #
#====================#

aln = pysam.AlignmentFile(f"{cramDir}/{chrom}_merged.cram", "rb", reference_filename = ref)

for rec in aln.pileup(chrom, truncate = True):
  mapq = rec.get_mapping_qualities()
  if (len(mapq) > 0):
    mean_mapq = statistics.mean(mapq)
    if (mean_mapq >= 50):
      print(rec.reference_name, rec.pos, rec.pos + 1, mean_mapq, sep = "\t")