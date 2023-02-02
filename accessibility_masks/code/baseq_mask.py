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


#============================#
# Generate base quality mask #
#============================#

aln = pysam.AlignmentFile(f"{cramDir}/{chrom}_merged.cram", "rb", reference_filename = ref)

for rec in aln.pileup(chrom, truncate = True):
  baseq = rec.get_query_qualities()
  if (len(baseq) > 0):
    q20_baseq = [i for i in baseq if i > 20]
    q20_prop = len(q20_baseq) / len(baseq)
    if (q20_prop >= 0.9):
      print(rec.reference_name, rec.pos, rec.pos + 1, q20_prop, sep = "\t")