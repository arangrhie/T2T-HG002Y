### THIS CODE WAS ORIGINALLY AUTHORED BY Wouter De Coster
### described in the following blog post: https://gigabaseorgigabyte.wordpress.com/2017/04/14/getting-the-edit-distance-from-a-bam-alignment-a-journey/
### stored in the following gist repository: https://gist.github.com/wdecoster/a0ca604306b5778a7167f705c597ee38
### the code below is presented in order to enhance the reproducibility 
### note that only very minor changes are present here compared to the original version

import sys
import os
import re
#import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats
import pysam


def processBam(bam):
	'''
	Processing function: calls pool of worker functions
	to extract from a bam file two definitions of the edit distances to the reference genome scaled by read length
	Returned in a pandas DataFrame
	'''
	samfile = pysam.AlignmentFile(bam, "rb")
	if not samfile.has_index():
		pysam.index(bam)
		samfile = pysam.AlignmentFile(bam, "rb")  # Need to reload the samfile after creating index
	chromosomes = samfile.references
	datadf = pd.DataFrame()
	output = extractFromBam([bam,"chrY_hg002"])

	#datadf["editDistancesNM"] = np.array([x for y in [elem[0] for elem in output] for x in y])
	#datadf["editDistancesMD"] = np.array([x for y in [elem[1] for elem in output] for x in y])

	with open('identities.txt', 'w') as f:
		f.write('\n'.join([ str((1-myelement)*100) for myelement in output[0] ]))

	return datadf


def extractFromBam(params):
	'''
	Worker function per chromosome
	loop over a bam file and create tuple with lists containing metrics:
	two definitions of the edit distances to the reference genome scaled by aligned read length
	'''
	bam, chromosome = params
	samfile = pysam.AlignmentFile(bam, "rb")
	editDistancesNM = []
	editDistancesMD = []
	for read in samfile.fetch(reference=chromosome, multiple_iterators=True):
		editDistancesNM.append(read.get_tag("NM")/read.query_alignment_length)
		editDistancesMD.append(
			(sum([len(item) for item in re.split('[0-9^]', read.get_tag("MD"))]) +  # Parse MD string to get mismatches/deletions
			sum([item[1] for item in read.cigartuples if item[0] == 1]))  # Parse cigar to get insertions
			/read.query_alignment_length)
	return (editDistancesNM, editDistancesMD)



df = processBam(sys.argv[1])

