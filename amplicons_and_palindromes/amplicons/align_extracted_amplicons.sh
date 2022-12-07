#!/bin/bash
#PBS -N amplicon_hg38.multialignment
#PBS -l select=1:ncpus=7:mem=100gb:scratch_local=100gb
#PBS -l walltime=12:00:00 
#PBS -m ae

set -e
set -x

# define a DATADIR variable: directory where the input files are taken from and where output will be copied to
DATADIR=/T2T/amplicon_extraction_HG002_v2.7/fa/hg38_Teitz/default 

echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" >> $DATADIR/jobs_info.txt

module load emboss-6.5.7 bbmap-38.42 bedtools-2.26.0 seqtk-1.0 clustal-wx-2.1

# test if scratch directory is set
# if scratch directory is not set, issue error message and exit
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

cp $DATADIR/*default.fa  $SCRATCHDIR || { echo >&2 "Error while copying input file(s)!"; exit 2; }

# move into scratch directory
cd $SCRATCHDIR 

#blue
emma -sequence blue.hg38.default.fa -outseq blue.hg38.default.aln -dendoutfile blue.hg38.default.emma &

#teal
emma -sequence teal.hg38.default.fa -outseq teal.hg38.default.aln -dendoutfile teal.hg38.default.emma &

#yellow
emma -sequence yellow.hg38.default.fa -outseq yellow.hg38.default.aln -dendoutfile yellow.hg38.default.emma &

#red
emma -sequence red.hg38.default.fa -outseq red.hg38.default.aln -dendoutfile red.hg38.default.emma &

#gray
emma -sequence gray.hg38.default.fa -outseq gray.hg38.default.aln -dendoutfile gray.hg38.default.emma &

#green
emma -sequence green.hg38.default.fa -outseq green.hg38.default.aln -dendoutfile green.hg38.default.emma &

wait

# move the output to user's DATADIR or exit in case of failure
cp *.default.aln $DATADIR/ || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }
cp *.default.emma $DATADIR/ || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }

echo "Done."

# clean the SCRATCH directory
clean_scratch














