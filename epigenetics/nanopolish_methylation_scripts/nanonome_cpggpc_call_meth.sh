#!/bin/bash

# PWH 10-19-2022
# The purpose of this script is to start methylation calling for the CHM13 v2.0 and v2.7 genomes. Note that these commands were actually run in parallel to speed up the process.

# cpggpc model - v2.7
./cpggpc_nanonome_CHM13v2.7_nanopolish_call_meth.sh 201119_HG002_nanonome_SREXL.fastq &> logs/all_cpggpc_nanonome_CHM13_v2.7_201119_HG002_nanonome_SREXL.log
./cpggpc_nanonome_CHM13v2.7_nanopolish_call_meth.sh 201126_HG002_nanonome_SRE.fastq &> logs/all_cpggpc_nanonome_CHM13_v2.7_201126_HG002_nanonome_SRE.log
./cpggpc_nanonome_CHM13v2.7_nanopolish_call_meth.sh 201201_HG002_nanonome_SREXL_4-5.fastq &> logs/all_cpggpc_nanonome_CHM13_v2.7_201201_HG002_nanonome_SREXL_4-5.log

# cpggpc model - v2.0
./cpggpc_nanonome_CHM13v2.0_nanopolish_call_meth.sh 201119_HG002_nanonome_SREXL.fastq &> logs/all_cpggpc_nanonome_CHM13_v2.0_201119_HG002_nanonome_SREXL.log
./cpggpc_nanonome_CHM13v2.0_nanopolish_call_meth.sh 201126_HG002_nanonome_SRE.fastq &> logs/all_cpggpc_nanonome_CHM13_v2.0_201126_HG002_nanonome_SRE.log
./cpggpc_nanonome_CHM13v2.0_nanopolish_call_meth.sh 201201_HG002_nanonome_SREXL_4-5.fastq &> logs/all_cpggpc_nanonome_CHM13_v2.0_201201_HG002_nanonome_SREXL_4-5.log
