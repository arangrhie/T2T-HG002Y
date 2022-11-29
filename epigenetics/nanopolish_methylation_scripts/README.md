# Scripts and commands used to call CpG and GpC methylation with [Nanopolish](https://github.com/jts/nanopolish/tree/nanonome)
## PWH 22-10-24

This directory contains shell scripts used to call CpG and GpC methylation on HG002 nanoNOMe data using [Nanopolish](https://github.com/jts/nanopolish/tree/nanonome).
- [align_hg002.sh](align_hg002.sh): A collection of commands used to download reference genomes, prepare genomes for [Winnowmap](https://github.com/marbl/Winnowmap) alignment, and perform [Winnowmap](https://github.com/marbl/Winnowmap) alignment
- [merge_hg002.sh](merge_hg002.sh): Contains command to merge [Winnowmap](https://github.com/marbl/Winnowmap) alignments
- [filter_hg002.sh](filter_hg002.sh): Contains commands filtering merged alignments so only primary alignments are retained
- [cpggpc_nanonome_CHM13v2.0_nanopolish_call_meth.sh](cpggpc_nanonome_CHM13v2.0_nanopolish_call_meth.sh) and [cpggpc_nanonome_CHM13v2.7_nanopolish_call_meth.sh](cpggpc_nanonome_CHM13v2.7_nanopolish_call_meth.sh): Contains the Nanopolish command used to call CpG and GpC methylation for each genome. These scripts were called for all FASTQ files separately using [nanonome_cpggpc_call_meth.sh](nanonome_cpggpc_call_meth.sh)
- [nanonome_cpggpc_call_meth.sh](nanonome_cpggpc_call_meth.sh): Commands used to run Nanopolish CpG/GpC methylation calling scripts
- [cpg_nanonome_CHM13v2.0_nanopolish_processing.sh](cpg_nanonome_CHM13v2.0_nanopolish_processing.sh) and [cpg_nanonome_CHM13v2.7_nanopolish_processing.sh](cpg_nanonome_CHM13v2.7_nanopolish_processing.sh): Contains all steps used to process Nanopolish methylation calls into aggregated CpG and GpC frequencies. These involve heavy use of [nanopore-methylation-utilities](https://github.com/timplab/nanopore-methylation-utilities) 
