# Processing of HG002 WGBS short-read, EM-seq short-read, and PacBio HiFi data
## PWH 22-10-24

This directory contains the [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline used to process HG002 short-read WGBS and EM-seq data from [Foox 2021](https://pubmed.ncbi.nlm.nih.gov/34872606/). The data are mainly processed using scripts from the popular methylation caller, [Bismark](https://github.com/FelixKrueger/Bismark).

Methylation calls from PacBio HiFi reads produced with [Primrose](https://github.com/PacificBiosciences/primrose) were also processed with a rule that uses [modbam2bed](https://github.com/epi2me-labs/modbam2bed) in this pipeline.

All rules for data processing can be found here: [/workflow/rules](https://github.com/timplab/pwh_projects/tree/master/hg002_methylation/wgbs_pacbio_processing/workflow/rules).
