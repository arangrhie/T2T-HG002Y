# HG002 ChrY Methylation Analysis
## Code associated with methylation analyses for the coming Rhie et al., HG002 ChrY manuscript
## Author: Paul W. Hook

This repo includes code analyzing the following HG002 data:
- CpG methylation from PacBio HiFi reads using Primrose: [wgbs_pacbio_processing](wgbs_pacbio_processing)
- CpG and GpC methylation calls from ONT reads from nanoNOMe experiments called with Nanopolish: [nanopolish_methylation_scripts](nanopolish_methylation_scripts)
- CpG methylation calls from ultra-long ONT reads called with Remora and Guppy 6.1.2: [remora_methylation_processing](remora_methylation_processing)
- CpG methylation calls from WGBS and EM-seq short read datasets: [wgbs_pacbio_processing](wgbs_pacbio_processing)

Comparison of the methylation identification technologies as well as the code used to produce figures can can be found in the [r_analysis](r_analysis) directory.

The code for production of summary methylation files (bigWigs and beds) can be found in the [r_analysis](r_analysis) directory.

Steps on how Nanopolish downloaded and installed are found here: [221019_nanopolish_nanonome_setup.pdf](221019_nanopolish_nanonome_setup.pdf)

Information on sequencing data resources can be found here: [hg002_chrY_methylation_resources.csv](hg002_chrY_methylation_resources.csv)
