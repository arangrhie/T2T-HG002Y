### Motivation

_The goal of this analysis is to aid in determining the boundaries between X-degenerate and ampliconic regions on the human Y chromosome by calculating the intrachromosomal similarity along the Y chromosome._

### Analysis

1. Split the newly assembled hardmasked Y chromosome into windows and map them back to the assembly with Winnowmap (Jain et al. 2022)
   ```sh
   ./split_Y_chromosome_into_windows_and_calculate_intrachromosomal_similarity.sh
   ```
2. Plot the interchromosomal similarity along the Y chromosome
   ```sh
   plot_intrachromosomal_identity.Rnw
   ```
