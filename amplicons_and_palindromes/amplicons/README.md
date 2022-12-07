### Motivation

_The goal of this analysis is to determine the coordinates of amplicons in the human Y chromosome assembly, using amplicon annotations as provided in Teitz et al. 2018._

### Analysis

1. Split the newly assembled hardmasked Y chromosome into windows and map them back to the assembly with Winnowmap 
   ```sh
   ./1_extract_amplicons.sh
   ```
2. Map the extracted amplicons to the Y chromosome using Winnowmap and determine their coordinates
   ```sh
   2_map_extracted_amplicons_to_T2T_Y.sh
   ```

Optionally, extracted amplicons can be used for a multi-alignment
   ```sh
   ./align_extracted_amplicons.sh
   ```