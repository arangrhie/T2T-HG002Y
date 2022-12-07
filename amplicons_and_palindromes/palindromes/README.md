### Motivation

_The goal of this analysis is to identify palindromes in the Y-chromosome assemblies by using dotplot analysis, followed by finding precise coordinates by an alignment._

### Analysis

1. Determine approximate boundaries of palindromes using dotplot (Gepard by Krumsiek et al. 2007; PMID: 17309896)

2. Extract sequences of palindrome arms using approximate coordinates
   ```sh
   ./1_extract_palindrome_arms_from_dotplot.sh
   ```
3. Align palindrome arms (with flanks) to each other using Stretcher (by EMBOSS) in order to determine the precise boundary
   ```sh
   ./2_align_palindrome_arms_that_include_flanks.sh
   ```
4. Extract the precise boundaries of palindromes
   ```sh
   ./3_extract_palindrome_arms_after_stitcher_adjustment.sh
   ```
5. Align palindrome arms in order to calculate arm-to-arm identity
   ```sh
   ./4_align_adjusted_palindrome_arms.sh
   ```
