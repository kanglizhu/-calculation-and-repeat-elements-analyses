CBS Analysis Pipeline: Alpha Calculation & Repeat Element Permutation Test

This repository contains pipelines for calculating Alpha values with confidence intervals for Copy Number Breakpoint Segments (CBSs) and performing stratified permutation tests to analyze the association between CBSs and Transposable Elements (TEs).

Requirements

Python Environment

python >= 3.7
pandas >= 1.5.0
numpy >= 1.23.0
matplotlib >= 3.6.0


R Environment

ggplot2
dplyr
stringr
scales
gridExtra
ggrepel
data.table
GenomicRanges
IRanges
tidyr
gridExtra


Part 1: Alpha Calculation

Calculate Alpha values and 95% Confidence Intervals (CI) for different types of CBSs.

1. Data Preparation

For each type of CBS, prepare your input file following the format of example_for_alpha_calculation.txt located in the data/alpha_calculation folder.

2. Run Calculation

Execute the bootstrap script to calculate Alpha and generate a preliminary plot:

bash
python3 bootstrap_alpha.py <input_file.txt> <output_plot.png>
Example:

bash
python3 bootstrap_alpha.py data/alpha_calculation/example_for_alpha_calculation.txt data/alpha_calculation/example_for_alpha_calculation.png
3. Visualization

After running the Python script:

Extract the Alpha value and its 95% CI for each CBS type from the output.
Compile these values into a summary file.
Use the R script plot_alpha_ci.R to generate the final figure:

bash
Rscript plot_alpha_ci.R

Part 2: Repeat Element Permutation Test

This module performs stratified permutation tests to evaluate the enrichment of Transposable Elements (TEs) in CBSs while controlling for genomic confounders.

1. Data Stratification Criteria

To ensure comparable lengths and genomic contexts, data is stratified based on the following features:

GC Content: Binned into Low (<0.3), Mid (0.3–0.5), and High (>0.5).
Mappability: Quantified using ENCODE-style mappability tracks, binned into Low (<0.7) and High (≥0.7).
Distance to Nearest TAD Boundary: Categorized as 'Near' (0–50 kb) or 'Far' (>50 kb).
Distance to Nearest Gene: Categorized as 'Promoter Proximal' (0–2 kb) or 'Intergenic' (>2 kb).
Note: The CBSs used in the "Supercool" file are derived from ENCODE data. An example input file is available in data/repeat_element/.
2. Run Permutation Test

Run the global stratified permutation test. The default script tests for CBSs with 40% overlapping TE, but this overlap criteria can be modified directly in the script.

bash
python3 run_stratified_te_global_permutation_0.4.py <input_file.txt>
Example:

bash
python3 run_stratified_te_global_permutation_0.4.py data/repeat_element/example_permutation_test.txt

3. Extract Metrics & Plot Robustness

After the permutation test completes:

Extract Data: Retrieve the P-value and Observed Fraction for each type of CBS from the output files.
Visualize: Use the R script plot_multi_group_robustness.R to plot the fractions across different groups.

bash
Rscript plot_multi_group_robustness.R

4. TE Type & Family Contribution Analysis

To illustrate the contribution of specific TE types and families to different CBS types, perform the permutation test using the R code:

bash
Rscript te_analysis_rstudio.R
This script generates a CSV file containing enrichment statistics.

5. Visualization of Enriched Terms

Use the generated CSV file to plot significantly enriched terms using the volcano plot script:

bash
Rscript Plot_TE_Log2FE_NewFormat.R
