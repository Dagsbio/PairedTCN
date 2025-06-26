# CompareTCN

Functions for analyzing and comparing Total Copy Number (TCN) profiles between paired and unpaired cancer sample cohorts.

## Overview

CompareTCN provides tools for:
- **Paired cohort analysis**: Compare matched samples (e.g., primary vs. metastatic tumors from the same patients)
- **Unpaired cohort analysis**: Compare independent cohorts of different sizes
- **Statistical testing**: Wilcoxon signed-rank test for paired data, Mann-Whitney U test for unpaired data
- **Visualization**: Heatmaps, volcano plots, Manhattan plots, and individual comparison plots
- **Gene annotation**: Integration with Biomart and OncoKB databases

## Features

### Paired Cohort Analysis
- Compares matched sample pairs (e.g., primary tumor vs. brain metastasis)
- Uses **Wilcoxon signed-rank test** for proper paired statistical testing
- Generates pairwise comparison plots showing TCN differences
- Creates summary statistics for each patient pair

### Unpaired Cohort Analysis  
- Compares two independent cohorts of potentially different sizes
- Uses **Mann-Whitney U test** for independent group comparisons
- Analyzes average TCN differences between cohorts
- Suitable for comparing different cancer types or treatment groups

## Visualizations

### Heatmaps
*Clustered heatmap showing TCN differences across genomic bins and sample comparisons*

<!-- ![Heatmap](plots/heatmap_example.png) -->

### Volcano Plots
*Volcano plot highlighting statistically significant genomic regions with gene annotations*

<!-- ![Volcano Plot](plots/volcano_plot_example.png) -->

### Manhattan Plots
*Manhattan plot displaying -log10(adjusted p-values) across chromosomes 1-22*

<!-- ![Manhattan Plot](plots/manhattan_plot_example.png) -->

### Pairwise Comparison Plots
*Individual comparison showing unique and shared TCN segments between paired samples*

<!-- ![Pairwise Comparison](plots/pairwise_comparison_example.png) -->

### Average TCN Difference Plot
*Genome-wide view of average TCN differences with oncogenic gene annotations*

<!-- ![Average TCN Difference](plots/average_tcn_diff_example.png) -->

## Installation

```r
# Install required packages
install.packages(c("tidyverse", "ggtext", "ggrepel"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "plyranges"))

# Development version of tidyheatmaps
devtools::install_github("stemangiola/tidyheatmaps")
