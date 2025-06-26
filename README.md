CompareTCN
A comprehensive R package for analyzing and comparing Total Copy Number (TCN) profiles between paired and unpaired cancer sample cohorts.

Overview
CompareTCN provides tools for:

Paired cohort analysis : Compare matched samples (e.g., primary vs. metastatic tumors from the same patients)
Unpaired cohort analysis : Compare independent cohorts of different sizes
Statistical testing : Wilcoxon signed-rank test for paired data, Mann-Whitney U test for unpaired data
Visualization : Heatmaps, volcano plots, Manhattan plots, and individual comparison plots
Gene annotation : Integration with Biomart and OncoKB databases
Features
Paired Cohort Analysis
Compares matched sample pairs (e.g., primary tumor vs. brain metastasis)
Uses Wilcoxon signed-rank test for proper paired statistical testing
Generates pairwise comparison plots showing TCN differences
Creates summary statistics for each patient pair
Unpaired Cohort Analysis
Compares two independent cohorts of potentially different sizes
Uses Mann-Whitney U test for independent group comparisons
Analyzes average TCN differences between cohorts
Suitable for comparing different cancer types or treatment groups
Visualizations
Heatmaps
Heatmap visualization showing TCN differences across the genome
Clustered heatmap showing TCN differences across genomic bins and sample comparisons

Volcano Plots
Volcano plot showing statistical significance vs effect size
Volcano plot highlighting statistically significant genomic regions with gene annotations

Manhattan Plots
Manhattan plot showing significance across chromosomes
Manhattan plot displaying -log10(adjusted p-values) across chromosomes 1-22

Pairwise Comparison Plots
Individual comparison showing TCN differences between paired samples
Example pairwise comparison showing unique and shared TCN segments

Average TCN Difference Plot
Bar plot showing average TCN differences across the genome
Genome-wide view of average TCN differences with oncogenic gene annotations

Installation
r
Save
Copy
1
2
3
4
5
6
7
8
9
10
# Install required packages
install.packages(c("tidyverse", "ggtext", "ggrepel"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("GenomicRanges", "plyranges"))

# Development version of tidyheatmaps
devtools::install_github("stemangiola/tidyheatmaps")
Usage
Paired Cohort Analysis
r
Save
Copy
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18
19
source("CompareTCN_functions.R")

# Define your paired cohorts
cohort_1 <- c("P-001-Primary", "P-002-Primary", "P-003-Primary")
cohort_2 <- c("P-001-Metastasis", "P-002-Metastasis", "P-003-Metastasis")

# Run paired analysis
generate_paired_cohort_analysis(
  cncf_path = "path/to/cncf_data.txt",
  biomart_gr_path = "path/to/biomart_annotations.rds",
  Cohort_1 = cohort_1,
  Cohort_2 = cohort_2,
  out_dir = "results/paired_analysis/",
  oncokb_cnas_path = "path/to/oncokb_cna_data.txt",
  bin_size = 1e6,
  tcn_diff_threshold = 1,
  exact_match = FALSE,
  prefer_hisens = TRUE
)
Unpaired Cohort Analysis
r
Save
Copy
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
# Define your independent cohorts (can be different sizes)
cohort_1_samples <- c("GBM-001", "GBM-002", "GBM-003", "GBM-004")
cohort_2_samples <- c("LGG-001", "LGG-002", "LGG-003")

# Run unpaired analysis
generate_unpaired_cohort_analysis(
  cncf_path = "path/to/cncf_data.txt",
  biomart_gr_path = "path/to/biomart_annotations.rds",
  cohort1_samples = cohort_1_samples,
  cohort2_samples = cohort_2_samples,
  out_dir = "results/unpaired_analysis/",
  oncokb_cnas_path = "path/to/oncokb_cna_data.txt",
  bin_size = 1e6,
  tcn_diff_threshold = 1,
  cohort1_name = "Glioblastoma",
  cohort2_name = "Low-grade Glioma"
)
Input Data Requirements
CNCF File
Tab-separated file with columns:

ID: Sample identifier
chrom: Chromosome number (1-22)
loc.start: Segment start position
loc.end: Segment end position
tcn: Total copy number
Biomart Annotations
RDS file containing a GRanges object with:

hgnc_symbol: Gene symbols
ensembl_gene_id: Ensembl gene IDs
Genomic coordinates
OncoKB CNA Data
Tab-separated file with columns:

HUGO_SYMBOL: Gene symbol
ONCOGENIC: Oncogenic classification
Output Files
Paired Analysis
heatmap.png: Clustered heatmap of TCN differences
heatmap_data.tsv: Matrix of TCN differences with gene annotations
volcano_plot.png: Statistical significance vs. effect size
manhattan_plot.png: Genome-wide significance plot
average_tcn_diff_barplot.png: Average differences across genome
diffPercent_lollipop.png: Percentage differences per comparison
sample_id_mapping.tsv: Sample ID resolution mapping
pairwise_plots/: Individual comparison plots and data
Unpaired Analysis
cohort_tcn_difference.png: TCN differences between cohorts
heatmap_cohort_comparison.png: Cohort comparison heatmap
combined_bin_data.csv: Detailed bin-level analysis
critical_regions.csv: Regions exceeding threshold
summary_statistics.csv: Cohort summary statistics
Statistical Methods
Paired Analysis
Wilcoxon signed-rank test : Tests if the median of paired differences is significantly different from zero
Paired t-test alternative : Non-parametric test suitable for copy number data
FDR correction : Benjamini-Hochberg correction for multiple testing
Unpaired Analysis
Mann-Whitney U test : Tests if two independent groups have different distributions
Wilcoxon rank-sum test : Non-parametric alternative to two-sample t-test
FDR correction : Benjamini-Hochberg correction for multiple testing
Parameters
Parameter
Description
Default
bin_size
Genomic bin size in base pairs
1e6
tcn_diff_threshold
Threshold for calling significant differences
1
exact_match
Use exact sample ID matching
TRUE
prefer_hisens
Prefer high-sensitivity samples when multiple matches
TRUE
p_threshold
P-value threshold for significance
0.05
effect_threshold
Effect size threshold
1

Advanced Features
Sample ID Resolution
Flexible sample matching with:

Exact string matching
Partial string matching (grep)
Priority for high-sensitivity samples
Detailed mapping output
Gene Annotation Pipeline
Genome-wide binning at specified resolution
Overlap detection with gene annotations
Integration with oncogenic gene databases
Multi-gene bin annotation support
Quality Control
Sample availability checking
Overlap validation
Missing data handling
Statistical power assessment
Examples
Basic Paired Analysis
r
Save
Copy
1
2
3
4
5
6
7
8
# Simple paired analysis with default settings
generate_paired_cohort_analysis(
  cncf_path = "data.cncf",
  biomart_gr_path = "annotations.rds", 
  Cohort_1 = primary_samples,
  Cohort_2 = metastasis_samples,
  out_dir = "results/"
)
Custom Unpaired Analysis
r
Save
Copy
1
2
3
4
5
6
7
8
9
10
11
12
# Unpaired analysis with custom parameters
generate_unpaired_cohort_analysis(
  cncf_path = "data.cncf",
  biomart_gr_path = "annotations.rds",
  cohort1_samples = treatment_group,
  cohort2_samples = control_group,
  out_dir = "treatment_comparison/",
  bin_size = 5e5,  # 500kb bins
  tcn_diff_threshold = 0.5,
  cohort1_name = "Treatment",
  cohort2_name = "Control"
)
Citation
If you use CompareTCN in your research, please cite:

Save
Copy
1
[Citation to be added]
License
This project is licensed under the MIT License - see the LICENSE file for details.

Contact
For questions and support, please open an issue on GitHub or contact [your-email@institution.edu ].
