CompareTCN
R functions to compare Total Copy Number (TCN) starting from copy number segment calls (such as FACETS cncf file).

Patient Analysis Function
This script provides a function to analyze copy number variation (CNV) data for specified patients, generate plots, and produce annotations. The script is organized as a group of functions, allowing straightforward reusability.

Authors
David Gómez Sánchez (gomezsd@mskcc.org )
Ramzi Homsi (homsir1@mskcc.org )
Usage
Function Name: generate_patient_analysis()
Parameters
cncf_path : Path to the CNV file in TSV format.
biomart_gr_path : Path to the biomart annotations in RDS format (GRanges).
patient_ids : Vector of patient IDs to process.
outDir : Directory to save the results.
oncokb_cnas_path : Path to OncoKB data containing CNA information in TSV format.
Example
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
generate_patient_analysis(
  cncf_path = "path/to/impact_facets_annotated.cncf.txt",
  biomart_gr_path = "path/to/useful_annotations/biomart_gr.rds",
  patient_ids = c("P-0070753", "P-0001876", "P-0000004", "P-0000964"),
  outDir = "output_directory/",
  oncokb_cnas_path = "path/to/data_CNA.oncokb.txt"
)
Output
Comparison Plots : Generates PNG images for patient comparisons.
Differential Segments Data : Saves TSV files containing differential segment data.
Annotated Genes : Produces TSV files with gene annotations overlapping differential segments.
Heatmap : Creates visual representation of CN differences across samples.
Average TCN Difference Barplot : Generates barplot summarizing gained/lost regions with OncoKB annotations.
Installation
Ensure the following R packages are installed before running the script:

r
Save
Copy
1
2
install.packages(c("tidyverse", "ggtext", "tidyheatmaps", "ggrepel"))
BiocManager::install(c("GenomicRanges", "plyranges"))
License
Include the license under which the project is distributed (e.g., MIT License, GPL-3, etc.).

Contact
For further questions or inquiries:

David Gómez Sánchez : gomezsd@mskcc.org
Ramzi Homsi : homsir1@mskcc.org
This README includes headings and lists that are formatted using markdown syntax, making it suitable for GitHub. If there are any specific requirements not addressed or further adjustments needed, please let me know!
