# CompareTCN

R functions to compare Total Copy Number (TCN) starting from copy number segment calls (such as FACETS cncf file).

## Description

This script provides functions to analyze copy number variation (CNV) data for specified patients, generate plots, and produce annotations. It is organized to allow straightforward reusability.

## Authors

- **David Gómez Sánchez** (gomezsd@mskcc.org)
- **Ramzi Homsi** (homsir1@mskcc.org)

## Usage

### Function Name: `generate_patient_analysis()`

#### Parameters

- **`cncf_path`**: Path to the CNV file in TSV format.
- **`biomart_gr_path`**: Path to the biomart annotations in RDS format (GRanges).
- **`patient_ids`**: Vector of patient IDs to process.
- **`outDir`**: Directory to save the results.
- **`oncokb_cnas_path`**: Path to OncoKB data containing CNA information in TSV format.

### Example

```r
generate_patient_analysis(
  cncf_path = "path/to/impact_facets_annotated.cncf.txt",
  biomart_gr_path = "path/to/useful_annotations/biomart_gr.rds",
  patient_ids = c("P-0070753", "P-0001876", "P-0000004", "P-0000964"),
  outDir = "output_directory/",
  oncokb_cnas_path = "path/to/data_CNA.oncokb.txt"
)
