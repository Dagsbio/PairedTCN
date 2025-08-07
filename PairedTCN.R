source("CompareTCN_functions.R")

library(tidyverse)

cncf_file = "impact_facets_annotated.cncf.txt"
biomart_file = "biomart_gr.rds"
oncokb_file = "data_CNA.oncokb.txt"
out_dir = "test"

generate_paired_cohort_analysis(tcn_diff_threshold = 0.5,
                                cncf_path = cncf_file,
                                biomart_gr_path = biomart_file,
                                Cohort_1_EXAMPLE = c("P-0000952-T03-IM5_P-0000952-N01-IM5_hisens","P-0000952-T01-IM3_P-0000952-N01-IM3_hisens","P-0000894-T01-IM3_P-0000894-N01-IM3_hisens"),
                                Cohort_2_EXAMPLE = c("P-0000964-T01-IM3_P-0000964-N01-IM3_hisens","P-0000964-T02-IM3_P-0000964-N01-IM3_hisens","P-0000894-T02-IM6_P-0000894-N02-IM6_hisens"),
                                out_dir = out_dir,
                                oncokb_cnas_path = oncokb_file,
                                exact_match = FALSE,    # Enable grep matching
                                prefer_hisens = TRUE    # Prefer hisens samples when multiple matches
)
