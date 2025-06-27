source("~/Desktop/projects/FACETS_group/CompareTCN_functions.R")

# Example data:
# ramzi <- read_csv("Desktop/projects/FACETS_group/primary_brain_met_pairs.csv")
# 
# primaries <- ramzi %>% group_by(patient_id) %>% slice_head %>% pull(primary)
# brain_mets <- ramzi %>% group_by(patient_id) %>% slice_head %>% pull(brain_met)

library(tidyverse)

cncf_file = "~/Desktop/bbdd/msk_impact/FACETS/impact_facets_annotated.cncf.txt"
biomart_file = "~/Desktop/bbdd/linux_9may25/scripts/useful_annotations/biomart_gr.rds"
oncokb_file = "~/Desktop/bbdd/msk_impact/clinical/data_CNA.oncokb.txt"
out_dir = "~/Desktop/projects/FACETS_group/test_ramzi3/"

generate_paired_cohort_analysis(tcn_diff_threshold = 0.5,
                                cncf_path = cncf_file,
                                biomart_gr_path = biomart_file,
                                Cohort_1 = brain_mets,
                                Cohort_2 = primaries,
                                out_dir = out_dir,
                                oncokb_cnas_path = oncokb_file,
                                exact_match = FALSE,    # Enable grep matching
                                prefer_hisens = TRUE    # Prefer hisens samples when multiple matches
)
