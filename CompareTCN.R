source("~/Desktop/projects/FACETS_group/CompareTCN_functions.R")

library(tidyverse)

cncf_file = "~/Desktop/bbdd/msk_impact/FACETS/impact_facets_annotated.cncf.txt"
biomart_file = "~/Desktop/bbdd/linux_9may25/scripts/useful_annotations/biomart_gr.rds"
oncokb_file = "~/Desktop/bbdd/msk_impact/clinical/data_CNA.oncokb.txt"
out_dir = "~/Desktop/projects/FACETS_group/newheatmaptest/"

ramzi <- read_csv("Desktop/projects/FACETS_group/primary_brain_met_pairs.csv")
primary <- ramzi %>% group_by(patient_id) %>% slice_head() %>% pull(primary)
brain_met <- ramzi %>% group_by(patient_id) %>% slice_head() %>% pull(brain_met)

# primary <- primary[1:40]
# brain_met <- brain_met[1:40]
generate_paired_cohort_analysis(tcn_diff_threshold = 0.5,
                                cncf_path = cncf_file,
                                biomart_gr_path = biomart_file,
                                # Cohort_1 = c("P-0000952-T03-IM5_P-0000952-N01-IM5_hisens","P-0000952-T01-IM3_P-0000952-N01-IM3_hisens","P-0000894-T01-IM3_P-0000894-N01-IM3_hisens"),
                                Cohort_1 = brain_met,
                                # Cohort_2 = c("P-0000964-T01-IM3_P-0000964-N01-IM3_hisens","P-0000964-T02-IM3_P-0000964-N01-IM3_hisens","P-0000894-T02-IM6_P-0000894-N02-IM6_hisens"),
                                Cohort_2 = primary,
                                out_dir = out_dir,
                                oncokb_cnas_path = oncokb_file,
                                exact_match = FALSE,    # Enable grep matching
                                prefer_hisens = TRUE    # Prefer hisens samples when multiple matches
)
