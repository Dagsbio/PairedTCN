#' Common setup and helper functions
initialize_analysis <- function(cncf_path, biomart_gr_path) {
  # Load required libraries
  suppressPackageStartupMessages({
    library(tidyverse)
    library(plyranges)
    library(ggtext)
    library(GenomicRanges)
    library(tidyheatmaps)
    library(ggrepel)
    library(viridis) 
  })
  
  # Define Chromosome lengths - only autosomes (no X/Y)
  ChrLengths <- c(
    249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
    159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
    115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
    59128983, 63025520, 48129895, 51304566
  )
  names(ChrLengths) <- as.character(1:22)
  
  # Load data
  cat("Loading data...\n")
  cncf <- read_tsv(cncf_path, show_col_types = FALSE) %>%
    filter(!(chrom %in% c(23, 24, "X", "Y")))
  
  biomart_gr <- read_rds(biomart_gr_path)
  
  return(list(
    ChrLengths = ChrLengths,
    cncf = cncf,
    biomart_gr = biomart_gr
  ))
}

create_Chromosome_colors <- function(chrom_levels) {
  chrom_colors <- rep(NA, length(chrom_levels))
  names(chrom_colors) <- chrom_levels
  for (i in seq_along(chrom_levels)) {
    chr <- chrom_levels[i]
    chr_num <- as.numeric(chr)
    chrom_colors[i] <- ifelse(chr_num %% 2 == 0, "grey", "black")
  }
  return(chrom_colors)
}

load_oncokb_data <- function(oncokb_cnas_path) {
  tryCatch({
    read_tsv(oncokb_cnas_path, show_col_types = FALSE) %>%
      filter(ONCOGENIC %in% c("Likely Oncogenic", "Oncogenic")) %>%
      select(gene = HUGO_SYMBOL, ONCOGENIC)
  }, error = function(e) {
    message("Error loading OncoKB data: ", e$message)
    tibble(gene = character(0), ONCOGENIC = character(0))
  })
}

create_bins <- function(ChrLengths, bin_size) {
  chrom_bins_df <- lapply(names(ChrLengths), function(chr) {
    max_end <- ChrLengths[chr]
    intervals <- seq(from = 1, to = max_end, by = bin_size)
    tibble(seqnames = chr, start = intervals[-length(intervals)], end = intervals[-1] - 1)
  }) %>% bind_rows()
  
  chrom_bins_df <- chrom_bins_df %>%
    mutate(bin_id = paste0(seqnames, ":", start, "-", end))
  
  return(chrom_bins_df)
}

annotate_heatmap_data <- function(heatmap_data_path, biomart_gr, bin_size = 1e6) {
  # Read the heatmap data
  heatmap_data <- read_tsv(heatmap_data_path, show_col_types = FALSE)
  
  # Convert heatmap bins to GRanges
  bins_gr <- GRanges(
    seqnames = str_extract(heatmap_data$bin, "^[^:]+"),
    ranges = IRanges(
      start = as.numeric(str_extract(heatmap_data$bin, "(?<=:)[0-9]+")),
      width = bin_size
    )
  )
  
  # Find overlaps with the biomart file
  overlapping_genes <- findOverlaps(bins_gr, biomart_gr)
  
  # Annotate genes to bins - collect ALL overlapping genes
  if (length(overlapping_genes) > 0) {
    annotated_genes <- tibble(
      bin = heatmap_data$bin[queryHits(overlapping_genes)],
      gene = biomart_gr$hgnc_symbol[subjectHits(overlapping_genes)]
    ) %>%
      filter(!is.na(gene) & gene != "") %>%
      group_by(bin) %>%
      summarize(gene = paste(unique(gene), collapse = ";"), .groups = "drop")
  } else {
    annotated_genes <- tibble(bin = character(0), gene = character(0))
  }
  
  # Add gene annotation to heatmap_data
  heatmap_data_annotated <- heatmap_data %>%
    left_join(annotated_genes, by = "bin")
  
  # Save updated heatmap data with gene annotations
  write_tsv(heatmap_data_annotated, heatmap_data_path)
  
  return(heatmap_data_annotated)
}

generate_heatmap <- function(heatmap_long, out_dir, title_suffix = "") {
  # Create colors for Chromosomes
  chrom_levels <- levels(heatmap_long$Chromosome)
  chrom_colors <- create_Chromosome_colors(chrom_levels)
  
  # Create simplified comparison names for display
  heatmap_long <- heatmap_long %>%
    mutate(
      Patient_id1 = str_extract(comparison, "^P-[0-9]+"),
      Patient_id2 = str_extract(comparison, "P-[0-9]+(?=.*vs.*P-[0-9]+)"),
      Comparison = str_replace(comparison, "^P-[0-9]+-T[^_]*_(.*?)_vs_P-[0-9]+-T[^_]*_(.*)", "\\1_vs_\\2")
    )
  
  # Generate the heatmap with clustering
  hm <- tidyheatmaps::tidy_heatmap(
    heatmap_long,
    fontsize = 12,
    rows = Comparison,
    columns = bin,
    values = tcn_diff,
    colors = c("#145afc", "white", "#ee4445"),
    color_legend_min = -3,
    color_legend_max = 3,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    show_colnames = FALSE,
    annotation_row = Patient_id,
    annotation_col = Chromosome,
    annotation_colors = list(
      Chromosome = chrom_colors
    )
  )
  
  # Extract cluster information from the heatmap object
  if ("ComplexHeatmap" %in% class(hm)) {
    # For ComplexHeatmap objects
    row_dend <- row_dend(hm)
    if (!is.null(row_dend)) {
      # Extract cluster assignments - you can adjust k (number of clusters) as needed
      k_clusters <- 3  # Adjust this number based on your needs
      cluster_assignments <- cutree(as.hclust(row_dend), k = k_clusters)
      
      # Create cluster mapping
      cluster_df <- tibble(
        Comparison = names(cluster_assignments),
        cluster = as.numeric(cluster_assignments)
      )
    } else {
      # Fallback if no dendrogram
      unique_comparisons <- unique(heatmap_long$Comparison)
      cluster_df <- tibble(
        Comparison = unique_comparisons,
        cluster = 1  # All in one cluster if no clustering
      )
    }
  } else {
    # Fallback for other heatmap types - perform manual clustering
    heatmap_wide_for_clustering <- heatmap_long %>%
      select(bin, Comparison, tcn_diff) %>%
      group_by(bin, Comparison) %>%
      summarise(tcn_diff = mean(tcn_diff, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = bin, values_from = tcn_diff, values_fill = 0)
    
    # Perform hierarchical clustering on rows (comparisons)
    if (nrow(heatmap_wide_for_clustering) > 1) {
      row_data <- as.matrix(heatmap_wide_for_clustering[, -1])
      rownames(row_data) <- heatmap_wide_for_clustering$Comparison
      
      # Calculate distance and perform clustering
      row_dist <- dist(row_data)
      row_hclust <- hclust(row_dist)
      
      # Cut tree to get clusters
      k_clusters <- min(3, nrow(heatmap_wide_for_clustering))  # Max 3 clusters or number of rows
      cluster_assignments <- cutree(row_hclust, k = k_clusters)
      
      cluster_df <- tibble(
        Comparison = names(cluster_assignments),
        cluster = as.numeric(cluster_assignments)
      )
    } else {
      cluster_df <- tibble(
        Comparison = heatmap_wide_for_clustering$Comparison,
        cluster = 1
      )
    }
  }
  
  # Add cluster information to heatmap_long
  heatmap_long_with_clusters <- heatmap_long %>%
    left_join(cluster_df, by = "Comparison")
  
  # If clusters weren't assigned properly, assign all to cluster 1
  if (any(is.na(heatmap_long_with_clusters$cluster))) {
    heatmap_long_with_clusters$cluster[is.na(heatmap_long_with_clusters$cluster)] <- 1
  }
  
  # Add percentage difference information - read from diffPercents if available
  diff_percents_path <- file.path(out_dir, "diffPercents.tsv")
  if (file.exists(diff_percents_path)) {
    diff_percents_data <- read_tsv(diff_percents_path, show_col_types = FALSE) %>%
      mutate(Comparison = str_replace(comparison, "^P-[0-9]+-T[^_]*_(.*?)_vs_P-[0-9]+-T[^_]*_(.*)", "\\1_vs_\\2"))
    
    heatmap_long_with_clusters <- heatmap_long_with_clusters %>%
      left_join(diff_percents_data %>% select(Comparison, diffPercent), by = "Comparison")
  } else {
    # If diffPercents not available, set to NA (will be handled later)
    heatmap_long_with_clusters$diffPercent <- NA
  }
  
  # Create cluster colors
  n_clusters <- length(unique(cluster_df$cluster))
  
  # Professional color palette for clusters
  if (n_clusters <= 8) {
    # Use ColorBrewer Set2 palette - professional and distinguishable
    professional_colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", 
                             "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")
    cluster_colors <- professional_colors[1:n_clusters]
  } else if (n_clusters <= 12) {
    # Use ColorBrewer Set3 for more clusters
    professional_colors <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", 
                             "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
                             "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")
    cluster_colors <- professional_colors[1:n_clusters]
  } else {
    # For many clusters, use a more sophisticated approach
    # Generate colors using viridis palette (professional and colorblind-friendly)
    cluster_colors <- viridis::viridis(n_clusters, option = "D", alpha = 0.8)
  }
  
  names(cluster_colors) <- paste0("Cluster_", sort(unique(cluster_df$cluster)))
  
  # Create professional gradient for percentage difference (0% to 100%)
  # Using a single color gradient - navy blue professional gradient
  diff_percent_colors <- colorRampPalette(c("#F8F9FA", "#1E3A8A"))(100)  # Light gray to navy blue
  names(diff_percent_colors) <- as.character(0:99)
  
  # Alternative professional gradients you can use instead:
  # Option 1: Warm professional (cream to dark brown)
  # diff_percent_colors <- colorRampPalette(c("#FEF7E0", "#92400E"))(100)
  
  # Option 2: Cool professional (light blue to dark teal)  
  # diff_percent_colors <- colorRampPalette(c("#F0F9FF", "#0F766E"))(100)
  
  # Option 3: Monochrome professional (light gray to charcoal)
  # diff_percent_colors <- colorRampPalette(c("#F9FAFB", "#374151"))(100)
  
  # Add cluster annotation to the comparison names
  cluster_df <- cluster_df %>%
    mutate(Cluster_label = paste0("Cluster_", cluster))
  
  heatmap_long_with_clusters <- heatmap_long_with_clusters %>%
    mutate(Cluster_label = paste0("Cluster_", cluster))
  
  # Prepare percentage difference for annotation
  heatmap_long_with_clusters <- heatmap_long_with_clusters %>%
    mutate(
      # Convert diffPercent to character for color mapping, handle NAs
      diffPercent_char = ifelse(is.na(diffPercent), "0", 
                                as.character(pmin(99, pmax(0, round(diffPercent)))))
    )
  
  # Regenerate heatmap with cluster and percentage annotations
  hm_with_clusters <- tidyheatmaps::tidy_heatmap(
    heatmap_long_with_clusters,
    fontsize = 8,  # Reduced overall font size
    fontsize_row = 8,  # Specific font size for row names
    rows = Comparison,
    columns = bin,
    values = tcn_diff,
    colors = c("#145afc", "white", "#ee4445"),
    color_legend_min = -3,
    color_legend_max = 3,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    show_colnames = FALSE,
    annotation_row = c(cluster_label, diffPercent_char),  # Add both annotations
    annotation_col = Chromosome,
    annotation_colors = list(
      Chromosome = chrom_colors,
      cluster_label = cluster_colors,
      diffPercent_char = diff_percent_colors
    ),
    annotation_legend_param = list(
      diffPercent_char = list(title = "% TCN\nDifference"),
      cluster_label = list(title = "Cluster"),
      Chromosome = list(title = "Chromosome")
    )
  )
  
  # Save the heatmap with clusters
  ggsave(
    file.path(out_dir, paste0("heatmap", title_suffix, ".png")),
    plot = hm_with_clusters,
    width = 20, height = 15, dpi = 300
  )
  
  # Save heatmap data in wide format WITH cluster information
  heatmap_wide_output <- heatmap_long_with_clusters %>%
    select(bin, comparison, tcn_diff, cluster) %>%  # Include cluster column
    # Handle duplicates by taking the mean
    group_by(bin, comparison, cluster) %>%
    summarise(tcn_diff = mean(tcn_diff, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = comparison, values_from = tcn_diff, values_fill = 0)
  
  # Also save a cluster mapping file with percentage difference
  cluster_mapping <- heatmap_long_with_clusters %>%
    select(comparison, Comparison, cluster, cluster_label, diffPercent) %>%
    distinct() %>%
    arrange(cluster, comparison)
  
  write_tsv(cluster_mapping, file.path(out_dir, paste0("cluster_mapping", title_suffix, ".tsv")))
  write_tsv(heatmap_wide_output, file.path(out_dir, paste0("heatmap_data", title_suffix, ".tsv")))
  
  # Print cluster summary
  cat("Cluster Summary:\n")
  cluster_summary <- cluster_mapping %>%
    count(cluster, name = "n_comparisons")
  print(cluster_summary)
  
  return(hm_with_clusters)
}

# New function for paired volcano plot with correct statistics
generate_volcano_plot_paired <- function(heatmap_data_path, out_dir, biomart_gr, p_threshold = 0.05, effect_threshold = tcn_diff_threshold) {
  # First annotate the heatmap data with all overlapping genes
  cat("Annotating bins with overlapping genes...\n")
  annotated_heatmap_data <- annotate_heatmap_data(heatmap_data_path, biomart_gr)
  
  # Prepare data for PAIRED statistical testing
  heatmap_long <- annotated_heatmap_data %>%
    select(-gene) %>%
    pivot_longer(-bin, names_to = "comparison", values_to = "tcn_diff")
  
  # For paired analysis, we test if the tcn_diff values are significantly different from 0
  # using Wilcoxon signed-rank test (non-parametric paired test)
  volcano_results <- heatmap_long %>%
    group_by(bin) %>%
    summarise(
      mean_tcn_diff = mean(tcn_diff, na.rm = TRUE),
      median_tcn_diff = median(tcn_diff, na.rm = TRUE),
      sd_tcn_diff = sd(tcn_diff, na.rm = TRUE),
      n_pairs = sum(!is.na(tcn_diff)),
      # Wilcoxon signed-rank test for paired data
      p_value = if(n_pairs >= 3) {
        tryCatch({
          # Test if paired differences are significantly different from 0
          wilcox.test(tcn_diff, mu = 0, alternative = "two.sided", 
                      exact = FALSE, paired = FALSE)$p.value
        }, error = function(e) {
          NA_real_
        })
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "fdr"),
      neg_log10_p = -log10(p_value),
      neg_log10_p_adj = -log10(p_adj),
      significance = case_when(
        is.na(p_adj) ~ "Not tested",
        p_adj < p_threshold & abs(mean_tcn_diff) > effect_threshold ~ "Significant",
        p_adj < p_threshold ~ "P < 0.05 only",
        abs(mean_tcn_diff) > effect_threshold ~ "Effect size only",
        TRUE ~ "Not significant"
      ),
      chromosome = str_extract(bin, "^[^:]+"),
      position = as.numeric(str_extract(bin, "(?<=:)[0-9]+"))
    ) %>%
    filter(!is.na(p_value), is.finite(neg_log10_p))
  
  # Add gene annotations
  gene_annotations <- annotated_heatmap_data %>% 
    select(bin, gene) %>% 
    distinct()
  
  volcano_results <- volcano_results %>%
    left_join(gene_annotations, by = "bin") %>%
    mutate(gene = ifelse(is.na(gene), "", gene))
  
  # Create volcano plot
  volcano_plot <- ggplot(volcano_results, aes(x = mean_tcn_diff, y = neg_log10_p_adj)) +
    geom_point(aes(color = significance), alpha = 0.7, size = 1.5) +
    scale_color_manual(
      values = c(
        "Significant" = "red",
        "P < 0.05 only" = "orange", 
        "Effect size only" = "blue",
        "Not significant" = "gray70",
        "Not tested" = "gray90"
      )
    ) +
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "red", alpha = 0.7) +
    geom_vline(xintercept = c(-effect_threshold, effect_threshold), linetype = "dashed", color = "blue", alpha = 0.7) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      text = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold"),
      legend.position = "right"
    ) +
    labs(
      x = "Mean TCN Difference",
      y = "-log10(Adjusted P-value)",
      title = "Volcano Plot: Paired TCN Differences",
      subtitle = paste0("Wilcoxon signed-rank test (paired), FDR correction, n=", nrow(volcano_results), " bins"),
      color = "Significance"
    ) +
    geom_label_repel(
      data = subset(volcano_results, significance == "Significant" & gene != ""),
      aes(label = gene),
      size = 3,
      max.overlaps = 10,
      point.padding = unit(0.5, "lines")
    )
  ggsave(file.path(out_dir, "volcano_plot.png"), volcano_plot, width = 12, height = 8, dpi = 300)
  
  # Save volcano results
  write_tsv(volcano_results, file.path(out_dir, "volcano_results.tsv"))
  
  # Create Manhattan plot by chromosome with proper numerical sorting
  manhattan_data <- volcano_results %>%
    mutate(chr_num = as.numeric(chromosome)) %>%
    filter(!is.na(chr_num)) %>%  # Remove any non-numeric chromosomes
    arrange(chr_num, position)
  
  # Calculate cumulative positions for Manhattan plot
  chr_lengths <- manhattan_data %>%
    group_by(chr_num) %>%
    summarise(
      max_pos = max(position, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(chr_num) %>%
    mutate(
      chr_start = cumsum(c(0, head(max_pos, -1))),
      chr_center = chr_start + max_pos / 2
    )
  
  # Add cumulative positions to manhattan_data
  manhattan_data <- manhattan_data %>%
    left_join(chr_lengths %>% select(chr_num, chr_start), by = "chr_num") %>%
    mutate(plot_pos = chr_start + position)
  
  # Create breaks and labels for chromosomes 1-22 in numerical order
  chr_breaks <- chr_lengths$chr_center
  chr_labels <- chr_lengths$chr_num
  
  manhattan_plot <- ggplot(manhattan_data, aes(x = plot_pos, y = neg_log10_p_adj)) +
    geom_point(aes(color = significance), alpha = 0.7, size = 1) +
    scale_color_manual(
      values = c(
        "Significant" = "red",
        "P < 0.05 only" = "orange",
        "Effect size only" = "blue",
        "Not significant" = "gray70",
        "Not tested" = "gray90"
      )
    ) +
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "red") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      text = element_text(size = 14),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      x = "Chromosome",
      y = "-log10(Adjusted P-value)",
      title = "Manhattan Plot: Significance Across Chromosomes"
    ) +
    scale_x_continuous(
      breaks = chr_breaks,
      labels = chr_labels
    ) +
    geom_label_repel(
      data = subset(manhattan_data, significance == "Significant" & gene != ""),
      aes(label = gene),
      size = 3,
      max.overlaps = 10,
      point.padding = unit(0.5, "lines")
    )
  ggsave(file.path(out_dir, "manhattan_plot.png"), manhattan_plot, width = 16, height = 6, dpi = 300)
  
  return(volcano_results)
}

#' Function to analyze paired cohorts
#' @param cncf_path Path to CNCF file
#' @param biomart_gr_path Path to biomart annotations
#' @param Cohort_1 Vector of sample IDs for first cohort
#' @param Cohort_2 Vector of sample IDs for second cohort (must be same length as Cohort_1)
#' @param out_dir Output directory
#' @param oncokb_cnas_path Path to OncoKB CNA data
#' @param bin_size Size of bins in base pairs (default: 1e6)
#' @param tcn_diff_threshold Threshold for critical regions (default: 1)
#' @param exact_match Logical, if TRUE uses exact matching, if FALSE uses grep matching (default: TRUE)
#' @param prefer_hisens Logical, if TRUE and multiple matches found, prefer samples with "hisens" in name (default: TRUE)
#' @return NULL

generate_paired_cohort_analysis <- function(cncf_path, biomart_gr_path, Cohort_1, Cohort_2,
                                            out_dir, oncokb_cnas_path, bin_size = 1e6, tcn_diff_threshold = 1,
                                            exact_match = TRUE, prefer_hisens = TRUE) {
  
  # Check that input vectors are of equal length
  if (length(Cohort_1) != length(Cohort_2)) {
    stop("Cohort_1 and Cohort_2 must be of equal length for paired analysis")
  }
  
  # Initialize analysis
  init_data <- initialize_analysis(cncf_path, biomart_gr_path)
  ChrLengths <- init_data$ChrLengths
  cncf <- init_data$cncf
  biomart_gr <- init_data$biomart_gr
  
  # Create output directory
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Helper function to find matching sample IDs
  find_sample_id <- function(pattern, available_ids, exact = TRUE, prefer_hisens = TRUE) {
    if (exact) {
      # Exact matching
      matches <- available_ids[available_ids == pattern]
    } else {
      # Grep matching
      matches <- available_ids[grepl(pattern, available_ids, fixed = TRUE)]
    }
    
    if (length(matches) == 0) {
      return(NULL)
    } else if (length(matches) == 1) {
      return(matches[1])
    } else {
      # Multiple matches found
      cat("Multiple matches found for '", pattern, "':\n")
      for (i in seq_along(matches)) {
        cat("  ", i, ": ", matches[i], "\n")
      }
      
      if (prefer_hisens) {
        hisens_matches <- matches[grepl("hisens", matches)]
        if (length(hisens_matches) == 1) {
          cat("  Selected hisens match: ", hisens_matches[1], "\n")
          return(hisens_matches[1])
        } else if (length(hisens_matches) > 1) {
          cat("  Multiple hisens matches found, using first one: ", hisens_matches[1], "\n")
          return(hisens_matches[1])
        }
      }
      
      # If no hisens preference or no hisens matches, use first match
      cat("  Using first match: ", matches[1], "\n")
      return(matches[1])
    }
  }
  
  # Get all available sample IDs
  available_ids <- unique(cncf$ID)
  
  # Resolve sample IDs and identify valid pairs
  cat("Resolving sample IDs...\n")
  valid_pairs <- c()
  sample_mapping_list <- list()
  
  for (i in seq_along(Cohort_1)) {
    # Resolve Cohort_1 sample
    resolved_1 <- find_sample_id(Cohort_1[i], available_ids, exact_match, prefer_hisens)
    
    # Resolve Cohort_2 sample
    resolved_2 <- find_sample_id(Cohort_2[i], available_ids, exact_match, prefer_hisens)
    
    # Check if both samples were found
    if (is.null(resolved_1) && is.null(resolved_2)) {
      warning("Pair ", i, ": Neither sample found - ", Cohort_1[i], " and ", Cohort_2[i], ". Skipping this pair.")
    } else if (is.null(resolved_1)) {
      warning("Pair ", i, ": Cohort_1 sample not found - ", Cohort_1[i], ". Skipping this pair.")
    } else if (is.null(resolved_2)) {
      warning("Pair ", i, ": Cohort_2 sample not found - ", Cohort_2[i], ". Skipping this pair.")
    } else {
      # Both samples found - add to valid pairs
      valid_pairs <- c(valid_pairs, i)
      sample_mapping_list[[length(sample_mapping_list) + 1]] <- tibble(
        original_pair = i,
        input_cohort_1 = Cohort_1[i],
        resolved_cohort_1 = resolved_1,
        input_cohort_2 = Cohort_2[i],
        resolved_cohort_2 = resolved_2
      )
      cat("Pair", i, ": ", Cohort_1[i], " -> ", resolved_1, " vs ", Cohort_2[i], " -> ", resolved_2, "\n")
    }
  }
  
  # Check if we have any valid pairs
  if (length(valid_pairs) == 0) {
    stop("No valid pairs found after sample ID resolution. Cannot proceed with analysis.")
  }
  
  cat("Found", length(valid_pairs), "valid pairs out of", length(Cohort_1), "input pairs.\n")
  
  # Create sample mapping table
  sample_mapping <- bind_rows(sample_mapping_list) %>%
    mutate(final_pair = row_number())
  
  write_tsv(sample_mapping, file.path(out_dir, "sample_id_mapping.tsv"))
  cat("Sample ID mapping saved to: sample_id_mapping.tsv\n")
  
  # Update cohorts to only include valid pairs
  valid_cohort_1 <- sapply(sample_mapping$resolved_cohort_1, function(x) x)
  valid_cohort_2 <- sapply(sample_mapping$resolved_cohort_2, function(x) x)
  
  # Helper function for baseline data
  getBaseline <- function(data) {
    data %>%
      as.data.frame() %>%
      rownames_to_column("Chromosome") %>%
      rename(. = 'end') %>%
      as_tibble() %>%
      mutate(seqnames = Chromosome, start = 1, tcn = 2)
  }
  
  # Store all differential segments and diffPercents
  all_diff_segments <- list()
  diff_percents <- tibble(comparison = character(), diffPercent = numeric(),
                          sample1 = character(), sample2 = character())
  pair_count <- 0
  
  # Process each valid pair of samples
  for (i in seq_along(valid_cohort_1)) {
    sample1_id <- valid_cohort_1[i]
    sample2_id <- valid_cohort_2[i]
    comparison_name <- paste0(sample1_id, "_vs_", sample2_id)
    cat(paste("Processing comparison:", comparison_name, "\n"))
    
    # Create pairwise plots directory (shared for all pairs)
    pairwise_plots_dir <- file.path(out_dir, "pairwise_plots")
    dir.create(pairwise_plots_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Extract sample data
    sample1 <- cncf %>%
      select(chrom, loc.start, loc.end, tcn, ID) %>%
      filter(ID == sample1_id)
    
    sample2 <- cncf %>%
      select(chrom, loc.start, loc.end, tcn, ID) %>%
      filter(ID == sample2_id)
    
    # Check if we found both samples (this should always pass now, but keeping as safety check)
    if (nrow(sample1) == 0 || nrow(sample2) == 0) {
      message("One or both samples not found in the data: ", sample1_id, ", ", sample2_id)
      next
    }
    
    # Convert to GRanges objects
    current_sample1_gr <- sample1 %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    current_sample2_gr <- sample2 %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    
    # Find overlapping segments
    overlapping_segments <- current_sample1_gr %>%
      join_overlap_intersect(current_sample2_gr) %>%
      mutate(segment_length = width(.))
    
    # Calculate total size of overlapping segments
    size <- sum(overlapping_segments$segment_length)
    
    # Skip if no overlapping segments
    if (size == 0) {
      message("No overlapping segments found between ", sample1_id, " and ", sample2_id)
      next
    }
    
    # Round TCN values
    overlapping_segments$tcn.x <- round(overlapping_segments$tcn.x)
    overlapping_segments$tcn.y <- round(overlapping_segments$tcn.y)
    
    # Separate segments with same and different TCN values
    sameSegs <- overlapping_segments %>%
      as_tibble() %>%
      group_by(seqnames, start, end) %>%
      filter(tcn.x == tcn.y)
    
    diffSegs <- overlapping_segments %>%
      as_tibble() %>%
      group_by(seqnames, start, end) %>%
      filter(tcn.x != tcn.y)
    
    # Calculate percentage of different segments
    diffPercent <- round((sum(diffSegs$segment_length) / size) * 100)
    
    # Store diffPercent
    diff_percents <- diff_percents %>%
      add_row(comparison = comparison_name, diffPercent = diffPercent,
              sample1 = sample1_id, sample2 = sample2_id)
    
    # Get baseline for plotting
    baseline <- ChrLengths %>% getBaseline()
    
    # Ensure seqnames are correctly defined as factors
    baseline <- baseline %>%
      mutate(seqnames = factor(seqnames, levels = as.character(1:22)))
    
    diffSegsToPlot <- diffSegs %>%
      mutate(seqnames = factor(seqnames, levels = as.character(1:22)))
    
    sameSegsToPlot <- sameSegs %>%
      mutate(seqnames = factor(seqnames, levels = as.character(1:22)))
    
    # Store diff segments for later use if there are any
    if (nrow(diffSegsToPlot) > 0) {
      Patient_id <- str_extract(sample1_id, "^P-[0-9]+")
      diffSegsToPlot$comparison <- comparison_name
      diffSegsToPlot$Patient_id <- Patient_id
      all_diff_segments[[comparison_name]] <- diffSegsToPlot
      pair_count <- pair_count + 1
    }
    
    # Create labels for plot
    label1 <- paste0("TCN unique to ", sample1_id)
    label2 <- paste0("TCN unique to ", sample2_id)
    
    # Define the color mapping
    color_values <- setNames(
      c("blue", "red", "green"),
      c(label1, label2, "Shared segments")
    )
    
    # Generate the plot
    p <- ggplot(baseline, aes(x = (start + end) / 2, y = tcn)) +
      facet_grid(~seqnames, scales = "free_x", space = "free_x", switch = "x") +
      theme_minimal() +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_text(size = 15, vjust = -0.5),
        axis.text.x = element_blank(),
        text = element_text(size = 15),
        strip.text = element_text(size = 10, face = "bold"),
        plot.background = element_rect(fill = "white", color = NA),
        plot.title = element_textbox_simple(size = 16, face = "bold")
      ) +
      scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 1)) +
      geom_segment(data = diffSegsToPlot, aes(x = start, xend = end, y = tcn.x, yend = tcn.x, color = label1), linewidth = 2) +
      geom_segment(data = diffSegsToPlot, aes(x = start, xend = end, y = tcn.y, yend = tcn.y, color = label2), linewidth = 2) +
      geom_segment(data = sameSegsToPlot, aes(x = start, xend = end, y = tcn.x, yend = tcn.x, color = "Shared segments"), linewidth = 2) +
      coord_cartesian(ylim = c(0, 10)) +
      labs(x = "Chromosome", y = "Copy Number", color = "") +
      scale_color_manual(values = color_values) +
      theme(legend.position = "right") +
      ggtitle(glue::glue("**{diffPercent}%** of total copy number difference"))
    
    # Save the plot
    output_plot_path <- file.path(pairwise_plots_dir, paste0(comparison_name, "_comparison.png"))
    ggsave(output_plot_path, p, width = 22, height = 6)
    # Save differential segments data
    diffSegsToPlot$diffPercent <- diffPercent
    diffSegsToPlot$sample1 <- sample1_id
    diffSegsToPlot$sample2 <- sample2_id
    output_diff_segments_path <- file.path(pairwise_plots_dir, paste0(comparison_name, "_diffSegments.tsv"))
    write_tsv(diffSegsToPlot, output_diff_segments_path)
    
    # Annotate genes in differential segments
    diff_segments <- GRanges(
      seqnames = diffSegsToPlot$seqnames,
      ranges = IRanges(start = diffSegsToPlot$start, end = diffSegsToPlot$end),
      strand = if("strand" %in% colnames(diffSegsToPlot)) diffSegsToPlot$strand else "*"
    )
    
    # Check if annotation file exists and has data
    has_annot_file <- !is.null(biomart_gr) && class(biomart_gr)[1] == "GRanges" && length(biomart_gr) > 0
    
    # Find overlaps and collect gene symbols for each segment
    if (has_annot_file && length(diff_segments) > 0) {
      tryCatch({
        annotated_genes <- diff_segments %>%
          join_overlap_left(biomart_gr) %>%
          as_tibble() %>%
          select(seqnames, start, end, hgnc_symbol, ensembl_gene_id)
        
        # Add sample info to annotated genes
        annotated_genes$sample1 <- sample1_id
        annotated_genes$sample2 <- sample2_id
        
        # Save annotated genes data
        output_gene_path <- file.path(pairwise_plots_dir, paste0(comparison_name, "_diffGenes.tsv"))
        write_tsv(annotated_genes, output_gene_path)
      }, error = function(e) {
        message("Error annotating genes: ", e$message)
      })
    }
  }
  # Save diffPercents to TSV
  write_tsv(diff_percents, file.path(out_dir, "diffPercents.tsv"))
  
  # Create boxplot of diffPercents
  if (nrow(diff_percents) > 0) {
    # Extract patient IDs for grouping
    diff_percents <- diff_percents %>%
      mutate(
        Patient_id1 = str_extract(sample1, "^P-[0-9]+"),
        Patient_id2 = str_extract(sample2, "^P-[0-9]+")
      ) %>% 
      mutate(Comparison =    str_replace(comparison, "^P-[0-9]+-T[^_]*_(.*?)_vs_P-[0-9]+-T[^_]*_(.*)", "\\1_vs_\\2"))
    
    mean_diff <- mean(diff_percents$diffPercent, na.rm = TRUE)
    
    # Create simple boxplot data
    boxplot_data <- diff_percents %>%
      mutate(group = "All Pairs")  # Create a single group for all pairs

    # # Create boxplot
    # boxplot_p <- ggplot(boxplot_data, aes(x = group, y = diffPercent)) +
    #   geom_boxplot(fill = "black", alpha = 0.7) +
    #   geom_jitter(color = "red", width = 0.2, alpha = 0.6, size = 2) +
    #   theme_minimal() +
    #   labs(
    #     title = "Distribution of % of TCN difference between pairs",
    #     subtitle = paste0("Mean difference: ", round(mean_diff, 2), "%"),
    #     x = "Pairwise comparisons",
    #     y = "TCN Difference Percentage (%)"
    #   ) +
    #   theme(plot.background = element_rect(fill = "white", color = NA)) +
    #   geom_label_repel(
    #     data = boxplot_data,
    #     aes(x = group, y = diffPercent, label = Comparison),
    #     size = 2,
    #     max.overlaps = 30,
    #     force = 20,
    #     direction = "y",
    #     nudge_y = 0.5,
    #     segment.color = "gray50",      # Color of the connecting line
    #     segment.size = 0.5,            # Thickness of the connecting line
    #     segment.alpha = 0.7,           # Transparency of the connecting line
    #     segment.linetype = "dashed"    # Line type (solid, dashed, dotted, etc.)
    #   )
    # ggsave(file.path(out_dir, "diffPercent_boxplot.png"), boxplot_p, width = 6, height = 10)
    
    # # Create a lollipop plot instead
    # Calculate median for vertical line
    median_diff <- median(diff_percents$diffPercent, na.rm = TRUE)
    
    # Add cluster information to boxplot_data
    # First, we need to get cluster info from the heatmap data if it exists
    cluster_info <- NULL
    heatmap_data_path <- file.path(out_dir, "heatmap_data.tsv")
    cluster_mapping_path <- file.path(out_dir, "cluster_mapping.tsv")
    
    # Try to read cluster mapping if it exists (will be created after heatmap generation)
    if (file.exists(cluster_mapping_path)) {
      cluster_info <- read_tsv(cluster_mapping_path, show_col_types = FALSE)
      
      # Match cluster info to boxplot_data
      boxplot_data <- boxplot_data %>%
        left_join(cluster_info %>% select(comparison, cluster), 
                  by = c("comparison" = "comparison")) %>%
        mutate(
          cluster_label = ifelse(is.na(cluster), "No Cluster", paste0("Cluster ", cluster)),
          Comparison_with_cluster = paste0(Comparison, " (", cluster_label, ")")
        )
    } else {
      # If no cluster info available, use regular labels
      boxplot_data <- boxplot_data %>%
        mutate(
          cluster_label = "No Cluster",
          Comparison_with_cluster = Comparison
        )
    }
    
    # Create a lollipop plot with median line and cluster labels
    lollipop_p <- ggplot(boxplot_data, aes(x = reorder(Comparison_with_cluster, diffPercent), y = diffPercent)) +
      geom_segment(aes(x = Comparison_with_cluster, xend = Comparison_with_cluster, y = 0, yend = diffPercent),
                   color = "gray50", size = 1) +
      geom_point(color = "red", size = 4) +
      # Add vertical line at median
      geom_hline(yintercept = median_diff, linetype = "dashed", color = "blue", size = 1, alpha = 0.7) +
      # Add median label
      annotate("text", x = length(unique(boxplot_data$Comparison_with_cluster)) * 0.1, 
               y = median_diff + max(boxplot_data$diffPercent) * 0.05, 
               label = paste0("Median: ", round(median_diff, 1), "%"), 
               color = "blue", size = 4, hjust = 0) +
      coord_flip() +
      theme_minimal() +
      labs(
        title = "TCN Difference Percentage by Comparison",
        subtitle = paste0("Mean: ", round(mean_diff, 2), "%, Median: ", round(median_diff, 2), "%"),
        x = "Pairwise Comparisons (with Cluster)",
        y = "TCN Difference Percentage (%)"
      ) +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        axis.text.y = element_text(size = 10),                    # Reduced font size for row names
        axis.text.x = element_text(size = 14),                    # X-axis text size
        axis.title.x = element_text(size = 18),                   # X-axis title
        axis.title.y = element_text(size = 18),                   # Y-axis title
        plot.title = element_text(size = 22, face = "bold"),      # Title size
        plot.subtitle = element_text(size = 20),                  # Increased subtitle size from 16 to 20
        text = element_text(size = 16),                           # Base text size
        panel.grid.major = element_line(color = "gray90"),        # Optional: lighter grid lines
        panel.grid.minor = element_blank()                        # Optional: remove minor grid lines
      )
    
    ggsave(file.path(out_dir, "diffPercent_lollipop.png"), lollipop_p, width = 20, height = 10, dpi = 300)
  }
  
  # Check if we processed any pairs
  if (pair_count == 0) {
    stop("No valid pairs were found. Cannot generate comparisons.")
  }
  
  # Check if we have any diffSegments data
  if (length(all_diff_segments) == 0) {
    stop("No differential segments were found for any comparisons.")
  }
  
  cat("Generating heatmap...\n")
  
  # Generate the heatmap (rest of the heatmap code remains the same but without sampling)
  # Combine all differential segments
  combined_diff_segments <- bind_rows(all_diff_segments)
  
  # Prepare Chromosome bins
  chrom_bins_df <- create_bins(ChrLengths, bin_size)
  
  # Convert combined_diff_segments to GRanges once
  diff_segments_gr <- GRanges(
    seqnames = combined_diff_segments$seqnames,
    ranges = IRanges(start = combined_diff_segments$start, end = combined_diff_segments$end),
    tcn_diff = combined_diff_segments$tcn.x - combined_diff_segments$tcn.y,
    comparison = combined_diff_segments$comparison,
    Patient_id = combined_diff_segments$Patient_id
    
  )
  
  # Convert chrom_bins to GRanges once
  bins_gr <- GRanges(
    seqnames = chrom_bins_df$seqnames,
    ranges = IRanges(start = chrom_bins_df$start, end = chrom_bins_df$end),
    bin_id = chrom_bins_df$bin_id
  )
  
  # Find all overlaps at once
  all_overlaps <- findOverlaps(diff_segments_gr, bins_gr)
  
  # Create a data frame of all overlaps with tcn_diff
  overlap_df <- tibble(
    segment_idx = queryHits(all_overlaps),
    bin_idx = subjectHits(all_overlaps),
    tcn_diff = diff_segments_gr$tcn_diff[queryHits(all_overlaps)],
    comparison = diff_segments_gr$comparison[queryHits(all_overlaps)],
    Patient_id = diff_segments_gr$Patient_id[queryHits(all_overlaps)],
    bin_id = bins_gr$bin_id[subjectHits(all_overlaps)]
  )
  
  # If multiple segments overlap the same bin for the same comparison, use the mean
  aggregated_overlaps <- overlap_df %>%
    group_by(comparison, bin_id) %>%
    summarize(tcn_diff = mean(tcn_diff), .groups = "drop")
  
  # Transform to wide format for heatmap
  heatmap_wide <- aggregated_overlaps %>%
    pivot_wider(names_from = comparison, values_from = tcn_diff, values_fill = 0)
  
  # Extract bin IDs
  bin_ids <- heatmap_wide$bin_id
  heatmap_matrix <- as.matrix(heatmap_wide[, -1])
  rownames(heatmap_matrix) <- bin_ids
  
  # NO SAMPLING - plot everything as requested
  cat("Creating heatmap with", nrow(heatmap_matrix), "bins and", ncol(heatmap_matrix), "comparisons\n")
  
  # Create Patient_id lookup table
  patient_lookup <- combined_diff_segments %>%
    select(comparison, Patient_id) %>%
    distinct()
  
  # Convert to long format for tidyheatmap
  heatmap_long <- heatmap_matrix %>%
    as.data.frame() %>%
    rownames_to_column("bin") %>%
    pivot_longer(-bin, names_to = "comparison", values_to = "tcn_diff") %>%
    # Create Patient_id directly from comparison name
    mutate(
      Patient_id = str_extract(comparison, "^P-[0-9]+"),
      pair = as.numeric(factor(comparison))  # Keep pair for backward compatibility
    ) %>%
    # Extract Chromosome and position information from bin ID
    mutate(
      Chromosome = str_extract(bin, "^[^:]+"),
      position = as.numeric(str_extract(bin, "(?<=:)[0-9]+"))
    ) %>%
    # Convert Chromosome to numeric for sorting
    mutate(
      chr_num = as.numeric(Chromosome),
      chr_parity = ifelse(chr_num %% 2 == 0, "even", "odd")
    ) %>%
    arrange(chr_num, position) %>%
    mutate(
      bin = factor(bin, levels = unique(bin)),
      Chromosome = factor(Chromosome, levels = as.character(1:22))
    )
  
  # Generate heatmap
  generate_heatmap(heatmap_long, out_dir)

  # Generate the avg_tcn plot (rest remains the same)
  cat("Generating average TCN plot...\n")
  
  # Aggregate the data to calculate the mean tcn_diff per bin
  avg_tcn_data <- heatmap_long %>%
    group_by(bin) %>%
    summarize(avg_tcn_diff = mean(tcn_diff, na.rm = TRUE), .groups = "drop")
  
  # Ensure bins are sorted naturally
  avg_tcn_data <- avg_tcn_data %>%
    mutate(
      Chromosome = str_extract(bin, "^[^:]+"),
      position = as.numeric(str_extract(bin, "(?<=:)[0-9]+")),
      chr_num = as.numeric(Chromosome)
    ) %>%
    arrange(chr_num, position) %>%
    mutate(bin = factor(bin, levels = unique(bin)),
           CNA_type = ifelse(avg_tcn_diff > 0, "Enrichment in 1",
                             ifelse(avg_tcn_diff < 0, "Enrichment in 2", "Neutral")))
  
  # Filter avg_tcn_data for regions where avg_tcn_diff crosses the threshold
  critical_regions <- avg_tcn_data %>%
    filter(avg_tcn_diff < -tcn_diff_threshold | avg_tcn_diff > tcn_diff_threshold)
  
  # Convert avg_tcn_data to GRanges for gene matching
  avg_tc_gr <- GRanges(
    seqnames = critical_regions$Chromosome,
    ranges = IRanges(start = critical_regions$position, width = bin_size),
    avg_tcn_diff = critical_regions$avg_tcn_diff
  )
  
  # Find overlaps with genes in biomart_gr on critical regions
  overlapping_genes <- findOverlaps(avg_tc_gr, biomart_gr)
  
  # Create an empty tibble if there are no overlaps
  if (length(overlapping_genes) == 0) {
    annotated_genes <- tibble(bin = character(0), gene = character(0))
  } else {
    annotated_genes <- tibble(
      bin = critical_regions$bin[queryHits(overlapping_genes)],
      gene = biomart_gr$hgnc_symbol[subjectHits(overlapping_genes)]
    )
  }
  
  # Read OncoKB data
  oncokb_cnas <- load_oncokb_data(oncokb_cnas_path)
  
  # Match oncogenic genes with critical regions
  if (nrow(annotated_genes) > 0 && nrow(oncokb_cnas) > 0) {
    # First filter for genes that are not NA and not empty
    filtered_annotated <- critical_regions %>%
      left_join(annotated_genes, by = "bin") %>%
      filter(!is.na(gene) & gene != "")
    
    # Then join with oncokb data
    if (nrow(filtered_annotated) > 0) {
      plotgenes <- filtered_annotated %>%
        inner_join(oncokb_cnas, by = "gene", relationship = "many-to-many") %>%
        distinct()
    } else {
      plotgenes <- tibble()
    }
  } else {
    plotgenes <- tibble()
  }
  
  # Add gene annotations to avg_tcn_data
  avg_tcn_data <- avg_tcn_data %>%
    left_join(plotgenes %>% select(bin, gene), by = "bin")
  
  # Merge genes within the same bin, handling NA values safely
  merged_avg_tcn_data <- avg_tcn_data %>%
    group_by(bin, avg_tcn_diff, Chromosome, position, chr_num, CNA_type) %>%
    summarise(
      gene = if (all(is.na(gene))) {
        NA_character_
      } else {
        paste(unique(na.omit(gene)), collapse = "\n")
      },
      .groups = "drop"
    ) %>%
    # Make empty strings NA
    mutate(gene = na_if(gene, ""))
  
  # Plot the avg_tc_data using ggplot2
  avg_tcn_plot <- ggplot(merged_avg_tcn_data, aes(x = bin, y = avg_tcn_diff, fill = CNA_type)) +
    geom_col(width = 1) +
    scale_fill_manual(
      values = c("Enrichment in 1" = "#FF7F7F", "Enrichment in 2" = "lightblue"),
      labels = c("Enrichment in 1" = "Average Enrichment in 1", "Enrichment in 2" = "Average Enrichment in 2")
    ) +
    facet_wrap(~ chr_num, scales = "free_x", nrow = 1, strip.position = "bottom") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.x = element_text(size = 20, vjust = -0.5),        # Increased from 15
      axis.title.y = element_text(size = 20),                      # Added explicit y-axis title size
      axis.text.y = element_text(size = 16),                       # Y-axis tick labels
      plot.title = element_text(size = 24, face = "bold"),         # Title size
      plot.subtitle = element_text(size = 18),                     # Subtitle size
      legend.text = element_text(size = 16),                       # Legend text
      legend.title = element_text(size = 18),                      # Legend title
      strip.text = element_text(size = 16, face = "bold"),         # Chromosome facet labels
      text = element_text(size = 18)                               # Base text size increased from 15
    ) +
    labs(fontsize=12,
      x = "Chromosome",
      y = "Average TCN Difference",
      title = "Average TCN Difference Between Paired Samples",
      subtitle = paste0("Comparing ", length(Cohort_1), " sample pairs with gene label threshold ", tcn_diff_threshold)
    )
  
  # Add gene labels if there are any non-NA genes
  gene_data <- merged_avg_tcn_data %>%
    filter(!is.na(gene) & gene != "")
  if (nrow(gene_data) > 0) {
    avg_tcn_plot <- avg_tcn_plot +
      geom_label_repel(
        data = gene_data,
        aes(x = bin, y = avg_tcn_diff, label = gene),
        size = 6,
        max.overlaps = Inf,
        force = 20,
        direction = "y",
        nudge_y = 0.5
      )
  }
  
  # Add points for genes if any
  point_data <- avg_tcn_data %>%
    filter(!is.na(gene))
  if (nrow(point_data) > 0) {
    avg_tcn_plot <- avg_tcn_plot +
      geom_point(data = point_data, aes(x = bin, y = avg_tcn_diff))
  }
  
  # Save the plot
  ggsave(
    file.path(out_dir, "average_tcn_diff_barplot.png"),
    plot = avg_tcn_plot,
    width = 40, height = 15, dpi = 300
  )
  # Generate volcano plot with CORRECTED paired statistics
  cat("Generating volcano plot with paired statistics...\n")
  volcano_results <- generate_volcano_plot_paired(
    heatmap_data_path = file.path(out_dir, "heatmap_data.tsv"),
    out_dir = out_dir,
    biomart_gr = biomart_gr,
    p_threshold = 0.05,
    effect_threshold = tcn_diff_threshold
  )
  
  # Print summary
  sig_bins <- sum(volcano_results$significance == "Significant", na.rm = TRUE)
  cat("Found", sig_bins, "significantly altered bins (FDR < 0.05 & |effect| > 1)\n")

  cat("Analysis completed successfully!\n")
  return(invisible(NULL))
}
