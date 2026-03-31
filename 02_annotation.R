# ============================================================================
# Script: STR Annotation and Feature Classification
# Description:
#   1. Load VEP annotation results (tab-separated) for RNA-derived pSTRs.
#   2. Load annotation priority table (Excel) defining group order and mapping.
#   3. For each STR locus, assign the highest-priority functional group based on
#      overlapping feature types (e.g., CDS > 5' UTR > 3' UTR > intron > intergenic).
#   4. Extract biotype information when available.
#   5. Output a table with locus, final group, gene ID, biotype, and original annotations.
# ============================================================================

rm(list = ls())

# ----------------------------- Libraries ------------------------------------
library(tidyr)
library(readxl)

# ----------------------------- Input files ----------------------------------
# Annotation priority table: columns: group (e.g., CDS), annotation (VEP terms), order (numeric)
annotation_order <- read_excel("annotation_order.xlsx")

# VEP annotation output for pSTRs (from RNA-seq)
raw_annotation <- read.table("pstr.rna.recode.annotated.vcf.gz", sep = "\t")

# ----------------------------- Helper logic --------------------------------
# For each locus, gather all VEP consequence terms, map to groups using annotation_order,
# then pick the group with the smallest order number (highest priority).
# 
# Note on special cases mentioned in original comments:
# - synonymous_variant within CDS is still considered CDS (no special exclusion here,
#   but if needed, the priority table can assign lower order to synonymous).
# - Noncoding_exon/intron and upstream/downstream are all categorized as intergenic
#   (ensured by annotation_order mapping).

feature_df <- data.frame()
total_loci <- length(unique(raw_annotation[, 1]))
progress <- 0

for (locus in unique(raw_annotation[, 1])) {
  progress <- progress + 1
  if (progress %% 100 == 0) {
    cat(sprintf("%.2f%% finished\n", progress / total_loci * 100))
  }
  
  # Get all VEP consequence terms for this locus (column 7, comma-separated)
  feature_type <- raw_annotation[raw_annotation[, 1] == locus, 7]
  feature_type <- unique(unlist(strsplit(feature_type, ",")))
  
  # Map each term to its group and order using annotation_order
  mapped <- annotation_order[annotation_order$annotation %in% feature_type, ]
  if (nrow(mapped) == 0) {
    # No matching annotation -> assign to "intergenic" or default?
    # Original code would have no rows; we set final_group = NA? Better to assign a default.
    # However, original code assumed at least one match. To be safe, we assign "intergenic" with high order.
    final_group <- "intergenic"
    final_order <- 999
    group_list <- "intergenic"
    order_list <- "999"
  } else {
    group_list <- unique(mapped$group)
    order_list <- unique(mapped$order)
    # Select the group with the smallest order (highest priority)
    best_idx <- which.min(mapped$order)
    final_group <- mapped$group[best_idx]
    final_order <- mapped$order[best_idx]
  }
  
  # Extract gene ID (column 5) and biotype (column 14, format "BIOTYPE=xxx")
  gene_id <- raw_annotation[raw_annotation[, 1] == locus, 5][1]  # take first if multiple
  extra_info <- raw_annotation[raw_annotation[, 1] == locus, 14][1]
  biotype <- NA
  if (grepl("BIOTYPE", extra_info)) {
    parts <- unlist(strsplit(extra_info, ";"))
    biotype_part <- parts[grepl("BIOTYPE", parts)]
    if (length(biotype_part) > 0) {
      biotype <- unlist(strsplit(biotype_part, "="))[2]
    }
  }
  
  # Store results
  temp <- data.frame(
    locus = locus,
    final_group = final_group,
    final_order = final_order,
    feature_type = paste(feature_type, collapse = ","),
    group_type = paste(if (exists("group_list")) group_list else "", collapse = ","),
    order_type = paste(if (exists("order_list")) order_list else "", collapse = ","),
    gene_id = gene_id,
    biotype = biotype,
    stringsAsFactors = FALSE
  )
  feature_df <- rbind(feature_df, temp)
}

# ----------------------------- Output ---------------------------------------
write.table(feature_df, file = "STR_annotation_group.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# End of script