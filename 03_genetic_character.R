# ============================================================================
# Script: Genetic characterization of polymorphic STRs (pSTRs)
# Description:
#   1. Load pSTR VCF, locus list, annotation group, and configuration.
#   2. Generate genetic characteristics data frame (period size, copy number, motif, etc.).
#   3. Analyze biotype distribution (bar plot).
#   4. Analyze pSTR composition across genomic features (stacked bar plots for total, mRNA, lncRNA).
#   5. Plot density of reference allele length and copy number by feature and period size.
#   6. Process motif patterns: collapse reverse complements and cyclic permutations, then visualize.
#   7. Generate supplementary plots (total density, pie charts for mRNA/lncRNA composition).
# ============================================================================

rm(list = ls())
library(readxl)
library(tidyr)
library(vcfR)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggsci)

# ----------------------------- Input data -----------------------------------
# pSTR VCF, locus list, and annotation group
pSTR_vcf <- read.vcfR("../pstr.rna.recode.vcf.gz")
pSTR_list <- read.table("../1_expression/locus_list.polymorphic.txt")
pSTR_group <- read.table("../2_annotation/STR_annotation_group.txt", header = TRUE)

# STR configuration (reference BED)
config <- read.table("../GRCh38.hipstr_reference.refine.bed")
config <- config[config[, 6] %in% pSTR_list[, 1], ]
config$reference_allele_length <- config[, 3] - config[, 2] + 1

# ----------------------------- (1) Genetic characteristics data frame ------
genetic_char <- data.frame(
  locus = pSTR_list[, 1],
  gene_id = pSTR_group$biotype[match(pSTR_list[, 1], pSTR_group$locus)],
  feature = pSTR_group$final_group[match(pSTR_list[, 1], pSTR_group$locus)],
  biotype = pSTR_group$biotype[match(pSTR_list[, 1], pSTR_group$locus)],
  period_size = config[match(pSTR_list[, 1], config[, 6]), 4],
  ncopy = config[match(pSTR_list[, 1], config[, 6]), 5],
  reference_allele_length = config$reference_allele_length[match(pSTR_list[, 1], config[, 6])],
  motif = config[match(pSTR_list[, 1], config[, 6]), 7]
)
write.table(genetic_char, file = "genetic_character_df.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# ----------------------------- (2) Biotype distribution ---------------------
biotype_vec <- genetic_char$biotype
biotype_vec[biotype_vec == "protein_coding"] <- "mRNA (protein_coding)"
biotype_vec[is.na(biotype_vec)] <- "NoAnnotation"
biotype_df <- as.data.frame(table(biotype_vec))
colnames(biotype_df) <- c("biotype", "counts")
biotype_df$percentage <- biotype_df$counts / sum(biotype_df$counts)
biotype_df$label <- paste0(biotype_df$counts, " (", round(biotype_df$percentage, 4) * 100, "%)")
biotype_df <- biotype_df[order(biotype_df$percentage, decreasing = TRUE), ]
biotype_df$biotype <- factor(biotype_df$biotype, levels = as.character(biotype_df$biotype))

biotype_barplot <- ggplot(biotype_df, aes(x = biotype, y = counts, fill = biotype)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = label), size = 3.5, vjust = -0.2) +
  theme_bw() +
  labs(title = "The Sourced Transcript Biotype of RNA-pSTRs",
       x = "Biotype", y = "pSTR Counts") +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 10, hjust = 1, angle = 45),
        axis.text.y = element_text(face = "bold", size = 10),
        legend.position = "none")
ggsave(biotype_barplot, file = "biotype_barplot.pdf", height = 6, width = 10)

# ----------------------------- (3) pSTR composition by feature ------------
# Helper: compute counts per period size per feature for given biotype
build_composition <- function(biotype_filter = "total") {
  if (biotype_filter == "total") {
    sub <- genetic_char
    features <- unique(sub$feature)
  } else {
    sub <- na.omit(genetic_char[genetic_char$biotype == biotype_filter, ])
    features <- unique(sub$feature)
  }
  
  counts_df <- data.frame()
  label_df <- data.frame()
  for (grp in features) {
    if (biotype_filter == "total") {
      ps <- sub[sub$feature == grp, "period_size"]
    } else {
      ps <- sub[sub$feature == grp, "period_size"]
    }
    tab <- as.data.frame(table(ps))
    tab$percentage <- tab$Freq / sum(tab$Freq)
    colnames(tab) <- c("period_size", "counts", "percentage")
    counts_df <- rbind(counts_df, data.frame(group = grp, tab))
    label_df <- rbind(label_df, data.frame(group = grp, group_sum = sum(tab$counts)))
  }
  label_df$percentage <- label_df$group_sum / sum(label_df$group_sum)
  label_df$label <- paste0(label_df$group_sum, " (", round(label_df$percentage, 4) * 100, "%)")
  return(list(counts = counts_df, labels = label_df))
}

# Total pSTRs
total_res <- build_composition("total")
total_counts <- total_res$counts
total_labels <- total_res$labels
total_counts$group <- factor(total_counts$group, levels = c("CDS", "5_UTR", "3_UTR", "intron", "intergenic"))
total_counts$period_size <- factor(total_counts$period_size, levels = 2:6,
                                   labels = c("Di", "Tri", "Tetra", "Penta", "Hexa"))
total_labels$group <- factor(total_labels$group, levels = c("CDS", "5_UTR", "3_UTR", "intron", "intergenic"))

pSTR_counts_barplot <- ggplot(total_counts, aes(x = group, y = percentage, fill = period_size)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  geom_text(data = total_labels, aes(x = group, y = 1, label = label), inherit.aes = FALSE, vjust = -0.2) +
  scale_fill_lancet() +
  scale_y_continuous(expand = c(0.01, 0.05)) +
  theme_bw() +
  labs(title = "The Composition of pSTRs", y = "Percentage", x = "Features") +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 10))
ggsave(pSTR_counts_barplot, file = "pSTR_counts_barplot.pdf", height = 6, width = 8)

# mRNA-derived pSTRs
mrna_res <- build_composition("protein_coding")
mrna_counts <- mrna_res$counts
mrna_labels <- mrna_res$labels
mrna_counts$group <- factor(mrna_counts$group, levels = c("CDS", "5_UTR", "3_UTR", "intron", "intergenic"))
mrna_counts$period_size <- factor(mrna_counts$period_size, levels = 2:6,
                                  labels = c("Di", "Tri", "Tetra", "Penta", "Hexa"))
mrna_labels$group <- factor(mrna_labels$group, levels = c("CDS", "5_UTR", "3_UTR", "intron", "intergenic"))

pSTR_counts_barplot_mRNA <- ggplot(mrna_counts, aes(x = group, y = percentage, fill = period_size)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  geom_text(data = mrna_labels, aes(x = group, y = 1, label = label), inherit.aes = FALSE, vjust = -0.2) +
  scale_fill_lancet() +
  scale_y_continuous(expand = c(0.01, 0.05)) +
  theme_bw() +
  labs(title = "The Composition of pSTRs Sourced from mRNA", y = "Percentage", x = "Features") +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 10))
ggsave(pSTR_counts_barplot_mRNA, file = "pSTR_counts_barplot_mRNA.pdf", height = 6, width = 8)

# lncRNA-derived pSTRs
lnc_res <- build_composition("lncRNA")
lnc_counts <- lnc_res$counts
lnc_labels <- lnc_res$labels
lnc_counts$group <- factor(lnc_counts$group, levels = c("CDS", "5_UTR", "3_UTR", "intron", "intergenic"))
lnc_counts$period_size <- factor(lnc_counts$period_size, levels = 2:6,
                                 labels = c("Di", "Tri", "Tetra", "Penta", "Hexa"))
lnc_labels$group <- factor(lnc_labels$group, levels = c("CDS", "5_UTR", "3_UTR", "intron", "intergenic"))

pSTR_counts_barplot_lncRNA <- ggplot(lnc_counts, aes(x = group, y = percentage, fill = period_size)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  geom_text(data = lnc_labels, aes(x = group, y = 1, label = label), inherit.aes = FALSE, vjust = -0.2) +
  scale_fill_lancet() +
  scale_y_continuous(expand = c(0.01, 0.05)) +
  theme_bw() +
  labs(title = "The Composition of pSTRs Sourced from lncRNA", y = "Percentage", x = "Features") +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 10))
ggsave(pSTR_counts_barplot_lncRNA, file = "pSTR_counts_barplot_lncRNA.pdf", height = 6, width = 8)

# ----------------------------- (4) Reference allele length distribution -----
# Combine with a "Total" category and filter to ≤100 bp
ral_df <- rbind(
  genetic_char[, c("feature", "period_size", "reference_allele_length")],
  data.frame(feature = "Total", genetic_char[, c("period_size", "reference_allele_length")])
)
ral_df$feature <- factor(ral_df$feature, levels = c("CDS", "5_UTR", "3_UTR", "intron", "intergenic", "Total"))
ral_df$period_size <- factor(ral_df$period_size, levels = 2:6, labels = c("Di", "Tri", "Tetra", "Penta", "Hexa"))
ral_df <- ral_df[ral_df$reference_allele_length <= 100, ]

ral_density <- ggplot(ral_df[ral_df$feature != "Total", ],
                      aes(x = reference_allele_length, color = period_size)) +
  geom_density(position = "identity", linewidth = 1) +
  scale_color_lancet() +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 0.2)) +
  facet_wrap(~feature, nrow = 1) +
  labs(title = "The Reference Allele Length of RNA-pSTRs",
       y = "Density", x = "Reference Allele Length (bp)") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 10))
ggsave(ral_density, file = "ral_distribution_density.pdf", height = 5, width = 15)

# ----------------------------- (5) Copy number distribution ----------------
ncopy_df <- rbind(
  genetic_char[, c("feature", "period_size", "ncopy")],
  data.frame(feature = "Total", genetic_char[, c("period_size", "ncopy")])
)
ncopy_df$feature <- factor(ncopy_df$feature, levels = c("CDS", "5_UTR", "3_UTR", "intron", "intergenic", "Total"))
ncopy_df$period_size <- factor(ncopy_df$period_size, levels = 2:6, labels = c("Di", "Tri", "Tetra", "Penta", "Hexa"))
ncopy_df <- ncopy_df[ncopy_df$ncopy <= 50, ]

ncopy_density <- ggplot(ncopy_df[ncopy_df$feature != "Total", ],
                        aes(x = ncopy, color = period_size)) +
  geom_density(position = "identity", linewidth = 1) +
  scale_color_lancet() +
  scale_x_continuous(limits = c(0, 50)) +
  scale_y_continuous(limits = c(0, 0.6)) +
  facet_wrap(~feature, nrow = 1) +
  labs(title = "The Reference Allele Copy Number of RNA-pSTRs",
       y = "Density", x = "Reference Allele Copy Number") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 10))
ggsave(ncopy_density, file = "ncopy_distribution_density.pdf", height = 5, width = 15)

# ----------------------------- (6) Motif pattern analysis -------------------
# Collapse reverse complements and circular permutations
motif_pattern <- data.frame(
  motif = unique(genetic_char$motif[!grepl("\\/", genetic_char$motif)]),
  motif_pattern = NA,
  stringsAsFactors = FALSE
)

for (p in 1:nrow(motif_pattern)) {
  if (is.na(motif_pattern$motif_pattern[p])) {
    chars <- unlist(strsplit(motif_pattern$motif[p], ""))
    # Reverse complement
    comp <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
    rev_comp <- rev(comp[chars])
    rev_comp_str <- paste(rev_comp, collapse = "")
    orig_str <- paste(chars, collapse = "")
    patterns <- c(orig_str, rev_comp_str)
    # Cyclic permutations
    n <- length(chars)
    if (n > 1) {
      for (i in 1:(n - 1)) {
        shifted <- c(chars[(i + 1):n], chars[1:i])
        shifted_str <- paste(shifted, collapse = "")
        shifted_rev_comp <- rev(comp[shifted])
        shifted_rev_str <- paste(shifted_rev_comp, collapse = "")
        patterns <- c(patterns, shifted_str, shifted_rev_str)
      }
    }
    patterns <- sort(unique(patterns))
    # Use the first as the canonical pattern
    motif_pattern$motif_pattern[motif_pattern$motif %in% patterns] <- patterns[1]
  }
}

motif_df <- data.frame(
  locus = genetic_char$locus,
  feature = genetic_char$feature,
  motif = genetic_char$motif
)
motif_df$pattern <- motif_pattern[match(motif_df$motif, motif_pattern$motif), 2]
motif_df <- na.omit(motif_df)
write.table(motif_df, file = "motif_df.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# Count patterns per feature
motif_counts <- data.frame()
for (grp in unique(motif_df$feature)) {
  pat_tab <- as.data.frame(table(motif_df[motif_df$feature == grp, "pattern"]))
  pat_tab$percentage <- pat_tab$Freq / sum(pat_tab$Freq)
  colnames(pat_tab) <- c("pattern", "counts", "percentage")
  motif_counts <- rbind(motif_counts, data.frame(group = grp, pat_tab))
}
motif_counts$pattern <- as.character(motif_counts$pattern)
motif_counts$period_size <- nchar(motif_counts$pattern)

# For visualization, keep top 10 patterns per period size (across all groups)
top_patterns <- c()
for (ps in 2:6) {
  ps_patterns <- motif_counts[motif_counts$period_size == ps, "pattern"]
  if (length(ps_patterns) <= 10) {
    top_patterns <- c(top_patterns, ps_patterns)
  } else {
    top_patterns <- c(top_patterns, ps_patterns[1:10])
  }
}
motif_counts_filtered <- motif_counts[motif_counts$pattern %in% top_patterns, ]
motif_counts_filtered$group <- factor(motif_counts_filtered$group,
                                      levels = c("CDS", "5_UTR", "3_UTR", "intron", "intergenic"))
motif_counts_filtered$period_size <- factor(motif_counts_filtered$period_size,
                                            levels = 2:6,
                                            labels = c("Di", "Tri", "Tetra", "Penta", "Hexa"))

motif_barplot <- ggplot(motif_counts_filtered, aes(x = reorder(pattern, -counts), y = counts, fill = period_size)) +
  geom_bar(stat = "identity") +
  facet_grid(group ~ period_size, scales = "free") +
  scale_fill_lancet() +
  labs(title = "The composition of motif patterns", y = "Counts", x = "Motif pattern") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = 10),
        legend.position = "none")
ggsave(motif_barplot, file = "motif_barplot.pdf", height = 12, width = 15)

# Subplots: CDS tri-nucleotide motifs
cds_tri <- motif_counts_filtered[motif_counts_filtered$group == "CDS" & motif_counts_filtered$period_size == "Tri", ]
cds_tri_plot <- ggplot(cds_tri, aes(x = reorder(pattern, -counts), y = counts, fill = period_size)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = "#2D6DB1") +
  labs(title = "Motif pattern composition of Tri-pSTR within CDS region",
       y = "Counts", x = "Motif pattern") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = 10),
        legend.position = "none")
ggsave(cds_tri_plot, file = "CDS_tri_motif_barplot.pdf", height = 6, width = 8)

# Intergenic tetra-nucleotide motifs
inter_tetra <- motif_counts_filtered[motif_counts_filtered$group == "intergenic" & motif_counts_filtered$period_size == "Tetra", ]
inter_tetra_plot <- ggplot(inter_tetra, aes(x = reorder(pattern, -counts), y = counts, fill = period_size)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = "#2D6DB1") +
  labs(title = "Motif pattern composition of Tetra-pSTR within intergenic region",
       y = "Counts", x = "Motif pattern") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = 10),
        legend.position = "none")
ggsave(inter_tetra_plot, file = "intergenic_tetra_motif_barplot.pdf", height = 6, width = 8)

# Intergenic penta-nucleotide motifs
inter_penta <- motif_counts_filtered[motif_counts_filtered$group == "intergenic" & motif_counts_filtered$period_size == "Penta", ]
inter_penta_plot <- ggplot(inter_penta, aes(x = reorder(pattern, -counts), y = counts, fill = period_size)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = "#2D6DB1") +
  labs(title = "Motif pattern composition of Penta-pSTR within intergenic region",
       y = "Counts", x = "Motif pattern") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = 10),
        legend.position = "none")
ggsave(inter_penta_plot, file = "intergenic_penta_motif_barplot.pdf", height = 6, width = 8)

# ----------------------------- Supplementary plots --------------------------
# (S1) Mean reference allele length per region per period size (heatmap-like table)
library(tidyr)
ral_mean <- data.frame()
ncopy_mean <- data.frame()
for (ps in c("Di", "Tri", "Tetra", "Penta", "Hexa")) {
  for (reg in c("CDS", "5_UTR", "3_UTR", "intron", "intergenic")) {
    ral_mean_val <- mean(ral_df$reference_allele_length[ral_df$feature == reg & ral_df$period_size == ps], na.rm = TRUE)
    ncopy_mean_val <- mean(ncopy_df$ncopy[ncopy_df$feature == reg & ncopy_df$period_size == ps], na.rm = TRUE)
    ral_mean <- rbind(ral_mean, data.frame(Period_size = ps, Region = reg, ral = ral_mean_val))
    ncopy_mean <- rbind(ncopy_mean, data.frame(Period_size = ps, Region = reg, ncopy = ncopy_mean_val))
  }
}
ral_mean_wide <- spread(ral_mean, key = "Region", value = "ral")
ncopy_mean_wide <- spread(ncopy_mean, key = "Region", value = "ncopy")
write.table(ral_mean_wide, file = "ral_mean_by_region.txt", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(ncopy_mean_wide, file = "ncopy_mean_by_region.txt", row.names = FALSE, quote = FALSE, sep = "\t")

# (S2) Total density plots (for main text or supplement)
ral_total_density <- ggplot(ral_df[ral_df$feature == "Total", ],
                            aes(x = reference_allele_length, fill = period_size)) +
  geom_density(position = "identity", alpha = 0.4, linewidth = 0.8) +
  scale_fill_lancet() +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 0.2)) +
  labs(title = "Reference Allele Length Distribution (Total)", y = "Density", x = "Length (bp)") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 10))
ggsave(ral_total_density, file = "ral_distribution_density2.pdf", height = 5, width = 15)

ncopy_total_density <- ggplot(ncopy_df[ncopy_df$feature == "Total", ],
                              aes(x = ncopy, fill = period_size)) +
  geom_density(position = "identity", alpha = 0.4, linewidth = 0.8) +
  scale_fill_lancet() +
  scale_x_continuous(limits = c(0, 50)) +
  scale_y_continuous(limits = c(0, 0.6)) +
  labs(title = "Reference Allele Copy Number Distribution (Total)", y = "Density", x = "Copy Number") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 10))
ggsave(ncopy_total_density, file = "ncopy_distribution_density2.pdf", height = 5, width = 15)

# (S3) Pie charts for mRNA and lncRNA composition
pie_df <- rbind(mrna_labels, lnc_labels)
pie_df$biotype <- c(rep("protein_coding", nrow(mrna_labels)), rep("lncRNA", nrow(lnc_labels)))
pie_df$label <- paste0(pie_df$group, " (", pie_df$group_sum, ")")
pie_df$biotype <- factor(pie_df$biotype, levels = c("protein_coding", "lncRNA"))

biotype_pieplot <- ggplot(pie_df, aes(x = "", y = percentage, fill = group)) +
  geom_col(color = "black") +
  geom_label_repel(aes(label = label), size = 3.5, fontface = "bold",
                   show.legend = FALSE, position = position_stack(vjust = 0.5)) +
  facet_wrap(~biotype) +
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_npg() +
  theme(strip.text = element_text(face = "bold", size = 13.5),
        legend.position = "none")
ggsave(biotype_pieplot, file = "biotype_pieplot.pdf", height = 5, width = 10)

# End of script