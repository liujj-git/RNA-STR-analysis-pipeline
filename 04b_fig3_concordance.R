# >>> RENAMED from fig3_revision_patch.R <<<
# Part of Section 2.4 (genotype-level concordance).
# RUN ORDER: source("04_concordance.R"); source("04b_fig3_concordance.R")  [same R session]
# Produces (final manuscript numbering): Figure 3 (A depth, B deviation,
# C substitution spectrum, D reproducibility) and renumbered concordance
# supplementary figures (S9-S12).
# ---------------------------------------------------------------------------
# =============================================================================
# fig3_revision_patch.R   (v7 — length-level Fig 3 split + A-to-I highlight)
#   * Fig 3 is now length-level only, 4 panels in citation order:
#       A depth, B deviation, C substitution spectrum, D reproducibility
#     (the old dropout panel moved to the allele-level Fig 4 / Fig S13).
#   * Fig 3C and the position supp now highlight the A-to-I editing PAIR
#     (A-to-G + T-to-C), not A-to-G alone.
#   * The old per-sample genotyped-count figure (S9) is dropped; the remaining
#     concordance supp figures are renumbered S10->S9, S11->S10, S12->S11,
#     S13->S12. Dropout-by-depth is produced by fig4_v3.R as Fig S13.
#
# Run AFTER concordance.R in the SAME session. Needs library(patchwork) and an
# Arial-capable device (cairo_pdf).
# =============================================================================
library(ggplot2); library(patchwork)

## ---- unified theme ----------------------------------------------------------
theme_jgg <- theme_classic(base_family = "Arial", base_size = 9) +
  theme(axis.title = element_text(colour = "black", size = 9),
        axis.text  = element_text(colour = "black", size = 8),
        axis.line  = element_line(colour = "black", linewidth = 0.4),
        axis.ticks = element_line(colour = "black", linewidth = 0.4),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(colour = "black", size = 8, margin = margin(2, 2, 2, 2)),
        legend.title = element_text(size = 8), legend.text = element_text(size = 7),
        legend.key.size = unit(3, "mm"),
        plot.tag = element_text(face = "bold", size = 11))

## ---- semantic concordance palette -------------------------------------------
conc_levels <- c("identical", "partial_identical", "non_identical")
conc_labels <- c("Identical", "Partial", "Non-identical")
conc_cols   <- c(identical = "#2D6DB1", partial_identical = "#E69F00", non_identical = "#D7263D")
conc_fill   <- scale_fill_manual(values = conc_cols, breaks = conc_levels,
                                 labels = conc_labels, name = "Concordance")
cmp_cols    <- c("RNA-DNA" = "#444444", "RNA-RNA" = "#0B7A75")

# ============================================================================
# A  read depth: discordant vs identical
# ============================================================================
depthDF <- rbind(
  data.frame(comparison = "RNA-DNA",
             grp = ifelse(ccd$concordance == 2, "Identical", "Discordant"), dp = ccd$DP_RNA),
  data.frame(comparison = "RNA-RNA",
             grp = ifelse(rrd$concordance == 2, "Identical", "Discordant"), dp = pmin(rrd$DP1, rrd$DP2)))
depthDF <- depthDF[is.finite(depthDF$dp), ]
depthDF$grp <- factor(depthDF$grp, levels = c("Discordant", "Identical"))
pDepth <- ggplot(depthDF, aes(grp, log2(dp + 1), fill = grp)) +
  geom_violin(colour = NA, alpha = 0.35, width = 0.9) +
  geom_boxplot(width = 0.28, outlier.size = 0.2, linewidth = 0.3, colour = "black") +
  facet_wrap(~comparison) +
  scale_fill_manual(values = c(Discordant = "#D7263D", Identical = "#2D6DB1"), guide = "none") +
  theme_jgg + labs(x = NULL, y = "log2(DP + 1)")

# ============================================================================
# B  identical-call rate vs deviation (line)
# ============================================================================
lineDF <- rbind(
  data.frame(comparison = "RNA-DNA", dev_bin = dev_tab$dev_bin,
             conc = as.character(dev_tab$concordance),    pct = dev_tab$pct),
  data.frame(comparison = "RNA-RNA", dev_bin = rr_dev_tab$dev_bin,
             conc = as.character(rr_dev_tab$concordance), pct = rr_dev_tab$pct))
lineDF <- lineDF[lineDF$conc == "identical", ]
lineDF$dev_bin <- factor(as.character(lineDF$dev_bin), levels = c("0", "1", "2", "3", "4", ">=5"))
pDev <- ggplot(lineDF, aes(dev_bin, pct, colour = comparison, linetype = comparison, group = comparison)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.8) +
  scale_colour_manual(values = cmp_cols, name = NULL) +
  scale_linetype_manual(values = c("RNA-DNA" = "solid", "RNA-RNA" = "longdash"), name = NULL) +
  scale_y_continuous(limits = c(50, 100)) +
  theme_jgg + theme(legend.position = c(0.98, 0.03), legend.justification = c(1, 0),
                    legend.background = element_rect(fill = "white", colour = "grey80")) +
  labs(x = "Max allele deviation from reference (copies)", y = "Identical calls (%)")

# ============================================================================
# C  substitution spectrum; A-to-I PAIR (A-to-G + T-to-C) highlighted (lollipop)
# ============================================================================
freq <- iav_changed_base_pattern_compared
freq <- freq[freq$region == "Total", c("pattern", "percentage")]
freq <- freq[order(-freq$percentage), ]
freq$pattern <- factor(as.character(freq$pattern), levels = as.character(freq$pattern))
.norm <- toupper(gsub("[^A-Za-z]", "", as.character(freq$pattern)))
freq$hl <- ifelse(.norm %in% c("ATOG", "TTOC"), "A-to-I (A>G/T>C)", "Other")
pSpec <- ggplot(freq, aes(pattern, 100 * percentage, colour = hl)) +
  geom_segment(aes(xend = pattern, y = 0, yend = 100 * percentage), linewidth = 0.9) +
  geom_point(aes(size = hl)) +
  scale_colour_manual(values = c("A-to-I (A>G/T>C)" = "#D7263D", "Other" = "grey60"), guide = "none") +
  scale_size_manual(values = c("A-to-I (A>G/T>C)" = 2.8, "Other" = 1.8), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  theme_jgg + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
  labs(x = "Substitution type", y = "% of substitutions")

# ============================================================================
# D  reproducible discordance per individual (15 pairs with matched DNA)
# ============================================================================
ind15 <- as.character(repeated_samples_compared$individual)
grid  <- expand.grid(individual = ind15,
                     group = c("partial_identical", "non_identical"), stringsAsFactors = FALSE)
cnt <- as.data.frame(table(individual = repeated_samples_compared_df$individual,
                           group      = repeated_samples_compared_df$group), stringsAsFactors = FALSE)
repDF <- merge(grid, cnt, all.x = TRUE); repDF$Freq[is.na(repDF$Freq)] <- 0
ord  <- names(sort(tapply(repDF$Freq, repDF$individual, sum)))
repDF$individual <- factor(repDF$individual, levels = ord)
repDF$group <- factor(repDF$group, levels = c("partial_identical", "non_identical"),
                      labels = c("Partial", "Non-identical"))
pRepro <- ggplot(repDF, aes(individual, Freq, fill = group)) +
  geom_col(position = position_dodge(0.78), colour = "black", width = 0.68, linewidth = 0.3) +
  scale_fill_manual(values = c("Partial" = conc_cols[["partial_identical"]],
                               "Non-identical" = conc_cols[["non_identical"]]),
                    name = "Reproducible discordance") +
  coord_flip() +
  theme_jgg + theme(axis.text.y = element_text(size = 7), legend.position = "top") +
  labs(x = "Individual", y = "Reproducible discordant loci (both replicates)")

# ============================================================================
# Assemble Fig 3 (length-level): tags A-D in citation order; cairo_pdf
# ============================================================================
fig3 <- (pDepth | pDev | pSpec) / pRepro +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11))
ggsave("Fig3_composite.pdf", fig3, width = 12, height = 9, device = cairo_pdf)
cat("Fig3_composite.pdf written (A depth / B deviation / C spectrum / D reproducibility).\n")

# ============================================================================
# Fig S9 (was S10) — concordance composition PER SAMPLE (A RNA-RNA, B RNA-DNA)
# ============================================================================
rr_per_id <- do.call(rbind, lapply(unique(concordance_repeated_df$individual), function(i) {
  s <- concordance_repeated_df[concordance_repeated_df$individual == i, ]
  data.frame(individual = i, identical = mean(s$concordance == 2),
             partial_identical = mean(s$concordance == 1), non_identical = mean(s$concordance == 0))
}))
to_long <- function(df, pcol) data.frame(
  individual = rep(df$individual, 3),
  conc = factor(rep(conc_levels, each = nrow(df)), levels = rev(conc_levels)),
  pct  = 100 * c(df$identical, df[[pcol]], df$non_identical))
order_ind <- function(df) df$individual[order(-df$identical)]
rr_long <- to_long(rr_per_id, "partial_identical"); rr_long$individual <- factor(rr_long$individual, levels = order_ind(rr_per_id))
rd_long <- to_long(per_ind,   "partial");           rd_long$individual <- factor(rd_long$individual, levels = order_ind(per_ind))
mk_comp <- function(d) ggplot(d, aes(individual, pct, fill = conc)) +
  geom_col(colour = "black", width = 0.85, linewidth = 0.25) + conc_fill +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  theme_jgg + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5)) +
  labs(x = "Individual", y = "% of loci")
ggsave("FigS9_concordance_composition_persample.pdf",
       (mk_comp(rr_long) / mk_comp(rd_long)) + plot_layout(guides = "collect") +
         plot_annotation(tag_levels = "A") &
         theme(plot.tag = element_text(face = "bold", size = 11), legend.position = "right"),
       width = 12, height = 8, device = cairo_pdf)

# ============================================================================
# Fig S10 (was S11) — concordance by region
# ============================================================================
recolour_region <- function(df) {
  df$group <- factor(as.character(df$group), levels = conc_levels)
  ggplot(df, aes(region, percentage, fill = group)) +
    geom_col(position = "fill", colour = "black", linewidth = 0.25) + conc_fill + theme_jgg +
    labs(x = "Region", y = "Proportion")
}
ggsave("FigS10_concordance_by_region.pdf",
       (recolour_region(concordance_region_repeated_df) | recolour_region(concordance_region_compared_df)) +
         plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
         theme(plot.tag = element_text(face = "bold", size = 11), legend.position = "right"),
       width = 10, height = 4.5, device = cairo_pdf)

# ============================================================================
# Fig S11 (was S12) — step size of partial-identical RNA-DNA calls
# ============================================================================
sx1 <- halfidentical_compared_df2[halfidentical_compared_df2$DPthreshold == 0, ]
sx1$group <- factor(as.character(sx1$group),
                    levels = c("negative_length_discordance", "sequence_discordance", "positive_length_discordance"),
                    labels = c("Length (RNA shorter)", "Same length (sequence)", "Length (RNA longer)"))
pSX1 <- ggplot(sx1, aes(dif_step, counts, fill = group)) +
  geom_col(colour = "black", width = 0.85, linewidth = 0.25) +
  scale_fill_manual(values = c("Length (RNA shorter)" = "#2D6DB1",
                               "Same length (sequence)" = "grey65",
                               "Length (RNA longer)" = "#D7263D"), name = NULL) +
  scale_x_continuous(breaks = seq(-10, 10, 2)) +
  theme_jgg +
  labs(x = "Allele-length difference, RNA minus DNA (repeat units)", y = "Genotype pairs")
ggsave("FigS11_stepsize_partial.pdf", pSX1, width = 7, height = 4.5, device = cairo_pdf)

# ============================================================================
# Fig S12 (was S13) — substitution position; DENSITY; A-to-I pair vs other
# ============================================================================
pos_src <- iav_changed_base_compared_pattern2
names(pos_src)[names(pos_src) == "Base_Variant_patterns"] <- "pattern"
.normp <- toupper(gsub("[^A-Za-z]", "", as.character(pos_src$pattern)))
pos_src$grp <- ifelse(.normp %in% c("ATOG", "TTOC"), "A-to-I (A>G/T>C)", "Other (10 types)")
ggsave("FigS12_subst_position_2curve.pdf",
       ggplot(pos_src, aes(dif_location_percentage, after_stat(density), colour = grp)) +
         geom_density(linewidth = 0.9) + scale_x_continuous(limits = c(0, 1)) +
         scale_colour_manual(values = c("A-to-I (A>G/T>C)" = "#D7263D", "Other (10 types)" = "grey55"), name = NULL) +
         theme_jgg +
         labs(x = "Relative position in allele (0 = 5', 1 = 3')", y = "Density"),
       width = 7, height = 4.5, device = cairo_pdf)

cat("All length-level figures written (Fig 3 + Fig S9-S12); dropout-by-depth is now Fig S13 via fig4_v3.R.\n")