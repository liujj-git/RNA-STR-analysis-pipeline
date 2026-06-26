# >>> RENAMED from fig4_revision_patch.R <<<
# NOTE: this script generates manuscript FIGURE 5 (Section 2.5), not "Fig 4";
# the old in-file comments predate the final numbering.
# RUN ORDER: source("05_polymorphism.R"); source("05b_fig5_allele_diversity.R")  [same R session]
# Produces (final manuscript numbering): Figure 5 (A major-allele frequency,
# B step-size, C number of alleles Na, D PIC) and Figure S14.
# ---------------------------------------------------------------------------
# =============================================================================
# fig4_revision_patch.R   (v1 — restyle + relayout of Fig 4 and Fig S14 for 2.4)
#
# Run AFTER polymorphism.R in the SAME R session. Reuses in-memory objects:
#   major_allele_df  ($group major_allele/alt_allele, $freq, $Ncopy, $ncopy_dif)
#   mar              ($Ncopy major copy, $reference_ncopy, $ncopy_dif, $PIC; |dif|<=10)
#   polymorphism_region ($group region, $period_size Di..Hexa, $Na, $PIC)
# Needs: library(patchwork) and an Arial-capable device (cairo_pdf).
#
# CHANGES vs the published Fig 4:
#   * dropped the major-vs-reference panel (old 4b): essentially one number
#     (95.76% match), better left in text -> 4 panels now (A,B,C,D)
#   * new layout: (A frequency | B step) / C Na (full width) / D PIC (full width)
#   * removed in-figure Mean/SD text; motif palette unified with Figure 2
#   * dual y-axes removed; theme_jgg + cairo_pdf throughout
#   * Fig S14 collapsed to ONE panel (PIC vs major-allele copy), label r (lower-case)
# =============================================================================
library(ggplot2); library(patchwork)

## ---- unified theme + motif palette (from genetic_character.R / Figure 2) -----
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
motif_cols  <- c(Di = "#2D6DB1", Tri = "#E69F00", Tetra = "#2E933C",
                 Penta = "#56B4E9", Hexa = "#B279A7")
region_labs <- c(CDS = "CDS", "5_UTR" = "5\u2032 UTR", "3_UTR" = "3\u2032 UTR",
                 intron = "intron", intergenic = "intergenic")

# =============================================================================
# VERIFICATION VALUES  (printed so the text numbers can be checked)
# =============================================================================
maj <- major_allele_df[major_allele_df$group == "major_allele", ]
alt <- major_allele_df[major_allele_df$group == "alt_allele", ]
alt_i <- alt[alt$ncopy_dif %% 1 == 0 & abs(alt$ncopy_dif) <= 10, ]   # integer diffs, +/-10

cat("\n================  2.4 VERIFICATION  ================\n")
cat(sprintf("Major-allele freq >= 0.5:             %.2f%% of loci\n", 100 * mean(maj$freq >= 0.5)))
cat(sprintf("Major allele == reference (dif==0):   %.2f%% (among |dif|<=10)\n", 100 * mean(mar$ncopy_dif == 0)))
cat(sprintf("Multiallelic (Na > 2):                %.2f%%   Biallelic (Na==2): %.2f%%\n",
            100 * mean(polymorphism_region$Na > 2), 100 * mean(polymorphism_region$Na == 2)))
cat(sprintf("Minor alleles: single-unit (|dif|=1): %.2f%% | shorter(<0): %.2f%% | longer(>0): %.2f%%\n",
            100 * mean(abs(alt_i$ncopy_dif) == 1),
            100 * mean(alt_i$ncopy_dif < 0), 100 * mean(alt_i$ncopy_dif > 0)))

cat("\n-- Mean Na by motif (di..hexa) --\n")
print(aggregate(Na ~ period_size, polymorphism_region, function(x) round(mean(x), 3)))
cat("\n-- Mean Na by region --\n")
print(aggregate(Na ~ group, polymorphism_region, function(x) round(mean(x), 3)))
cat("\n-- Mean PIC by motif (di..hexa) --\n")
print(aggregate(PIC ~ period_size, polymorphism_region, function(x) round(mean(x), 3)))
cat("\n-- Mean PIC by region --\n")
print(aggregate(PIC ~ group, polymorphism_region, function(x) round(mean(x), 3)))

ct_major <- cor.test(mar$Ncopy,           mar$PIC)
ct_ref   <- cor.test(mar$reference_ncopy, mar$PIC)
cat("\n-- PIC vs copy-number correlation (Pearson) --\n")
cat(sprintf("  PIC ~ major-allele copy:  r = %.4f  (P = %.3g)\n", ct_major$estimate, ct_major$p.value))
cat(sprintf("  PIC ~ reference copy:     r = %.4f  (P = %.3g)\n", ct_ref$estimate,   ct_ref$p.value))
cat("===================================================\n\n")

# =============================================================================
# Fig 4 panels
# =============================================================================
# A  major-allele frequency distribution (log-scaled count)
pA <- ggplot(maj, aes(freq)) +
  geom_histogram(bins = 20, fill = "#2D6DB1", colour = "black", linewidth = 0.2) +
  scale_y_log10() +
  scale_x_continuous(limits = c(0, 1)) +
  theme_jgg + labs(x = "Major-allele frequency", y = "Loci (log scale)")

# B  step-size of minor (alt) alleles relative to the major allele (central range)
stepDF <- as.data.frame(table(dif = alt_i$ncopy_dif[abs(alt_i$ncopy_dif) <= 6]))
stepDF$dif <- as.integer(as.character(stepDF$dif))
stepDF$pct <- 100 * stepDF$Freq / nrow(alt_i)        # % of all minor alleles
pB <- ggplot(stepDF, aes(factor(dif), pct)) +
  geom_col(fill = "#2D6DB1", colour = "black", linewidth = 0.2, width = 0.85) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_jgg +
  labs(x = "Minor-allele size minus major (repeat units)", y = "% of minor alleles")

# C  Na by region x motif (full width, no text annotations)
pC <- ggplot(polymorphism_region, aes(group, Na, fill = period_size)) +
  geom_boxplot(position = position_dodge(0.8), outlier.size = 0.2,
               linewidth = 0.3, colour = "black") +
  scale_fill_manual(values = motif_cols, name = "Motif") +
  scale_x_discrete(labels = region_labs) +
  coord_cartesian(ylim = c(2, 12)) +     # trim extreme outliers for readability
  theme_jgg + labs(x = NULL, y = "Number of alleles (Na)")

# D  PIC by region x motif (full width)
pD <- ggplot(polymorphism_region, aes(group, PIC, fill = period_size)) +
  geom_boxplot(position = position_dodge(0.8), outlier.size = 0.2,
               linewidth = 0.3, colour = "black") +
  scale_fill_manual(values = motif_cols, name = "Motif") +
  scale_x_discrete(labels = region_labs) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.03))) +
  theme_jgg + labs(x = "Region", y = "PIC")

# ---- assemble: A|B on top, C and D each full width; collect the motif legend ----
fig4 <- (pA | pB) / pC / pD +
  plot_layout(heights = c(0.9, 1, 1), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11), legend.position = "right")
ggsave("Fig4_composite.pdf", fig4, width = 11, height = 11, device = cairo_pdf)
cat("Fig4_composite.pdf written (A frequency | B step / C Na / D PIC).\n")

# =============================================================================
# Fig S14 — PIC vs copy number, single panel (major-allele copy)
# =============================================================================
r_lab <- paste0("italic(r) == ", formatC(ct_major$estimate, format = "f", digits = 2))
pS14 <- ggplot(mar, aes(Ncopy, PIC)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#2D6DB1", fill = "#A6BDDB") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.15, vjust = 1.6,
           label = r_lab, parse = TRUE, size = 3) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.12, vjust = 3.3,
           label = "italic(P) < 0.001", parse = TRUE, size = 3) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_jgg +
  labs(x = "Copy number of major allele", y = "PIC")
ggsave("FigS14_PIC_vs_copynumber.pdf", pS14, width = 6, height = 4, device = cairo_pdf)
cat("FigS14_PIC_vs_copynumber.pdf written (single panel, r lower-case).\n")
