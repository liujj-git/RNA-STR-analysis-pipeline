# >>> RENAMED from figS1_S2_S4_patch.R <<<
# Part of Section 2.1 supplementary figures.
# RUN ORDER: source("01_expression.R"); source("01b_figS1_S2_S4.R")
# Produces (final manuscript numbering): Figure S1, Figure S2, Figure S4.
# (Figure S3 = per-locus DP pheatmaps, produced directly by 01_expression.R.)
# ---------------------------------------------------------------------------
# =============================================================================
# figS1_S2_S4_patch.R  —  unified-style (theme_jgg) rebuild of the Section 2.1
# supplementary figures, with MANUSCRIPT-FAITHFUL numbering.
#
# Source AFTER expression.R (uses objects: fig1, dp_bin_df).
#   source("expression.R"); source("figS1_S2_S4_patch.R")
#
# Renumbering note — expression.R's in-script comments use an OLD scheme that no
# longer matches the locked text. These are the FINAL numbers:
#   Fig S1  DP & genotyping rate vs degree of polymorphism   -> FigureS1.pdf  (NEW: GR panel added)
#   Fig S2  mean DP by motif length x reference-allele length -> FigureS2.pdf  (was "Fig 1c"/heatmap_summary)
#   Fig S4  GR & DP vs reference-allele length (correlations) -> FigureS4.pdf  (was "Fig S3a/S3b"/corplots)
# Fig S3 (per-locus DP pheatmaps) is left as pheatmap output — title only, no patch.
# Fig S5 is built in annotation.R, already in this style.
#
# Style matches Figure 2 and Figure S5: theme_jgg (Arial, thin black axes,
# grey92 strips), cairo_pdf vector export.
# =============================================================================

stopifnot(exists("fig1"), exists("dp_bin_df"))   # must run after expression.R

library(ggplot2)
library(patchwork)
library(ggpubr)      # stat_cor (Fig S4)

## ---- canonical paper theme (identical to genetic_character.R / annotation.R) -
theme_jgg <- theme_classic(base_family = "Arial", base_size = 9) +
  theme(
    axis.title      = element_text(colour = "black", size = 9),
    axis.text       = element_text(colour = "black", size = 8),
    axis.line       = element_line(colour = "black", linewidth = 0.4),
    axis.ticks      = element_line(colour = "black", linewidth = 0.4),
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text      = element_text(colour = "black", size = 8, margin = margin(2, 2, 2, 2)),
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 7),
    legend.key.size = unit(3, "mm"),
    plot.subtitle   = element_text(size = 7, hjust = 0, margin = margin(b = 2)),
    plot.tag        = element_text(face = "bold", size = 11))

grp_labeller <- as_labeller(c(polymorphic = "pSTR", monomorphic = "mSTR"))

# =============================================================================
# Figure S1 : read depth and genotyping rate vs degree of polymorphism
#   A = log2 DP by heterozygosity bin ; B = genotyping rate by heterozygosity bin
#   (the GR panel is new; the DP panel mirrors dp_quality_by_polymorphism_bin)
# =============================================================================
poly_bin <- cut(fig1$Heterozygosity,
                breaks = c(-Inf, 0, 0.1, 0.2, 0.3, 0.5, Inf),
                labels = c("0", "(0,0.1]", "(0.1,0.2]", "(0.2,0.3]", "(0.3,0.5]", ">0.5"))
s1 <- data.frame(Ho_bin = poly_bin,
                 logDP  = log2(fig1$meanDP + 1),
                 GR     = fig1$genotyping_rate)

# Spearman correlation on the continuous heterozygosity (the trend the text cites)
rho_dp <- suppressWarnings(cor(fig1$Heterozygosity, fig1$meanDP,          method = "spearman"))
rho_gr <- suppressWarnings(cor(fig1$Heterozygosity, fig1$genotyping_rate, method = "spearman"))

# sequential blues so the ordinal bins read low -> high (ColorBrewer Blues, 6-class)
blues6 <- c("#EFF3FF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C")

pS1a <- ggplot(s1, aes(Ho_bin, logDP, fill = Ho_bin)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.3, width = 0.7) +
  scale_fill_manual(values = blues6, guide = "none") +
  labs(x = NULL, y = expression(log[2] * "(mean DP + 1)"),
       subtitle = bquote("Spearman " * italic(rho) == .(round(rho_dp, 2)))) +
  theme_jgg +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

pS1b <- ggplot(s1, aes(Ho_bin, GR, fill = Ho_bin)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.3, width = 0.7) +
  scale_fill_manual(values = blues6, guide = "none") +
  labs(x = "Heterozygosity", y = "Genotyping rate",
       subtitle = bquote("Spearman " * italic(rho) == .(round(rho_gr, 2)))) +
  theme_jgg +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

figS1 <- pS1a / pS1b + plot_annotation(tag_levels = "A")
ggsave("FigureS1.pdf", figS1, width = 5, height = 5.6, device = cairo_pdf)

# =============================================================================
# Figure S2 : mean DP by motif length x reference-allele length (pSTR | mSTR)
#   dp_bin_df is built in expression.R; here we only restyle.
# =============================================================================
figS2 <- ggplot(dp_bin_df, aes(motif_length, reference_allele_length, fill = meanDP)) +
  facet_wrap(~ group, nrow = 2, labeller = grp_labeller) +
  geom_tile(colour = "grey85", linewidth = 0.3) +
  geom_text(aes(label = label), colour = "black", size = 2.4) +
  scale_fill_gradient(low = "white", high = "firebrick3",
                      na.value = "grey92", name = "Mean DP") +
  coord_fixed() +
  labs(x = "Motif length", y = "Reference-allele length (bp)") +
  theme_jgg +
  theme(legend.position = "right")
ggsave("FigureS2.pdf", figS2, width = 4.8, height = 6.2, device = cairo_pdf)

# =============================================================================
# Figure S4 : genotyping rate and DP vs reference-allele length
#   A = genotyping rate ; B = log2 DP ; each faceted pSTR | mSTR with Pearson r
# =============================================================================
pS4a <- ggplot(fig1, aes(reference_allele_length, genotyping_rate)) +
  geom_smooth(method = "lm", formula = y ~ x,
              colour = "#2D6DB1", fill = "grey80", linewidth = 0.6) +
  facet_grid(~ group_f, scales = "free", labeller = grp_labeller) +
  stat_cor(method = "pearson", size = 2.6,
           label.x.npc = "left", label.y.npc = "top") +
  labs(x = NULL, y = "Genotyping rate") +
  theme_jgg

pS4b <- ggplot(fig1, aes(reference_allele_length, log2(meanDP + 1))) +
  geom_smooth(method = "lm", formula = y ~ x,
              colour = "#2D6DB1", fill = "grey80", linewidth = 0.6) +
  facet_grid(~ group_f, scales = "free", labeller = grp_labeller) +
  stat_cor(method = "pearson", size = 2.6,
           label.x.npc = "left", label.y.npc = "top") +
  labs(x = "Reference-allele length (bp)", y = expression(log[2] * "(mean DP + 1)")) +
  theme_jgg

figS4 <- pS4a / pS4b + plot_annotation(tag_levels = "A")
ggsave("FigureS4.pdf", figS4, width = 6, height = 5.6, device = cairo_pdf)

cat("Wrote FigureS1.pdf, FigureS2.pdf, FigureS4.pdf (theme_jgg, cairo_pdf).\n")
cat(sprintf("Fig S1 annotations: Spearman rho DP = %.2f , GR = %.2f\n", rho_dp, rho_gr))
