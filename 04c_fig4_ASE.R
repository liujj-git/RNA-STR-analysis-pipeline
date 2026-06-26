# >>> RENAMED from fig4_v3.R  (fig4_v2.R was the discarded richer draft) <<<
# Part of Section 2.4 (allele-level / allele-specific expression).
# RUN ORDER: source("04_concordance.R") -> ase analysis -> source("04c_fig4_ASE.R")
# Produces (final manuscript numbering): Figure 4 (A allelic dropout vs
# allele-length difference; B minor-allele read fraction by RNA depth) and
# Figure S13. Reads the ASE tables written during the Section 2.4 run.
# ---------------------------------------------------------------------------
# =============================================================================
# fig4_v3.R  — trimmed allele-level figures from ase_patch.R outputs.
#   MAIN  Fig 4  : A dropout vs allele-length difference  | B minor-fraction by depth
#   SUPP  Fig S13: A dropout vs RNA depth                 | B minor-fraction distribution
# Reads ase_dropout_by_copydiff.txt, ase_dropout_by_depth.txt, ase_balance.txt
# (written by ase_patch.R). Run after ase_patch.R, same working dir.
# =============================================================================
library(ggplot2); library(patchwork)

theme_jgg <- theme_classic(base_family = "Arial", base_size = 9) +
  theme(axis.title = element_text(colour = "black", size = 9),
        axis.text  = element_text(colour = "black", size = 8),
        axis.line  = element_line(colour = "black", linewidth = 0.4),
        axis.ticks = element_line(colour = "black", linewidth = 0.4),
        plot.tag = element_text(face = "bold", size = 11))
dp_lv <- c("<10", "10-19", "20-29", ">=30")

drc <- read.table("ase_dropout_by_copydiff.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
drd <- read.table("ase_dropout_by_depth.txt",    header = TRUE, sep = "\t", stringsAsFactors = FALSE)
bal <- read.table("ase_balance.txt",             header = TRUE, sep = "\t", stringsAsFactors = FALSE)
bal$RNA_DP_bin <- factor(as.character(bal$RNA_DP_bin), levels = dp_lv)

# ---- MAIN A: dropout vs copy-number difference (lollipop; "0" -> "<1") -------
drc$dCopy_bin <- as.character(drc$dCopy_bin); drc$dCopy_bin[drc$dCopy_bin == "0"] <- "<1"
drc$dCopy_bin <- factor(drc$dCopy_bin, levels = c("<1", "1", "2", "3", "4", ">=5"))
pCopy <- ggplot(drc, aes(dCopy_bin, 100 * dropout_rate)) +
  geom_segment(aes(xend = dCopy_bin, y = 0, yend = 100 * dropout_rate), colour = "grey65", linewidth = 0.6) +
  geom_point(colour = "#D7263D", size = 2.6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10)), limits = c(0, NA)) +
  theme_jgg + labs(x = "Copy-number difference between DNA alleles", y = "Allelic dropout (%)")

# ---- MAIN B: minor-allele read fraction by RNA depth (violin + box) ----------
pDepthBal <- ggplot(bal, aes(RNA_DP_bin, minor_frac, fill = RNA_DP_bin)) +
  geom_violin(colour = NA, alpha = 0.35, width = 0.9, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.20, outlier.size = 0.12, linewidth = 0.3, colour = "black") +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.4) +
  scale_fill_manual(values = c("#9ECAE1", "#6BAED6", "#3182BD", "#08519C"), guide = "none") +
  theme_jgg + labs(x = "RNA read depth", y = "Minor-allele read fraction")

# ---- SUPP A: dropout vs RNA read depth (line) --------------------------------
drd$RNA_DP_bin <- factor(as.character(drd$RNA_DP_bin), levels = dp_lv)
pDepthDrop <- ggplot(drd, aes(RNA_DP_bin, 100 * dropout_rate, group = 1)) +
  geom_line(colour = "#D7263D", linewidth = 0.9) + geom_point(colour = "#D7263D", size = 2.4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10)), limits = c(0, NA)) +
  theme_jgg + labs(x = "RNA read depth", y = "Allelic dropout (%)")

# ---- SUPP B: minor-allele read fraction distribution (histogram) -------------
pDist <- ggplot(bal, aes(minor_frac)) +
  geom_histogram(binwidth = 0.02, fill = "#2D6DB1", colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_jgg + labs(x = "Minor-allele read fraction", y = "Heterozygous calls")

fig4 <- (pCopy | pDepthBal) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11))
ggsave("Fig4_composite.pdf", fig4, width = 9, height = 4.2, device = cairo_pdf)

figS13 <- (pDepthDrop | pDist) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11))
ggsave("FigS13_allele_depth.pdf", figS13, width = 9, height = 4.2, device = cairo_pdf)

cat("Fig4_composite.pdf (A dropout/copy, B minor-frac/depth) and",
    "FigS13_allele_depth.pdf (A dropout/depth, B minor-frac distribution) written.\n")
