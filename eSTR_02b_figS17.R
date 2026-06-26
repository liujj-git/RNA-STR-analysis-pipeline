# >>> RENAMED from eSTR_S18_patch.R <<<
# NOTE: generates manuscript FIGURE S17 (top-3 eSTR example scatters), not "S18";
# the old filename predates the final numbering.
# RUN ORDER: source("eSTR_02_mapping.R"); source("eSTR_02b_figS17.R")  [same R session]
# ---------------------------------------------------------------------------
# =============================================================================
# eSTR_S18_patch.R — combine the top-3 eSTR example scatters into one Fig S18
#
# Run AFTER eSTR_02_mapping.R in the SAME R session. Reuses in-memory objects:
#   results (ordered by p), ann_intron, d_sub, expr_mat, cov_sub
#   (re-derives locus_col / gene_col to be safe).
# Needs: library(patchwork) and an Arial-capable device (cairo_pdf).
#
# Produces FigS18_top_eSTRs.pdf: panels A/B/C = the three strongest eSTRs,
# each STR length (bp) vs inverse-normal expression, points coloured by cohort,
# with the linear fit and beta / FDR / n. theme_jgg, cairo_pdf.
# =============================================================================
library(ggplot2); library(patchwork)

theme_jgg <- theme_classic(base_family = "Arial", base_size = 9) +
  theme(axis.title = element_text(colour = "black", size = 9),
        axis.text  = element_text(colour = "black", size = 8),
        axis.line  = element_line(colour = "black", linewidth = 0.4),
        axis.ticks = element_line(colour = "black", linewidth = 0.4),
        legend.title = element_text(size = 8), legend.text = element_text(size = 7),
        legend.key.size = unit(3, "mm"),
        plot.title    = element_text(face = "italic", size = 9),
        plot.subtitle = element_text(size = 7, colour = "grey30"),
        plot.tag      = element_text(face = "bold", size = 11))

# ENSG -> gene symbol for the CURRENT top-3 (verify the order matches results!)
sym_map     <- c(ENSG00000156110 = "ADK", ENSG00000188647 = "PTAR1", ENSG00000144118 = "RALB")
cohort_cols <- c(R = "#264653", VB = "#E76F51")

locus_col <- if ("locus" %in% colnames(ann_intron)) "locus" else
             if ("ID" %in% colnames(ann_intron)) "ID" else
             if ("locus_id" %in% colnames(ann_intron)) "locus_id" else "STR_id"
gene_col  <- "gene_id"

top_n  <- min(3L, nrow(results))
panels <- vector("list", top_n)
for (i in seq_len(top_n)) {
  locus_i <- results$locus[i]; gene_i <- results$gene_id[i]; ensg_i <- results$gene_ensg[i]
  idx <- which(ann_intron[[locus_col]] == locus_i & ann_intron[[gene_col]] == gene_i)[1]
  if (is.na(idx)) next
  g <- ensg_i
  if (is.na(g) || !nzchar(g)) g <- ann_intron$gene_id_strip[idx]
  if (is.na(g) || !(g %in% rownames(expr_mat))) next
  d <- d_sub[idx, ]; y <- expr_mat[g, ]
  ok <- !is.na(d) & !is.na(y)
  if (sum(ok) < 2L) next
  df_plot <- data.frame(dosage = d[ok], expr = y[ok], cohort = cov_sub$cohort[ok])
  sym <- if (!is.na(ensg_i) && ensg_i %in% names(sym_map)) sym_map[[ensg_i]] else g
  sub <- sprintf("\u03b2 = %+.3f per bp    FDR = %.1e    n = %d",
                 results$beta_bp[i], results$fdr[i], sum(ok))
  sub <- gsub("e-", "e\u2212", sub)               # tidy minus in the exponent
  panels[[i]] <- ggplot(df_plot, aes(dosage, expr)) +
    geom_point(aes(colour = cohort), alpha = 0.65, size = 1.3) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "black", linewidth = 0.5) +
    scale_colour_manual(values = cohort_cols, name = "Cohort") +
    labs(title = sym, subtitle = sub, x = "STR length (bp)", y = "Inverse-normal expression") +
    theme_jgg
}
panels <- panels[!vapply(panels, is.null, logical(1))]

figS18 <- wrap_plots(panels, nrow = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11), legend.position = "bottom")
ggsave("FigS18_top_eSTRs.pdf", figS18, width = 12, height = 4.2, device = cairo_pdf)
cat("FigS18_top_eSTRs.pdf written (top-3 eSTR scatters: A/B/C).\n")
