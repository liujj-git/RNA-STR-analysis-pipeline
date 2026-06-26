# =============================================================================
# expression.R  —  clean, manuscript-faithful re-run (revision build)
# Manuscript: Sections 2.1-2.2 ; Figures 1a-1e, S1, S2, S3a-S3b ; Tables S1, S2
#
# Runs BEFORE annotation.R and uses NO VEP terms -> unaffected by the intergenic
# re-coding of scripts 3-6.
#
# Reproduces the reported counts UNCHANGED:
#   796,550 targeted | 44,328 genotyped (GR>=0.5) | 24,013 pSTR | 20,315 mSTR
#   538,033 (67.55%) ungenotyped
#
# Two things relative to the very first script:
#   (1) Fig 1d/1e significance test is now an explicit pSTR-vs-mSTR comparison
#       (two-sided Wilcoxon) within each motif class; conclusion unchanged.
#   (2) Section 4 writes EVERY revision-relevant value to a named file so each
#       item in the writing list is traceable:
#         revision_values.txt        master reference (headline + thresholds + r)
#         chromosome_counts.txt       per-chromosome pSTR/mSTR counts  (Fig S1, 2.1)
#         fig1c_dp_by_bin.txt         binned mean DP behind Figure 1c
#         fig1de_test_pvalues.txt     exact Wilcoxon p for Figures 1d/1e
#         mono_multiallele_loci.txt   the monomorphic loci carrying >=2 alleles
# =============================================================================

rm(list = ls())

library(dplyr)      # left_join
library(readxl)
library(vcfR)
library(ggplot2)
library(ggpubr)     # stat_compare_means, stat_cor
library(pheatmap)
library(ggsci)      # lancet palette

## ---- shared plot theme ------------------------------------------------------
base_theme <- theme_bw() +
  theme(axis.title   = element_text(face = "bold", colour = "black", size = 15),
        axis.title.x = element_text(face = "bold", colour = "black", size = 12),
        axis.title.y = element_text(face = "bold", colour = "black", size = 12),
        axis.text.x  = element_text(face = "bold", colour = "black", size = 10),
        axis.text.y  = element_text(face = "bold", colour = "black", size = 10))

period_levels <- 2:6
period_labels <- c("Di", "Tri", "Tetra", "Penta", "Hexa")
group_levels  <- c("polymorphic", "monomorphic")
fct_period <- function(x) factor(x, levels = period_levels, labels = period_labels)
fct_group  <- function(x) factor(x, levels = group_levels)

# =============================================================================
# 0. Inputs
# =============================================================================
unrelated_individuals <- read_excel("../VBsampleinfo.xlsx", sheet = "unrelated_individuals")
colnames(unrelated_individuals) <- unrelated_individuals[1, ]
unrelated_individuals <- unrelated_individuals[-1, ]

sample_gender <- read_excel("../VBsampleinfo.xlsx", sheet = "total")
colnames(sample_gender) <- sample_gender[1, ]
sample_gender <- sample_gender[-1, ]
unrelated_individuals$gender <-
  sample_gender$gender[match(unrelated_individuals$ID, sample_gender$ID)]

n_total <- nrow(unrelated_individuals)
n_male  <- sum(unrelated_individuals$gender == "male")

# config cols: 1 chr | 2 start | 3 end | 4 period | 5 copy_number | 6 locus | 7 motif
config <- read.table("../GRCh38.hipstr_reference.refine.bed")
config$reference_allele_length <- config[, 3] - config[, 2] + 1   # 1-based inclusive

wrong_locus <- read.table("removedlocus_with_wrongallele.txt", fill = TRUE)

raw_vcf <- read.vcfR("hipstr.filtered.rna.vcf.gz")
GT <- extract.gt(raw_vcf)
GT <- GT[!rownames(GT) %in% wrong_locus[, 1],
         colnames(GT) %in% unrelated_individuals$ID]
DP <- extract.gt(raw_vcf, element = "DP")
DP <- DP[, colnames(DP) %in% unrelated_individuals$ID]
rm(raw_vcf)

# =============================================================================
# 1. Genotyping rate + pSTR / mSTR classification
# =============================================================================
# Y loci: denominator = males only.  X loci: full cohort (males hemizygous; a
# simplification kept to reproduce the reported X = 595 pSTR / 626 mSTR).
gr_df <- data.frame(counts = apply(GT, 1, function(x) sum(!is.na(x))))
gr_df$chr   <- config[match(rownames(gr_df), config[, 6]), 1]
gr_df$denom <- ifelse(gr_df$chr == "Y", n_male, n_total)
gr_df$genotyping_rate <- gr_df$counts / gr_df$denom

gr_threshold <- 0.5
locus_with_sufficient_gr <- rownames(gr_df)[gr_df$genotyping_rate >= gr_threshold]

## Table S2: loci never expressed (GR == 0) ----------------------------------
expressed_locus    <- rownames(gr_df)[gr_df$genotyping_rate != 0]
loci_no_expression <- config[!config[, 6] %in% expressed_locus, -ncol(config)]
colnames(loci_no_expression) <-
  c("Chr", "Start", "End", "Motif_length", "Copy_Number", "Locus", "Motif_Sequence")
write.table(loci_no_expression, "Loci_withnoexpression.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

GT <- GT[rownames(GT) %in% locus_with_sufficient_gr, ]

## ---- polymorphic vs monomorphic + observed heterozygosity ------------------
# polymorphic (pSTR) = >1 DISTINCT GENOTYPE observed across the cohort;
# monomorphic (mSTR) = a single genotype shared by all genotyped individuals.
# QC (Section 4) reports the only locus where this differs from an allele-level
# definition. Keeping this definition preserves 24,013 / 20,315.
is_het <- function(gt) length(unique(strsplit(gt, "\\|")[[1]])) > 1
genotype_summary <- function(x) {
  g   <- na.omit(x)
  het <- vapply(g, is_het, logical(1))
  data.frame(genotype_types = length(unique(g)),
             genotype_counts = length(g),
             het_counts      = sum(het),
             Heterozygosity  = sum(het) / length(g))
}
genotype_types_df <- do.call(rbind, apply(GT, 1, genotype_summary))
rownames(genotype_types_df) <- rownames(GT)

polymorphism_group <- data.frame(
  locus = rownames(genotype_types_df),
  group = ifelse(genotype_types_df$genotype_types != 1, "polymorphic", "monomorphic"))

# =============================================================================
# 2. STR expression (locus-level mean DP)
# =============================================================================
# DP = mean read depth across the FULL cohort, non-genotyped samples = 0.
# This is locus- and cohort-level (NOT allele-resolved) and blends expression
# with detectability; see Methods note and the new allelic-balance analysis.
DP <- DP[rownames(DP) %in% locus_with_sufficient_gr, ]
stopifnot(nrow(DP) == nrow(GT))
DP[is.na(DP)] <- 0
dp_num <- apply(DP, 2, as.numeric)
rownames(dp_num) <- rownames(DP)
DP <- data.frame(locus = rownames(DP), dp_num, check.names = FALSE)

expr_df <- data.frame(
  locus  = DP$locus,
  meanDP = rowMeans(DP[, -1]),
  group  = polymorphism_group$group[match(DP$locus, polymorphism_group$locus)])
expr_df$period_size <- config[match(expr_df$locus, config[, 6]), 4]
expr_df$reference_allele_length <-
  config[match(expr_df$locus, config[, 6]), "reference_allele_length"]
coords <- config[match(expr_df$locus, config[, 6]), c(1, 2, 3, 6)]
colnames(coords) <- c("chr", "start", "end", "locus")
expr_df <- left_join(expr_df, coords, by = "locus")
expr_df <- expr_df[order(expr_df$meanDP, decreasing = TRUE), ]

write.table(expr_df, "locus_list.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(expr_df$locus[expr_df$group == "polymorphic"], "locus_list.polymorphic.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

# =============================================================================
# 3. Figures
# =============================================================================
fig1 <- polymorphism_group
fig1$genotyping_rate <- gr_df$genotyping_rate[match(fig1$locus, rownames(gr_df))]
fig1$Heterozygosity  <- genotype_types_df$Heterozygosity[match(fig1$locus, rownames(genotype_types_df))]
fig1$meanDP          <- expr_df$meanDP[match(fig1$locus, expr_df$locus)]
fig1$period_size     <- config[match(fig1$locus, config[, 6]), 4]
fig1$reference_allele_length <- config[match(fig1$locus, config[, 6]), "reference_allele_length"]
fig1$period_f <- fct_period(fig1$period_size)
fig1$group_f  <- fct_group(fig1$group)

## ---- Figure S1 : chromosomal distribution (ideogram input) -----------------
idio_block <- function(loci, colour, shape) {
  m <- match(loci, config[, 6])
  data.frame(group = "mapping", font_size = ".", color = colour, locus = ".",
             shape = shape, chr = paste0("chr", config[m, 1]),
             start = config[m, 2], end = config[m, 3], strand = ".")
}
idio <- rbind(
  idio_block(polymorphism_group$locus[polymorphism_group$group == "polymorphic"],  "144,0,33", "filled_box"),
  idio_block(polymorphism_group$locus[polymorphism_group$group == "monomorphic"], "0,47,167", "filled_circle"))
write.table(idio, "idio.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

## ---- Figure 1a : cumulative distribution of heterozygosity (pSTR) ----------
grid_h <- seq(0.1, 1, 0.1)
ho_cfl_df <- do.call(rbind, lapply(period_levels, function(ms) {
  sub <- fig1[fig1$period_size == ms & fig1$group == "polymorphic", ]
  rbind(data.frame(motif_length = ms, heterozygosity = 0, cumulative_fraction = 0),
        data.frame(motif_length = ms, heterozygosity = grid_h,
                   cumulative_fraction = vapply(grid_h, function(p)
                     1 - sum(sub$Heterozygosity > p) / nrow(sub), numeric(1))))
}))
ho_cfl_df$motif_length <- factor(ho_cfl_df$motif_length, levels = period_levels)
ho_cfl <- ggplot(ho_cfl_df, aes(heterozygosity, cumulative_fraction,
                                group = motif_length, colour = motif_length)) +
  geom_line(linewidth = 1) + scale_color_lancet() + base_theme +
  labs(x = "Heterozygosity", y = "Cumulative Fraction")
ggsave("ho_cfl.pdf", ho_cfl, width = 8, height = 5)

## ---- Figure 1b : cumulative distribution of genotyping rate ----------------
grid_g <- seq(0.5, 1, 0.1)
gr_cfl_df <- do.call(rbind, lapply(period_levels, function(ms) {
  do.call(rbind, lapply(group_levels, function(grp) {
    sub <- fig1[fig1$period_size == ms & fig1$group == grp, ]
    data.frame(motif_length = ms, group = grp, genotyping_rate = grid_g,
               cumulative_fraction = vapply(grid_g, function(p)
                 1 - sum(sub$genotyping_rate > p) / nrow(sub), numeric(1)))
  }))
}))
gr_cfl_df$motif_length <- factor(gr_cfl_df$motif_length, levels = period_levels)
gr_cfl_df$group <- fct_group(gr_cfl_df$group)
gr_cfl <- ggplot(gr_cfl_df, aes(genotyping_rate, cumulative_fraction)) +
  geom_line(aes(colour = motif_length, linetype = group), linewidth = 1) +
  scale_color_lancet() + base_theme +
  labs(x = "Genotyping Rate", y = "Cumulative Fraction")
ggsave("gr_cfl.pdf", gr_cfl, width = 8, height = 5)

## ---- Figure 1c : mean DP by motif length x reference-allele-length bin ------
brks <- c(0, 20, 50, 100)
bin_labels <- c(sprintf("(%d,%d]", brks[-length(brks)], brks[-1]), sprintf(">%d", brks[length(brks)]))
dp_bin_df <- do.call(rbind, lapply(period_levels, function(ms) {
  do.call(rbind, lapply(group_levels, function(grp) {
    sub <- fig1[fig1$period_size == ms & fig1$group == grp, ]
    do.call(rbind, lapply(seq_along(brks), function(k) {
      if (k < length(brks)) {
        sel <- sub$reference_allele_length > brks[k] & sub$reference_allele_length <= brks[k + 1]
        lab <- sprintf("(%d,%d]", brks[k], brks[k + 1])
      } else {
        sel <- sub$reference_allele_length > brks[k]; lab <- sprintf(">%d", brks[k])
      }
      data.frame(group = grp, motif_length = period_labels[ms - 1],
                 reference_allele_length = lab,
                 n_loci = sum(sel), meanDP = mean(sub$meanDP[sel]))
    }))
  }))
}))
# data behind Figure 1c (so each Fig-1c statement is checkable)
write.table(dp_bin_df, "fig1c_dp_by_bin.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

dp_bin_df$label  <- ifelse(is.na(dp_bin_df$meanDP), "NA", round(dp_bin_df$meanDP, 2))
dp_bin_df$motif_length <- factor(dp_bin_df$motif_length, levels = period_labels)
dp_bin_df$group  <- fct_group(dp_bin_df$group)
dp_bin_df$reference_allele_length <- factor(dp_bin_df$reference_allele_length, levels = bin_labels)
heatmap_summary <- ggplot(dp_bin_df, aes(motif_length, reference_allele_length, fill = meanDP)) +
  facet_wrap(~group, nrow = 2) +
  geom_tile(colour = "gray") +
  geom_text(aes(label = label), colour = "black", size = 3) +
  scale_fill_gradient2(na.value = "white", low = "#003366", mid = "#ffffff", high = "#990033") +
  coord_fixed() + base_theme +
  labs(x = "Motif Length", y = "Reference Allele Length")
ggsave("heatmap_summary.pdf", heatmap_summary, width = 8, height = 5)

## ---- Figure 1d : genotyping rate, pSTR vs mSTR by motif length -------------
genotyping_rate_boxplot <- ggplot(fig1, aes(group_f, genotyping_rate, fill = group_f)) +
  geom_boxplot(position = position_dodge(0.9)) +
  facet_grid(~period_f) +
  stat_compare_means(method = "wilcox.test", comparisons = list(group_levels), label = "p.format") +
  scale_fill_lancet() + base_theme +
  labs(title = "Genotyping Rate (pSTR vs mSTR)", x = NULL, y = "Genotyping Rate") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
ggsave("genotyping_rate_boxplot.pdf", genotyping_rate_boxplot, width = 15, height = 5)

## ---- Figure 1e : log2 DP, pSTR vs mSTR by motif length ---------------------
expression_boxplot <- ggplot(fig1, aes(group_f, log2(meanDP + 1), fill = group_f)) +
  geom_boxplot(position = position_dodge(0.9)) +
  facet_grid(~period_f) +
  stat_compare_means(method = "wilcox.test", comparisons = list(group_levels), label = "p.format") +
  scale_fill_lancet() + base_theme +
  labs(title = "log2 mean DP (pSTR vs mSTR)", x = NULL, y = "log2(mean DP + 1)") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
ggsave("expression_boxplot.pdf", expression_boxplot, width = 15, height = 5)

## ---- Figure S3a/S3b : GR / log2 DP vs reference allele length ---------------
# n is large (~24k / 20k): report r as the effect; p is near-zero by construction.
genotyping_rate_corplot <- ggplot(fig1, aes(reference_allele_length, genotyping_rate)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#2D6DB1", fill = "#756bb1") +
  facet_grid(~group_f, scales = "free") + stat_cor(method = "pearson") + base_theme +
  labs(x = "Reference Allele Length", y = "Genotyping Rate")
ggsave("genotyping_rate_corplot.pdf", genotyping_rate_corplot, width = 10, height = 5)

expression_corplot <- ggplot(fig1, aes(reference_allele_length, log2(meanDP + 1))) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#2D6DB1", fill = "#756bb1") +
  facet_grid(~group_f, scales = "free") + stat_cor(method = "pearson") + base_theme +
  labs(x = "Reference Allele Length", y = "log2(mean DP + 1)")
ggsave("expression_corplot.pdf", expression_corplot, width = 10, height = 5)

## ---- Figure S2 : DP heatmap (inter- vs intra-locus variation) --------------
dp_mat <- function(grp) {
  loci <- expr_df$locus[expr_df$group == grp]
  log2(t(as.matrix(DP[match(loci, DP$locus), -1])) + 1)
}
pheatmap(dp_mat("polymorphic"),
         show_rownames = TRUE, show_colnames = FALSE,
         cluster_rows = FALSE, cluster_cols = FALSE,
         fontsize_row = 7.5, angle_row = 90,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         filename = "pheatmap_polymorphic.pdf", width = 15, height = 15)

# Optional monomorphic panel: text says "a heatmap" (singular) and pSTRs are the
# focus, so off by default. Set TRUE only if Figure S2 was a two-panel figure.
make_mono_heatmap <- TRUE   # Fig S2 is a two-panel (pSTR + mSTR) figure
if (make_mono_heatmap) {
  pheatmap(dp_mat("monomorphic"),
           show_rownames = TRUE, show_colnames = FALSE,
           cluster_rows = FALSE, cluster_cols = FALSE,
           fontsize_row = 7.5, angle_row = 90,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           filename = "pheatmap_monomorphic.pdf", width = 15, height = 15)
}

# =============================================================================
# 4. Revision values  ->  written to files for the writing list
# =============================================================================
## per-chromosome pSTR / mSTR counts  (Section 2.1; Figure S1) ----------------
chrom <- config[match(polymorphism_group$locus, config[, 6]), 1]
chr_counts <- as.data.frame.matrix(table(chrom, polymorphism_group$group))
chr_counts <- data.frame(chromosome = rownames(chr_counts), chr_counts,
                         check.names = FALSE, row.names = NULL)
chr_counts <- chr_counts[, c("chromosome", "polymorphic", "monomorphic")]
colnames(chr_counts) <- c("chromosome", "pSTR", "mSTR")
chr_order  <- c(1:22, "X", "Y")
chr_counts <- chr_counts[match(chr_order, chr_counts$chromosome), ]
chr_counts <- chr_counts[!is.na(chr_counts$chromosome), ]
write.table(chr_counts, "chromosome_counts.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

## exact Wilcoxon p behind Figures 1d / 1e  (overall + per motif) -------------
fig1de_p <- do.call(rbind, lapply(c("Overall", period_labels), function(lbl) {
  s <- if (lbl == "Overall") fig1 else fig1[fig1$period_size == match(lbl, period_labels) + 1, ]
  data.frame(set = lbl,
             p_genotyping_rate = suppressWarnings(wilcox.test(genotyping_rate ~ group, data = s)$p.value),
             p_log2DP          = suppressWarnings(wilcox.test(log2(meanDP + 1) ~ group, data = s)$p.value))
}))
write.table(fig1de_p, "fig1de_test_pvalues.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

## QC(a): monomorphic loci that actually carry >=2 alleles --------------------
mono_loci <- polymorphism_group$locus[polymorphism_group$group == "monomorphic"]
mono_n_alleles <- vapply(mono_loci, function(L)
  length(unique(unlist(strsplit(na.omit(GT[L, ]), "\\|")))), integer(1))
mono_multi <- data.frame(locus = mono_loci[mono_n_alleles >= 2])
write.table(mono_multi, "mono_multiallele_loci.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

## QC(b): heterozygosity threshold '>' vs '>=' 0.1 ----------------------------
H <- fig1$Heterozygosity[fig1$group == "polymorphic"]
het_thr <- do.call(rbind, c(
  list(data.frame(group = "pSTR_total", H_gt0.1 = 100 * mean(H > 0.1), H_ge0.1 = 100 * mean(H >= 0.1))),
  lapply(period_levels, function(ms) {
    h <- fig1$Heterozygosity[fig1$group == "polymorphic" & fig1$period_size == ms]
    data.frame(group = period_labels[ms - 1], H_gt0.1 = 100 * mean(h > 0.1), H_ge0.1 = 100 * mean(h >= 0.1))
  })))

## QC(c): DP ~ reference-length correlation, with-zeros vs expressed-only ------
expressed_meanDP <- rowSums(DP[, -1]) / pmax(rowSums(DP[, -1] > 0), 1)
chk <- data.frame(len = fig1$reference_allele_length, grp = fig1$group,
                  dp0 = log2(fig1$meanDP + 1),
                  dpE = log2(expressed_meanDP[match(fig1$locus, DP$locus)] + 1))
dp_len_r <- do.call(rbind, lapply(group_levels, function(grp)
  data.frame(group = grp,
             r_DP_with0     = cor(chk$len[chk$grp == grp], chk$dp0[chk$grp == grp]),
             r_DP_expressed = cor(chk$len[chk$grp == grp], chk$dpE[chk$grp == grp]))))

## ---- master reference file --------------------------------------------------
ung <- nrow(loci_no_expression)
rv <- character(0); add <- function(...) rv[[length(rv) + 1]] <<- sprintf(...)
add("# expression.R revision values  (auto-generated %s)", as.character(Sys.Date()))
add("")
add("## [HEADLINE COUNTS]  (Section 2.1)")
add("targeted_loci\t%d\t(GRCh38.hipstr_reference.refine.bed, nrow)", nrow(config))
add("genotyped_GR>=0.5\t%d\t(locus_list.txt, nrow)", nrow(polymorphism_group))
add("pSTR\t%d\t(locus_list.txt, group==polymorphic)", sum(polymorphism_group$group == "polymorphic"))
add("mSTR\t%d\t(locus_list.txt, group==monomorphic)", sum(polymorphism_group$group == "monomorphic"))
add("ungenotyped_GR0\t%d\t(%.2f%%)\t(Loci_withnoexpression.txt, nrow)", ung, 100 * ung / nrow(config))
add("# per-chromosome ranges & X/Y -> chromosome_counts.txt")
add("")
add("## [HETEROZYGOSITY THRESHOLD]  (Section 2.2, Fig 1a) -- unify manuscript text to '>= 0.1'")
add("group\tH>0.1(pct)\tH>=0.1(pct)")
for (i in seq_len(nrow(het_thr)))
  add("%s\t%.2f\t%.2f", het_thr$group[i], het_thr$H_gt0.1[i], het_thr$H_ge0.1[i])
add("")
add("## [DP ~ REFERENCE LENGTH r]  (Section 2.2, Fig S3b) -- report r, soften causal language")
add("group\tr_DP_with0\tr_DP_expressed")
for (i in seq_len(nrow(dp_len_r)))
  add("%s\t%+.3f\t%+.3f", dp_len_r$group[i], dp_len_r$r_DP_with0[i], dp_len_r$r_DP_expressed[i])
add("")
add("## [FIG 1d/1e TEST]  exact two-sided Wilcoxon p -> fig1de_test_pvalues.txt")
add("")
add("## [DEFINITION QC]  (Section 2.1 Methods)")
add("monomorphic_loci_with_>=2_alleles\t%d / %d\t(mono_multiallele_loci.txt)",
    sum(mono_n_alleles >= 2), length(mono_loci))
add("")
add("## [R2b: DP & GENOTYPING RATE BY DEGREE OF POLYMORPHISM]  (Ho bins, not just mSTR/pSTR) -> dp_quality_by_polymorphism_bin.txt")
poly_bin <- cut(fig1$Heterozygosity, c(-Inf, 0, 0.1, 0.2, 0.3, 0.5, Inf),
                labels = c("Ho=0", "(0,0.1]", "(0.1,0.2]", "(0.2,0.3]", "(0.3,0.5]", ">0.5"))
poly_strat <- do.call(rbind, lapply(levels(poly_bin), function(L) {
  s <- fig1[poly_bin == L, ]
  if (!nrow(s)) return(NULL)
  data.frame(Ho_bin = L, n_loci = nrow(s),
             mean_genotyping_rate = mean(s$genotyping_rate),
             mean_DP = mean(s$meanDP), median_DP = median(s$meanDP))
}))
write.table(poly_strat, "dp_quality_by_polymorphism_bin.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
sp_dp <- suppressWarnings(cor.test(fig1$Heterozygosity, fig1$meanDP,          method = "spearman"))
sp_gr <- suppressWarnings(cor.test(fig1$Heterozygosity, fig1$genotyping_rate, method = "spearman"))
kw_dp <- kruskal.test(meanDP ~ Ho_bin, data = data.frame(meanDP = fig1$meanDP, Ho_bin = poly_bin))
kw_gr <- kruskal.test(genotyping_rate ~ Ho_bin, data = data.frame(genotyping_rate = fig1$genotyping_rate, Ho_bin = poly_bin))
ggsave("dp_quality_by_polymorphism_bin.pdf",
       ggplot(data.frame(Ho_bin = poly_bin, meanDP = fig1$meanDP), aes(Ho_bin, log2(meanDP + 1), fill = Ho_bin)) +
         geom_boxplot(outlier.size = 0.3) + theme_bw() +
         theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1, face = "bold")) +
         labs(title = "Read depth by degree of polymorphism", x = "Heterozygosity bin", y = "log2(meanDP + 1)"),
       width = 9, height = 5)
add("Ho_bin\tn\tmean_GR\tmean_DP\tmedian_DP")
for (i in seq_len(nrow(poly_strat)))
  add("%s\t%d\t%.3f\t%.2f\t%.2f", poly_strat$Ho_bin[i], poly_strat$n_loci[i],
      poly_strat$mean_genotyping_rate[i], poly_strat$mean_DP[i], poly_strat$median_DP[i])
add("Spearman Ho~DP\trho=%.3f\tp=%.2e", unname(sp_dp$estimate), sp_dp$p.value)
add("Spearman Ho~GR\trho=%.3f\tp=%.2e", unname(sp_gr$estimate), sp_gr$p.value)
add("Kruskal-Wallis DP across bins\tp=%.2e", kw_dp$p.value)
add("Kruskal-Wallis GR across bins\tp=%.2e", kw_gr$p.value)
writeLines(rv, "revision_values.txt")

cat("Revision files written: revision_values.txt, chromosome_counts.txt,",
    "fig1c_dp_by_bin.txt, fig1de_test_pvalues.txt, mono_multiallele_loci.txt,",
    "dp_quality_by_polymorphism_bin.txt\n")



####重做图：
# =============================================================================
# Figure 1 (JGG layout): A = genotyping rate, B = log2 DP, pSTR vs mSTR
# Two-panel composite, vector PDF, journal styling.
# Requires: patchwork. install.packages("patchwork") if absent.
# =============================================================================
library(patchwork)
library(ggplot2)

## ---- JGG-style theme --------------------------------------------------------
# White background, no grid, thin black axes, sans-serif, sizes tuned for a
# ~half-to-full-width figure at final print size (>= 7 pt).
theme_jgg <- theme_classic(base_family = "Arial", base_size = 9) +
  theme(
    axis.title   = element_text(colour = "black", size = 9),
    axis.text    = element_text(colour = "black", size = 8),
    axis.line    = element_line(colour = "black", linewidth = 0.4),
    axis.ticks   = element_line(colour = "black", linewidth = 0.4),
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text   = element_text(colour = "black", size = 8, margin = margin(2,0,2,0)),
    legend.position = "none",
    plot.margin  = margin(4, 6, 2, 4),
    plot.tag     = element_text(face = "bold", size = 11))

## Fixed pSTR/mSTR colours, reused everywhere (colour-blind-safe blue/orange).
grp_cols <- c(polymorphic = "#2D6DB1", monomorphic = "#E69F00")

## Format the Wilcoxon p-value label as P < 2e-16 style (italic P, superscript).
## stat_compare_means prints ggpubr's default; we override to a clean label by
## precomputing per-facet p and placing it with geom_text.
pval_lab <- function(p) {
  ifelse(p < 2.2e-16, "italic(P) < 2 %*% 10^-16",
         sprintf("italic(P) == %.1e", p))     # fallback; rarely triggered here
}
pdat <- do.call(rbind, lapply(period_labels, function(lbl) {
  s <- fig1[fig1$period_f == lbl, ]
  data.frame(period_f = factor(lbl, levels = period_labels),
             p_gr = wilcox.test(genotyping_rate ~ group, s)$p.value,
             p_dp = wilcox.test(log2(meanDP + 1) ~ group, s)$p.value)
}))
pdat$lab_gr <- pval_lab(pdat$p_gr)
pdat$lab_dp <- pval_lab(pdat$p_dp)

y_gr <- max(fig1$genotyping_rate) * 1.02
y_dp <- max(log2(fig1$meanDP + 1)) * 1.02

## ---- Panel A : genotyping rate ---------------------------------------------
pA <- ggplot(fig1, aes(group_f, genotyping_rate, fill = group_f)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.3, width = 0.65) +
  facet_grid(~ period_f) +
  geom_text(data = pdat, aes(x = 1.5, y = y_gr, label = lab_gr),
            parse = TRUE, size = 2.4, inherit.aes = FALSE) +
  scale_fill_manual(values = grp_cols) +
  scale_x_discrete(labels = c(polymorphic = "pSTR", monomorphic = "mSTR")) +
  labs(x = NULL, y = "Genotyping rate") +
  theme_jgg +
  theme(axis.text.x = element_text(angle = 0))

## ---- Panel B : log2 DP ------------------------------------------------------
pB <- ggplot(fig1, aes(group_f, log2(meanDP + 1), fill = group_f)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.3, width = 0.65) +
  facet_grid(~ period_f) +
  geom_text(data = pdat, aes(x = 1.5, y = y_dp, label = lab_dp),
            parse = TRUE, size = 2.4, inherit.aes = FALSE) +
  scale_fill_manual(values = grp_cols) +
  scale_x_discrete(labels = c(polymorphic = "pSTR", monomorphic = "mSTR")) +
  labs(x = NULL, y = expression(log[2]*"(mean DP + 1)")) +
  theme_jgg +
  theme(axis.text.x = element_text(angle = 0))

## ---- Compose A over B -------------------------------------------------------
fig1_AB <- pA / pB + plot_annotation(tag_levels = "A")

## Export as vector PDF. Width ~ 180 mm (double column) -> 7.1 in; height to taste.
ggsave("Figure1.pdf", fig1_AB, width = 7.1, height = 6.2, device = cairo_pdf)


source("figS1_S2_S4_patch.R")



## ---- Figure S3 : DP heatmap (inter- vs intra-locus variation), A/B panels ---
## pheatmap returns a grid object ($gtable), not a ggplot. We capture each panel
## with silent=TRUE, draw ONE shared legend, and assemble with cowplot.
library(cowplot)
library(grid)        # for editing/Extracting gtable grobs
library(pheatmap)
library(ggplot2)
dp_mat <- function(grp) {
  loci <- expr_df$locus[expr_df$group == grp]
  log2(t(as.matrix(DP[match(loci, DP$locus), -1])) + 1)
}

## ---- ONE toggle: layout direction ------------------------------------------
## "v" = stacked vertically (A on top, B below)  -> better for wide heatmaps
## "h" = side by side (A left, B right)
S3_LAYOUT <- "v"

## ---- shared colour scale across BOTH panels (required for fair comparison) --
.rng  <- range(c(dp_mat("polymorphic"), dp_mat("monomorphic")), na.rm = TRUE)
.brks <- seq(.rng[1], .rng[2], length.out = 51)
.cols <- colorRampPalette(c("navy", "white", "firebrick3"))(50)

## build one pheatmap gtable; legend kept only for the panel we harvest it from
ph_gtable <- function(grp, with_legend) {
  pheatmap(dp_mat(grp),
           show_rownames = TRUE, show_colnames = FALSE,
           cluster_rows = FALSE, cluster_cols = FALSE,
           fontsize_row = 7.5, angle_row = 90,
           color = .cols, breaks = .brks,
           legend = with_legend, silent = TRUE)$gtable
}

## panels WITHOUT legends (so neither carries its own colour bar)
g_poly <- ph_gtable("polymorphic", with_legend = FALSE)
g_mono <- ph_gtable("monomorphic", with_legend = FALSE)

## a full pheatmap WITH legend, from which we extract just the legend grob
g_leg_full <- ph_gtable("polymorphic", with_legend = TRUE)
leg_grob   <- gtable::gtable_filter(g_leg_full, "legend")   # the colour bar only

## assemble the two heatmaps, tagged A / B
panels <- plot_grid(
  g_poly, g_mono,
  labels = c("A", "B"), label_size = 18,
  ncol = if (S3_LAYOUT == "h") 2 else 1,
  rel_widths  = if (S3_LAYOUT == "h") c(1, 1) else 1,
  rel_heights = if (S3_LAYOUT == "h") 1 else c(1, 1)
)

## attach the single shared legend on the right
figS3 <- plot_grid(panels, leg_grob, ncol = 2, rel_widths = c(1, 0.08))

## size: wide heatmaps -> tall when stacked, wide when side by side
.w <- if (S3_LAYOUT == "h") 26 else 15
.h <- if (S3_LAYOUT == "h") 14 else 24
ggsave("figS3_dp_heatmap_AB.pdf", figS3, width = .w, height = .h, limitsize = FALSE)

library(ragg)   # install.packages("ragg") if needed
agg_png("figS3_dp_heatmap_AB.png",
        width = .w, height = .h, units = "in", res = 200)
print(figS3)      # 或 grid::grid.draw(figS3)
dev.off()