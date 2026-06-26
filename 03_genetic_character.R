# =============================================================================
# genetic_character.R  —  clean, revision build
# Manuscript: Section 2.3 ; Figures 2a-2d, S4, S5, S6
# Inputs: pSTR vcf, locus_list.polymorphic.txt (script 1),
#         STR_annotation_group.txt (script 2, with reconciled biotype +
#         gene_id + intergenic_subtype), config.
#
# Changes vs original:
#   * FIX: gene_id was set to biotype by mistake -> now reads the real gene_id
#          column produced by the revised annotation.R.
#   * FIX: facet_wrap(~biotype, nrow=) empty argument -> nrow = 1.
#   * CLEAN: the `ral_distribution_density2 <- ral_distribution_density <- ...`
#          double assignment; the vestigial "NC" factor level; the per-iteration
#          print() in the motif loop.
#   * Reconciled biotype flows in automatically -> Figure S4 headline split is
#     unchanged (mRNA 92.08% / lncRNA 3.91%); Figure S5 is now consistent
#     (lncRNA cannot fall in CDS/UTR). Re-run regenerates S5.
#   * Charts not referenced in the text (intergenic tetra/penta motif panels,
#     biotype pie) are gated behind make_optional_figs.
#   * Section 4 writes traceability files for the §2.3 numbers.
# =============================================================================

rm(list = ls())

library(readxl)
library(tidyr)
library(vcfR)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggsci)

make_optional_figs <- FALSE   # TRUE to also draw intergenic motif panels + biotype pie

base_theme <- theme_bw() +
  theme(axis.title   = element_text(face = "bold", colour = "black", size = 15),
        axis.title.x = element_text(face = "bold", colour = "black", size = 12),
        axis.title.y = element_text(face = "bold", colour = "black", size = 12),
        axis.text.x  = element_text(face = "bold", colour = "black", size = 10),
        axis.text.y  = element_text(face = "bold", colour = "black", size = 10))

region_levels <- c("CDS", "5_UTR", "3_UTR", "intron", "intergenic")
period_levels <- 2:6
period_labels <- c("Di", "Tri", "Tetra", "Penta", "Hexa")
fct_region <- function(x) factor(x, levels = region_levels)
fct_period <- function(x) factor(x, levels = period_levels, labels = period_labels)

# =============================================================================
# 0. Inputs
# =============================================================================
pSTR_list  <- read.table("../1_expression/locus_list.polymorphic.txt")
pSTR_group <- read.table("../2_annotation/STR_annotation_group.txt", header = TRUE,
                         sep = "\t", quote = "", stringsAsFactors = FALSE)
config <- read.table("../GRCh38.hipstr_reference.refine.bed")
config <- config[config[, 6] %in% pSTR_list[, 1], ]
config$reference_allele_length <- config[, 3] - config[, 2] + 1

m <- match(pSTR_list[, 1], pSTR_group$locus)
mc <- match(pSTR_list[, 1], config[, 6])

# =============================================================================
# 1. genetic_character_df  (gene_id FIX)
# =============================================================================
genetic_character_df <- data.frame(
  locus                   = pSTR_list[, 1],
  gene_id                 = pSTR_group$gene_id[m],      # FIX: was pSTR_group$biotype
  feature                 = pSTR_group$final_group[m],
  biotype                 = pSTR_group$biotype[m],      # reconciled in annotation.R
  intergenic_subtype      = pSTR_group$intergenic_subtype[m],
  period_size             = config[mc, 4],
  ncopy                   = config[mc, 5],
  reference_allele_length = config$reference_allele_length[mc],
  motif                   = config[mc, 7],
  stringsAsFactors = FALSE)
write.table(genetic_character_df, "genetic_character_df.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# =============================================================================
# 2. Biotype composition  (Figure S4)
# =============================================================================
biotype_df <- genetic_character_df$biotype
biotype_df[biotype_df == "protein_coding"] <- "mRNA (protein_coding)"
biotype_df[is.na(biotype_df)] <- "NoAnnotation"
biotype_df <- as.data.frame(table(biotype_df))
colnames(biotype_df) <- c("biotype", "counts")
biotype_df$percentage <- biotype_df$counts / sum(biotype_df$counts)
biotype_df$label <- paste0(biotype_df$counts, " (", round(biotype_df$percentage, 4) * 100, "%)")
biotype_df <- biotype_df[order(biotype_df$percentage, decreasing = TRUE), ]
biotype_df$biotype <- factor(biotype_df$biotype, levels = as.character(biotype_df$biotype))
biotype_barplot <- ggplot(biotype_df, aes(biotype, counts, fill = biotype)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = label), size = 3.5, vjust = -0.2) +
  base_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(title = "Sourced transcript biotype of RNA-pSTRs", x = "Biotype", y = "pSTR Counts")
ggsave("biotype_barplot.pdf", biotype_barplot, height = 6, width = 10)

# =============================================================================
# 3. pSTR composition by region x motif  (Figure 2a; Figure S5 by biotype)
# =============================================================================
gc <- genetic_character_df
gc$region <- fct_region(gc$feature)
gc$period <- fct_period(gc$period_size)

build_counts <- function(sub) {
  regions <- intersect(region_levels, as.character(unique(sub$region)))
  do.call(rbind, lapply(regions, function(rg) {
    pc <- as.data.frame(table(period = droplevels(sub$period[sub$region == rg])))
    data.frame(region = rg, period = pc$period, count = pc$Freq,
               pct_within_region = pc$Freq / sum(pc$Freq))
  }))
}
build_label <- function(sub) {
  regions <- intersect(region_levels, as.character(unique(sub$region)))
  rs <- vapply(regions, function(rg) sum(sub$region == rg), numeric(1))
  data.frame(region = regions, region_sum = rs, pct_within_biotype = rs / sum(rs))
}

scopes <- list(total = gc,
               protein_coding = gc[gc$biotype %in% "protein_coding", ],
               lncRNA = gc[gc$biotype %in% "lncRNA", ])
counts_df <- do.call(rbind, Map(function(b, s) data.frame(biotype = b, build_counts(s)),
                                names(scopes), scopes))
label_df  <- do.call(rbind, Map(function(b, s) data.frame(biotype = b, build_label(s)),
                                names(scopes), scopes))
rownames(counts_df) <- NULL; rownames(label_df) <- NULL
counts_df$region <- fct_region(counts_df$region)
label_df$region  <- fct_region(label_df$region)
label_df$label   <- paste0(label_df$region_sum, " (", round(label_df$pct_within_biotype, 4) * 100, "%)")

composition_plot <- function(bt, ttl) {
  ggplot(counts_df[counts_df$biotype == bt, ], aes(region, pct_within_region, fill = period)) +
    geom_bar(stat = "identity", position = "fill", colour = "black") +
    geom_text(data = label_df[label_df$biotype == bt, ],
              aes(region, 1, label = label), inherit.aes = FALSE, vjust = -0.2) +
    scale_fill_lancet() + scale_y_continuous(expand = c(0.01, 0.05)) +
    base_theme + labs(title = ttl, x = "Features", y = "Percentage")
}
ggsave("pSTR_counts_barplot.pdf",        composition_plot("total",          "Composition of pSTRs"),                  height = 6, width = 8)
ggsave("pSTR_counts_barplot_mRNA.pdf",   composition_plot("protein_coding", "Composition of pSTRs sourced from mRNA"),  height = 6, width = 8)
ggsave("pSTR_counts_barplot_lncRNA.pdf", composition_plot("lncRNA",         "Composition of pSTRs sourced from lncRNA"), height = 6, width = 8)

# =============================================================================
# 4. Reference allele length & copy number distributions
#    per-region (Figure S6) + overall (Figures 2b/2c)
# =============================================================================
ral_df <- rbind(gc[, c("feature", "period_size", "reference_allele_length")],
                data.frame(feature = "Total", gc[, c("period_size", "reference_allele_length")]))
ral_df$feature <- factor(ral_df$feature, levels = c(region_levels, "Total"))
ral_df$period_size <- fct_period(ral_df$period_size)
ral_df <- ral_df[ral_df$reference_allele_length <= 100, ]   # display threshold

ncopy_df <- rbind(gc[, c("feature", "period_size", "ncopy")],
                  data.frame(feature = "Total", gc[, c("period_size", "ncopy")]))
ncopy_df$feature <- factor(ncopy_df$feature, levels = c(region_levels, "Total"))
ncopy_df$period_size <- fct_period(ncopy_df$period_size)
ncopy_df <- ncopy_df[ncopy_df$ncopy <= 50, ]                # display threshold

## per-region (Figure S6) ----------------------------------------------------
ral_distribution_density <- ggplot(ral_df[ral_df$feature != "Total", ],
                                   aes(reference_allele_length, colour = period_size)) +
  geom_density(linewidth = 1) + scale_color_lancet() +
  scale_x_continuous(limits = c(0, 100)) + scale_y_continuous(limits = c(0, 0.2)) +
  facet_wrap(~feature, nrow = 1) + base_theme +
  labs(title = "Reference allele length of RNA-pSTRs", x = "Reference Allele Length", y = "Density")
ggsave("ral_distribution_density.pdf", ral_distribution_density, height = 5, width = 15)

ncopy_distribution_density <- ggplot(ncopy_df[ncopy_df$feature != "Total", ],
                                     aes(ncopy, colour = period_size)) +
  geom_density(linewidth = 1) + scale_color_lancet() +
  scale_x_continuous(limits = c(0, 50)) + scale_y_continuous(limits = c(0, 0.6)) +
  facet_wrap(~feature, nrow = 1) + base_theme +
  labs(title = "Reference allele copy number of RNA-pSTRs", x = "Reference Allele Copy Number", y = "Density")
ggsave("ncopy_distribution_density.pdf", ncopy_distribution_density, height = 5, width = 15)

## overall (Figures 2b/2c) ---------------------------------------------------
ral_distribution_total <- ggplot(ral_df[ral_df$feature == "Total", ],
                                 aes(reference_allele_length, fill = period_size)) +
  geom_density(alpha = 0.4, linewidth = 0.8) + scale_fill_lancet() +
  scale_x_continuous(limits = c(0, 100)) + scale_y_continuous(limits = c(0, 0.2)) +
  base_theme + labs(title = "Reference allele length of RNA-pSTRs", x = "Reference Allele Length", y = "Density")
ggsave("ral_distribution_density2.pdf", ral_distribution_total, height = 5, width = 15)

ncopy_distribution_total <- ggplot(ncopy_df[ncopy_df$feature == "Total", ],
                                   aes(ncopy, fill = period_size)) +
  geom_density(alpha = 0.4, linewidth = 0.8) + scale_fill_lancet() +
  scale_x_continuous(limits = c(0, 50)) + scale_y_continuous(limits = c(0, 0.6)) +
  base_theme + labs(title = "Reference allele copy number of RNA-pSTRs", x = "Reference Allele Copy Number", y = "Density")
ggsave("ncopy_distribution_density2.pdf", ncopy_distribution_total, height = 5, width = 15)

# =============================================================================
# 5. Motif patterns  (Figure 2d + CDS-tri call-out)
# =============================================================================
# canonical pattern = first (sorted) of all rotations of the motif and its
# reverse complement; multi-motif ("/") entries are dropped.
motif_pattern <- data.frame(
  motif = unique(genetic_character_df$motif[!grepl("/", genetic_character_df$motif)]),
  motif_pattern = NA_character_, stringsAsFactors = FALSE)
comp <- c(A = "T", T = "A", C = "G", G = "C")
for (p in seq_len(nrow(motif_pattern))) {
  if (!is.na(motif_pattern$motif_pattern[p])) next
  raw <- strsplit(motif_pattern$motif[p], "")[[1]]
  cp  <- rev(comp[raw])
  reps <- character(0)
  for (i in seq_along(raw)) {
    rr <- c(raw[i:length(raw)], raw[seq_len(i - 1)])
    cc <- c(cp[i:length(cp)],  cp[seq_len(i - 1)])
    reps <- c(reps, paste(rr, collapse = ""), paste(cc, collapse = ""))
  }
  reps <- sort(unique(reps))
  motif_pattern$motif_pattern[motif_pattern$motif %in% reps] <- reps[1]
}
motif_df <- data.frame(locus = genetic_character_df$locus,
                       feature = genetic_character_df$feature,
                       motif = genetic_character_df$motif, stringsAsFactors = FALSE)
motif_df$pattern <- motif_pattern$motif_pattern[match(motif_df$motif, motif_pattern$motif)]
motif_df <- na.omit(motif_df)
motif_df$period_size <- nchar(motif_df$pattern)
write.table(motif_df, "motif_df.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# counts per region (within-region percentage), top 10 patterns per motif length
motif_df2 <- do.call(rbind, lapply(unique(motif_df$feature), function(g) {
  t <- as.data.frame(table(pattern = motif_df$pattern[motif_df$feature == g]))
  data.frame(group = g, pattern = as.character(t$pattern), counts = t$Freq,
             percentage = t$Freq / sum(t$Freq), stringsAsFactors = FALSE)
}))
motif_df2$period_size <- nchar(motif_df2$pattern)
keep_patterns <- unlist(lapply(period_levels, function(i) {
  pp <- unique(motif_df2$pattern[motif_df2$period_size == i])
  pp[seq_len(min(10, length(pp)))]
}))
motif_df2 <- motif_df2[motif_df2$pattern %in% keep_patterns, ]
motif_df2$group <- fct_region(motif_df2$group)
motif_df2$period_size <- fct_period(motif_df2$period_size)

motif_barplot <- ggplot(motif_df2, aes(reorder(pattern, counts, decreasing = TRUE), counts, fill = period_size)) +
  geom_bar(stat = "identity") + facet_grid(group ~ period_size, scale = "free") +
  scale_fill_lancet() + base_theme +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none") +
  labs(title = "Composition of motif patterns", x = "Motif pattern", y = "Counts")
ggsave("motif_barplot.pdf", motif_barplot, height = 12, width = 15)

CDS_tri_motif_barplot <- ggplot(motif_df2[motif_df2$group == "CDS" & motif_df2$period_size == "Tri", ],
                                aes(reorder(pattern, counts, decreasing = TRUE), counts)) +
  geom_bar(stat = "identity", fill = "#2D6DB1") + base_theme + theme(legend.position = "none") +
  labs(title = "Motif pattern composition of Tri-pSTR within CDS", x = "Motif pattern", y = "Counts")
ggsave("CDS_tri_motif_barplot.pdf", CDS_tri_motif_barplot, height = 6, width = 8)

# =============================================================================
# 6. Revision values  ->  written to files for the writing list
# =============================================================================
## region x biotype composition table (backs Fig 2a / S5) --------------------
write.table(counts_df, "region_by_biotype.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

## mean ref allele length & copy number by region x motif (backs Fig S6) -----
lr <- expand.grid(period = period_labels, region = region_levels, stringsAsFactors = FALSE)
lr$mean_ref_length <- mapply(function(p, r)
  mean(ral_df$reference_allele_length[ral_df$feature == r & ral_df$period_size == p]), lr$period, lr$region)
lr$mean_copy_number <- mapply(function(p, r)
  mean(ncopy_df$ncopy[ncopy_df$feature == r & ncopy_df$period_size == p]), lr$period, lr$region)
write.table(lr, "allele_length_repeat_by_region.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

## motif pattern composition per region x motif (backs Fig 2d call-outs) -----
mt <- as.data.frame(table(region = motif_df$feature, period_size = motif_df$period_size, pattern = motif_df$pattern))
mt <- mt[mt$Freq > 0, ]
mt$pct_within_region_period <- ave(mt$Freq, mt$region, mt$period_size, FUN = function(x) x / sum(x))
mt$pct_within_region        <- ave(mt$Freq, mt$region, FUN = function(x) x / sum(x))  # manuscript denominator
mt <- mt[order(mt$region, mt$period_size, -mt$Freq), ]
write.table(mt, "motif_pattern_counts.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

## master reference file ------------------------------------------------------
gcf <- genetic_character_df
u_f <- gcf$ncopy[gcf$ncopy <= 50]                                      # repeat-unit display filter
b_f <- gcf$reference_allele_length[gcf$reference_allele_length <= 100]  # length display filter
pct_intron <- function(bt) {
  s <- label_df[label_df$biotype == bt & label_df$region == "intron", "pct_within_biotype"]
  if (length(s)) 100 * s else NA
}
trihexa_frac <- function(rg) {
  d <- counts_df[counts_df$biotype == "total" & counts_df$region == rg, ]
  100 * sum(d$count[d$period %in% c("Tri", "Hexa")]) / sum(d$count)
}
motif_pct <- function(rg, per, pat) {   # within region x motif-length (manuscript denominator)
  v <- mt$pct_within_region_period[mt$region == rg & mt$period_size == per & mt$pattern == pat]
  if (length(v)) 100 * v else 0
}

rv <- character(0); add <- function(...) rv[[length(rv) + 1]] <<- sprintf(...)
add("# genetic_character.R revision values  (auto-generated %s)", as.character(Sys.Date()))
add("# pSTRs: %d", nrow(genetic_character_df))
add("")
add("## [BIOTYPE]  (Figure S4) -> headline should reproduce 92.08%% mRNA / 3.91%% lncRNA")
for (i in seq_len(nrow(biotype_df)))
  add("%s\t%d\t%.2f%%", as.character(biotype_df$biotype[i]), biotype_df$counts[i], 100 * biotype_df$percentage[i])
add("")
add("## [WITHIN-BIOTYPE REGION]  (Figure S5) -- re-check vs 85.09%% / 74.15%% (may move slightly)")
add("mRNA intronic\t%.2f%%", pct_intron("protein_coding"))
add("lncRNA intronic\t%.2f%%", pct_intron("lncRNA"))
add("# full region x biotype table -> region_by_biotype.txt")
add("")
add("## [TRI+HEXA ENRICHMENT]  (Section 2.3) -- vs 80.52%% (CDS), 69.15%% (5'UTR)")
add("CDS tri+hexa\t%.2f%%", trihexa_frac("CDS"))
add("5_UTR tri+hexa\t%.2f%%", trihexa_frac("5_UTR"))
add("")
add("## [ALLELE LENGTH / REPEAT]  (Figures 2b/2c)")
add("# display-filtered (units<=50; length<=100bp) -- reproduces manuscript 3-49/7.07 ; 12-99/23.19 ; 97.61%%")
add("repeat_units\trange %g-%g\tmean %.2f", min(u_f), max(u_f), mean(u_f))
add("ref_length_bp\trange %g-%g\tmean %.2f", min(b_f), max(b_f), mean(b_f))
add("ref_length<50bp\t%.2f%%", 100 * mean(b_f < 50))
for (p in period_levels) {
  pu <- gcf$ncopy[gcf$period_size == p]
  pb <- gcf$reference_allele_length[gcf$period_size == p]
  add("%s\tfiltered_units %.2f / filtered_bp %.2f\tfull_units %.2f / full_bp %.2f",
      period_labels[p - 1], mean(pu[pu <= 50]), mean(pb[pb <= 100]), mean(pu), mean(pb))
}
add("# FULL unfiltered (true tail; recommended to report these):")
add("repeat_units(full)\trange %g-%g\tmean %.2f", min(gcf$ncopy), max(gcf$ncopy), mean(gcf$ncopy))
add("ref_length_bp(full)\trange %g-%g\tmean %.2f", min(gcf$reference_allele_length), max(gcf$reference_allele_length), mean(gcf$reference_allele_length))
add("ref_length<50bp(full)\t%.2f%%", 100 * mean(gcf$reference_allele_length < 50))
add("# per region x motif means (display-filtered) -> allele_length_repeat_by_region.txt")
add("")
add("## [MOTIF]  (Figure 2d) -- vs 306 patterns; CDS-tri AGC 31.35%%/AGG 29.45%%; 5'UTR-tri CCG 65.26%%")
add("distinct_patterns\t%d", length(unique(motif_df$pattern)))
add("CDS Tri AGC\t%.2f%%", motif_pct("CDS", 3, "AGC"))
add("CDS Tri AGG\t%.2f%%", motif_pct("CDS", 3, "AGG"))
add("5_UTR Tri CCG\t%.2f%%", motif_pct("5_UTR", 3, "CCG"))
add("# full region x motif x pattern table -> motif_pattern_counts.txt")
writeLines(rv, "revision_values_genetic_character.txt")

# =============================================================================
# 7. Optional figures (gated; not referenced in the manuscript text)
# =============================================================================
if (make_optional_figs) {
  for (per in c("Tetra", "Penta")) {
    g <- ggplot(motif_df2[motif_df2$group == "intergenic" & motif_df2$period_size == per, ],
                aes(reorder(pattern, counts, decreasing = TRUE), counts)) +
      geom_bar(stat = "identity", fill = "#2D6DB1") + base_theme + theme(legend.position = "none") +
      labs(title = sprintf("Motif pattern composition of %s-pSTR within intergenic", per),
           x = "Motif pattern", y = "Counts")
    ggsave(sprintf("intergenic_%s_motif_barplot.pdf", tolower(per)), g, height = 6, width = 8)
  }
  pie_df <- label_df[label_df$biotype != "total", ]
  pie_df$label <- paste0(pie_df$region, " (", pie_df$region_sum, ")")
  pie_df$biotype <- factor(pie_df$biotype, levels = c("protein_coding", "lncRNA"))
  biotype_pieplot <- ggplot(pie_df, aes("", pct_within_biotype, fill = region)) +
    geom_col(colour = "black") +
    geom_label_repel(aes(label = label), size = 3.5, fontface = "bold", colour = "black",
                     show.legend = FALSE, position = position_stack(vjust = 0.5)) +
    facet_wrap(~biotype, nrow = 1) + coord_polar(theta = "y") + theme_void() + scale_fill_npg() +
    theme(strip.text = element_text(face = "bold", colour = "black", size = 13.5), legend.position = "none")
  ggsave("biotype_pieplot.pdf", biotype_pieplot, height = 5, width = 10)
}

cat("genetic_character.R done. gene_id fixed; reconciled biotype used; revision files written.\n")



# ============================================================================
# VERIFICATION: Fig S8 "significantly longer/shorter" — return console output
# ============================================================================
cat("\n========== Fig S8 significance (genetic_character.R) ==========\n")
gcf <- genetic_character_df
func_reg <- c("CDS","5_UTR"); bg_reg <- c("intron","intergenic")
for (per in period_labels) {
  ps <- which(period_labels==per)+1
  for (metric in c("reference_allele_length","ncopy")) {
    a <- gcf[[metric]][gcf$period_size==ps & gcf$feature %in% func_reg]
    b <- gcf[[metric]][gcf$period_size==ps & gcf$feature %in% bg_reg]
    if (length(a)>2 && length(b)>2) {
      p <- wilcox.test(a,b)$p.value
      cat(sprintf("%-5s %-22s funct(med=%.1f,n=%d) vs bg(med=%.1f,n=%d)  P=%.2e\n",
                  per, metric, median(a), length(a), median(b), length(b), p))
    }
  }
}
# ============================================================================
# Figure 2 (JGG layout): A region composition | B length | C copy number | D motifs
# Reuses data frames already built above: counts_df, label_df, ral_df, ncopy_df, motif_df2
# Requires: patchwork
# ============================================================================
library(patchwork)
library(ggplot2)
theme_jgg <- theme_classic(base_family = "Arial", base_size = 9) +
  theme(axis.title = element_text(colour="black", size=9),
        axis.text  = element_text(colour="black", size=8),
        axis.line  = element_line(colour="black", linewidth=0.4),
        axis.ticks = element_line(colour="black", linewidth=0.4),
        strip.background = element_rect(fill="grey92", colour=NA),
        strip.text = element_text(colour="black", size=8, margin=margin(2,2,2,2)),
        legend.title = element_text(size=8), legend.text = element_text(size=7),
        legend.key.size = unit(3,"mm"),
        plot.tag = element_text(face="bold", size=11))

# ============================================================================
# Shared palette (harmonised with Fig 1: Di=Fig1-blue, Tri=Fig1-orange)
# ============================================================================
motif_cols <- c(Di="#2D6DB1", Tri="#E69F00", Tetra="#2E933C",
                Penta="#56B4E9", Hexa="#B279A7")

region_lab <- c(CDS="CDS", "5_UTR"="5\u2032 UTR", "3_UTR"="3\u2032 UTR",
                intron="intron", intergenic="intergenic")

# Build a "region (count, pct)" two-line axis label for Fig 2A and S7
make_region_n_lab <- function(lab_df_sub) {
  setNames(sprintf("%s\n(%d, %.2f%%)",
                   region_lab[as.character(lab_df_sub$region)],
                   lab_df_sub$region_sum,
                   100*lab_df_sub$pct_within_biotype),
           as.character(lab_df_sub$region))
}

# ---- Fig 2A: region composition, n moved to x-axis (no crowded top labels) ----
dA <- counts_df[counts_df$biotype=="total", ]
lA <- label_df[label_df$biotype=="total", ]
axis_lab_A <- make_region_n_lab(lA)
pA <- ggplot(dA, aes(region, pct_within_region, fill=period)) +
  geom_bar(stat="identity", position="fill", colour="black", linewidth=0.25) +
  scale_fill_manual(values=motif_cols, name="Motif", drop=FALSE) +
  scale_x_discrete(labels=axis_lab_A) +
  scale_y_continuous(expand=expansion(mult=c(0.01,0.02))) +
  labs(x=NULL, y="Proportion") + theme_jgg +
  theme(axis.text.x=element_text(angle=30, hjust=1, size=7))

pB <- ggplot(ral_df[ral_df$feature=="Total", ],
             aes(reference_allele_length, fill=period_size)) +
  geom_density(alpha=0.3, linewidth=0.3) +
  scale_fill_manual(values=motif_cols, name="Motif", drop=FALSE) +
  coord_cartesian(xlim=c(0,100), ylim=c(0,0.2)) +
  labs(x="Reference-allele length (bp)", y="Density") + theme_jgg

pC <- ggplot(ncopy_df[ncopy_df$feature=="Total", ],
             aes(ncopy, fill=period_size)) +
  geom_density(alpha=0.3, linewidth=0.3) +
  scale_fill_manual(values=motif_cols, name="Motif", drop=FALSE) +
  coord_cartesian(xlim=c(0,50), ylim=c(0,0.6)) +
  labs(x="Reference-allele copy number", y="Density") + theme_jgg

dD <- motif_df2
dD$group <- factor(dD$group, levels=region_levels, labels=region_lab[region_levels])
pD <- ggplot(dD, aes(reorder(pattern, counts, decreasing=TRUE), counts, fill=period_size)) +
  geom_bar(stat="identity") +
  facet_grid(group ~ period_size, scales="free") +
  scale_fill_manual(values=motif_cols, name="Motif", drop=FALSE) +
  labs(x="Motif pattern", y="Counts") + theme_jgg +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=6),
        legend.position="none")   # D shares the collected legend

# compose; collect ONE legend on the right
top <- pA + pB + pC + plot_layout(widths=c(1,1,1))
fig2 <- (top / pD) + plot_layout(heights=c(1,1.6), guides="collect") +
  plot_annotation(tag_levels="A") &
  theme(legend.position="right")
ggsave("Figure2.pdf", fig2, width=10, height=10, device=cairo_pdf)

# ============================================================================
# Figure 2 (JGG layout): A region composition | B length | C copy number | D motifs
# Reuses data frames already built above: counts_df, label_df, ral_df, ncopy_df, motif_df2
# Requires: patchwork
# ============================================================================
library(patchwork)

theme_jgg <- theme_classic(base_family = "Arial", base_size = 9) +
  theme(axis.title = element_text(colour="black", size=9),
        axis.text  = element_text(colour="black", size=8),
        axis.line  = element_line(colour="black", linewidth=0.4),
        axis.ticks = element_line(colour="black", linewidth=0.4),
        strip.background = element_rect(fill="grey92", colour=NA),
        strip.text = element_text(colour="black", size=8, margin=margin(2,2,2,2)),
        legend.title = element_text(size=8), legend.text = element_text(size=7),
        legend.key.size = unit(3,"mm"),
        plot.tag = element_text(face="bold", size=11))

# Fixed motif colours (colour-blind-safe), reused in Fig 2 AND S6/S7/S8
motif_cols <- c(Di="#2D6DB1", Tri="#D7263D", Tetra="#2E933C", Penta="#3FA7D6", Hexa="#7B5EA7")

region_lab <- c(CDS="CDS", "5_UTR"="5\u2032 UTR", "3_UTR"="3\u2032 UTR",
                intron="intron", intergenic="intergenic")  # prime symbols

## --- A: region composition (stacked, motif-coloured) ---
dA <- counts_df[counts_df$biotype=="total", ]
lA <- label_df[label_df$biotype=="total", ]
pA <- ggplot(dA, aes(region, pct_within_region, fill=period)) +
  geom_bar(stat="identity", position="fill", colour="black", linewidth=0.25) +
  geom_text(data=lA, aes(region, 1.02, label=label), inherit.aes=FALSE,
            vjust=0, size=2.1) +
  scale_fill_manual(values=motif_cols, name="Motif") +
  scale_x_discrete(labels=region_lab) +
  scale_y_continuous(expand=expansion(mult=c(0.01,0.08))) +
  labs(x=NULL, y="Proportion") + theme_jgg +
  theme(axis.text.x=element_text(angle=30, hjust=1))

## --- B: reference-allele length (overall) ---
pB <- ggplot(ral_df[ral_df$feature=="Total", ],
             aes(reference_allele_length, fill=period_size)) +
  geom_density(alpha=0.45, linewidth=0.3) +
  scale_fill_manual(values=motif_cols, name="Motif") +
  coord_cartesian(xlim=c(0,100), ylim=c(0,0.2)) +
  labs(x="Reference-allele length (bp)", y="Density") + theme_jgg

## --- C: reference-allele copy number (overall) ---
pC <- ggplot(ncopy_df[ncopy_df$feature=="Total", ],
             aes(ncopy, fill=period_size)) +
  geom_density(alpha=0.45, linewidth=0.3) +
  scale_fill_manual(values=motif_cols, name="Motif") +
  coord_cartesian(xlim=c(0,50), ylim=c(0,0.6)) +
  labs(x="Reference-allele copy number", y="Density") + theme_jgg

## --- D: motif patterns (region x motif facet) ---
dD <- motif_df2
dD$group <- factor(dD$group, levels=region_levels, labels=region_lab[region_levels])
pD <- ggplot(dD, aes(reorder(pattern, counts, decreasing=TRUE), counts, fill=period_size)) +
  geom_bar(stat="identity") +
  facet_grid(group ~ period_size, scales="free") +
  scale_fill_manual(values=motif_cols, guide="none") +
  labs(x="Motif pattern", y="Counts") + theme_jgg +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=6))

## --- compose: (A B C) top row, D full-width bottom ---
top <- pA + pB + pC + plot_layout(widths=c(1,1,1), guides="collect")
fig2 <- (top / pD) + plot_layout(heights=c(1,1.6)) +
  plot_annotation(tag_levels="A") &
  theme(legend.position="right")

ggsave("Figure2.pdf", fig2, width=7.1, height=8.0, device=cairo_pdf)
# ============================================================================
# S5–S8 JGG restyle. Put S5 block in annotation.R, S6/S7/S8 in genetic_character.R.
# Shares theme_jgg + motif_cols defined above (re-define if running separately).
# ============================================================================

## ---- Fig S6: biotype barplot (genetic_character.R) ----
# keep all biotypes; clean labels; log10 y so the small classes are visible
bt6 <- biotype_df
bt6$biotype <- factor(bt6$biotype, levels=bt6$biotype[order(bt6$counts, decreasing=TRUE)])
figS6 <- ggplot(bt6, aes(biotype, counts)) +
  geom_col(fill="#2D6DB1", colour="black", linewidth=0.3, width=0.7) +
  geom_text(aes(label=label), vjust=-0.4, size=2.2) +
  scale_y_continuous(expand=expansion(mult=c(0,0.12))) +
  labs(x="Transcript biotype", y="pSTR count") + theme_jgg +
  theme(axis.text.x=element_text(angle=40, hjust=1, size=7))
ggsave("FigureS6.pdf", figS6, width=7, height=4.2, device=cairo_pdf)

# ---- Fig S7: stacked bars by biotype, n in axis labels per facet ----
d7 <- counts_df[counts_df$biotype %in% c("protein_coding","lncRNA"), ]
l7 <- label_df[label_df$biotype %in% c("protein_coding","lncRNA"), ]
bt_lab <- c(protein_coding="mRNA (protein-coding)", lncRNA="lncRNA")
figS7 <- ggplot(d7, aes(region, pct_within_region, fill=period)) +
  geom_bar(stat="identity", position="fill", colour="black", linewidth=0.25) +
  facet_wrap(~biotype, nrow=1, labeller=as_labeller(bt_lab), scales="free_x") +
  scale_fill_manual(values=motif_cols, name="Motif", drop=FALSE) +
  scale_x_discrete(labels=region_lab) +
  scale_y_continuous(expand=expansion(mult=c(0.01,0.02))) +
  labs(x=NULL, y="Proportion") + theme_jgg +
  theme(axis.text.x=element_text(angle=30, hjust=1, size=7))
ggsave("FigureS7.pdf", figS7, width=7, height=4, device=cairo_pdf)

## ---- Fig S8: length & copy-number density by region x motif (genetic_character.R) ----
d8a <- ral_df[ral_df$feature!="Total", ]
d8a$feature <- factor(d8a$feature, levels=region_levels, labels=region_lab[region_levels])
s8a <- ggplot(d8a, aes(reference_allele_length, colour=period_size)) +
  geom_density(linewidth=0.6) + facet_wrap(~feature, nrow=1) +
  scale_colour_manual(values=motif_cols, name="Motif") +
  coord_cartesian(xlim=c(0,100), ylim=c(0,0.2)) +
  labs(x="Reference-allele length (bp)", y="Density") + theme_jgg

d8b <- ncopy_df[ncopy_df$feature!="Total", ]
d8b$feature <- factor(d8b$feature, levels=region_levels, labels=region_lab[region_levels])
s8b <- ggplot(d8b, aes(ncopy, colour=period_size)) +
  geom_density(linewidth=0.6) + facet_wrap(~feature, nrow=1) +
  scale_colour_manual(values=motif_cols, name="Motif") +
  coord_cartesian(xlim=c(0,50), ylim=c(0,0.6)) +
  labs(x="Reference-allele copy number", y="Density") + theme_jgg

figS8 <- (s8a / s8b) + plot_layout(guides="collect") +
  plot_annotation(tag_levels="A") & theme(legend.position="right")
ggsave("FigureS8.pdf", figS8, width=9, height=6, device=cairo_pdf)