# ============================================================================
# Script: STR Genotyping Rate, Polymorphism Classification, and Expression Analysis
# Description: 
#   1. Load sample info, STR configuration, and HipSTR VCF output.
#   2. Compute genotyping rate per locus, filter by rate >= 0.5.
#   3. Classify loci as polymorphic (multiple alleles) or monomorphic (single allele).
#   4. Visualize genotyping rate, heterozygosity, and expression (depth of coverage, DP).
#   5. Output locus lists for downstream analysis.
# ============================================================================

rm(list = ls())

# ----------------------------- Libraries ------------------------------------
library(dplyr)
library(tidyr)
library(readxl)
library(vcfR)
library(ggplot2)
library(ggpubr)
library(EnvStats)
library(pheatmap)
library(Hmisc)
library(ggsci)

# ----------------------------- Input data -----------------------------------

# Sample information
unrelated_individuals <- read_excel("../VBsampleinfo.xlsx", sheet = "unrelated_individuals")
colnames(unrelated_individuals) <- unrelated_individuals[1, ]
unrelated_individuals <- unrelated_individuals[-1, ]

sample_gender <- read_excel("../VBsampleinfo.xlsx", sheet = "total")
colnames(sample_gender) <- sample_gender[1, ]
sample_gender <- sample_gender[-1, ]
unrelated_individuals$gender <- sample_gender$gender[match(unrelated_individuals$ID, sample_gender$ID)]

# STR configuration file (reference BED)
config <- read.table("../GRCh38.hipstr_reference.refine.bed")
config$reference_allele_length <- config[, 3] - config[, 2] + 1

# Loci to remove (wrong alleles)
wrong_locus <- read.table("removedlocus_with_wrongallele.txt", fill = TRUE)

# HipSTR VCF (RNA-seq)
raw_vcf <- read.vcfR("hipstr.filtered.rna.vcf.gz")
GT <- extract.gt(raw_vcf)   # Genotype matrix
GT <- GT[!rownames(GT) %in% wrong_locus[, 1], colnames(GT) %in% unrelated_individuals$ID]

DP <- extract.gt(raw_vcf, element = "DP")   # Depth of coverage matrix
DP <- DP[, colnames(DP) %in% unrelated_individuals$ID]
rm(raw_vcf)

# ----------------------------- Helper function ------------------------------
# Count number of samples with non-NA genotype for a locus
calculate_genotype_counts <- function(x) {
  length(na.omit(x))
}

# ----------------------------- Step 1: Genotyping rate ----------------------

genotyping_rate_df <- as.data.frame(apply(GT, 1, calculate_genotype_counts))
colnames(genotyping_rate_df)[1] <- "counts"
genotyping_rate_df$chr <- config[match(rownames(genotyping_rate_df), config[, 6]), 1]

# Genotyping rate: for Y chromosome, denominator = number of males; otherwise all samples
genotyping_rate_df$genotyping_rate[genotyping_rate_df$chr == "Y"] <- 
  genotyping_rate_df$counts[genotyping_rate_df$chr == "Y"] / nrow(unrelated_individuals[unrelated_individuals$gender == "male", ])
genotyping_rate_df$genotyping_rate[genotyping_rate_df$chr != "Y"] <- 
  genotyping_rate_df$counts[genotyping_rate_df$chr != "Y"] / nrow(unrelated_individuals)

# Filter loci with genotyping rate >= 0.5
gr_threshold <- 0.5
loci_sufficient_gr <- rownames(genotyping_rate_df)[genotyping_rate_df$genotyping_rate >= gr_threshold]

# Export loci with zero expression (genotyping rate = 0)
zero_expr_loci <- config[!config[, 6] %in% rownames(genotyping_rate_df[genotyping_rate_df$genotyping_rate != 0, ]), -ncol(config)]
colnames(zero_expr_loci) <- c("Chr", "Start", "End", "Motif_length", "Copy_Number", "Locus", "Motif_Sequence")
write.table(zero_expr_loci, file = "Loci_withnoexpression.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# Apply filter to GT and DP matrices
GT <- GT[rownames(GT) %in% loci_sufficient_gr, ]
DP <- DP[rownames(DP) %in% loci_sufficient_gr, ]

# ----------------------------- Step 2: Polymorphism classification ----------
# For each locus, determine number of distinct alleles and heterozygosity
calculate_genotype_types <- function(x) {
  genotype <- na.omit(x)
  unique_genotypes <- unique(genotype)
  
  # Determine homozygous/heterozygous status for each unique genotype
  het_status <- sapply(unique_genotypes, function(gt) {
    alleles <- unlist(strsplit(gt, "\\|"))
    ifelse(length(unique(alleles)) == 1, "homo", "het")
  })
  
  # For all called genotypes (including repeated)
  all_het_status <- sapply(genotype, function(gt) {
    alleles <- unlist(strsplit(gt, "\\|"))
    ifelse(length(unique(alleles)) == 1, "homo", "het")
  })
  
  data.frame(
    genotype_types = length(unique_genotypes),
    het_types = sum(het_status == "het"),
    genotype_counts = length(genotype),
    het_counts = sum(all_het_status == "het")
  )
}

genotype_types_df <- do.call(rbind, apply(GT, 1, calculate_genotype_types))
genotype_types_df$Heterozygosity <- genotype_types_df$het_counts / genotype_types_df$genotype_counts

# Polymorphic: more than one distinct genotype; Monomorphic: only one genotype
polymorphism_group <- rbind(
  data.frame(group = "polymorphic", locus = rownames(genotype_types_df)[genotype_types_df$genotype_types != 1]),
  data.frame(group = "monomorphic",  locus = rownames(genotype_types_df)[genotype_types_df$genotype_types == 1])
)

# ----------------------------- Step 3: Visualization ------------------------

# Prepare merged data frame
temp <- polymorphism_group
temp$genotyping_rate <- genotyping_rate_df$genotyping_rate[match(temp$locus, rownames(genotyping_rate_df))]
temp$period_size <- config[match(temp$locus, config[, 6]), 4]
temp$reference_allele_length <- config[match(temp$locus, config[, 6]), "reference_allele_length"]
temp$group <- factor(temp$group, levels = c("polymorphic", "monomorphic"))
temp$period_size <- factor(temp$period_size, levels = 2:6,
                           labels = c("Di", "Tri", "Tetra", "Penta", "Hexa"))

# (3.1) Boxplot: genotyping rate by period size and group
my_comparison <- list(c("polymorphic", "monomorphic"))
gr_boxplot <- ggplot(temp, aes(x = period_size, y = genotyping_rate, fill = group)) +
  geom_boxplot(position = position_dodge(0.9)) +
  facet_grid(~period_size, scale = "free") +
  stat_compare_means(method = "kruskal.test", label = "p.format") +
  scale_fill_lancet() +
  theme_bw() +
  labs(title = "Distribution of Genotyping Rate", x = "Groups", y = "Genotyping Rate") +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 10))

ggsave(gr_boxplot, file = "genotyping_rate_boxplot.pdf", width = 15, height = 5)

# (3.2) Correlation: genotyping rate vs reference allele length
gr_cor <- ggplot(temp, aes(x = reference_allele_length, y = genotyping_rate)) +
  geom_smooth(method = "lm", formula = y ~ x, color = '#2D6DB1', fill = "#756bb1") +
  facet_grid(~group, scale = "free") +
  stat_cor() +
  theme_bw() +
  labs(title = "Correlation between Genotyping Rate and Reference Allele Length",
       x = "Reference Length (bp)", y = "Genotyping Rate") +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 10))

ggsave(gr_cor, file = "genotyping_rate_corplot.pdf", width = 10, height = 5)

# ----------------------------- Step 4: Expression (DP) analysis -------------

# Convert DP matrix to numeric
DP_numeric <- as.data.frame(sapply(as.data.frame(DP), as.numeric))
DP_numeric[is.na(DP_numeric)] <- 0
DP_df <- data.frame(locus = rownames(DP_numeric), DP_numeric)

# Mean DP per locus
meanDP_df <- data.frame(
  locus = DP_df$locus,
  meanDP = rowMeans(DP_df[, -1])
)
meanDP_df$group <- polymorphism_group$group[match(meanDP_df$locus, polymorphism_group$locus)]
meanDP_df <- meanDP_df[order(meanDP_df$meanDP, decreasing = TRUE), ]
meanDP_df$period_size <- config[match(meanDP_df$locus, config[, 6]), 4]
meanDP_df$reference_allele_length <- config[match(meanDP_df$locus, config[, 6]), "reference_allele_length"]
meanDP_df$group <- factor(meanDP_df$group, levels = c("polymorphic", "monomorphic"))
meanDP_df$period_size <- factor(meanDP_df$period_size, levels = 2:6,
                                labels = c("Di", "Tri", "Tetra", "Penta", "Hexa"))

# Export locus list with coordinates
coord <- config[match(meanDP_df$locus, config$V6), c(1, 2, 3, 6)]
colnames(coord) <- c("chr", "start", "end", "locus")
meanDP_df <- left_join(meanDP_df, coord, by = "locus")
write.table(meanDP_df, file = "locus_list.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(meanDP_df$locus[meanDP_df$group == "polymorphic"], file = "locus_list.polymorphic.txt",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

# Heatmap: polymorphic loci
poly_loci <- meanDP_df$locus[meanDP_df$group == "polymorphic"]
poly_mat <- t(DP_numeric[rownames(DP_numeric) %in% poly_loci, ])
poly_mat_log <- log(poly_mat + 1, 2)
pheatmap_poly <- pheatmap(poly_mat_log,
                          show_rownames = TRUE, show_colnames = FALSE,
                          cluster_rows = FALSE, cluster_cols = FALSE,
                          fontsize_row = 7.5,
                          angle_row = 90,
                          color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
ggsave(pheatmap_poly, file = "pheatmap_polymorphic.pdf", width = 15, height = 15)

# Heatmap: monomorphic loci
mono_loci <- meanDP_df$locus[meanDP_df$group == "monomorphic"]
mono_mat <- t(DP_numeric[rownames(DP_numeric) %in% mono_loci, ])
mono_mat_log <- log(mono_mat + 1, 2)
pheatmap_mono <- pheatmap(mono_mat_log,
                          show_rownames = TRUE, show_colnames = FALSE,
                          cluster_rows = FALSE, cluster_cols = FALSE,
                          fontsize_row = 7.5,
                          angle_row = 90,
                          color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
ggsave(pheatmap_mono, file = "pheatmap_monomorphic.pdf", width = 15, height = 15)

# Boxplot: log2(meanDP+1) by period size and group
dp_boxplot <- ggplot(meanDP_df, aes(x = period_size, y = log(meanDP + 1, 2), fill = group)) +
  geom_boxplot(position = position_dodge(0.9)) +
  facet_grid(~period_size, scale = "free") +
  stat_compare_means(method = "kruskal.test", label = "p.format") +
  scale_fill_lancet() +
  theme_bw() +
  labs(title = "Distribution of log2-Average Call Depth", x = "Groups", y = "log2-Average Call Depth") +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 10))

ggsave(dp_boxplot, file = "expression_boxplot.pdf", width = 15, height = 5)

# Correlation: log2(meanDP+1) vs reference allele length
dp_cor <- ggplot(meanDP_df, aes(x = reference_allele_length, y = log(meanDP + 1, 2))) +
  geom_smooth(method = "lm", formula = y ~ x, color = '#2D6DB1', fill = "#756bb1") +
  facet_grid(~group, scale = "free") +
  stat_cor() +
  theme_bw() +
  labs(title = "Correlation between log2-Average Call Depth and Reference Allele Length",
       x = "Reference Length (bp)", y = "log2-Average Call Depth") +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 10))

ggsave(dp_cor, file = "expression_corplot.pdf", width = 10, height = 5)

# ----------------------------- Step 5: Ideogram file -------------------------
# Generate input for Ideogram.js visualization
idio <- rbind(
  data.frame(group = "mapping", font_size = ".", color = "144,0,33", locus = ".",
             shape = "filled_box", chr = paste0("chr", config[match(polymorphism_group$locus[polymorphism_group$group == "polymorphic"], config[, 6]), 1]),
             start = config[match(polymorphism_group$locus[polymorphism_group$group == "polymorphic"], config[, 6]), 2],
             end = config[match(polymorphism_group$locus[polymorphism_group$group == "polymorphic"], config[, 6]), 3],
             strand = "."),
  data.frame(group = "mapping", font_size = ".", color = "0,47,167", locus = ".",
             shape = "filled_circle", chr = paste0("chr", config[match(polymorphism_group$locus[polymorphism_group$group == "monomorphic"], config[, 6]), 1]),
             start = config[match(polymorphism_group$locus[polymorphism_group$group == "monomorphic"], config[, 6]), 2],
             end = config[match(polymorphism_group$locus[polymorphism_group$group == "monomorphic"], config[, 6]), 3],
             strand = ".")
)
write.table(idio, file = "idio.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

# ----------------------------- Supplementary Figures -------------------------

# (S1) Cumulative fraction of genotyping rate stratified by motif length and group
temp1 <- polymorphism_group
temp1$genotyping_rate <- genotyping_rate_df$genotyping_rate[match(temp1$locus, rownames(genotyping_rate_df))]
temp1$heterozygosity <- genotype_types_df$Heterozygosity[match(temp1$locus, rownames(genotype_types_df))]
temp1$reference_allele_length <- config$reference_allele_length[match(temp1$locus, config[, 6])]
temp1$motif_length <- config[match(temp1$locus, config[, 6]), 4]
temp1$meanDP <- meanDP_df$meanDP[match(temp1$locus, meanDP_df$locus)]

cumulative_df <- data.frame()
for (ml in 2:6) {
  for (grp in c("polymorphic", "monomorphic")) {
    sub <- temp1[temp1$motif_length == ml & temp1$group == grp, ]
    for (p in seq(0.5, 1, 0.1)) {
      cumulative_df <- rbind(cumulative_df,
                             data.frame(motif_length = ml, group = grp,
                                        genotyping_rate = p,
                                        cumulative_fraction = 1 - nrow(sub[sub$genotyping_rate > p, ]) / nrow(sub)))
    }
  }
}
cumulative_df$motif_length <- factor(cumulative_df$motif_length, levels = 2:6)
cumulative_df$group <- factor(cumulative_df$group, levels = c("polymorphic", "monomorphic"))

gr_cfl <- ggplot(cumulative_df, aes(x = genotyping_rate, y = cumulative_fraction)) +
  geom_line(aes(col = motif_length, linetype = group), size = 1) +
  scale_color_lancet() +
  theme_bw() +
  labs(x = "Genotyping Rate", y = "Cumulative Fraction") +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 10))
ggsave(gr_cfl, file = "gr_cfl.pdf", width = 8, height = 5)

# (S2) Cumulative fraction of heterozygosity for polymorphic loci
het_cumul <- data.frame()
for (ml in 2:6) {
  sub <- temp1[temp1$motif_length == ml & temp1$group == "polymorphic", ]
  for (p in seq(0.1, 1, 0.1)) {
    het_cumul <- rbind(het_cumul,
                       data.frame(motif_length = ml, heterozygosity = p,
                                  cumulative_fraction = 1 - nrow(sub[sub$heterozygosity > p, ]) / nrow(sub)))
  }
}
het_cumul$motif_length <- factor(het_cumul$motif_length, levels = 2:6)

ho_cfl <- ggplot(het_cumul, aes(x = heterozygosity, y = cumulative_fraction, group = motif_length, col = motif_length)) +
  geom_line(size = 1) +
  scale_color_lancet() +
  theme_bw() +
  labs(x = "Heterozygosity", y = "Cumulative Fraction") +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 10))
ggsave(ho_cfl, file = "ho_cfl.pdf", width = 8, height = 5)

# (S3) Heatmap of mean DP by reference allele length bins and motif length
breaks <- c(0, 20, 50, 100)
heat_df <- data.frame()
for (ml in 2:6) {
  for (grp in c("polymorphic", "monomorphic")) {
    sub <- temp1[temp1$motif_length == ml & temp1$group == grp, ]
    for (i in 1:length(breaks)) {
      if (i != length(breaks)) {
        mean_val <- mean(sub$meanDP[sub$reference_allele_length > breaks[i] & sub$reference_allele_length <= breaks[i+1]], na.rm = TRUE)
        bin_label <- paste0("(", breaks[i], ",", breaks[i+1], "]")
      } else {
        mean_val <- mean(sub$meanDP[sub$reference_allele_length > breaks[i]], na.rm = TRUE)
        bin_label <- paste0(">", breaks[i])
      }
      heat_df <- rbind(heat_df,
                       data.frame(group = grp, motif_length = ml,
                                  ref_bin = bin_label, meanDP = mean_val))
    }
  }
}
heat_df$label <- ifelse(!is.na(heat_df$meanDP), round(heat_df$meanDP, 2), "NA")
heat_df$motif_length <- factor(heat_df$motif_length, levels = 2:6)
heat_df$ref_bin <- factor(heat_df$ref_bin, levels = c("(0,20]", "(20,50]", "(50,100]", ">100"))

heatmap_summary <- ggplot(heat_df, aes(x = motif_length, y = ref_bin, fill = meanDP)) +
  facet_wrap(~group, nrow = 2) +
  geom_raster() +
  geom_tile(col = "gray") +
  geom_text(aes(label = label), color = "black", size = 3) +
  scale_fill_gradient2(na.value = "white", low = "#003366", high = "#990033", mid = "#ffffff") +
  scale_y_discrete(position = "left") +
  coord_fixed() +
  theme_bw() +
  labs(x = "Motif Length", y = "Reference Allele Length") +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 10))
ggsave(heatmap_summary, file = "heatmap_summary.pdf", width = 8, height = 5)

# End of script