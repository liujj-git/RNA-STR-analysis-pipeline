# =============================================================================
# disease.R  —  clean, revision build
# Manuscript: Section 2.7 (disease-associated STRs in the blood transcriptome)
#
# Computation kept VERBATIM. Changes vs original:
#   * The "filter_0,0" column relied on data.frame(check.names=TRUE) silently
#     renaming the comma to a dot ("filter_0.0"). Replaced by an explicit
#     `filter_zero` logical (TRUE = drop "0,0" genotypes) — same logic, no rename.
#   * Shared theme; bare object-name echoes removed.
#   * Section "REVISION VALUES" writes traceability files for §2.7.
# Inputs are ExpansionHunter-specific and read from the disease working dir.
# =============================================================================

rm(list = ls())

library(rjson)
library(jsonlite)
library(vcfR)
library(tidyr)
library(readxl)
library(ggplot2)
library(ggsci)
library(dplyr)

base_theme <- theme_bw() +
  theme(axis.title   = element_text(face = "bold", colour = "black", size = 12),
        axis.title.x = element_text(face = "bold", colour = "black", size = 10),
        axis.title.y = element_text(face = "bold", colour = "black", size = 10),
        axis.text.x  = element_text(face = "bold", colour = "black", size = 8),
        axis.text.y  = element_text(face = "bold", colour = "black", size = 8))

# =============================================================================
# Inputs
# =============================================================================
config_hipstr  <- read.table("../GRCh38.hipstr_reference.refine.bed")
config_disease <- fromJSON("eh.hg38.variant_catalog.disease.json")

disease_position <- read.table("disease.locusposition.txt")
colnames(disease_position) <- c("locus", "subtype", "position")
disease_position <- separate(disease_position, col = "position", into = c("chr", "position"), sep = ":")
disease_position <- separate(disease_position, col = "position", into = c("start", "end"), sep = "-")

locus_intersection <- read.table("locus_intersection.txt")

disease_threshold <- read.table("../../disease_threshold.txt")
temp1 <- disease_threshold[disease_threshold$V1 %in% config_disease$LocusId, ]
colnames(temp1) <- c("locus", "NormalMax", "PathogenicMin")
write.table(temp1, "disease_threshold.config.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# REPCN genotypes (ExpansionHunter)
disease_rna <- read.vcfR("disease.rna.vcf")
disease_rna_gt <- as.data.frame(extract.gt(disease_rna, element = "REPCN"))
disease_dna <- read.vcfR("disease.dna.vcf")
disease_dna_gt <- as.data.frame(extract.gt(disease_dna, element = "REPCN"))
colnames(disease_dna_gt) <- gsub(".deduped", "", colnames(disease_dna_gt))
colnames(disease_dna_gt)[grepl("VB", colnames(disease_dna_gt))] <-
  paste(colnames(disease_dna_gt)[grepl("VB", colnames(disease_dna_gt))], "_dna", sep = "")   # sample-name fix
disease_rna_gt <- data.frame("locus" = rownames(disease_rna_gt), disease_rna_gt)
disease_rna_gt <- gather(disease_rna_gt, key = "sample", value = "genotype", -locus)
disease_dna_gt <- data.frame("locus" = rownames(disease_dna_gt), disease_dna_gt)
disease_dna_gt <- gather(disease_dna_gt, key = "sample", value = "genotype", -locus)
disease_gt <- rbind(data.frame("group" = "dna", disease_dna_gt[, c(2, 1, 3)]),
                    data.frame("group" = "rna", disease_rna_gt[, c(2, 1, 3)]))
disease_gt <- na.omit(disease_gt)
disease_gt$locus <- disease_position$subtype[match(disease_gt$locus, paste(disease_position$chr, disease_position$start, sep = "_"))]
disease_gt$genotype <- gsub("\\/", ",", disease_gt$genotype)
disease_gt <- separate(disease_gt, col = "genotype", into = c("A1", "A2"), sep = ",")
disease_gt[as.numeric(disease_gt$A1) > as.numeric(disease_gt$A2), c("A1", "A2")] <-
  disease_gt[as.numeric(disease_gt$A1) > as.numeric(disease_gt$A2), c("A2", "A1")]
disease_gt <- unite(disease_gt, "genotype", c("A1", "A2"), sep = ",")

# sample info
sample_info <- as.data.frame(read_xlsx("../VBsampleinfo.xlsx", sheet = "concordance"))
colnames(sample_info) <- sample_info[1, ]; sample_info <- sample_info[-1, ]
repeated_samples <- as.data.frame(read_xlsx("../VBsampleinfo.xlsx", sheet = "repeatability"))
colnames(repeated_samples) <- repeated_samples[1, ]; repeated_samples <- repeated_samples[-1, ]
unrelated_individuals <- read_excel("../VBsampleinfo.xlsx", sheet = "unrelated_individuals")
colnames(unrelated_individuals) <- unrelated_individuals[1, ]; unrelated_individuals <- unrelated_individuals[-1, ]

# =============================================================================
# Genotyping rate of disease STRs
# =============================================================================
rate_by_locus <- function(gt_long, samples, n_samples) {
  by_sub <- data.frame()
  for (i in disease_position$subtype) {
    t1 <- gt_long[gt_long$sample %in% samples & gt_long$locus == i, ]
    by_sub <- rbind(by_sub,
                    data.frame(filter_zero = FALSE, subtype = i, genotype_counts = nrow(t1)),
                    data.frame(filter_zero = TRUE,  subtype = i, genotype_counts = nrow(t1[t1$genotype != "0,0", ])))
  }
  by_sub$locus <- disease_position$locus[match(by_sub$subtype, disease_position$subtype)]
  out <- data.frame()
  for (i in unique(by_sub$locus)) {
    a <- by_sub[by_sub$locus == i & !by_sub$filter_zero, ]
    f <- by_sub[by_sub$locus == i &  by_sub$filter_zero, ]
    out <- rbind(out,
                 data.frame(group = "all",      locus = i, genotyping_rate = (sum(a$genotype_counts) / nrow(a)) / n_samples),
                 data.frame(group = "filtered", locus = i, genotyping_rate = (sum(f$genotype_counts) / nrow(f)) / n_samples))
  }
  out
}

## unrelated RNA samples
genotyping_rate_df <- rate_by_locus(disease_gt, unrelated_individuals$ID, nrow(unrelated_individuals))
keep <- genotyping_rate_df$locus[genotyping_rate_df$group == "filtered" & genotyping_rate_df$genotyping_rate >= 0.1]
genotyping_rate_df <- genotyping_rate_df[genotyping_rate_df$locus %in% keep, ]
ord <- genotyping_rate_df[genotyping_rate_df$group == "all", ]
ord <- ord[order(ord$genotyping_rate, decreasing = TRUE), ]
genotyping_rate_df$label <- paste(round(genotyping_rate_df$genotyping_rate, 4) * 100, "%", sep = "")
genotyping_rate_df$locus <- factor(genotyping_rate_df$locus, levels = ord$locus)
genotyping_rate_barplot <- ggplot(genotyping_rate_df, aes(locus, genotyping_rate, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single", width = 0.8), width = 0.8) +
  geom_text(aes(label = label), position = position_dodge(width = 0.8), size = 3, vjust = -0.2) +
  scale_fill_rickandmorty() + base_theme + labs(x = "Locus", y = "Percentage")
ggsave("genotyping_rate_barplot.pdf", genotyping_rate_barplot, height = 5, width = 15)

## matched DNA samples
genotyping_rate_df2 <- rate_by_locus(disease_gt, unique(sample_info$DNA), length(unique(sample_info$DNA)))
keep2 <- genotyping_rate_df2$locus[genotyping_rate_df2$group == "filtered" & genotyping_rate_df2$genotyping_rate >= 0.1]
genotyping_rate_df2 <- genotyping_rate_df2[genotyping_rate_df2$locus %in% keep2, ]
ord2 <- genotyping_rate_df2[genotyping_rate_df2$group == "all", ]
ord2 <- ord2[order(ord2$genotyping_rate, decreasing = TRUE), ]
genotyping_rate_df2$label <- paste(round(genotyping_rate_df2$genotyping_rate, 4) * 100, "%", sep = "")
genotyping_rate_df2$locus <- factor(genotyping_rate_df2$locus, levels = ord2$locus)
genotyping_rate_barplot2 <- ggplot(genotyping_rate_df2, aes(locus, genotyping_rate, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single", width = 0.8), width = 0.8) +
  geom_text(aes(label = label), position = position_dodge(width = 0.8), size = 1.5, vjust = -0.2) +
  scale_fill_manual(values = c("#0073C2FF", "#EFC000FF")) + base_theme +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Genotyping rate of disease-associated STRs in blood genome", x = "Locus", y = "Percentage")
ggsave("genotyping_rate_barplot2.pdf", genotyping_rate_barplot2, height = 5, width = 15)

# =============================================================================
# Concordance for disease STRs  (RNA-RNA replicates, RNA-DNA pairs)
# =============================================================================
kept_loci <- unique(genotyping_rate_df$locus)

# helper: per-pair genotypes for the subtypes of each kept locus
pair_genotypes <- function(s1, s2) {
  out <- data.frame()
  for (lc in kept_loci) {
    subs <- disease_position$subtype[disease_position$locus == lc]
    blk <- data.frame()
    for (sub in subs) {
      g <- disease_gt[disease_gt$locus == sub, ]
      blk <- rbind(blk, data.frame(GT1 = g$genotype[match(s1, g$sample)],
                                   GT2 = g$genotype[match(s2, g$sample)]))
    }
    out <- rbind(out, data.frame(locus = lc, sample1 = s1, sample2 = s2, blk))
  }
  na.omit(out)
}

# helper: identical / partial / non-identical counts per locus, all + filtered
concordance_counts <- function(pairs) {
  res <- data.frame()
  for (lc in unique(pairs$locus)) {
    d <- pairs[pairs$locus == lc, ]
    same <- d[d$GT1 == d$GT2, ]
    res <- rbind(res,
                 data.frame(locus = lc, group = "identical", filter_zero = FALSE, counts = nrow(same)),
                 data.frame(locus = lc, group = "identical", filter_zero = TRUE,
                            counts = nrow(same[same$GT1 != "0,0" & same$GT2 != "0,0", ])))
    diff <- d[d$GT1 != d$GT2, ]
    if (nrow(diff)) {
      diff$inter <- vapply(seq_len(nrow(diff)), function(p)
        length(intersect(strsplit(diff$GT1[p], ",")[[1]], strsplit(diff$GT2[p], ",")[[1]])), integer(1))
      res <- rbind(res,
                   data.frame(locus = lc, group = "partial_identical", filter_zero = FALSE, counts = nrow(diff[diff$inter == 1, ])),
                   data.frame(locus = lc, group = "partial_identical", filter_zero = TRUE,
                              counts = nrow(diff[diff$inter == 1 & diff$GT1 != "0,0" & diff$GT2 != "0,0", ])),
                   data.frame(locus = lc, group = "non_identical", filter_zero = FALSE, counts = nrow(diff[diff$inter == 0, ])),
                   data.frame(locus = lc, group = "non_identical", filter_zero = TRUE,
                              counts = nrow(diff[diff$inter == 0 & diff$GT1 != "0,0" & diff$GT2 != "0,0", ])))
    } else {
      for (gp in c("partial_identical", "non_identical"))
        res <- rbind(res,
                     data.frame(locus = lc, group = gp, filter_zero = FALSE, counts = 0),
                     data.frame(locus = lc, group = gp, filter_zero = TRUE, counts = 0))
    }
  }
  res
}

repeated_df2 <- concordance_counts(pair_genotypes(repeated_samples$R1, repeated_samples$R2))
compared_df2 <- concordance_counts(pair_genotypes(sample_info$RNA, sample_info$DNA))

# plotting helper (filter to loci with >=10 shared pairs)
concordance_plot <- function(cc, zero, title, palette_startrek = TRUE) {
  d <- cc[cc$filter_zero == zero, ]
  for (lc in unique(d$locus)) d$sum[d$locus == lc] <- sum(d$counts[d$locus == lc])
  d <- d[d$sum >= 10, ]
  d$percentage <- d$counts / d$sum
  o <- d[d$group == "identical", ]; o <- o[order(o$percentage, decreasing = TRUE), ]
  d$locus <- factor(d$locus, levels = o$locus)
  d$group <- factor(d$group, levels = c("non_identical", "partial_identical", "identical"))
  g <- ggplot(d, aes(locus, counts, fill = group)) +
    geom_bar(stat = "identity", position = position_dodge(preserve = "single", width = 0.8), width = 0.8) +
    geom_text(aes(label = counts), position = position_dodge(width = 0.8), size = 3, vjust = -0.2) +
    base_theme + labs(title = title, x = "Locus", y = "Shared Genotype Counts")
  if (palette_startrek) g + scale_fill_startrek() else g + scale_fill_brewer(palette = "Set2")
}
ggsave("repeated_all_concordance_barplot.pdf",       concordance_plot(repeated_df2, FALSE, NULL, TRUE),  width = 8,  height = 5)
ggsave("repeated_filtered_concordance_barplot.pdf",  concordance_plot(repeated_df2, TRUE,
                                                                      "Genotype (filtered) concordance of disease-associated STRs between repeated RNA samples", FALSE), width = 10, height = 5)
ggsave("compared_all_concordance_barplot.pdf",       concordance_plot(compared_df2, FALSE, NULL, TRUE),  width = 8,  height = 5)
ggsave("compared_filtered_concordance_barplot.pdf",  concordance_plot(compared_df2, TRUE,
                                                                      "Genotype (filtered) concordance of disease-associated STRs between RNA-DNA sample pairs", FALSE), width = 10, height = 5)

# =============================================================================
# Disease alleles (DNA samples; alleles >= PathogenicMin)
# =============================================================================
disease_allele_df <- data.frame()
for (i in unique(disease_position$locus)) {
  t1 <- disease_gt[disease_gt$sample %in% unique(sample_info$DNA) &
                     disease_gt$locus %in% disease_position$subtype[disease_position$locus == i], ]
  t2 <- as.data.frame(table(unlist(strsplit(t1$genotype, ","))))
  if (!nrow(t2)) next
  t2[, 1] <- as.numeric(as.character(t2[, 1]))
  t2$percentage <- t2$Freq / sum(t2$Freq)
  t2 <- t2[order(t2$percentage, decreasing = TRUE), ]
  t3 <- t2[t2[, 1] >= disease_threshold[disease_threshold[, 1] == i, 3], ]
  if (nrow(t3))
    disease_allele_df <- rbind(disease_allele_df,
                               data.frame(locus = i, allele = t3[, 1], diseaseallele_counts = t3[, 2],
                                          totalallele_counts = sum(t2[, 2]), allele_percentage = t3[, 3]))
}
if (nrow(disease_allele_df))
  disease_allele_df$threshold <- disease_threshold[match(disease_allele_df$locus, disease_threshold[, 1]), 3]
write.table(disease_allele_df, "disease_allele_df.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# =============================================================================
# Genotyping rate: transcriptome vs genome (all disease catalog loci)
# =============================================================================
t_rate <- genotyping_rate_df[genotyping_rate_df$group == "all", ]
g_rate <- genotyping_rate_df2[genotyping_rate_df2$group == "all", ]
genotyping_rate_df3 <- rbind(data.frame(group = "transcriptome", t_rate[, c("locus", "genotyping_rate")]),
                             data.frame(group = "genome",        g_rate[, c("locus", "genotyping_rate")]))
genotyping_rate_df3 <- spread(genotyping_rate_df3, key = "group", value = "genotyping_rate")
genotyping_rate_df3[is.na(genotyping_rate_df3)] <- 0
genotyping_rate_df3 <- left_join(data.frame(locus = config_disease$LocusId), genotyping_rate_df3)
genotyping_rate_df3 <- gather(genotyping_rate_df3, key = "group", value = "genotyping_rate", -locus)
ord3 <- genotyping_rate_df3[genotyping_rate_df3$group == "transcriptome", ]
ord3 <- ord3[order(ord3$genotyping_rate, decreasing = TRUE), ]
genotyping_rate_df3$locus <- factor(genotyping_rate_df3$locus, levels = ord3$locus)
genotyping_rate_df3$label <- paste(round(genotyping_rate_df3$genotyping_rate, 4) * 100, "%", sep = "")
genotyping_rate_barplot3 <- ggplot(genotyping_rate_df3, aes(locus, genotyping_rate, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single", width = 0.8), width = 0.8) +
  geom_text(aes(label = label), position = position_dodge(width = 0.8), size = 2, vjust = -0.2) +
  scale_fill_rickandmorty() + base_theme + theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Locus", y = "Genotyping Rate")
ggsave("genotyping_rate_barplot3.pdf", genotyping_rate_barplot3, height = 5, width = 15)
genotyping_rate_df3_wide <- spread(genotyping_rate_df3[, c("locus", "group", "genotyping_rate")], key = "group", value = "genotyping_rate")
write.table(genotyping_rate_df3_wide, "disease_genotyping_rate.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# =============================================================================
# REVISION VALUES  ->  written to files for the writing list
# =============================================================================
# pooled concordance across kept disease loci (sum over loci), all + filtered
pool <- function(cc, zero) {
  d <- cc[cc$filter_zero == zero, ]
  tot <- sum(d$counts)
  if (tot == 0) return(c(identical = NA, partial = NA, non = NA, n = 0))
  c(identical = 100 * sum(d$counts[d$group == "identical"]) / tot,
    partial   = 100 * sum(d$counts[d$group == "partial_identical"]) / tot,
    non       = 100 * sum(d$counts[d$group == "non_identical"]) / tot,
    n = tot)
}
rr_all <- pool(repeated_df2, FALSE); rr_flt <- pool(repeated_df2, TRUE)
rd_all <- pool(compared_df2, FALSE); rd_flt <- pool(compared_df2, TRUE)

n_catalog <- length(unique(config_disease$LocusId))
n_rna_geno <- length(unique(genotyping_rate_df$locus))     # loci passing >=10% in RNA
n_dna_geno <- length(unique(genotyping_rate_df2$locus))

# --- remaining §2.7 manuscript numbers (read-outs over existing objects) -----
n_DNA <- length(unique(sample_info$DNA))                   # confirm == 55 (manuscript)
n_RNA <- nrow(unrelated_individuals)                       # transcriptome denominator
mean_dna_rate <- 100 * mean(g_rate$genotyping_rate)        # gDNA mean rate (99.55%)
mean_rna_rate <- 100 * mean(t_rate$genotyping_rate)        # transcriptome mean over passing loci
n_rna_gt50    <- sum(t_rate$genotyping_rate > 0.5)         # loci >50% in RNA (6)

# pathogenic-range alleles in DNA (Table S7)
n_path_alleles <- nrow(disease_allele_df)
n_path_loci    <- length(unique(disease_allele_df$locus))
mean_path_freq <- if (n_path_alleles > 0) 100 * mean(disease_allele_df$allele_percentage) else NA

# per-locus maximum RNA-RNA identical rate (Fig S13; manuscript "< 80%")
max_rr_identical <- function(cc, zero) {
  d <- cc[cc$filter_zero == zero, ]
  for (lc in unique(d$locus)) d$sum[d$locus == lc] <- sum(d$counts[d$locus == lc])
  d <- d[d$sum >= 10, ]
  ident <- d[d$group == "identical", ]
  if (nrow(ident) == 0) return(NA)
  max(100 * ident$counts / ident$sum)
}
rr_max_all <- max_rr_identical(repeated_df2, FALSE)
rr_max_flt <- max_rr_identical(repeated_df2, TRUE)

rv <- character(0); add <- function(...) rv[[length(rv) + 1]] <<- sprintf(...)
add("# disease.R revision values  (auto-generated %s)", as.character(Sys.Date()))
add("")
add("## [GENOTYPING RATE]  disease-associated STRs")
add("catalog loci\t%d", n_catalog)
add("DNA samples (denominator; confirm == 55)\t%d", n_DNA)
add("RNA unrelated samples (denominator)\t%d", n_RNA)
add("loci genotyped >=10%% in blood transcriptome (RNA)\t%d", n_rna_geno)
add("loci genotyped >50%% in blood transcriptome (RNA)\t%d", n_rna_gt50)
add("loci genotyped >=10%% in blood genome (DNA)\t%d", n_dna_geno)
add("mean genome (DNA) rate over its loci\t%.2f%%", mean_dna_rate)
add("mean transcriptome (RNA) rate over its loci\t%.2f%%", mean_rna_rate)
add("# per-locus transcriptome vs genome -> disease_genotyping_rate.txt")
add("")
add("## [CONCORDANCE of disease STRs]  (loci with >=10 shared pairs; identical/partial/non)")
add("RNA-RNA all\tidentical %.1f%%\tpartial %.1f%%\tnon %.1f%%\tn=%d", rr_all["identical"], rr_all["partial"], rr_all["non"], as.integer(rr_all["n"]))
add("RNA-RNA filtered(0,0 dropped)\tidentical %.1f%%\tpartial %.1f%%\tnon %.1f%%\tn=%d", rr_flt["identical"], rr_flt["partial"], rr_flt["non"], as.integer(rr_flt["n"]))
add("RNA-DNA all\tidentical %.1f%%\tpartial %.1f%%\tnon %.1f%%\tn=%d", rd_all["identical"], rd_all["partial"], rd_all["non"], as.integer(rd_all["n"]))
add("RNA-DNA filtered(0,0 dropped)\tidentical %.1f%%\tpartial %.1f%%\tnon %.1f%%\tn=%d", rd_flt["identical"], rd_flt["partial"], rd_flt["non"], as.integer(rd_flt["n"]))
add("per-locus MAX RNA-RNA identical (Fig S13; manuscript '<80%%')\tall %.1f%%\tfiltered %.1f%%", rr_max_all, rr_max_flt)
add("")
add("## [DISEASE ALLELES]  (DNA; alleles >= PathogenicMin)  -> disease_allele_df.txt")
add("pathogenic-range alleles\t%d", n_path_alleles)
add("loci carrying a pathogenic-range allele\t%d", n_path_loci)
add("mean pathogenic-allele frequency\t%.2f%%", mean_path_freq)
if (n_path_alleles > 0) for (i in seq_len(n_path_alleles))
  add("  %s\tallele=%g\tthreshold=%g\tfreq=%.2f%%", disease_allele_df$locus[i], disease_allele_df$allele[i],
      disease_allele_df$threshold[i], 100 * disease_allele_df$allele_percentage[i])
writeLines(rv, "revision_values_disease.txt")

cat("disease.R done. filter_zero cleanup applied; revision files written.\n")