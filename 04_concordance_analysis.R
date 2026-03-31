# ============================================================================
# Script: Concordance analysis for RNA replicates and RNA-DNA comparisons
# Description:
#   1. Load STR configuration, sample info, VCF (GT/DP), and pSTR annotation.
#   2. Build allele sequence ladder (with isoallele handling).
#   3. RNA technical replicate analysis:
#      - Genotype counts, shared loci, concordance (identical/partial/non-identical)
#      - Step difference distribution for partial-identical calls
#      - Sequence discordance (isoallele variants): PHP count, base change patterns, positional distribution
#   4. RNA-DNA matched pair analysis:
#      - Similar steps with DNA DP filtering (>=3)
#      - Concordance, step differences, isoallele analysis
#   5. Regional concordance and reference allele length effects
#   6. Supplementary: Venn diagrams of discordant loci across replicates
# ============================================================================

rm(list = ls())
library(readxl)
library(vcfR)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(EnvStats)
library(data.table)
library(ggsci)
library(ggvenn)

# ----------------------------- Input data -----------------------------------
# STR configuration (reference BED)
config <- read.table("../GRCh38.hipstr_reference.refine.bed")

# Sample info for replicates and concordance pairs
repeated_samples <- as.data.frame(read_xlsx("../VBsampleinfo.xlsx", sheet = "repeatability"))
colnames(repeated_samples) <- repeated_samples[1, ]
repeated_samples <- repeated_samples[-1, ]

sample_info <- as.data.frame(read_xlsx("../VBsampleinfo.xlsx", sheet = "concordance"))
colnames(sample_info) <- sample_info[1, ]
sample_info <- sample_info[-1, ]

# VCF (genotypes and depth)
raw_vcf <- read.vcfR("../pstr.recode.vcf.gz")
GT <- extract.gt(raw_vcf, element = "GT")
DP <- extract.gt(raw_vcf, element = "DP")
DP <- as.data.frame(DP)
DP <- data.frame(locus = rownames(DP), gather(DP, key = "sample", value = "DP"))
DP <- na.omit(DP)

# pSTR annotation group
pSTR_group <- read.table("../2_annotation/STR_annotation_group.txt", header = TRUE)

# ----------------------------- Allele sequence ladder ------------------------
# Build allele sequence mapping with isoallele resolution
allele_ladder <- data.frame(
  locus = raw_vcf@fix[, 3],
  Ref = raw_vcf@fix[, 4],
  Alt = raw_vcf@fix[, 5]
)
allele_ladder <- na.omit(gather(allele_ladder, key = "allele_group", value = "sequence", -locus))

# Split multi-allelic Alt entries
temp1 <- allele_ladder[!grepl(",", allele_ladder$sequence), ]
temp2 <- allele_ladder[grepl(",", allele_ladder$sequence), ]

if (nrow(temp2) > 0) {
  temp3 <- do.call(rbind, apply(temp2, 1, function(x) {
    alt_allele <- unlist(strsplit(x[3], ","))
    data.frame(locus = x[1],
               allele_group = 1:length(alt_allele),
               sequence = alt_allele,
               stringsAsFactors = FALSE)
  }))
  temp3 <- temp3[nchar(temp3$sequence) != 0, ]
  temp1 <- rbind(temp1, temp3)
}

temp1$allele_group[temp1$allele_group == "Ref"] <- 0
temp1$allele_group[temp1$allele_group == "Alt"] <- 1

# Compute copy number and handle isoalleles (same length, different sequence)
temp1$period_size <- as.numeric(config[match(temp1$locus, config[, 6]), 4])
ncopy_int <- nchar(temp1$sequence) %/% temp1$period_size
remainder <- nchar(temp1$sequence) %% temp1$period_size
remainder[remainder != 0] <- paste0(".", remainder[remainder != 0])
remainder[remainder == 0] <- ""
temp1$ncopy <- paste0(ncopy_int, remainder)
temp1$ncopy_label <- paste(temp1$locus, temp1$ncopy, sep = "@")

dup_labels <- as.data.frame(table(temp1$ncopy_label))
dup_labels <- dup_labels[dup_labels$Freq != 1, 1]
setDT(temp1)
temp1[!ncopy_label %in% dup_labels, ncopy_label := ncopy]
temp1[ncopy_label %in% dup_labels, ncopy_label := paste(ncopy, seq_len(.N), sep = "_"), by = ncopy_label]
temp1$unique_label <- paste(temp1$locus, temp1$allele_group, sep = "@")
allele_ladder <- temp1[, c("locus", "period_size", "allele_group", "sequence", "ncopy", "ncopy_label")]
allele_ladder <- allele_ladder[order(allele_ladder$locus), ]
write.table(allele_ladder, file = "allele_sequence.concordance.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

allele_sequence <- allele_ladder
allele_sequence$ncopy <- as.numeric(allele_sequence$ncopy)
rm(raw_vcf, allele_ladder)

# ----------------------------- RNA technical replicates ----------------------
# Genotype counts per replicate sample
sample_GTcounts <- data.frame()
for (i in seq_along(repeated_samples$individual)) {
  temp1 <- data.frame(individual = repeated_samples$individual[i],
                      group = "duplication1",
                      sample = repeated_samples$R1[i],
                      GTcounts = length(na.omit(GT[, colnames(GT) == repeated_samples$R1[i]])))
  temp2 <- data.frame(individual = repeated_samples$individual[i],
                      group = "duplication2",
                      sample = repeated_samples$R2[i],
                      GTcounts = length(na.omit(GT[, colnames(GT) == repeated_samples$R2[i]])))
  sample_GTcounts <- rbind(sample_GTcounts, temp1, temp2)
}
sample_GTcounts <- sample_GTcounts[order(sample_GTcounts$GTcounts, decreasing = TRUE), ]
sample_GTcounts$individual <- factor(sample_GTcounts$individual, levels = unique(sample_GTcounts$individual))
sample_GTcounts$group <- factor(sample_GTcounts$group, levels = c("duplication1", "duplication2"))

repeatedsample_GTcounts_barplot <- ggplot(sample_GTcounts, aes(x = individual, y = GTcounts)) +
  geom_bar(aes(fill = group), stat = "identity", position = position_dodge(0.8), width = 0.7) +
  scale_fill_rickandmorty() +
  theme_bw() +
  labs(title = "Genotype Counts of Repeated RNA Samples", x = "Individuals", y = "Genotype Counts") +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face = "bold", size = 10, angle = 90),
        axis.text.y = element_text(face = "bold", size = 10))
ggsave(repeatedsample_GTcounts_barplot, file = "repeatedsample_GTcounts_barplot.pdf", height = 6, width = 12)

# Shared loci counts per individual
sharedlocus_counts_df <- data.frame()
for (i in repeated_samples$individual) {
  s1 <- as.character(repeated_samples$R1[repeated_samples$individual == i])
  s2 <- as.character(repeated_samples$R2[repeated_samples$individual == i])
  temp <- GT[, match(c(s2, s1), colnames(GT))]
  temp_df <- data.frame(individual = i,
                        duplication1 = s1,
                        duplication2 = s2,
                        duplication1_locus_counts = length(na.omit(temp[, 1])),
                        duplication2_locus_counts = length(na.omit(temp[, 2])),
                        shared_locus_counts = nrow(na.omit(temp)))
  sharedlocus_counts_df <- rbind(sharedlocus_counts_df, temp_df)
}

# Concordance per locus between replicates
concordance_repeated_df <- data.frame()
for (i in repeated_samples$individual) {
  s1 <- repeated_samples$R1[repeated_samples$individual == i]
  s2 <- repeated_samples$R2[repeated_samples$individual == i]
  cat(paste(s2, s1, collapse = "--"), "\n")
  temp <- as.data.frame(na.omit(GT[, match(c(s2, s1), colnames(GT))]))
  temp <- gather(data.frame(locus = rownames(temp), temp), key = "sample", value = "GT", -locus)
  temp <- separate(temp, col = "GT", into = c("v1", "v2"), sep = "\\|")
  # Sort alleles to avoid "0,1" vs "1,0" misclassification
  idx <- as.numeric(temp$v1) > as.numeric(temp$v2)
  temp[idx, c("v1", "v2")] <- temp[idx, c("v2", "v1")]
  temp <- unite(temp, c("v1", "v2"), col = "GT", sep = ",")
  temp <- spread(temp, key = "sample", value = "GT")
  rownames(temp) <- temp$locus
  temp <- temp[, -1]
  
  # Identical calls
  identical_loci <- temp[temp[, 1] == temp[, 2], ]
  if (nrow(identical_loci) > 0) {
    df_ident <- data.frame(individual = i, duplication1 = s1, duplication2 = s2,
                           locus = rownames(identical_loci),
                           GT1 = identical_loci[, 1], GT2 = identical_loci[, 2],
                           concordance = 2)  # 2 = identical
    concordance_repeated_df <- rbind(concordance_repeated_df, df_ident)
  }
  # Non-identical
  diff_loci <- temp[temp[, 1] != temp[, 2], ]
  if (nrow(diff_loci) > 0) {
    df_diff <- data.frame()
    for (p in 1:nrow(diff_loci)) {
      gt1 <- sort(as.numeric(unlist(strsplit(diff_loci[p, 1], ","))))
      gt2 <- sort(as.numeric(unlist(strsplit(diff_loci[p, 2], ","))))
      shared <- length(intersect(gt1, gt2))
      conc <- ifelse(shared == 1, 1, 0)   # 1 = partial-identical, 0 = non-identical
      df_diff <- rbind(df_diff, data.frame(individual = i, duplication1 = s1, duplication2 = s2,
                                           locus = rownames(diff_loci)[p],
                                           GT1 = paste(gt1, collapse = ","),
                                           GT2 = paste(gt2, collapse = ","),
                                           concordance = conc))
    }
    concordance_repeated_df <- rbind(concordance_repeated_df, df_diff)
  }
}
# Add DP information
concordance_repeated_df$DP1 <- as.numeric(DP$DP[match(paste(concordance_repeated_df$duplication1, concordance_repeated_df$locus, sep = "@"),
                                                      paste(DP$sample, DP$locus, sep = "@"))])
concordance_repeated_df$DP2 <- as.numeric(DP$DP[match(paste(concordance_repeated_df$duplication2, concordance_repeated_df$locus, sep = "@"),
                                                      paste(DP$sample, DP$locus, sep = "@"))])
write.table(concordance_repeated_df, file = "concordance_repeated_df.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# Concordance composition at different DP thresholds (0,10,20,30)
conc_rep_sample <- data.frame()
conc_rep_sample_counts <- data.frame()
for (th in c(0, 10, 20, 30)) {
  sub <- concordance_repeated_df[concordance_repeated_df$DP1 >= th & concordance_repeated_df$DP2 >= th, ]
  if (nrow(sub) > 0) {
    for (ind in unique(sub$individual)) {
      tab <- as.data.frame(table(sub[sub$individual == ind, "concordance"]))
      tab$percentage <- tab$Freq / sum(tab$Freq)
      conc_rep_sample <- rbind(conc_rep_sample, data.frame(DPthreshold = th, individual = ind,
                                                           group = as.character(tab$Var1),
                                                           percentage = tab$percentage,
                                                           counts = tab$Freq))
      conc_rep_sample_counts <- rbind(conc_rep_sample_counts, data.frame(DPthreshold = th, individual = ind,
                                                                         counts = sum(tab$Freq)))
    }
  }
}
# Order individuals by total counts at DP=0
order_ind <- conc_rep_sample_counts[conc_rep_sample_counts$DPthreshold == 0, ]
order_ind <- order_ind[order(order_ind$counts, decreasing = FALSE), ]
conc_rep_sample$individual <- factor(conc_rep_sample$individual, levels = order_ind$individual)
conc_rep_sample$DPthreshold <- factor(conc_rep_sample$DPthreshold, levels = c(0, 10, 20, 30))
conc_rep_sample$group <- factor(conc_rep_sample$group, levels = c(0, 1, 2),
                                labels = c("non_identical", "partial_identical", "identical"))

# Plot for DP=0 (main)
p_rep0 <- ggplot(conc_rep_sample[conc_rep_sample$DPthreshold == 0, ],
                 aes(x = individual, y = counts, fill = group)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  scale_fill_startrek() +
  theme_bw() +
  labs(title = "The Concordance between Repeated RNA Sample Pairs", x = "Individual", y = "Loci Counts") +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 8))
ggsave(p_rep0, file = "concordance_repeated_sample_plot.pdf", width = 12, height = 12)

# Plot for DP>0
p_rep_other <- ggplot(conc_rep_sample[conc_rep_sample$DPthreshold != 0, ],
                      aes(x = individual, y = counts, fill = group)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~DPthreshold) +
  coord_flip() +
  scale_fill_startrek() +
  theme_bw() +
  labs(title = "The Concordance between Repeated RNA Sample Pairs under a Sequence of DP Thresholds",
       x = "Individual", y = "Loci Counts") +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.y = element_text(face = "bold", size = 8))
ggsave(p_rep_other, file = "concordance_repeated_sample_plot2.pdf", width = 15, height = 6)

# Step difference analysis for partial-identical calls
halfidentical_repeated <- concordance_repeated_df[concordance_repeated_df$concordance == 1,
                                                  c("individual", "duplication1", "duplication2",
                                                    "locus", "GT1", "GT2", "DP1", "DP2")]
halfidentical_repeated <- separate(halfidentical_repeated, col = "GT1", into = c("GT1_A1", "GT1_A2"), sep = ",")
halfidentical_repeated <- separate(halfidentical_repeated, col = "GT2", into = c("GT2_A1", "GT2_A2"), sep = ",")

for (i in 1:nrow(halfidentical_repeated)) {
  cat(round(i / nrow(halfidentical_repeated) * 100, 2), "%\n")
  gt1 <- halfidentical_repeated[i, c("GT1_A1", "GT1_A2")]
  gt2 <- halfidentical_repeated[i, c("GT2_A1", "GT2_A2")]
  common <- intersect(gt1, gt2)
  halfidentical_repeated$dif_GT1[i] <- ifelse(length(gt1[gt1 != common]) > 0, gt1[gt1 != common], common)
  halfidentical_repeated$dif_GT2[i] <- ifelse(length(gt2[gt2 != common]) > 0, gt2[gt2 != common], common)
}
# Map to copy numbers
idx1 <- match(paste(halfidentical_repeated$locus, halfidentical_repeated$dif_GT1, sep = "@"),
              paste(allele_sequence$locus, allele_sequence$allele_group, sep = "@"))
idx2 <- match(paste(halfidentical_repeated$locus, halfidentical_repeated$dif_GT2, sep = "@"),
              paste(allele_sequence$locus, allele_sequence$allele_group, sep = "@"))
halfidentical_repeated$ncopy1 <- as.numeric(allele_sequence$ncopy[idx1])
halfidentical_repeated$ncopy2 <- as.numeric(allele_sequence$ncopy[idx2])
halfidentical_repeated$dif_step <- abs(halfidentical_repeated$ncopy1 - halfidentical_repeated$ncopy2)
halfidentical_repeated$sequence1 <- allele_sequence$sequence[idx1]
halfidentical_repeated$sequence2 <- allele_sequence$sequence[idx2]
halfidentical_repeated <- na.omit(halfidentical_repeated)

# Step difference distribution at different DP thresholds
step_rep_summary <- data.frame()
for (th in c(0, 10, 20, 30)) {
  sub <- halfidentical_repeated[halfidentical_repeated$DP1 >= th & halfidentical_repeated$DP2 >= th, ]
  if (nrow(sub) > 0) {
    sub_step <- sub$dif_step[sub$dif_step %% 1 == 0 & sub$dif_step <= 10]
    tab <- as.data.frame(table(sub_step))
    step_rep_summary <- rbind(step_rep_summary, data.frame(DPthreshold = th, dif_step = as.numeric(as.character(tab$Var1)),
                                                           counts = tab$Freq))
  }
}
step_rep_summary$DPthreshold <- factor(step_rep_summary$DPthreshold, levels = c(0, 10, 20, 30))
step_rep_summary$group <- ifelse(step_rep_summary$dif_step == 0, "sequence_discordance",
                                 ifelse(step_rep_summary$dif_step > 0, "positive_length_discordance", "negative_length_discordance"))

# Plot step difference (DP=0)
p_step_rep0 <- ggplot(step_rep_summary[step_rep_summary$DPthreshold == 0, ],
                      aes(x = dif_step, y = counts, fill = group)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = counts), size = 3.5, vjust = -0.2) +
  scale_fill_manual(values = c("black", "#DC1623", "#2D6DB1")) +
  scale_x_continuous(limits = c(-0.5, 10.5), breaks = 0:10) +
  theme_bw() +
  labs(title = "Step Difference among Partial-identical Genotype Pairs (RNA-RNA)", x = "Step Difference", y = "Counts") +
  theme(legend.position = "none")
ggsave(p_step_rep0, file = "stepdif_halfidentical_plot.pdf", width = 10, height = 6)

# For DP>0
p_step_rep_other <- ggplot(step_rep_summary[step_rep_summary$DPthreshold != 0, ],
                           aes(x = dif_step, y = counts, fill = group)) +
  geom_bar(stat = "identity") +
  facet_wrap(~DPthreshold, nrow = 3) +
  scale_fill_manual(values = c("black", "#DC1623", "#2D6DB1")) +
  scale_x_continuous(limits = c(-0.5, 10.5), breaks = 0:10) +
  theme_bw() +
  labs(title = "Step Difference among Partial-identical Genotype Pairs (RNA-RNA) with DP filters", x = "Step Difference", y = "Counts") +
  theme(legend.position = "none")
ggsave(p_step_rep_other, file = "stepdif_halfidentical_plot2.pdf", width = 10, height = 6)

# Isoallele variant (IAV) analysis for RNA-RNA
iav_locus <- as.data.frame(table(halfidentical_repeated$locus))
colnames(iav_locus) <- c("locus", "halfidentical_counts")
tmp <- as.data.frame(table(halfidentical_repeated$locus[halfidentical_repeated$dif_step == 0]))
iav_locus$iav_counts <- tmp$Freq[match(iav_locus$locus, tmp$Var1)]
iav_locus$iav_counts[is.na(iav_locus$iav_counts)] <- 0
iav_locus$percentage <- iav_locus$iav_counts / iav_locus$halfidentical_counts

# Detailed sequence differences for IAVs
iav_sequences <- halfidentical_repeated[halfidentical_repeated$dif_step == 0,
                                        c("locus", "ncopy1", "ncopy2", "sequence1", "sequence2")]
iav_diff_df <- data.frame()
for (i in 1:nrow(iav_sequences)) {
  cat(round(i / nrow(iav_sequences) * 100, 2), "%\n")
  s1 <- unlist(strsplit(iav_sequences$sequence1[i], ""))
  s2 <- unlist(strsplit(iav_sequences$sequence2[i], ""))
  mismatches <- which(s1 != s2)
  if (length(mismatches) > 0) {
    for (pos in mismatches) {
      base_pair <- paste(sort(c(s1[pos], s2[pos])), collapse = "_")
      iav_diff_df <- rbind(iav_diff_df, data.frame(locus = iav_sequences$locus[i],
                                                   allele_length = nchar(iav_sequences$sequence1[i]),
                                                   dif_location = pos,
                                                   base1 = s1[pos],
                                                   base2 = s2[pos],
                                                   changed_base_pattern = base_pair,
                                                   dif_location_percentage = pos / nchar(iav_sequences$sequence1[i])))
    }
  }
}
iav_diff_df$changed_base_pattern <- gsub("_", "|", iav_diff_df$changed_base_pattern)

# PHP count per IAV locus-allele
iav_php_per_locus <- aggregate(dif_location ~ paste(locus, allele_length, sep = "_"), data = iav_diff_df, function(x) length(unique(x)))
colnames(iav_php_per_locus) <- c("locus_len", "PHP_count")
p_php <- ggplot(iav_php_per_locus, aes(x = PHP_count)) +
  geom_histogram(fill = "#2D6DB1", bins = max(iav_php_per_locus$PHP_count), alpha = 0.6) +
  scale_x_continuous(breaks = 1:max(iav_php_per_locus$PHP_count)) +
  labs(y = "Genotype Pair Counts", x = "Base Variant Counts") +
  theme_bw()
ggsave(p_php, file = "PHPcounts_histogram.pdf", width = 6, height = 3)

# Base change pattern pie charts per region (for RNA-RNA)
# This part replicates the original pie chart generation; we'll simplify but keep essential.
iav_changed_base_pattern <- data.frame()
for (reg in c(unique(pSTR_group$final_group), "Total")) {
  if (reg != "Total") {
    locus_sub <- pSTR_group$locus[pSTR_group$final_group == reg]
    sub <- iav_diff_df[iav_diff_df$locus %in% locus_sub, ]
  } else {
    sub <- iav_diff_df
  }
  # Keep only single-PHP events (for fair comparison)
  php_counts <- aggregate(dif_location ~ locus + allele_length, data = sub, function(x) length(unique(x)))
  single_php_loci <- paste(php_counts$locus, php_counts$allele_length, sep = "_")[php_counts$dif_location == 1]
  sub_single <- sub[paste(sub$locus, sub$allele_length, sep = "_") %in% single_php_loci, ]
  if (nrow(sub_single) > 0) {
    tab <- as.data.frame(table(sub_single$changed_base_pattern))
    colnames(tab) <- c("pattern", "counts")
    tab$percentage <- tab$counts / sum(tab$counts)
    tab$label <- paste0(tab$counts, " (", round(tab$percentage * 100, 2), "%)")
    tab$region <- reg
    iav_changed_base_pattern <- rbind(iav_changed_base_pattern, tab)
  }
}
iav_changed_base_pattern$region <- factor(iav_changed_base_pattern$region,
                                          levels = c("CDS", "5_UTR", "3_UTR", "intron", "intergenic", "Total"))
p_pie_rep <- ggplot(iav_changed_base_pattern, aes(x = "", y = percentage, fill = pattern)) +
  geom_col(color = "black") +
  geom_label(aes(label = label), size = 2, position = position_stack(vjust = 0.5), show.legend = FALSE) +
  facet_wrap(~region, nrow = 2) +
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_rickandmorty() +
  labs(fill = "Base Variant Patterns")
ggsave(p_pie_rep, file = "iav_changed_base_pattern_pieplot.pdf", width = 15, height = 12)

# Base variant pattern density along allele position (RNA-RNA)
p_density_rep <- ggplot(iav_diff_df[paste(iav_diff_df$locus, iav_diff_df$allele_length, sep = "_") %in% single_php_loci, ],
                        aes(x = dif_location_percentage, after_stat(count), color = changed_base_pattern)) +
  geom_density(linewidth = 1) +
  scale_color_rickandmorty() +
  scale_x_continuous(limits = c(0, 1)) +
  labs(title = "Base Variant Pattern across Allele Positions (RNA-RNA)", y = "Density", x = "Allele Position") +
  theme_bw()
ggsave(p_density_rep, file = "iav_changed_base_pattern_regions.pdf", width = 8, height = 5)

