# =============================================================================
# polymorphism.R  —  clean, revision build
# Manuscript: Section 2.5 (forensic informativeness of pSTRs)
#
# Forensic-parameter math (He, Ho, MP, DP, PE2, PE3, PIC, Ae) is kept VERBATIM
# from the original — all of these were verified correct (PE2/PE3 still to be
# checked against the cited reference). Changes vs original:
#   * FIX: empty-genotype loci returned "" (character) for the numeric columns,
#          which coerced whole numeric columns to character on rbind. Now they
#          return NA_real_ (and allele_freq returns NULL), so types are clean.
#   * FIX: reference_ncopy_PIC_corplot mapped x = Ncopy (major-allele copy),
#          duplicating the major-allele plot; it now maps x = reference_ncopy.
#   * Shared theme; progress print() calls removed.
#   * Section "REVISION VALUES" writes traceability files for §2.5.
# =============================================================================

rm(list = ls())

library(tidyr)
library(vcfR)
library(data.table)
library(dplyr)
library(readxl)
library(ggplot2)
library(ggpubr)
library(scales)
library(EnvStats)
library(ggsci)

base_theme <- theme_bw() +
  theme(axis.title   = element_text(face = "bold", colour = "black", size = 15),
        axis.title.x = element_text(face = "bold", colour = "black", size = 12),
        axis.title.y = element_text(face = "bold", colour = "black", size = 12),
        axis.text.x  = element_text(face = "bold", colour = "black", size = 10),
        axis.text.y  = element_text(face = "bold", colour = "black", size = 10))

# =============================================================================
# Allele sequence ladder + copy-number genotype table
# =============================================================================
raw_vcf <- read.vcfR("../pstr.rna.recode.vcf.gz")
config  <- read.table("../GRCh38.hipstr_reference.refine.bed")

allele_ladder <- data.frame("locus" = raw_vcf@fix[, 3],
                            "Ref" = raw_vcf@fix[, 4],
                            "Alt" = raw_vcf@fix[, 5])
allele_ladder <- na.omit(gather(allele_ladder, key = "allele_group", value = "sequence", -locus))
temp1 <- allele_ladder[!grepl(",", allele_ladder$sequence), ]
temp2 <- allele_ladder[grepl(",", allele_ladder$sequence), ]
temp3 <- do.call(rbind, apply(temp2, 1, function(x) {
  alt_allele <- unlist(strsplit(x[3], ","))
  data.frame("locus" = x[1], "allele_group" = 1:length(alt_allele), "sequence" = alt_allele, row.names = NULL)
}))
temp3 <- temp3[nchar(temp3$sequence) != 0, ]   # individual malformed alt (Human_STR_904666)
temp1 <- rbind(temp1, temp3)
temp1$allele_group[temp1$allele_group == "Ref"] <- 0
temp1$allele_group[temp1$allele_group == "Alt"] <- 1
temp1$period_size <- as.numeric(config[match(temp1$locus, config[, 6]), 4])
temp2 <- nchar(temp1[, 3]) %/% temp1[, 4]
temp3 <- nchar(temp1[, 3]) %% temp1[, 4]
temp3[temp3 != 0] <- paste(".", temp3[temp3 != 0], sep = "")
temp3[temp3 == 0] <- ""
temp1$ncopy <- paste(temp2, temp3, sep = "")
temp1$ncopy_label <- paste(temp1$locus, temp1$ncopy, sep = "@")
temp2 <- as.data.frame(table(temp1$ncopy_label))
temp2[, 1] <- as.character(temp2[, 1])
temp2 <- temp2[temp2$Freq != 1, 1]                       # same-length, different-sequence alleles
setDT(temp1)
temp1[!temp1$ncopy_label %in% temp2, ncopy_label := ncopy]
temp1[temp1$ncopy_label %in% temp2, ncopy_label := { paste(ncopy, seq_along(ncopy), sep = "_") }, by = ncopy_label]
temp1$unique_label <- paste(temp1$locus, temp1$allele_group, sep = "@")
allele_ladder <- temp1

# GT table (allele-index) -> copy-number labels
gt_df <- as.data.frame(extract.gt(raw_vcf, element = "GT"))
gt_df <- data.frame("locus" = rownames(gt_df), gt_df)
gt_df <- na.omit(gather(gt_df, key = "sample", value = "GT", -locus))
gt_df <- separate(gt_df, col = "GT", into = c("allele1", "allele2"), sep = "\\|")
gt_df[as.numeric(gt_df$allele1) > as.numeric(gt_df$allele2), c("allele1", "allele2")] <-
  gt_df[as.numeric(gt_df$allele1) > as.numeric(gt_df$allele2), c("allele2", "allele1")]   # keep lower index first
gt_df$A1_label <- paste(gt_df$locus, gt_df$allele1, sep = "@")
gt_df$A2_label <- paste(gt_df$locus, gt_df$allele2, sep = "@")
gt_df$ncopy1 <- unlist(allele_ladder[match(gt_df$A1_label, allele_ladder$unique_label), "ncopy_label"])
gt_df$ncopy2 <- unlist(allele_ladder[match(gt_df$A2_label, allele_ladder$unique_label), "ncopy_label"])
gt_df <- gt_df[, c("locus", "sample", "ncopy1", "ncopy2")]
colnames(gt_df) <- c("locus", "sample", "allele1", "allele2")
gt_df <- unite(gt_df, "GT", c("allele1", "allele2"), sep = ",")
gt_df <- spread(gt_df, key = "sample", value = "GT")

allele_ladder <- allele_ladder[, c("locus", "period_size", "allele_group", "sequence", "ncopy", "ncopy_label")]
allele_ladder <- allele_ladder[order(allele_ladder$locus), ]
write.table(allele_ladder, "allele_sequence.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(allele_ladder[!grepl("STR", allele_ladder$locus), ],
            "allele_sequence.reported.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(gt_df, "GT.copynumber.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
rm(raw_vcf)

# =============================================================================
# Forensic-parameter functions  (math kept verbatim)
# =============================================================================
identify_het <- function(x) {
  allele <- unlist(strsplit(x[1], ","))
  if (length(unique(allele)) == 1) "Homo" else "Het"
}
MPvalue <- function(freq = "") {
  if (length(freq) == 1) return(0)
  Frequency <- c()
  for (i in 1:length(freq)) {
    Frequency <- c(Frequency, freq[i] * freq[i])
    for (j in (1 + i):length(freq)) {
      if (j > length(freq)) break
      Frequency <- c(Frequency, 2 * freq[i] * freq[j])
    }
  }
  sum(Frequency * Frequency)
}
DP_value <- function(x) {
  if (length(x) == 1) return(0)
  p1 <- 0; p2 <- 0
  for (i in seq_along(x)) { p1 <- p1 + (x[i])^2; p2 <- p2 + (x[i])^4 }
  1 - 2 * p1^2 + p2
}
PE2_value <- function(x) {
  if (length(x) == 1) return(0)
  p1 <- 0; p2 <- 0
  for (i in seq_along(x)) {
    p1 <- p1 + ((x[i])^2) * ((1 - x[i])^2)
    if (i != length(x)) for (j in (i + 1):length(x)) p2 <- p2 + 2 * x[i] * x[j] * ((1 - x[i] - x[j])^2)
  }
  p1 + p2
}
PE3_value <- function(x) {
  if (length(x) == 1) return(0)
  p1 <- 0; p2 <- 0
  for (i in seq_along(x)) {
    p1 <- p1 + x[i] * ((1 - x[i])^2)
    if (i != length(x)) for (j in (i + 1):length(x)) p2 <- p2 + 0.5 * (x[i]^2) * (x[j]^2) * (4 - 3 * x[i] - 3 * x[j])
  }
  p1 - p2
}
PIC_value <- function(x) {
  if (length(x) == 1) return(0)
  p1 <- 0; p2 <- 0
  for (i in seq_along(x)) { p1 <- p1 + (x[i])^2; p2 <- p2 + (x[i])^4 }
  1 - p1 - p1^2 + p2
}
Ae_value <- function(x) {
  if (length(x) == 1) return(0)
  p1 <- 0
  for (i in seq_along(x)) p1 <- p1 + (x[i])^2
  1 / p1
}

forensic_parameters_calculator <- function(x) {
  locus <- x[1]
  GT <- na.omit(x[-1])
  if (length(GT) == 0)
    return(data.frame(locus = locus, N = 0L, Na = 0L, commonNa = 0L,
                      Ho = NA_real_, He = NA_real_, MP = NA_real_, DP = NA_real_,
                      PE2 = NA_real_, PE3 = NA_real_, PIC = NA_real_, Ae = NA_real_, row.names = NULL))
  alleles <- unlist(strsplit(GT, ","))
  allele_number <- length(alleles)                       # N (allele observations)
  allele_type   <- length(unique(alleles))               # Na
  allele_table  <- as.data.frame(table(alleles))
  allele_freq   <- allele_table[, 2] / sum(allele_table[, 2])
  common_allele_type <- length(allele_freq[allele_freq > 0.01])
  GT_table  <- as.data.frame(table(GT))
  Het_grepl <- apply(GT_table, 1, identify_het)
  Observe_Het <- round(sum(GT_table[which(Het_grepl == "Het"), 2]) / sum(GT_table[, 2]), 4)
  Expect_Het  <- round(allele_number / (allele_number - 1) * (1 - sum((allele_table[, 2] / sum(allele_table[, 2]))^2)), 4)
  data.frame(locus = locus, N = allele_number, Na = allele_type, commonNa = common_allele_type,
             Ho = Observe_Het, He = Expect_Het,
             MP = round(MPvalue(allele_freq), 4), DP = round(DP_value(allele_freq), 4),
             PE2 = round(PE2_value(allele_freq), 4), PE3 = round(PE3_value(allele_freq), 4),
             PIC = round(PIC_value(allele_freq), 4), Ae = round(Ae_value(allele_freq), 4), row.names = NULL)
}
allele_freq_calculator <- function(x) {
  locus <- x[1]
  GT <- na.omit(x[-1])
  if (length(GT) == 0) return(NULL)                      # no alleles to report
  alleles <- unlist(strsplit(GT, ","))
  allele_table <- as.data.frame(table(alleles))
  data.frame(locus = locus, allele = as.character(allele_table[, 1]),
             frequency = round(allele_table[, 2] / sum(allele_table[, 2]), 3), row.names = NULL)
}

# =============================================================================
# Per-locus forensic parameters + allele frequencies (unrelated individuals)
# =============================================================================
unrelated_individuals <- read_excel("../VBsampleinfo.xlsx", sheet = "unrelated_individuals")
colnames(unrelated_individuals) <- unrelated_individuals[1, ]
unrelated_individuals <- unrelated_individuals[-1, ]
individual_gt <- gt_df[, colnames(gt_df) %in% c("locus", unrelated_individuals$ID)]

fp_df <- do.call(rbind, apply(individual_gt, 1, forensic_parameters_calculator))
write.table(fp_df, "forensic_parameters.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
allele_freq_df <- do.call(rbind, apply(individual_gt, 1, allele_freq_calculator))
write.table(allele_freq_df, "allele_freq.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# reported (named marker) loci
pSTR_group <- read.table("../2_annotation/STR_annotation_group.txt", header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
reported_locus <- allele_freq_df[!grepl("Human", allele_freq_df$locus), ]
reported_locus$sequence <- allele_ladder$sequence[match(paste(reported_locus$locus, reported_locus$allele, sep = "@"),
                                                        paste(allele_ladder$locus, allele_ladder$ncopy_label, sep = "@"))]
reported_locus$group <- pSTR_group$final_group[match(reported_locus$locus, pSTR_group$locus)]
reported_locus <- reported_locus[, c("locus", "group", "allele", "sequence", "frequency")]
write.table(reported_locus, "reported_locus.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(fp_df[fp_df$locus %in% reported_locus$locus, ], "forensic_parameters.reported_locus.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# =============================================================================
# Na, commonNa, PIC by region x motif
# =============================================================================
polymorphism_region <- pSTR_group[, c("locus", "final_group")]
colnames(polymorphism_region)[2] <- "group"
polymorphism_region$period_size <- config[match(polymorphism_region$locus, config[, 6]), 4]
polymorphism_region$Na       <- fp_df$Na[match(polymorphism_region$locus, fp_df$locus)]
polymorphism_region$commonNa <- fp_df$commonNa[match(polymorphism_region$locus, fp_df$locus)]
polymorphism_region$PIC      <- fp_df$PIC[match(polymorphism_region$locus, fp_df$locus)]
polymorphism_region$group <- factor(polymorphism_region$group, levels = c("CDS", "5_UTR", "3_UTR", "intron", "intergenic"))
polymorphism_region$period_size <- factor(polymorphism_region$period_size, levels = 2:6, labels = c("Di", "Tri", "Tetra", "Penta", "Hexa"))
polymorphism_region <- na.omit(polymorphism_region)

region_na_boxplot <- ggplot(polymorphism_region, aes(group, Na, fill = period_size)) +
  geom_boxplot(position = position_dodge(0.9)) + scale_fill_lancet() +
  stat_mean_sd_text(digits = 3, label.size = 0.6) + base_theme +
  labs(title = "The Na of pSTRs", x = "Region", y = "Na")
ggsave("region_na_boxplot.pdf", region_na_boxplot, width = 12, height = 6)

region_pic_boxplot <- ggplot(polymorphism_region, aes(group, PIC, fill = period_size)) +
  geom_violin(trim = TRUE) + geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  scale_fill_lancet() + stat_mean_sd_text(digits = 3, label.size = 0.6) + base_theme +
  labs(title = "The PIC of pSTRs", x = "Region", y = "PIC")
ggsave("region_pic_boxplot.pdf", region_pic_boxplot, width = 12, height = 6)

# =============================================================================
# Major allele frequency + copy-number difference analyses
# =============================================================================
allele_freq_df$period_size <- config[match(allele_freq_df$locus, config[, 6]), 4]
temp4 <- c(); temp5 <- c()
for (i in 1:nrow(allele_freq_df)) {
  temp1 <- allele_freq_df$allele[i]
  if (grepl("_", temp1)) temp1 <- as.numeric(unlist(strsplit(temp1, "_"))[1])
  temp4 <- c(temp4, temp1)
  if (grepl("\\.", temp1)) {
    temp1 <- as.character(temp1)
    temp2 <- as.numeric(unlist(strsplit(temp1, "\\.")))
    temp3 <- temp2[1] * allele_freq_df$period_size[i] + temp2[2]
  } else {
    temp1 <- as.numeric(temp1)
    temp3 <- temp1[1] * allele_freq_df$period_size[i]
  }
  temp5 <- c(temp5, temp3)
}
allele_freq_df$Ncopy <- temp4
allele_freq_df$allele_length <- temp5

major_allele_df <- data.frame()
for (i in unique(allele_freq_df$locus)) {
  temp1 <- allele_freq_df[allele_freq_df$locus == i, ]
  temp2 <- data.frame()
  for (p in unique(temp1$Ncopy)) {
    temp2 <- rbind(temp2, data.frame("locus" = i, "Ncopy" = p, "freq" = sum(temp1$frequency[temp1$Ncopy == p])))
  }
  temp2 <- temp2[order(temp2$freq, decreasing = TRUE), ]
  temp2$Ncopy <- as.numeric(temp2$Ncopy)
  temp2$group <- NA
  temp2$group[1] <- "major_allele"
  temp2$group[is.na(temp2$group)] <- "alt_allele"
  temp2$ncopy_dif <- 0
  if (nrow(temp2[temp2$group == "alt_allele", ]) != 0)
    temp2$ncopy_dif[temp2$group == "alt_allele"] <- temp2$Ncopy[temp2$group == "alt_allele"] - temp2$Ncopy[temp2$group == "major_allele"]
  major_allele_df <- rbind(major_allele_df, temp2)
}
major_allele_df$freq[major_allele_df$freq > 1] <- 1

major_allele_freq_histogram <- ggplot(major_allele_df[major_allele_df$group == "major_allele", ], aes(freq)) +
  geom_histogram(colour = "black", fill = "gray", alpha = 0.3, position = "identity", bins = 20) +
  scale_y_continuous(trans = log10_trans()) + base_theme +
  labs(title = "The Frequency Distribution of Major Allele", x = "Frequency", y = "Counts")
ggsave("major_allele_freq_histogram.pdf", major_allele_freq_histogram, width = 6, height = 4)

# major vs other alleles: copy-number difference (integer diffs within +/-10)
# Counts are UNWEIGHTED distinct alt alleles (one per locus x alt-Ncopy). This matches
# the published Fig 4c (single 78.96% / shorter 67.86% / longer 32.14%); occurrence-
# weighting (round(freq*N)) gives a different split and zeroes out rare alleles. See
# revision_values [MINOR-ALLELE STEP DISTRIBUTION] for all three weightings.
temp1 <- major_allele_df[major_allele_df$ncopy_dif %% 1 == 0 & major_allele_df$ncopy_dif <= 10 &
                           major_allele_df$ncopy_dif >= -10 & major_allele_df$group == "alt_allele", ]
temp2 <- data.frame()
for (i in unique(temp1$ncopy_dif))
  temp2 <- rbind(temp2, data.frame("ncopy_dif" = i, "counts" = nrow(temp1[temp1$ncopy_dif == i, ])))
temp2$percentage <- temp2$counts / sum(temp2$counts)
temp2$group[temp2$ncopy_dif > 0] <- "positive"
temp2$group[temp2$ncopy_dif < 0] <- "negative"
temp2$ncopy_dif <- factor(temp2$ncopy_dif, levels = c(seq(-10, -1, 1), seq(1, 10, 1)),
                          labels = as.character(c(seq(-10, -1, 1), seq(1, 10, 1))))
major_alt_ncopy_dif_barplot <- ggplot(temp2, aes(ncopy_dif)) +
  geom_bar(aes(y = counts), stat = "identity", colour = "black", alpha = 0.3) +
  geom_line(aes(y = percentage * max(counts) * 2, group = 1), linewidth = 0.8) +
  geom_point(aes(y = percentage * max(counts) * 2, group = 1)) +
  scale_y_continuous(name = "Allele Counts", sec.axis = sec_axis(trans = ~ . / (max(temp2$counts) * 2), name = "Percentage")) +
  base_theme + theme(legend.position = "none") +
  labs(title = "The Copy Number Difference between Major Allele and Other Alleles", x = "Copy Number Difference")
ggsave("major_alt_ncopy_dif_barplot.pdf", major_alt_ncopy_dif_barplot, width = 6, height = 4)

# major allele vs reference: copy-number difference (integer diffs within +/-10)
mar <- major_allele_df[major_allele_df$group == "major_allele", ]
mar$N <- fp_df[match(mar$locus, fp_df$locus), "N"]
mar$reference_length <- config[match(mar$locus, config[, 6]), 3] - config[match(mar$locus, config[, 6]), 2] + 1
mar$period_size <- config[match(mar$locus, config[, 6]), 4]
mar$reference_ncopy <- as.numeric(paste(mar$reference_length %/% mar$period_size, mar$reference_length %% mar$period_size, sep = "."))
mar$ncopy_dif <- mar$Ncopy - mar$reference_ncopy
mar <- mar[mar$ncopy_dif >= -10 & mar$ncopy_dif <= 10 & mar$ncopy_dif %% 1 == 0, ]
mar$PIC <- fp_df[match(mar$locus, fp_df$locus), "PIC"]
temp2 <- data.frame()
for (i in unique(mar$ncopy_dif))
  temp2 <- rbind(temp2, data.frame("ncopy_dif" = i, "counts" = nrow(mar[mar$ncopy_dif == i, ])))
temp2$percentage <- temp2$counts / sum(temp2$counts)
temp2$group[temp2$ncopy_dif == 0] <- "identical"
temp2$group[temp2$ncopy_dif > 0]  <- "positive"
temp2$group[temp2$ncopy_dif < 0]  <- "negative"
temp2$ncopy_dif <- factor(temp2$ncopy_dif, levels = seq(-10, 10, 1), labels = as.character(seq(-10, 10, 1)))
major_reference_ncopy_dif_barplot <- ggplot(temp2, aes(ncopy_dif, fill = group)) +
  geom_bar(aes(y = counts), stat = "identity", colour = "black", alpha = 0.3) +
  geom_line(aes(y = percentage * max(counts), group = 1), linewidth = 0.8) +
  geom_point(aes(y = percentage * max(counts), group = 1)) +
  scale_fill_manual(values = c("gray", "gray", "gray")) +
  scale_y_continuous(name = "Counts", sec.axis = sec_axis(trans = ~ . / (max(temp2$counts)), name = "Percentage")) +
  base_theme + theme(legend.position = "none") +
  labs(title = "The Copy Number Difference Between Major Allele and Reference Allele", x = "Copy Number Difference", y = "Counts")
ggsave("major_reference_ncopy_dif_barplot.pdf", major_reference_ncopy_dif_barplot, width = 6, height = 4)

# PIC vs copy number (major allele; reference allele)
major_ncopy_PIC_corplot <- ggplot(mar, aes(Ncopy, PIC)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "#2D6DB1", fill = "#756bb1") + stat_cor() + base_theme +
  labs(title = "Correlation between PIC and Copy Number of Major Allele", x = "Copy Number of Major Allele", y = "PIC")
ggsave("major_ncopy_PIC_corplot.pdf", major_ncopy_PIC_corplot, width = 7, height = 3.5)

reference_ncopy_PIC_corplot <- ggplot(mar, aes(reference_ncopy, PIC)) +   # FIX: was Ncopy
  geom_smooth(method = "lm", formula = y ~ x, colour = "#2D6DB1", fill = "#756bb1") + stat_cor() + base_theme +
  labs(title = "Correlation between PIC and Copy Number of Reference Allele", x = "Copy Number of Reference Allele", y = "PIC")
ggsave("reference_ncopy_PIC_corplot.pdf", reference_ncopy_PIC_corplot, width = 7, height = 3.5)

# =============================================================================
# REVISION VALUES  ->  written to files for the writing list
# =============================================================================
poly <- fp_df[!is.na(fp_df$PIC) & fp_df$Na >= 2, ]      # polymorphic loci (Na >= 2)
num_cols <- c("Na", "commonNa", "Ho", "He", "MP", "DP", "PE2", "PE3", "PIC", "Ae")

summ <- do.call(rbind, lapply(num_cols, function(cn) {
  v <- as.numeric(poly[[cn]])
  data.frame(parameter = cn, mean = round(mean(v, na.rm = TRUE), 4), median = round(median(v, na.rm = TRUE), 4),
             min = round(min(v, na.rm = TRUE), 4), max = round(max(v, na.rm = TRUE), 4))
}))
write.table(summ, "forensic_summary.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# cumulative panel power for the reported (named marker) loci
rep_fp <- fp_df[fp_df$locus %in% reported_locus$locus & !is.na(fp_df$MP), ]
comb_PM  <- prod(rep_fp$MP)
comb_PD  <- 1 - comb_PM
comb_PE2 <- 1 - prod(1 - rep_fp$PE2)
comb_PE3 <- 1 - prod(1 - rep_fp$PE3)

mar_cc  <- mar[!is.na(mar$PIC) & !is.na(mar$Ncopy) & !is.na(mar$reference_ncopy), ]
ct_major <- suppressWarnings(cor.test(mar_cc$Ncopy, mar_cc$PIC))
ct_ref   <- suppressWarnings(cor.test(mar_cc$reference_ncopy, mar_cc$PIC))
ref_dist <- as.data.frame(table(factor(sign(mar$ncopy_dif), levels = c(-1, 0, 1),
                                       labels = c("shorter", "equal", "longer"))))
ref_dist$pct <- 100 * ref_dist$Freq / sum(ref_dist$Freq)

# --- remaining §2.5 manuscript numbers (read-outs over existing objects) -----
maj_all      <- major_allele_df[major_allele_df$group == "major_allele", ]   # one row per locus
n_loci_total <- nrow(maj_all)                                                # 24,013
maj_ge05     <- sum(maj_all$freq >= 0.5)                                      # Fig 4a (98.78%)
maj_lt05     <- sum(maj_all$freq <  0.5)                                      # the TRUE "<0.5" count
n_filter_dropped <- n_loci_total - nrow(mar)                                 # 24013-23827 (= the real "186")
n_equal_all  <- sum(mar$ncopy_dif == 0)                                       # major == reference
pct_equal_all <- 100 * n_equal_all / n_loci_total                            # over ALL pSTR (95.76%)

# minor (length-differing) alleles, occurrence-weighted, integer diffs +/-10 (Fig 4c)
alt <- major_allele_df[major_allele_df$group == "alt_allele" &
                         major_allele_df$ncopy_dif %% 1 == 0 & abs(major_allele_df$ncopy_dif) <= 10, ]
alt$eN <- round(alt$freq * as.numeric(fp_df$N[match(alt$locus, fp_df$locus)]))
tot_eN <- sum(alt$eN)
pct_minor_single  <- 100 * sum(alt$eN[abs(alt$ncopy_dif) == 1]) / tot_eN
pct_minor_shorter <- 100 * sum(alt$eN[alt$ncopy_dif < 0]) / tot_eN
pct_minor_longer  <- 100 * sum(alt$eN[alt$ncopy_dif > 0]) / tot_eN
# alternative weightings (to identify which one Fig 4c actually used):
#   unweighted = each distinct alt allele counts once (keeps rare alleles that
#   round(freq*N) would zero out); freq_weighted = summed allele frequency.
n_alt        <- nrow(alt)
us_single    <- 100 * sum(abs(alt$ncopy_dif) == 1) / n_alt
us_shorter   <- 100 * sum(alt$ncopy_dif < 0) / n_alt
us_longer    <- 100 * sum(alt$ncopy_dif > 0) / n_alt
fw           <- sum(alt$freq)
fw_single    <- 100 * sum(alt$freq[abs(alt$ncopy_dif) == 1]) / fw
fw_shorter   <- 100 * sum(alt$freq[alt$ncopy_dif < 0]) / fw
fw_longer    <- 100 * sum(alt$freq[alt$ncopy_dif > 0]) / fw

na_v <- as.numeric(fp_df$Na)
n_multi <- sum(na_v >  2, na.rm = TRUE)                                       # Table S4: Na > 2
n_bi    <- sum(na_v == 2, na.rm = TRUE)                                       # Table S4: Na == 2
n_common_poly <- sum(as.numeric(fp_df$commonNa) >= 2, na.rm = TRUE)           # common-allele polymorphic

hp <- fp_df[na_v >= 3 & as.numeric(fp_df$PIC) >= 0.5, ]                       # highly polymorphic
hp <- hp[!is.na(hp$locus), ]
hp_motif <- table(factor(config[match(hp$locus, config[, 6]), 4], levels = 2:6,
                         labels = c("Di", "Tri", "Tetra", "Penta", "Hexa")))

pr <- polymorphism_region
pr$PIC <- as.numeric(pr$PIC); pr$Na <- as.numeric(pr$Na)
region_pic <- aggregate(PIC ~ group,       pr, function(x) mean(x, na.rm = TRUE))   # Fig 4e
na_motif   <- aggregate(Na  ~ period_size, pr, function(x) mean(x, na.rm = TRUE))   # Fig 4d

# PIC ~ copy number restricted to introns (Figure S11 scope check vs all-pSTR)
mar$region <- pSTR_group$final_group[match(mar$locus, pSTR_group$locus)]
mi <- mar[!is.na(mar$PIC) & mar$region == "intron", ]
ct_maj_in <- suppressWarnings(cor.test(as.numeric(mi$Ncopy),          as.numeric(mi$PIC)))
ct_ref_in <- suppressWarnings(cor.test(as.numeric(mi$reference_ncopy), as.numeric(mi$PIC)))

rv <- character(0); add <- function(...) rv[[length(rv) + 1]] <<- sprintf(...)
add("# polymorphism.R revision values  (auto-generated %s)", as.character(Sys.Date()))
add("# loci with forensic parameters: %d  |  genotyped (N>0): %d  |  polymorphic (Na>=2): %d",
    nrow(fp_df), sum(fp_df$N > 0, na.rm = TRUE), nrow(poly))
add("")
add("## [FORENSIC SUMMARY over polymorphic loci]  -> forensic_summary.txt")
for (i in seq_len(nrow(summ)))
  add("%s\tmean %.4f\tmedian %.4f\trange %.4f-%.4f", summ$parameter[i], summ$mean[i], summ$median[i], summ$min[i], summ$max[i])
add("")
add("## [REPORTED MARKER LOCI]  n=%d  -> forensic_parameters.reported_locus.txt", nrow(rep_fp))
add("combined PM (prod MP)\t%.3e", comb_PM)
add("combined PD (1-prod MP)\t%.6f", comb_PD)
add("combined PE2 (1-prod(1-PE2))\t%.6f", comb_PE2)
add("combined PE3 (1-prod(1-PE3))\t%.6f", comb_PE3)
add("")
add("## [PIC vs COPY NUMBER]  (report r; large N -> de-emphasize p)")
add("PIC ~ major-allele copy\tr=%.3f\tp=%.2e\tn=%d", unname(ct_major$estimate), ct_major$p.value, nrow(mar_cc))
add("PIC ~ reference copy (corrected)\tr=%.3f\tp=%.2e\tn=%d", unname(ct_ref$estimate), ct_ref$p.value, nrow(mar_cc))
add("")
add("## [MAJOR vs REFERENCE copy difference]")
for (i in seq_len(nrow(ref_dist)))
  add("%s\t%d\t%.2f%% (of %d in +/-10 integer set)", as.character(ref_dist$Var1[i]), ref_dist$Freq[i], ref_dist$pct[i], nrow(mar))
add("major == reference over ALL pSTR\t%d / %d\t%.2f%%", n_equal_all, n_loci_total, pct_equal_all)
add("")
add("## [MAJOR-ALLELE FREQUENCY]  (Fig 4a/4b)")
add("major-allele freq >= 0.5\t%d / %d\t%.2f%%", maj_ge05, n_loci_total, 100 * maj_ge05 / n_loci_total)
add("major-allele freq <  0.5 (TRUE count)\t%d\t%.2f%%", maj_lt05, 100 * maj_lt05 / n_loci_total)
add("loci dropped by +/-10 integer filter (= the manuscript's '186')\t%d", n_filter_dropped)
add("")
add("## [MINOR-ALLELE STEP DISTRIBUTION]  (Fig 4c; 3 weightings -> match to manuscript 78.96/67.86/32.14)")
add("occurrence-weighted (round freq*N)\tsingle %.2f%%\tshorter %.2f%%\tlonger %.2f%%", pct_minor_single, pct_minor_shorter, pct_minor_longer)
add("unweighted (distinct alt alleles)\tsingle %.2f%%\tshorter %.2f%%\tlonger %.2f%%", us_single, us_shorter, us_longer)
add("frequency-weighted (sum freq)\tsingle %.2f%%\tshorter %.2f%%\tlonger %.2f%%", fw_single, fw_shorter, fw_longer)
add("")
add("## [ALLELE-NUMBER CLASSES]  (Table S4)")
add("multiallelic (Na>2)\t%d\t%.2f%%", n_multi, 100 * n_multi / nrow(fp_df))
add("biallelic (Na==2)\t%d\t%.2f%%", n_bi, 100 * n_bi / nrow(fp_df))
add("polymorphic on common alleles (commonNa>=2)\t%d", n_common_poly)
add("highly polymorphic (Na>=3 & PIC>=0.5)\t%d", nrow(hp))
for (m in names(hp_motif)) add("  %s\t%d", m, as.integer(hp_motif[[m]]))
add("")
add("## [MEAN PIC BY REGION]  (Fig 4e)")
for (i in seq_len(nrow(region_pic))) add("%s\t%.4f", as.character(region_pic$group[i]), region_pic$PIC[i])
add("")
add("## [MEAN Na BY MOTIF]  (Fig 4d)")
for (i in seq_len(nrow(na_motif))) add("%s\t%.4f", as.character(na_motif$period_size[i]), na_motif$Na[i])
add("")
add("## [PIC ~ COPY NUMBER, INTRON-ONLY]  (Fig S11 scope; cf. all-pSTR above)")
add("intron PIC ~ major-allele copy\tr=%.3f\tp=%.2e\tn=%d", unname(ct_maj_in$estimate), ct_maj_in$p.value, nrow(mi))
add("intron PIC ~ reference copy\tr=%.3f\tp=%.2e\tn=%d", unname(ct_ref_in$estimate), ct_ref_in$p.value, nrow(mi))
writeLines(rv, "revision_values_polymorphism.txt")

cat("polymorphism.R done. empty-GT and reference-copy fixes applied; revision files written.\n")

source("fig4_revision_patch.R")
