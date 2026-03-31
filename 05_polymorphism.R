# ============================================================================
# Script: Polymorphism analysis and forensic parameter calculation for pSTRs
# Description:
#   1. Build allele sequence ladder from VCF and convert genotypes to copy number format.
#   2. Compute forensic parameters (Na, Ho, He, MP, DP, PE2, PE3, PIC, Ae) per locus.
#   3. Generate allele frequency tables.
#   4. Produce regional plots for Na and PIC.
#   5. Analyze major allele frequency and copy number differences between major/reference alleles.
# ============================================================================

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

# ----------------------------- Input data -----------------------------------
raw_vcf <- read.vcfR("../pstr.rna.recode.vcf.gz")
config <- read.table("../GRCh38.hipstr_reference.refine.bed")

# ----------------------------- Allele sequence ladder ------------------------
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
allele_ladder <- temp1

# ----------------------------- Convert GT to copy number format -------------
gt_df <- as.data.frame(extract.gt(raw_vcf, element = "GT"))
gt_df <- data.frame(locus = rownames(gt_df), gt_df)
gt_df <- na.omit(gather(gt_df, key = "sample", value = "GT", -locus))
gt_df <- separate(gt_df, col = "GT", into = c("allele1", "allele2"), sep = "\\|")
# Sort alleles to keep lower index first
idx <- as.numeric(gt_df$allele1) > as.numeric(gt_df$allele2)
gt_df[idx, c("allele1", "allele2")] <- gt_df[idx, c("allele2", "allele1")]

gt_df$A1_label <- paste(gt_df$locus, gt_df$allele1, sep = "@")
gt_df$A2_label <- paste(gt_df$locus, gt_df$allele2, sep = "@")

gt_df$ncopy1 <- unlist(allele_ladder[match(gt_df$A1_label, allele_ladder$unique_label), "ncopy_label"])
gt_df$ncopy2 <- unlist(allele_ladder[match(gt_df$A2_label, allele_ladder$unique_label), "ncopy_label"])

gt_df <- gt_df[, c("locus", "sample", "ncopy1", "ncopy2")]
colnames(gt_df) <- c("locus", "sample", "allele1", "allele2")
gt_df <- unite(gt_df, "GT", c("allele1", "allele2"), sep = ",")
gt_df <- spread(gt_df, key = "sample", value = "GT")

# Output allele sequence and GT copy number tables
allele_ladder_out <- allele_ladder[, c("locus", "period_size", "allele_group", "sequence", "ncopy", "ncopy_label")]
allele_ladder_out <- allele_ladder_out[order(allele_ladder_out$locus), ]
write.table(allele_ladder_out, file = "allele_sequence.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(allele_ladder_out[!grepl("STR", allele_ladder_out$locus), ],
            file = "allele_sequence.reported.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(gt_df, file = "GT.copynumber.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
rm(raw_vcf)

# ----------------------------- Forensic parameter functions ------------------
identify_het <- function(x) {
  alleles <- unlist(strsplit(x[1], ","))
  if (length(unique(alleles)) == 1) return("Homo") else return("Het")
}

MPvalue <- function(freq = "") {
  if (length(freq) == 1) return(0)
  Frequency <- c()
  for (i in seq_along(freq)) {
    Frequency <- c(Frequency, freq[i]^2)
    if (i < length(freq)) {
      for (j in (i + 1):length(freq)) {
        Frequency <- c(Frequency, 2 * freq[i] * freq[j])
      }
    }
  }
  mp <- sum(Frequency^2)
  return(mp)
}

DP_value <- function(x) {
  if (length(x) == 1) return(0)
  p1 <- sum(x^2)
  p2 <- sum(x^4)
  DP <- 1 - 2 * p1^2 + p2
  return(DP)
}

PE2_value <- function(x) {
  if (length(x) == 1) return(0)
  p1 <- sum((x^2) * ((1 - x)^2))
  p2 <- 0
  if (length(x) > 1) {
    for (i in seq_along(x)) {
      if (i < length(x)) {
        for (j in (i + 1):length(x)) {
          p2 <- p2 + 2 * x[i] * x[j] * ((1 - x[i] - x[j])^2)
        }
      }
    }
  }
  PE2 <- p1 + p2
  return(PE2)
}

PE3_value <- function(x) {
  if (length(x) == 1) return(0)
  p1 <- sum(x * ((1 - x)^2))
  p2 <- 0
  if (length(x) > 1) {
    for (i in seq_along(x)) {
      if (i < length(x)) {
        for (j in (i + 1):length(x)) {
          p2 <- p2 + 0.5 * (x[i]^2) * (x[j]^2) * (4 - 3 * x[i] - 3 * x[j])
        }
      }
    }
  }
  PE3 <- p1 - p2
  return(PE3)
}

PIC_value <- function(x) {
  if (length(x) == 1) return(0)
  p1 <- sum(x^2)
  p2 <- sum(x^4)
  PIC <- 1 - p1 - p1^2 + p2
  return(PIC)
}

Ae_value <- function(x) {
  if (length(x) == 1) return(0)
  p1 <- sum(x^2)
  Ae <- 1 / p1
  return(Ae)
}

forensic_parameters_calculator <- function(x) {
  locus <- x[1]
  GT <- na.omit(x[-1])
  if (length(GT) == 0) {
    return(data.frame(locus = locus, N = 0, Na = 0, commonNa = 0,
                      Ho = "", He = "", MP = "", DP = "", PE2 = "", PE3 = "", PIC = "", Ae = "",
                      stringsAsFactors = FALSE))
  }
  alleles <- unlist(strsplit(GT, ","))
  allele_number <- length(alleles)          # N
  allele_type <- length(unique(alleles))    # Na
  allele_table <- as.data.frame(table(alleles))
  allele_freq <- allele_table[, 2] / sum(allele_table[, 2])
  common_allele_type <- sum(allele_freq > 0.01)  # commonNa
  
  GT_table <- as.data.frame(table(GT))
  het_status <- apply(GT_table, 1, identify_het)
  Observe_Het <- round(sum(GT_table[het_status == "Het", 2]) / sum(GT_table[, 2]), 4)
  
  # Expected heterozygosity (unbiased)
  Expect_Het <- round(allele_number / (allele_number - 1) * (1 - sum((allele_table[, 2] / sum(allele_table[, 2]))^2)), 4)
  
  DP_val <- round(DP_value(allele_freq), 4)
  PE2_val <- round(PE2_value(allele_freq), 4)
  PE3_val <- round(PE3_value(allele_freq), 4)
  PIC_val <- round(PIC_value(allele_freq), 4)
  Ae_val <- round(Ae_value(allele_freq), 4)
  MP_val <- round(MPvalue(allele_freq), 4)
  
  return(data.frame(locus = locus, N = allele_number, Na = allele_type, commonNa = common_allele_type,
                    Ho = Observe_Het, He = Expect_Het, MP = MP_val, DP = DP_val,
                    PE2 = PE2_val, PE3 = PE3_val, PIC = PIC_val, Ae = Ae_val,
                    stringsAsFactors = FALSE))
}

allele_freq_calculator <- function(x) {
  locus <- x[1]
  GT <- na.omit(x[-1])
  if (length(GT) == 0) {
    return(data.frame(locus = locus, allele = "", frequency = "", stringsAsFactors = FALSE))
  }
  alleles <- unlist(strsplit(GT, ","))
  allele_table <- as.data.frame(table(alleles))
  allele_freq <- allele_table[, 2] / sum(allele_table[, 2])
  allele_df <- data.frame(locus = locus,
                          allele = as.character(allele_table[, 1]),
                          frequency = round(allele_freq, 3),
                          stringsAsFactors = FALSE)
  return(allele_df)
}

# ----------------------------- Compute parameters for unrelated individuals --
unrelated_individuals <- read_excel("../VBsampleinfo.xlsx", sheet = "unrelated_individuals")
colnames(unrelated_individuals) <- unrelated_individuals[1, ]
unrelated_individuals <- unrelated_individuals[-1, ]

individual_gt <- gt_df[, colnames(gt_df) %in% c("locus", unrelated_individuals$ID)]

fp_df <- do.call(rbind, apply(individual_gt, 1, forensic_parameters_calculator))
write.table(fp_df, file = "forensic_parameters.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

allele_freq_df <- do.call(rbind, apply(individual_gt, 1, allele_freq_calculator))
write.table(allele_freq_df, file = "allele_freq.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# ----------------------------- Reported loci (non-Human_STR) -----------------
pSTR_group <- read.table("../2_annotation/STR_annotation_group.txt", header = TRUE)

reported_locus <- allele_freq_df[!grepl("Human", allele_freq_df$locus), ]
reported_locus$sequence <- allele_ladder_out$sequence[match(paste(reported_locus$locus, reported_locus$allele, sep = "@"),
                                                            paste(allele_ladder_out$locus, allele_ladder_out$ncopy_label, sep = "@"))]
reported_locus$group <- pSTR_group$final_group[match(reported_locus$locus, pSTR_group$locus)]
reported_locus <- reported_locus[, c("locus", "group", "allele", "sequence", "frequency")]
write.table(reported_locus, file = "reported_locus.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(fp_df[fp_df$locus %in% reported_locus$locus, ],
            file = "forensic_parameters.reported_locus.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# ----------------------------- Downstream: Na and PIC by region -------------
polymorphism_region <- pSTR_group[, c("locus", "final_group")]
colnames(polymorphism_region)[2] <- "group"
polymorphism_region$period_size <- config[match(polymorphism_region$locus, config[, 6]), 4]
polymorphism_region$Na <- fp_df$Na[match(polymorphism_region$locus, fp_df$locus)]
polymorphism_region$commonNa <- fp_df$commonNa[match(polymorphism_region$locus, fp_df$locus)]
polymorphism_region$PIC <- fp_df$PIC[match(polymorphism_region$locus, fp_df$locus)]
polymorphism_region$group <- factor(polymorphism_region$group,
                                    levels = c("CDS", "5_UTR", "3_UTR", "intron", "intergenic"))
polymorphism_region$period_size <- factor(polymorphism_region$period_size, levels = 2:6,
                                          labels = c("Di", "Tri", "Tetra", "Penta", "Hexa"))
polymorphism_region <- na.omit(polymorphism_region)

region_na_boxplot <- ggplot(polymorphism_region, aes(x = group, y = Na, fill = period_size)) +
  geom_boxplot(position = position_dodge(0.9)) +
  scale_fill_lancet() +
  stat_mean_sd_text(digits = 3, label.size = 0.6) +
  theme_bw() +
  labs(title = "The Na of pSTRs", x = "Region", y = "Na") +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 10))
ggsave(region_na_boxplot, file = "region_na_boxplot.pdf", width = 12, height = 6)

region_pic_boxplot <- ggplot(polymorphism_region, aes(x = group, y = PIC, fill = period_size)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  scale_fill_lancet() +
  stat_mean_sd_text(digits = 3, label.size = 0.6) +
  theme_bw() +
  labs(title = "The PIC of pSTRs", x = "Region", y = "PIC") +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 10))
ggsave(region_pic_boxplot, file = "region_pic_boxplot.pdf", width = 12, height = 6)

# ----------------------------- Major allele frequency analysis ---------------
allele_freq_df$period_size <- config[match(allele_freq_df$locus, config[, 6]), 4]

# Compute Ncopy (repeat count) and allele length
ncopy_vec <- numeric()
len_vec <- numeric()
for (i in 1:nrow(allele_freq_df)) {
  allele_str <- allele_freq_df$allele[i]
  if (grepl("_", allele_str)) {
    ncopy_val <- as.numeric(unlist(strsplit(allele_str, "_"))[1])
  } else {
    ncopy_val <- as.numeric(allele_str)
  }
  ncopy_vec <- c(ncopy_vec, ncopy_val)
  
  if (grepl("\\.", allele_str)) {
    parts <- as.numeric(unlist(strsplit(as.character(allele_str), "\\.")))
    len_val <- parts[1] * allele_freq_df$period_size[i] + parts[2]
  } else {
    len_val <- ncopy_val * allele_freq_df$period_size[i]
  }
  len_vec <- c(len_vec, len_val)
}
allele_freq_df$Ncopy <- ncopy_vec
allele_freq_df$allele_length <- len_vec

# Aggregate frequencies by copy number per locus
major_allele_df <- data.frame()
for (l in unique(allele_freq_df$locus)) {
  sub <- allele_freq_df[allele_freq_df$locus == l, ]
  agg <- data.frame()
  for (cn in unique(sub$Ncopy)) {
    agg <- rbind(agg, data.frame(locus = l, Ncopy = cn,
                                 freq = sum(sub$frequency[sub$Ncopy == cn])))
  }
  agg <- agg[order(agg$freq, decreasing = TRUE), ]
  agg$Ncopy <- as.numeric(agg$Ncopy)
  agg$group <- "alt_allele"
  agg$group[1] <- "major_allele"
  agg$ncopy_dif <- 0
  if (nrow(agg[agg$group == "alt_allele", ]) > 0) {
    agg$ncopy_dif[agg$group == "alt_allele"] <- agg$Ncopy[agg$group == "alt_allele"] - agg$Ncopy[agg$group == "major_allele"]
  }
  major_allele_df <- rbind(major_allele_df, agg)
}
major_allele_df$freq[major_allele_df$freq > 1] <- 1

# Histogram of major allele frequency
p_major_freq <- ggplot(major_allele_df[major_allele_df$group == "major_allele", ], aes(x = freq)) +
  geom_histogram(color = "black", fill = "gray", alpha = 0.3, bins = 20) +
  scale_y_continuous(trans = log10_trans()) +
  labs(title = "The Frequency Distribution of Major Allele", x = "Frequency", y = "Counts") +
  theme_bw()
ggsave(p_major_freq, file = "major_allele_freq_histogram.pdf", width = 6, height = 4)

# Copy number difference between major allele and other alleles (within ±10, integer)
alt_sub <- major_allele_df[major_allele_df$ncopy_dif %% 1 == 0 &
                             abs(major_allele_df$ncopy_dif) <= 10 &
                             major_allele_df$group == "alt_allele", ]
alt_sub$N <- fp_df$N[match(alt_sub$locus, fp_df$locus)]
alt_sub$eN <- alt_sub$freq * alt_sub$N

dif_summary <- data.frame()
for (dif in unique(alt_sub$ncopy_dif)) {
  sub2 <- alt_sub[alt_sub$ncopy_dif == dif, ]
  cnt <- sum(round(sub2$freq * sub2$N, 0))
  dif_summary <- rbind(dif_summary, data.frame(ncopy_dif = dif, counts = cnt))
}
dif_summary$percentage <- dif_summary$counts / sum(dif_summary$counts)
dif_summary$group <- ifelse(dif_summary$ncopy_dif > 0, "positive", "negative")
dif_summary$ncopy_dif <- factor(dif_summary$ncopy_dif,
                                levels = c(seq(-10, -1, 1), seq(1, 10, 1)),
                                labels = as.character(c(seq(-10, -1, 1), seq(1, 10, 1))))

p_dif_alt <- ggplot(dif_summary, aes(x = ncopy_dif)) +
  geom_bar(aes(y = counts), stat = "identity", color = "black", alpha = 0.3) +
  geom_line(aes(y = percentage * max(counts) * 2, group = 1), linewidth = 0.8) +
  geom_point(aes(y = percentage * max(counts) * 2, group = 1)) +
  scale_y_continuous(name = "Counts", sec.axis = sec_axis(trans = ~ . / (max(dif_summary$counts) * 2), name = "Percentage")) +
  theme_bw() +
  labs(title = "The Copy Number Difference between Major Allele and Other Alleles", x = "Copy Number Difference") +
  theme(legend.position = "none")
ggsave(p_dif_alt, file = "major_alt_ncopy_dif_barplot.pdf", width = 6, height = 4)

# Comparison between major allele and reference allele
major_sub <- major_allele_df[major_allele_df$group == "major_allele", ]
major_sub$N <- fp_df$N[match(major_sub$locus, fp_df$locus)]
major_sub$reference_length <- config[match(major_sub$locus, config[, 6]), 3] - config[match(major_sub$locus, config[, 6]), 2] + 1
major_sub$period_size <- config[match(major_sub$locus, config[, 6]), 4]
major_sub$reference_ncopy <- as.numeric(paste(major_sub$reference_length %/% major_sub$period_size,
                                              major_sub$reference_length %% major_sub$period_size, sep = "."))
major_sub$ncopy_dif <- major_sub$Ncopy - major_sub$reference_ncopy
major_sub <- major_sub[abs(major_sub$ncopy_dif) <= 10 & major_sub$ncopy_dif %% 1 == 0, ]

dif_ref_summary <- data.frame()
for (dif in unique(major_sub$ncopy_dif)) {
  cnt <- nrow(major_sub[major_sub$ncopy_dif == dif, ])
  dif_ref_summary <- rbind(dif_ref_summary, data.frame(ncopy_dif = dif, counts = cnt))
}
dif_ref_summary$percentage <- dif_ref_summary$counts / sum(dif_ref_summary$counts)
dif_ref_summary$group <- ifelse(dif_ref_summary$ncopy_dif == 0, "identical",
                                ifelse(dif_ref_summary$ncopy_dif > 0, "positive", "negative"))
dif_ref_summary$ncopy_dif <- factor(dif_ref_summary$ncopy_dif, levels = seq(-10, 10, 1),
                                    labels = as.character(seq(-10, 10, 1)))

p_dif_ref <- ggplot(dif_ref_summary, aes(x = ncopy_dif, fill = group)) +
  geom_bar(aes(y = counts), stat = "identity", color = "black", alpha = 0.3) +
  geom_line(aes(y = percentage * max(counts), group = 1), linewidth = 0.8) +
  geom_point(aes(y = percentage * max(counts), group = 1)) +
  scale_fill_manual(values = c("gray", "gray", "gray")) +
  scale_y_continuous(name = "Counts", sec.axis = sec_axis(trans = ~ . / max(dif_ref_summary$counts), name = "Percentage")) +
  theme_bw() +
  labs(title = "The Copy Number Difference Between Major Allele and Reference Allele", x = "Copy Number Difference") +
  theme(legend.position = "none")
ggsave(p_dif_ref, file = "major_reference_ncopy_dif_barplot.pdf", width = 6, height = 4)

# Correlation between PIC and copy number (major and reference)
major_sub$PIC <- fp_df$PIC[match(major_sub$locus, fp_df$locus)]

p_cor_major <- ggplot(major_sub, aes(x = Ncopy, y = PIC)) +
  geom_smooth(method = "lm", formula = y ~ x, color = '#2D6DB1', fill = "#756bb1") +
  stat_cor() +
  theme_bw() +
  labs(title = "The Correlation between PIC and Copy Number of Major Allele", x = "Copy Number of Major Allele", y = "PIC") +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 10))
ggsave(p_cor_major, file = "major_ncopy_PIC_corplot.pdf", width = 7, height = 3.5)

p_cor_ref <- ggplot(major_sub, aes(x = reference_ncopy, y = PIC)) +
  geom_smooth(method = "lm", formula = y ~ x, color = '#2D6DB1', fill = "#756bb1") +
  stat_cor() +
  theme_bw() +
  labs(title = "The Correlation between PIC and Copy Number of Reference Allele", x = "Copy Number of Reference Allele", y = "PIC") +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 10))
ggsave(p_cor_ref, file = "reference_ncopy_PIC_corplot.pdf", width = 7, height = 3.5)

# End of script