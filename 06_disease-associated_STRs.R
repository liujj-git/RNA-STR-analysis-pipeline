# ============================================================================
# Script: Analysis of disease-associated STRs
# Description:
#   1. Load configuration and VCFs for disease STRs (ExpansionHunter output).
#   2. Compute genotyping rates in RNA and DNA samples.
#   3. Evaluate concordance between technical RNA replicates and RNA-DNA pairs.
#   4. Identify potentially pathogenic alleles exceeding disease thresholds.
#   5. Generate summary plots for genotyping rates and concordance.
# ============================================================================

rm(list = ls())
library(rjson)
library(jsonlite)
library(vcfR)
library(tidyr)
library(readxl)
library(ggplot2)
library(ggsci)
library(dplyr)

# ----------------------------- Input data -----------------------------------
# HipSTR reference config (for period sizes etc., though not heavily used)
config_hipstr <- read.table("../GRCh38.hipstr_reference.refine.bed")

# Disease STR catalog from ExpansionHunter (JSON)
config_disease <- fromJSON("eh.hg38.variant_catalog.disease.json")

# Disease locus positions
disease_position <- read.table("disease.locusposition.txt")
colnames(disease_position) <- c("locus", "subtype", "position")
disease_position <- separate(disease_position, col = "position", into = c("chr", "position"), sep = ":")
disease_position <- separate(disease_position, col = "position", into = c("start", "end"), sep = "-")

# Locus intersection (optional, from original)
locus_intersection <- read.table("locus_intersection.txt")

# Disease threshold (normal max, pathogenic min)
disease_threshold <- read.table("../../disease_threshold.txt")
temp1 <- disease_threshold[disease_threshold$V1 %in% config_disease$LocusId, ]
colnames(temp1) <- c("locus", "NormalMax", "PathogenicMin")
write.table(temp1, file = "disease_threshold.config.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# Genotype VCFs from ExpansionHunter (REPCN field)
disease_rna <- read.vcfR("disease.rna.vcf")
disease_rna_gt <- as.data.frame(extract.gt(disease_rna, element = "REPCN"))

disease_dna <- read.vcfR("disease.dna.vcf")
disease_dna_gt <- as.data.frame(extract.gt(disease_dna, element = "REPCN"))
colnames(disease_dna_gt) <- gsub(".deduped", "", colnames(disease_dna_gt))
colnames(disease_dna_gt)[grepl("VB", colnames(disease_dna_gt))] <- paste(colnames(disease_dna_gt)[grepl("VB", colnames(disease_dna_gt))], "_dna", sep = "")

# Reshape to long format
disease_rna_gt <- data.frame(locus = rownames(disease_rna_gt), disease_rna_gt)
disease_rna_gt <- gather(disease_rna_gt, key = "sample", value = "genotype", -locus)

disease_dna_gt <- data.frame(locus = rownames(disease_dna_gt), disease_dna_gt)
disease_dna_gt <- gather(disease_dna_gt, key = "sample", value = "genotype", -locus)

disease_gt <- rbind(data.frame(group = "dna", sample = disease_dna_gt$sample, locus = disease_dna_gt$locus, genotype = disease_dna_gt$genotype),
                    data.frame(group = "rna", sample = disease_rna_gt$sample, locus = disease_rna_gt$locus, genotype = disease_rna_gt$genotype))
disease_gt <- na.omit(disease_gt)

# Map locus coordinates to subtype names
disease_gt$locus <- disease_position$subtype[match(disease_gt$locus, paste(disease_position$chr, disease_position$start, sep = "_"))]

# Clean genotype: replace "/" with "," and sort alleles
disease_gt$genotype <- gsub("\\/", ",", disease_gt$genotype)
disease_gt <- separate(disease_gt, col = "genotype", into = c("A1", "A2"), sep = ",")
idx <- as.numeric(disease_gt$A1) > as.numeric(disease_gt$A2)
disease_gt[idx, c("A1", "A2")] <- disease_gt[idx, c("A2", "A1")]
disease_gt <- unite(disease_gt, "genotype", c("A1", "A2"), sep = ",")

# Sample information
sample_info <- as.data.frame(read_xlsx("../VBsampleinfo.xlsx", sheet = "concordance"))
colnames(sample_info) <- sample_info[1, ]
sample_info <- sample_info[-1, ]

repeated_samples <- as.data.frame(read_xlsx("../VBsampleinfo.xlsx", sheet = "repeatability"))
colnames(repeated_samples) <- repeated_samples[1, ]
repeated_samples <- repeated_samples[-1, ]

unrelated_individuals <- read_excel("../VBsampleinfo.xlsx", sheet = "unrelated_individuals")
colnames(unrelated_individuals) <- unrelated_individuals[1, ]
unrelated_individuals <- unrelated_individuals[-1, ]

# ----------------------------- Genotyping rate (RNA unrelated samples) -------
genotyping_rate_df <- data.frame()
for (sub in disease_position$subtype) {
  sub_gt <- disease_gt[disease_gt$sample %in% unrelated_individuals$ID & disease_gt$locus == sub, ]
  cnt_all <- nrow(sub_gt)
  cnt_filtered <- nrow(sub_gt[sub_gt$genotype != "0,0", ])
  genotyping_rate_df <- rbind(genotyping_rate_df,
                              data.frame(filter_0_0 = FALSE, subtype = sub, genotype_counts = cnt_all),
                              data.frame(filter_0_0 = TRUE,  subtype = sub, genotype_counts = cnt_filtered))
}
genotyping_rate_df$locus <- disease_position$locus[match(genotyping_rate_df$subtype, disease_position$subtype)]

# Aggregate by locus
rate_summary <- data.frame()
for (l in unique(genotyping_rate_df$locus)) {
  sub_all <- genotyping_rate_df[genotyping_rate_df$locus == l & genotyping_rate_df$filter_0_0 == FALSE, ]
  sub_filt <- genotyping_rate_df[genotyping_rate_df$locus == l & genotyping_rate_df$filter_0_0 == TRUE, ]
  rate_all <- (sum(sub_all$genotype_counts) / nrow(sub_all)) / nrow(unrelated_individuals)
  rate_filt <- (sum(sub_filt$genotype_counts) / nrow(sub_filt)) / nrow(unrelated_individuals)
  rate_summary <- rbind(rate_summary,
                        data.frame(group = "all", locus = l, genotyping_rate = rate_all),
                        data.frame(group = "filtered", locus = l, genotyping_rate = rate_filt))
}
# Keep loci with filtered rate >= 0.1
keep_loci <- rate_summary$locus[rate_summary$group == "filtered" & rate_summary$genotyping_rate >= 0.1]
rate_summary <- rate_summary[rate_summary$locus %in% keep_loci, ]

# Order by all-rate
order_loci <- rate_summary[rate_summary$group == "all", ]
order_loci <- order_loci[order(order_loci$genotyping_rate, decreasing = TRUE), ]
rate_summary$locus <- factor(rate_summary$locus, levels = order_loci$locus)
rate_summary$label <- paste0(round(rate_summary$genotyping_rate * 100, 2), "%")

p_rate_rna <- ggplot(rate_summary, aes(x = locus, y = genotyping_rate, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  geom_text(aes(label = label), position = position_dodge(width = 0.8), size = 3, vjust = -0.2) +
  scale_fill_rickandmorty() +
  theme_bw() +
  labs(x = "Locus", y = "Genotyping Rate") +
  theme(axis.text.x = element_text(face = "bold", size = 8),
        axis.title = element_text(face = "bold", size = 10))
ggsave(p_rate_rna, file = "genotyping_rate_barplot.pdf", height = 5, width = 15)

# ----------------------------- Genotyping rate (DNA compared samples) -------
genotyping_rate_df2 <- data.frame()
for (sub in disease_position$subtype) {
  sub_gt <- disease_gt[disease_gt$sample %in% unique(sample_info$DNA) & disease_gt$locus == sub, ]
  cnt_all <- nrow(sub_gt)
  cnt_filtered <- nrow(sub_gt[sub_gt$genotype != "0,0", ])
  genotyping_rate_df2 <- rbind(genotyping_rate_df2,
                               data.frame(filter_0_0 = FALSE, subtype = sub, genotype_counts = cnt_all),
                               data.frame(filter_0_0 = TRUE,  subtype = sub, genotype_counts = cnt_filtered))
}
genotyping_rate_df2$locus <- disease_position$locus[match(genotyping_rate_df2$subtype, disease_position$subtype)]

rate_summary2 <- data.frame()
for (l in unique(genotyping_rate_df2$locus)) {
  sub_all <- genotyping_rate_df2[genotyping_rate_df2$locus == l & genotyping_rate_df2$filter_0_0 == FALSE, ]
  sub_filt <- genotyping_rate_df2[genotyping_rate_df2$locus == l & genotyping_rate_df2$filter_0_0 == TRUE, ]
  rate_all <- (sum(sub_all$genotype_counts) / nrow(sub_all)) / length(unique(sample_info$DNA))
  rate_filt <- (sum(sub_filt$genotype_counts) / nrow(sub_filt)) / length(unique(sample_info$DNA))
  rate_summary2 <- rbind(rate_summary2,
                         data.frame(group = "all", locus = l, genotyping_rate = rate_all),
                         data.frame(group = "filtered", locus = l, genotyping_rate = rate_filt))
}
keep_loci2 <- rate_summary2$locus[rate_summary2$group == "filtered" & rate_summary2$genotyping_rate >= 0.1]
rate_summary2 <- rate_summary2[rate_summary2$locus %in% keep_loci2, ]
order_loci2 <- rate_summary2[rate_summary2$group == "all", ]
order_loci2 <- order_loci2[order(order_loci2$genotyping_rate, decreasing = TRUE), ]
rate_summary2$locus <- factor(rate_summary2$locus, levels = order_loci2$locus)
rate_summary2$label <- paste0(round(rate_summary2$genotyping_rate * 100, 2), "%")

p_rate_dna <- ggplot(rate_summary2, aes(x = locus, y = genotyping_rate, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  geom_text(aes(label = label), position = position_dodge(width = 0.8), size = 2, vjust = -0.2) +
  scale_fill_manual(values = c("#0073C2FF", "#EFC000FF")) +
  theme_bw() +
  labs(title = "Genotyping rate of disease-associated STRs in blood genome", x = "Locus", y = "Genotyping Rate") +
  theme(axis.text.x = element_text(face = "bold", size = 8, angle = 90),
        axis.title = element_text(face = "bold", size = 10))
ggsave(p_rate_dna, file = "genotyping_rate_barplot2.pdf", height = 5, width = 15)

# ----------------------------- Concordance in RNA replicates -----------------
# Build repeated sample concordance table
repeated_df <- data.frame()
for (l in unique(rate_summary$locus)) {
  subtypes <- disease_position$subtype[disease_position$locus == l]
  sub_gt_list <- list()
  for (st in subtypes) {
    sub_gt <- disease_gt[disease_gt$locus == st, ]
    gt1 <- sub_gt$genotype[match(repeated_samples$R1, sub_gt$sample)]
    gt2 <- sub_gt$genotype[match(repeated_samples$R2, sub_gt$sample)]
    sub_gt_list <- rbind(sub_gt_list, data.frame(GT1 = gt1, GT2 = gt2))
  }
  temp <- data.frame(locus = l, sample1 = repeated_samples$R1, sample2 = repeated_samples$R2,
                     GT1 = sub_gt_list$GT1, GT2 = sub_gt_list$GT2)
  repeated_df <- rbind(repeated_df, temp)
}
repeated_df <- na.omit(repeated_df)

# Classify concordance
repeated_conc <- data.frame()
for (l in unique(repeated_df$locus)) {
  sub <- repeated_df[repeated_df$locus == l, ]
  identical_all <- sub[sub$GT1 == sub$GT2, ]
  identical_filt <- identical_all[identical_all$GT1 != "0,0" & identical_all$GT2 != "0,0", ]
  cnt_ident_all <- nrow(identical_all)
  cnt_ident_filt <- nrow(identical_filt)
  
  non_identical <- sub[sub$GT1 != sub$GT2, ]
  if (nrow(non_identical) > 0) {
    inter_allele <- sapply(1:nrow(non_identical), function(i) {
      a1 <- unlist(strsplit(non_identical$GT1[i], ","))
      a2 <- unlist(strsplit(non_identical$GT2[i], ","))
      length(intersect(a1, a2))
    })
    non_identical$inter <- inter_allele
    partial_all <- non_identical[non_identical$inter == 1, ]
    partial_filt <- partial_all[partial_all$GT1 != "0,0" & partial_all$GT2 != "0,0", ]
    non_all <- non_identical[non_identical$inter == 0, ]
    non_filt <- non_all[non_all$GT1 != "0,0" & non_all$GT2 != "0,0", ]
    cnt_partial_all <- nrow(partial_all)
    cnt_partial_filt <- nrow(partial_filt)
    cnt_non_all <- nrow(non_all)
    cnt_non_filt <- nrow(non_filt)
  } else {
    cnt_partial_all <- cnt_partial_filt <- cnt_non_all <- cnt_non_filt <- 0
  }
  repeated_conc <- rbind(repeated_conc,
                         data.frame(locus = l, group = "identical", filter_0_0 = FALSE, counts = cnt_ident_all),
                         data.frame(locus = l, group = "identical", filter_0_0 = TRUE,  counts = cnt_ident_filt),
                         data.frame(locus = l, group = "partial_identical", filter_0_0 = FALSE, counts = cnt_partial_all),
                         data.frame(locus = l, group = "partial_identical", filter_0_0 = TRUE,  counts = cnt_partial_filt),
                         data.frame(locus = l, group = "non_identical", filter_0_0 = FALSE, counts = cnt_non_all),
                         data.frame(locus = l, group = "non_identical", filter_0_0 = TRUE,  counts = cnt_non_filt))
}

# Plot RNA-RNA concordance (all genotypes)
temp_rep_all <- repeated_conc[repeated_conc$filter_0_0 == FALSE, ]
for (l in unique(temp_rep_all$locus)) {
  temp_rep_all$total[temp_rep_all$locus == l] <- sum(temp_rep_all$counts[temp_rep_all$locus == l])
}
temp_rep_all <- temp_rep_all[temp_rep_all$total >= 10, ]  # at least 10 shared pairs
temp_rep_all$group <- factor(temp_rep_all$group, levels = c("non_identical", "partial_identical", "identical"))
# Order by identical percentage
ident_pct <- temp_rep_all[temp_rep_all$group == "identical", ]
ident_pct <- ident_pct[order(ident_pct$counts / ident_pct$total, decreasing = TRUE), ]
temp_rep_all$locus <- factor(temp_rep_all$locus, levels = ident_pct$locus)

p_rep_all <- ggplot(temp_rep_all, aes(x = locus, y = counts, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  geom_text(aes(label = counts), position = position_dodge(width = 0.8), size = 3, vjust = -0.2) +
  scale_fill_startrek() +
  theme_bw() +
  labs(x = "Locus", y = "Shared Genotype Counts") +
  theme(axis.text.x = element_text(face = "bold", size = 8),
        axis.title = element_text(face = "bold", size = 10))
ggsave(p_rep_all, file = "repeated_all_concordance_barplot.pdf", width = 8, height = 5)

# Filtered (exclude "0,0") version
temp_rep_filt <- repeated_conc[repeated_conc$filter_0_0 == TRUE, ]
for (l in unique(temp_rep_filt$locus)) {
  temp_rep_filt$total[temp_rep_filt$locus == l] <- sum(temp_rep_filt$counts[temp_rep_filt$locus == l])
}
temp_rep_filt <- temp_rep_filt[temp_rep_filt$total >= 10, ]
temp_rep_filt$group <- factor(temp_rep_filt$group, levels = c("non_identical", "partial_identical", "identical"))
ident_pct_filt <- temp_rep_filt[temp_rep_filt$group == "identical", ]
ident_pct_filt <- ident_pct_filt[order(ident_pct_filt$counts / ident_pct_filt$total, decreasing = TRUE), ]
temp_rep_filt$locus <- factor(temp_rep_filt$locus, levels = ident_pct_filt$locus)

p_rep_filt <- ggplot(temp_rep_filt, aes(x = locus, y = counts, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  geom_text(aes(label = counts), position = position_dodge(width = 0.8), size = 3, vjust = -0.2) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  labs(title = "Genotype (filtered) concordance of disease-associated STRs between repeated RNA samples",
       x = "Locus", y = "Shared Genotype Counts (Filtered)") +
  theme(axis.text.x = element_text(face = "bold", size = 8),
        axis.title = element_text(face = "bold", size = 10))
ggsave(p_rep_filt, file = "repeated_filtered_concordance_barplot.pdf", width = 10, height = 5)

# ----------------------------- Concordance in RNA-DNA pairs -----------------
compared_df <- data.frame()
for (l in unique(rate_summary$locus)) {
  subtypes <- disease_position$subtype[disease_position$locus == l]
  sub_gt_list <- list()
  for (st in subtypes) {
    sub_gt <- disease_gt[disease_gt$locus == st, ]
    gt_rna <- sub_gt$genotype[match(sample_info$RNA, sub_gt$sample)]
    gt_dna <- sub_gt$genotype[match(sample_info$DNA, sub_gt$sample)]
    sub_gt_list <- rbind(sub_gt_list, data.frame(GT1 = gt_rna, GT2 = gt_dna))
  }
  temp <- data.frame(locus = l, sample1 = sample_info$RNA, sample2 = sample_info$DNA,
                     GT1 = sub_gt_list$GT1, GT2 = sub_gt_list$GT2)
  compared_df <- rbind(compared_df, temp)
}
compared_df <- na.omit(compared_df)

# Classify
compared_conc <- data.frame()
for (l in unique(compared_df$locus)) {
  sub <- compared_df[compared_df$locus == l, ]
  identical_all <- sub[sub$GT1 == sub$GT2, ]
  identical_filt <- identical_all[identical_all$GT1 != "0,0" & identical_all$GT2 != "0,0", ]
  cnt_ident_all <- nrow(identical_all)
  cnt_ident_filt <- nrow(identical_filt)
  
  non_identical <- sub[sub$GT1 != sub$GT2, ]
  if (nrow(non_identical) > 0) {
    inter_allele <- sapply(1:nrow(non_identical), function(i) {
      a1 <- unlist(strsplit(non_identical$GT1[i], ","))
      a2 <- unlist(strsplit(non_identical$GT2[i], ","))
      length(intersect(a1, a2))
    })
    non_identical$inter <- inter_allele
    partial_all <- non_identical[non_identical$inter == 1, ]
    partial_filt <- partial_all[partial_all$GT1 != "0,0" & partial_all$GT2 != "0,0", ]
    non_all <- non_identical[non_identical$inter == 0, ]
    non_filt <- non_all[non_all$GT1 != "0,0" & non_all$GT2 != "0,0", ]
    cnt_partial_all <- nrow(partial_all)
    cnt_partial_filt <- nrow(partial_filt)
    cnt_non_all <- nrow(non_all)
    cnt_non_filt <- nrow(non_filt)
  } else {
    cnt_partial_all <- cnt_partial_filt <- cnt_non_all <- cnt_non_filt <- 0
  }
  compared_conc <- rbind(compared_conc,
                         data.frame(locus = l, group = "identical", filter_0_0 = FALSE, counts = cnt_ident_all),
                         data.frame(locus = l, group = "identical", filter_0_0 = TRUE,  counts = cnt_ident_filt),
                         data.frame(locus = l, group = "partial_identical", filter_0_0 = FALSE, counts = cnt_partial_all),
                         data.frame(locus = l, group = "partial_identical", filter_0_0 = TRUE,  counts = cnt_partial_filt),
                         data.frame(locus = l, group = "non_identical", filter_0_0 = FALSE, counts = cnt_non_all),
                         data.frame(locus = l, group = "non_identical", filter_0_0 = TRUE,  counts = cnt_non_filt))
}

# Plot RNA-DNA concordance (all)
temp_comp_all <- compared_conc[compared_conc$filter_0_0 == FALSE, ]
for (l in unique(temp_comp_all$locus)) {
  temp_comp_all$total[temp_comp_all$locus == l] <- sum(temp_comp_all$counts[temp_comp_all$locus == l])
}
temp_comp_all <- temp_comp_all[temp_comp_all$total >= 10, ]
temp_comp_all$group <- factor(temp_comp_all$group, levels = c("non_identical", "partial_identical", "identical"))
ident_pct_comp <- temp_comp_all[temp_comp_all$group == "identical", ]
ident_pct_comp <- ident_pct_comp[order(ident_pct_comp$counts / ident_pct_comp$total, decreasing = TRUE), ]
temp_comp_all$locus <- factor(temp_comp_all$locus, levels = ident_pct_comp$locus)

p_comp_all <- ggplot(temp_comp_all, aes(x = locus, y = counts, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  geom_text(aes(label = counts), position = position_dodge(width = 0.8), size = 3, vjust = -0.2) +
  scale_fill_startrek() +
  theme_bw() +
  labs(x = "Locus", y = "Shared Genotype Counts (All)") +
  theme(axis.text.x = element_text(face = "bold", size = 8),
        axis.title = element_text(face = "bold", size = 10))
ggsave(p_comp_all, file = "compared_all_concordance_barplot.pdf", width = 8, height = 5)

# Filtered version
temp_comp_filt <- compared_conc[compared_conc$filter_0_0 == TRUE, ]
for (l in unique(temp_comp_filt$locus)) {
  temp_comp_filt$total[temp_comp_filt$locus == l] <- sum(temp_comp_filt$counts[temp_comp_filt$locus == l])
}
temp_comp_filt <- temp_comp_filt[temp_comp_filt$total >= 10, ]
temp_comp_filt$group <- factor(temp_comp_filt$group, levels = c("non_identical", "partial_identical", "identical"))
ident_pct_comp_filt <- temp_comp_filt[temp_comp_filt$group == "identical", ]
ident_pct_comp_filt <- ident_pct_comp_filt[order(ident_pct_comp_filt$counts / ident_pct_comp_filt$total, decreasing = TRUE), ]
temp_comp_filt$locus <- factor(temp_comp_filt$locus, levels = ident_pct_comp_filt$locus)

p_comp_filt <- ggplot(temp_comp_filt, aes(x = locus, y = counts, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  geom_text(aes(label = counts), position = position_dodge(width = 0.8), size = 3, vjust = -0.2) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  labs(title = "Genotype (filtered) concordance of disease-associated STRs between RNA-DNA sample pairs",
       x = "Locus", y = "Shared Genotype Counts (Filtered)") +
  theme(axis.text.x = element_text(face = "bold", size = 8),
        axis.title = element_text(face = "bold", size = 10))
ggsave(p_comp_filt, file = "compared_filtered_concordance_barplot.pdf", width = 10, height = 5)

# ----------------------------- Pathogenic allele detection -------------------
disease_allele_df <- data.frame()
for (l in unique(disease_position$locus)) {
  subtypes <- disease_position$subtype[disease_position$locus == l]
  sub_gt <- disease_gt[disease_gt$sample %in% unique(sample_info$DNA) &
                         disease_gt$locus %in% subtypes, ]
  if (nrow(sub_gt) == 0) next
  alleles <- unlist(strsplit(sub_gt$genotype, ","))
  allele_counts <- as.data.frame(table(alleles))
  allele_counts$Var1 <- as.numeric(as.character(allele_counts$Var1))
  allele_counts$freq <- allele_counts$Freq / sum(allele_counts$Freq)
  allele_counts <- allele_counts[order(allele_counts$freq, decreasing = TRUE), ]
  threshold_val <- disease_threshold[disease_threshold$V1 == l, 3]  # PathogenicMin
  pathogenic <- allele_counts[allele_counts$Var1 >= threshold_val, ]
  if (nrow(pathogenic) > 0) {
    disease_allele_df <- rbind(disease_allele_df,
                               data.frame(locus = l,
                                          allele = pathogenic$Var1,
                                          disease_allele_counts = pathogenic$Freq,
                                          total_allele_counts = sum(allele_counts$Freq),
                                          allele_percentage = pathogenic$freq,
                                          threshold = threshold_val))
  }
}
write.table(disease_allele_df, file = "disease_allele_summary.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# ----------------------------- Supplementary combined rate plot ------------
# Combine RNA and DNA genotyping rates for the same set of loci
rate_combined <- merge(rate_summary[rate_summary$group == "all", c("locus", "genotyping_rate")],
                       rate_summary2[rate_summary2$group == "all", c("locus", "genotyping_rate")],
                       by = "locus", all = TRUE)
colnames(rate_combined) <- c("locus", "transcriptome", "genome")
rate_combined[is.na(rate_combined)] <- 0
rate_combined_long <- gather(rate_combined, key = "group", value = "genotyping_rate", -locus)
order_long <- rate_combined_long[rate_combined_long$group == "transcriptome", ]
order_long <- order_long[order(order_long$genotyping_rate, decreasing = TRUE), ]
rate_combined_long$locus <- factor(rate_combined_long$locus, levels = order_long$locus)
rate_combined_long$label <- paste0(round(rate_combined_long$genotyping_rate * 100, 2), "%")

p_rate_comb <- ggplot(rate_combined_long, aes(x = locus, y = genotyping_rate, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  geom_text(aes(label = label), position = position_dodge(width = 0.8), size = 2, vjust = -0.2) +
  scale_fill_rickandmorty() +
  theme_bw() +
  labs(x = "Locus", y = "Genotyping Rate") +
  theme(axis.text.x = element_text(face = "bold", size = 8, angle = 90),
        axis.title = element_text(face = "bold", size = 10))
ggsave(p_rate_comb, file = "genotyping_rate_barplot3.pdf", height = 5, width = 15)

# End of script