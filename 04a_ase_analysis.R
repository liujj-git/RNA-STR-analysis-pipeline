# >>> RENAMED from ase_patch.R <<<
# Part of Section 2.4 (allele-specific expression, allele-level analysis).
# Standalone ASE analysis on the relaxed-allele-bias VCF (RNA keeps allele
# imbalance = the ASE signal; DNA uses the fully-filtered VCF as a presence mask).
# RUN ORDER: run after 04_concordance.R; writes the ASE tables
#   (ase_dropout_by_copydiff.txt, ase_dropout_by_depth.txt, ase_balance.txt)
#   that 04c_fig4_ASE.R then reads to draw Figure 4 / Figure S13.
# So the Section 2.4 sequence is:
#   04_concordance.R  ->  04a_ase_analysis.R  ->  04b_fig3_concordance.R
#                                             ->  04c_fig4_ASE.R
# ---------------------------------------------------------------------------
# =============================================================================
# ase_patch.R  — standalone allele-level (ASE) analysis on the relaxed-min-ab VCF
# -----------------------------------------------------------------------------
# Produces the NEW numbers and the NEW main figure (Fig 4) for the allele-level
# section (§2.4), using pstr.ase.recode.vcf.gz (min-ab NOT applied, 24,013-locus
# subset). Design:
#   * RNA side : min-ab relaxed  -> allele bias retained = the ASE signal
#   * DNA side : full filtering  -> high-accuracy germline reference
#
# DNA full-filtering is reproduced EXACTLY and index-safely: we read the original
# fully-filtered VCF (pstr.recode.vcf.gz) only as a BOOLEAN PRESENCE MASK and drop
# any DNA call that the catalogue removed. All GT VALUES, DP and PDP/AB/DAB come
# from the single relaxed VCF, so RNA and DNA share ONE allele-index space — no
# cross-file index mismatch even if vcftools --recode pruned ALT differently.
# (min-ab marks calls no-call, it does not change retained genotypes, so a call
#  present in both files has identical alleles.)
#
# The data logic (allele_sequence, concordance_compared_df, dropout, imbalance)
# is reused VERBATIM from concordance.R. Only THREE things are new/changed and are
# flagged inline below:
#   [NEW 1] input VCF = the relaxed file
#   [NEW 2] DNA presence mask from the full file
#   [NEW 3] the RNA DP>=3 filter is applied NUMERICALLY (see note) so the
#           dropout/imbalance-by-depth bins are populated and correct.
#
# Outputs: revision_values_ase.txt (numbers for the §2.4 placeholders),
#          Fig4_composite.pdf, and the ase_*.txt data tables.
# Run as its OWN session; paths are relative to the concordance working dir.
# =============================================================================
rm(list = ls())
library(vcfR); library(readxl); library(tidyr); library(data.table)
library(ggplot2); library(patchwork)

ASE_VCF  <- "../pstr.ase.recode.vcf.gz"   # [NEW 1] relaxed (no min-ab), 24,013 subset
FULL_VCF <- "../pstr.recode.vcf.gz"       # original fully-filtered catalogue VCF (DNA presence mask)

# -----------------------------------------------------------------------------
# NOTE on the DP>=3 filter.
# In concordance.R the DP table column DP$DP is CHARACTER at the point of the
# "DP$DP >= 3" filter (it is only coerced with as.numeric() afterwards). A
# character ">= 3" is a STRING comparison, so it keeps DP {3-9, 30-99, 300-999,
# ...} and silently DROPS {10-29, 100-299, ...}. That is almost certainly
# unintended. For the allele-level analysis the read-depth axis is central
# (dropout/imbalance vs depth), so reproducing the string behaviour would leave
# the 10-19 and 20-29 depth bins EMPTY. This patch therefore coerces DP to
# numeric BEFORE the filter, i.e. a clean "at least 3 reads" QC. See the message
# accompanying this file for what this means for §2.3.
# -----------------------------------------------------------------------------

## ---- inputs (same sheets/paths as concordance.R) ----------------------------
config <- read.table("../GRCh38.hipstr_reference.refine.bed")
sample_info <- as.data.frame(read_xlsx("../VBsampleinfo.xlsx", sheet = "concordance"))
colnames(sample_info) <- sample_info[1, ]; sample_info <- sample_info[-1, ]

## ---- read the relaxed VCF: GT / DP / PDP / AB / DAB -------------------------
raw_vcf <- read.vcfR(ASE_VCF)
GT  <- extract.gt(raw_vcf, element = "GT")
DPm <- extract.gt(raw_vcf, element = "DP")
.ase_have <- all(c("PDP","AB","DAB") %in% strsplit(raw_vcf@gt[1, "FORMAT"], ":")[[1]])
if (!.ase_have) stop("Relaxed VCF lacks PDP/AB/DAB; point ASE_VCF at the HipSTR-format VCF.")
PDP_mat <- extract.gt(raw_vcf, element = "PDP")
AB_mat  <- extract.gt(raw_vcf, element = "AB")
DAB_mat <- extract.gt(raw_vcf, element = "DAB")

rna_cols <- intersect(as.character(sample_info$RNA), colnames(GT))
dna_cols <- intersect(as.character(sample_info$DNA), colnames(GT))
stopifnot(setequal(rna_cols, unique(as.character(sample_info$RNA))),
          setequal(dna_cols, unique(as.character(sample_info$DNA))))

## ---- [NEW 2] DNA full-filter: keep only DNA calls present in the full VCF ----
## (full VCF used purely as a boolean presence mask; no index/allele transfer)
full_vcf <- read.vcfR(FULL_VCF)
GT_full  <- extract.gt(full_vcf, element = "GT")
ridx <- match(rownames(GT), rownames(GT_full))      # relaxed loci -> rows in full VCF
for (cc in dna_cols) {
  present_full <- !is.na(GT_full[ridx, cc])         # was this DNA call kept by the catalogue?
  GT[!present_full, cc] <- NA                        # drop DNA calls the full-filter VCF removed
}
rm(full_vcf, GT_full)
## RNA columns are left untouched (relaxed) -> exactly the ASE signal we want.

## ---- DP long table (concordance.R lines 44-46) ------------------------------
DPm <- as.data.frame(DPm)
DP  <- data.frame("locus" = rownames(DPm), gather(DPm, key = "sample", value = "DP"))
DP  <- na.omit(DP)
DP$DP <- as.numeric(DP$DP)                            # [NEW 3] numeric -> clean DP>=3 (see note above)

## ---- allele_sequence  (VERBATIM concordance.R lines 51-95, on the ASE VCF) ---
allele_ladder <- data.frame("locus" = raw_vcf@fix[, 3],
                            "Ref"   = raw_vcf@fix[, 4],
                            "Alt"   = raw_vcf@fix[, 5])
allele_ladder <- na.omit(gather(allele_ladder, key = "allele_group", value = "sequence", -locus))
temp1 <- allele_ladder[!grepl(",", allele_ladder$sequence), ]
temp2 <- allele_ladder[grepl(",", allele_ladder$sequence), ]
temp3 <- do.call(rbind, apply(temp2, 1, function(x) {
  alt_allele <- unlist(strsplit(x[3], ","))
  data.frame("locus" = x[1], "allele_group" = 1:length(alt_allele),
             "sequence" = alt_allele, row.names = NULL)
}))
temp3 <- temp3[nchar(temp3$sequence) != 0, ]
temp1 <- rbind(temp1, temp3)
temp1$allele_group[temp1$allele_group == "Ref"] <- 0
temp1$allele_group[temp1$allele_group == "Alt"] <- 1
temp1$period_size <- as.numeric(config[match(temp1$locus, config[, 6]), 4])
temp2 <- nchar(temp1[, 3]) %/% temp1[, 4]
temp3 <- nchar(temp1[, 3]) %%  temp1[, 4]
temp3[temp3 != 0] <- paste(".", temp3[temp3 != 0], sep = "")
temp3[temp3 == 0] <- ""
temp1$ncopy <- paste(temp2, temp3, sep = "")
temp1$ncopy_label <- paste(temp1$locus, temp1$ncopy, sep = "@")
temp2 <- as.data.frame(table(temp1$ncopy_label)); temp2[, 1] <- as.character(temp2[, 1])
temp2 <- temp2[temp2$Freq != 1, 1]
setDT(temp1)
temp1[!temp1$ncopy_label %in% temp2, ncopy_label := ncopy]
temp1[temp1$ncopy_label %in% temp2, ncopy_label := { paste(ncopy, seq_along(ncopy), sep = "_") }, by = ncopy_label]
temp1$unique_label <- paste(temp1$locus, temp1$allele_group, sep = "@")
allele_sequence <- as.data.frame(temp1)[, c("locus","period_size","allele_group","sequence","ncopy","ncopy_label")]
allele_sequence <- allele_sequence[order(allele_sequence$locus), ]
allele_sequence$ncopy <- as.numeric(allele_sequence$ncopy)
rm(raw_vcf)

## ---- concordance_compared_df  (VERBATIM concordance.R lines 513-562) ---------
## RNA GT come from the relaxed VCF; DNA GT are the full-filter calls (masked above)
concordance_compared_df <- data.frame()
for (i in 1:nrow(sample_info)) {
  temp1 <- as.character(sample_info$RNA[i])
  temp2 <- as.character(sample_info$DNA[i])
  temp3 <- as.data.frame(na.omit(GT[rownames(GT) %in% DP$locus[DP$sample == temp1 & DP$DP >= 3],
                                    match(c(temp1, temp2), colnames(GT))]))
  temp3 <- gather(data.frame("locus" = rownames(temp3), temp3), key = "sample", value = "GT", -"locus")
  temp3 <- separate(temp3, col = "GT", into = c("v1","v2"), sep = "\\|")
  temp3[as.numeric(temp3$v1) > as.numeric(temp3$v2), c("v1","v2")] <-
    temp3[as.numeric(temp3$v1) > as.numeric(temp3$v2), c("v2","v1")]
  temp3 <- unite(temp3, c("v1","v2"), col = "GT", sep = ",")
  temp3 <- spread(temp3, key = "sample", value = "GT")
  rownames(temp3) <- temp3$locus; temp3 <- temp3[, -1]
  temp4 <- temp3[temp3[, 1] == temp3[, 2], ]
  df1 <- data.frame("individual" = temp1, "RNA" = temp1, "DNA" = temp2,
                    "locus" = rownames(temp4),
                    "GT_RNA" = temp4[, colnames(temp4) == temp1],
                    "GT_DNA" = temp4[, colnames(temp4) == temp2], "concordance" = 2)
  concordance_compared_df <- rbind(concordance_compared_df, df1)
  temp3 <- temp3[temp3[, 1] != temp3[, 2], ]
  if (nrow(temp3) != 0) {
    df2 <- data.frame()
    for (p in 1:nrow(temp3)) {
      GT1 <- sort(as.numeric(unlist(strsplit(temp3[p, colnames(temp3) == temp1], ","))))
      GT2 <- sort(as.numeric(unlist(strsplit(temp3[p, colnames(temp3) == temp2], ","))))
      temp4 <- ifelse(length(intersect(GT1, GT2)) == 1, 1, 0)
      df2 <- rbind(df2, data.frame("individual" = temp1, "RNA" = temp1, "DNA" = temp2,
                                   "locus" = rownames(temp3)[p],
                                   "GT_RNA" = paste(GT1, collapse = ","),
                                   "GT_DNA" = paste(GT2, collapse = ","), "concordance" = temp4))
    }
    concordance_compared_df <- rbind(concordance_compared_df, df2)
  }
}
concordance_compared_df$DP_RNA <- as.numeric(DP$DP[match(paste(concordance_compared_df$RNA, concordance_compared_df$locus, sep = "@"),
                                                         paste(DP$sample, DP$locus, sep = "@"))])
concordance_compared_df$DP_DNA <- as.numeric(DP$DP[match(paste(concordance_compared_df$DNA, concordance_compared_df$locus, sep = "@"),
                                                         paste(DP$sample, DP$locus, sep = "@"))])

## ===========================================================================
## ASE setup + dropout + imbalance  (VERBATIM concordance.R lines 1131-1226,
## data only; the original ase_theme PDFs are replaced by Fig 4 below)
## ===========================================================================
AB_SIG <- log10(0.05)   # ~ -1.30 ; AB below this = significant allele bias (p < 0.05)
ccd <- concordance_compared_df
al  <- allele_sequence
ncopy_lu <- setNames(as.numeric(al$ncopy), paste(al$locus, al$allele_group, sep = "@"))
ref_lu   <- setNames(as.numeric(al$ncopy[al$allele_group == 0]), al$locus[al$allele_group == 0])
ase_idx     <- function(gt, k) vapply(strsplit(gt, ","), function(x) x[k], character(1))
ase_gf      <- function(mat, locus, sample) mat[cbind(match(locus, rownames(mat)), match(sample, colnames(mat)))]
ase_dpbin   <- function(x) cut(x, c(-Inf, 10, 20, 30, Inf), c("<10","10-19","20-29",">=30"), right = FALSE)
ase_copybin <- function(x) cut(round(abs(x)), c(-Inf, 0.5, 1.5, 2.5, 3.5, 4.5, Inf), c("0","1","2","3","4",">=5"))
ccd$rna_i1 <- ase_idx(ccd$GT_RNA, 1); ccd$rna_i2 <- ase_idx(ccd$GT_RNA, 2)
ccd$dna_i1 <- ase_idx(ccd$GT_DNA, 1); ccd$dna_i2 <- ase_idx(ccd$GT_DNA, 2)
ccd$rna_c1 <- ncopy_lu[paste(ccd$locus, ccd$rna_i1, sep = "@")]
ccd$rna_c2 <- ncopy_lu[paste(ccd$locus, ccd$rna_i2, sep = "@")]
ccd$dna_c1 <- ncopy_lu[paste(ccd$locus, ccd$dna_i1, sep = "@")]
ccd$dna_c2 <- ncopy_lu[paste(ccd$locus, ccd$dna_i2, sep = "@")]
ccd$ref_c  <- ref_lu[ccd$locus]
ccd <- ccd[complete.cases(ccd[, c("rna_c1","rna_c2","dna_c1","dna_c2","ref_c")]), ]

## (2a) allelic dropout: DNA length-hets; dropout = RNA called homozygous -------
lh <- ccd[ccd$dna_i1 != ccd$dna_i2 & ccd$dna_c1 != ccd$dna_c2, ]
lh$rna_hom    <- lh$rna_i1 == lh$rna_i2
lh$dCopy_bin  <- ase_copybin(abs(lh$dna_c1 - lh$dna_c2))
lh$RNA_DP_bin <- ase_dpbin(lh$DP_RNA)
dropout_overall <- mean(lh$rna_hom)
drop_by_copy <- aggregate(rna_hom ~ dCopy_bin, lh, function(x) c(n = length(x), rate = mean(x)))
drop_by_copy <- data.frame(dCopy_bin = drop_by_copy$dCopy_bin,
                           n = as.integer(drop_by_copy$rna_hom[, "n"]), dropout_rate = drop_by_copy$rna_hom[, "rate"])
drop_by_dp <- aggregate(rna_hom ~ RNA_DP_bin, lh, function(x) c(n = length(x), rate = mean(x)))
drop_by_dp <- data.frame(RNA_DP_bin = drop_by_dp$RNA_DP_bin,
                         n = as.integer(drop_by_dp$rna_hom[, "n"]), dropout_rate = drop_by_dp$rna_hom[, "rate"])
write.table(drop_by_copy, "ase_dropout_by_copydiff.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(drop_by_dp,   "ase_dropout_by_depth.txt",    col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

## (2b) allelic balance: RNA het length-het calls; balance from PDP, AB ---------
bal <- ccd[ccd$rna_i1 != ccd$rna_i2 & ccd$rna_c1 != ccd$rna_c2, ]
bal$PDP <- ase_gf(PDP_mat, bal$locus, bal$RNA)
bal$AB  <- suppressWarnings(as.numeric(ase_gf(AB_mat,  bal$locus, bal$RNA)))
bal$DAB <- suppressWarnings(as.numeric(ase_gf(DAB_mat, bal$locus, bal$RNA)))
ase_minor <- function(s) {
  v <- suppressWarnings(as.numeric(strsplit(s, "\\|")[[1]]))
  if (length(v) != 2 || any(is.na(v)) || sum(v) == 0) return(NA_real_)
  min(v) / sum(v)
}
bal$minor_frac <- vapply(bal$PDP, ase_minor, numeric(1))
bal <- bal[!is.na(bal$minor_frac), ]
bal$dCopy_bin  <- ase_copybin(abs(bal$rna_c1 - bal$rna_c2))
bal$RNA_DP_bin <- ase_dpbin(bal$DP_RNA)
bal$AB_sig <- !is.na(bal$AB) & bal$AB < AB_SIG
write.table(bal[, c("individual","RNA","locus","GT_RNA","DP_RNA","PDP","minor_frac","AB","DAB","dCopy_bin","RNA_DP_bin","AB_sig")],
            "ase_balance.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
bal_summary <- rbind(
  data.frame(stratum = "all", level = "all", n = nrow(bal),
             median_minor_frac = median(bal$minor_frac), pct_AB_significant = 100 * mean(bal$AB_sig)),
  do.call(rbind, lapply(levels(bal$dCopy_bin), function(L) { s <- bal[bal$dCopy_bin == L, ]; if (!nrow(s)) return(NULL)
    data.frame(stratum = "copy_diff", level = L, n = nrow(s), median_minor_frac = median(s$minor_frac), pct_AB_significant = 100 * mean(s$AB_sig)) })),
  do.call(rbind, lapply(levels(bal$RNA_DP_bin), function(L) { s <- bal[bal$RNA_DP_bin == L, ]; if (!nrow(s)) return(NULL)
    data.frame(stratum = "RNA_DP", level = L, n = nrow(s), median_minor_frac = median(s$minor_frac), pct_AB_significant = 100 * mean(s$AB_sig)) })))
write.table(bal_summary, "ase_balance_summary.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

## ===========================================================================
## Numbers for the §2.4 manuscript placeholders
## ===========================================================================
fmtp <- function(x) sprintf("%.2f%%", 100 * x)
rv <- c(
  sprintf("# allele-level (ASE) values on the relaxed-min-ab VCF: %s", ASE_VCF),
  sprintf("# generated %s", as.character(Sys.time())),
  sprintf("RNA-DNA matched genotype pairs (n)\t%d", nrow(ccd)),
  "",
  "## [DROPOUT]  scope = DNA length-heterozygotes;  dropout = RNA called homozygous",
  sprintf("DNA length-het genotypes (n)\t%d", nrow(lh)),
  sprintf("overall dropout rate\t%s", fmtp(dropout_overall)),
  "  by copy-number difference between DNA alleles:")
for (i in seq_len(nrow(drop_by_copy)))
  rv <- c(rv, sprintf("    dCopy=%s\tn=%d\tdropout=%s",
                      as.character(drop_by_copy$dCopy_bin[i]), drop_by_copy$n[i], fmtp(drop_by_copy$dropout_rate[i])))
rv <- c(rv, "  by RNA read depth:")
for (i in seq_len(nrow(drop_by_dp)))
  rv <- c(rv, sprintf("    RNA_DP=%s\tn=%d\tdropout=%s",
                      as.character(drop_by_dp$RNA_DP_bin[i]), drop_by_dp$n[i], fmtp(drop_by_dp$dropout_rate[i])))
rv <- c(rv, "",
        "## [IMBALANCE]  scope = RNA heterozygous length-het calls",
        sprintf("RNA het length-het calls (n)\t%d", nrow(bal)),
        sprintf("median minor-allele read fraction\t%.3f", median(bal$minor_frac)),
        sprintf("significant allele bias overall (AB < %.2f, p<0.05)\t%s", AB_SIG, fmtp(mean(bal$AB_sig))),
        "  by RNA read depth:")
for (L in levels(bal$RNA_DP_bin)) { s <- bal[bal$RNA_DP_bin == L, ]; if (nrow(s))
  rv <- c(rv, sprintf("    RNA_DP=%s\tn=%d\tsignificant=%s", L, nrow(s), fmtp(mean(s$AB_sig)))) }
writeLines(rv, "revision_values_ase.txt")
cat("\n", paste(rv, collapse = "\n"), "\n\n", sep = "")
cat("revision_values_ase.txt written.\n")

## ===========================================================================
## Figure 4  (allele-level; theme_jgg, cairo_pdf)
##   A dropout vs copy-difference   B dropout vs RNA depth
##   C minor-allele read fraction   D % significant allele bias vs RNA depth
## ===========================================================================
theme_jgg <- theme_classic(base_family = "Arial", base_size = 9) +
  theme(axis.title = element_text(colour = "black", size = 9),
        axis.text  = element_text(colour = "black", size = 8),
        axis.line  = element_line(colour = "black", linewidth = 0.4),
        axis.ticks = element_line(colour = "black", linewidth = 0.4),
        plot.tag = element_text(face = "bold", size = 11))

# A: dropout vs copy-number difference  ("0" bin relabelled "<1")
dc <- drop_by_copy; dc$dCopy_bin <- as.character(dc$dCopy_bin)
dc$dCopy_bin[dc$dCopy_bin == "0"] <- "<1"
dc$dCopy_bin <- factor(dc$dCopy_bin, levels = c("<1","1","2","3","4",">=5"))
pDrop <- ggplot(dc, aes(dCopy_bin, 100 * dropout_rate)) +
  geom_col(fill = "#D7263D", colour = "black", width = 0.72, linewidth = 0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_jgg + labs(x = "Copy-number difference between DNA alleles", y = "RNA homozygous-call rate (%)")

# B: dropout vs RNA read depth
dd <- drop_by_dp; dd$RNA_DP_bin <- factor(as.character(dd$RNA_DP_bin), levels = c("<10","10-19","20-29",">=30"))
pDropDepth <- ggplot(dd, aes(RNA_DP_bin, 100 * dropout_rate)) +
  geom_col(fill = "#D7263D", colour = "black", width = 0.6, linewidth = 0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_jgg + labs(x = "RNA read depth", y = "RNA homozygous-call rate (%)")

# C: minor-allele read fraction distribution
pBalDist <- ggplot(bal, aes(minor_frac)) +
  geom_histogram(binwidth = 0.02, fill = "#2D6DB1", colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_jgg + labs(x = "Minor-allele read fraction", y = "Heterozygous calls")

# D: % significant allele bias vs RNA read depth
sigdp <- bal_summary[bal_summary$stratum == "RNA_DP", ]
sigdp$level <- factor(sigdp$level, levels = c("<10","10-19","20-29",">=30"))
pBalSig <- ggplot(sigdp, aes(level, pct_AB_significant)) +
  geom_col(fill = "#E69F00", colour = "black", width = 0.6, linewidth = 0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  theme_jgg + labs(x = "RNA read depth", y = "Calls with significant allele bias (%)")

fig4 <- (pDrop | pDropDepth) / (pBalDist | pBalSig) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11))
ggsave("Fig4_composite.pdf", fig4, width = 9, height = 8, device = cairo_pdf)
cat("Fig4_composite.pdf written (A dropout/copy, B dropout/depth, C imbalance dist, D significant bias/depth).\n")
