# =============================================================================
# eSTR_01_expression_prep.R   (revised: uses fread, exports TSS, manifest)
#
# Purpose: Build a normalized gene expression matrix for the eSTR analysis,
#          merging the two separately-quantified cohorts:
#            * VB cohort  : 45 BAMs  / 40 individuals (5 with technical reps)
#            * R  cohort  : 147 BAMs / 127 individuals (20 with technical reps)
#          into a single 167-individual matrix, with covariates.
#
# Inputs  (paths relative to 8_eSTR/):
#   featurecounts_raw.45VB.txt    -- VB cohort featureCounts gene-level matrix
#   featurecounts_raw.147VB.txt   -- R  cohort featureCounts gene-level matrix
#   ../VBsampleinfo.xlsx          -- master sample sheet (4 sheets)
#
# Outputs (8_eSTR/):
#   expression_for_eSTR.tsv              -- gene x 167 ID, INT-transformed
#   eSTR_covariates.tsv                  -- ID, sex, cohort, PC1..PCk
#   eSTR_sample_info.tsv                 -- per-ID QC after rep folding
#   gene_info_eSTR.tsv                   -- Geneid, chr, start, end, strand,
#                                           length, TSS (kept genes only)
#   eSTR_prep_manifest.tsv               -- input paths + MD5 + run params
#   revision_values_eSTR_prep.txt        -- QC summary (for writing list)
#   expression_pca_screeplot.pdf
#   expression_pca_pc1pc2_by_cohort.pdf
#
# Key decisions:
#   * Replicate folding: raw counts SUMMED across an individual's BAMs (treats
#     reps as deeper sequencing). Merged column name follows
#     unrelated_individuals$ID (e.g., RF002_1, VBF02A_1, RF004), so columns
#     align with the existing STR-genotype VCF used in expression.R/concordance.R.
#   * Gene filter: TPM > 0.1 in >=50% of samples in BOTH cohorts.
#   * Normalization: edgeR TMM -> log2(CPM+1) -> per-gene inverse-normal
#     transform (GTEx convention).
#   * Hidden factors: PCA on the INT matrix; PCs with >1% variance retained
#     (capped at 20). The eSTR model in Phase 2 will fit:
#         expr_INT ~ STR_dosage + sex + PC1..PCk
#   * TSS for cis-window in Phase 2: TSS = start if strand=='+' else end
#     (gene-level boundary; sufficient for +/-100 kb mapping).
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(edgeR)
  library(tools)
})

# ---- 0. paths ---------------------------------------------------------------
fc_vb_path <- "featurecounts_raw.45VB.txt"
fc_r_path  <- "featurecounts_raw.147VB.txt"
info_path  <- "../VBsampleinfo.xlsx"

out_expr     <- "expression_for_eSTR.tsv"
out_cov      <- "eSTR_covariates.tsv"
out_info     <- "eSTR_sample_info.tsv"
out_gene     <- "gene_info_eSTR.tsv"
out_manifest <- "eSTR_prep_manifest.tsv"
out_rv       <- "revision_values_eSTR_prep.txt"
out_scree    <- "expression_pca_screeplot.pdf"
out_pcplot   <- "expression_pca_pc1pc2_by_cohort.pdf"

stopifnot(file.exists(fc_vb_path), file.exists(fc_r_path), file.exists(info_path))

t0 <- Sys.time()
rv <- character()
add <- function(...) { rv[length(rv)+1L] <<- sprintf(...) }

# ---- helpers ----------------------------------------------------------------
# Read a featureCounts gene-level table. Skips the leading "# Program ..."
# comment line via skip="Geneid". Column headers in the raw file are full BAM
# paths; we shorten to the directory containing each BAM (canonical sample ID).
read_fc <- function(path) {
  d <- fread(path, sep = "\t", header = TRUE, skip = "Geneid",
             check.names = FALSE, data.table = FALSE)
  meta_cols <- c("Geneid","Chr","Start","End","Strand","Length")
  stopifnot(all(meta_cols %in% colnames(d)))
  sample_cols <- setdiff(colnames(d), meta_cols)
  short_ids <- vapply(sample_cols, function(p) {
    parts <- strsplit(p, "/", fixed = TRUE)[[1]]
    if (length(parts) >= 2L) parts[length(parts) - 1L] else parts[length(parts)]
  }, character(1))
  list(meta = d[, meta_cols, drop = FALSE],
       counts = as.matrix(d[, sample_cols, drop = FALSE]),
       sample_ids = unname(short_ids))
}

strip_rep <- function(x) sub("_[12]$", "", x)

INT <- function(x) {
  ok <- !is.na(x)
  r  <- rank(x, ties.method = "average", na.last = "keep")
  out <- rep(NA_real_, length(x))
  out[ok] <- qnorm((r[ok] - 0.5) / sum(ok))
  out
}

# ---- 1. read both featureCounts & align gene sets ---------------------------
vb <- read_fc(fc_vb_path); colnames(vb$counts) <- vb$sample_ids
r  <- read_fc(fc_r_path);  colnames(r$counts)  <- r$sample_ids
cat(sprintf("VB cohort: %d genes x %d BAMs\n", nrow(vb$counts), ncol(vb$counts)))
cat(sprintf("R  cohort: %d genes x %d BAMs\n", nrow(r$counts),  ncol(r$counts)))

# tolerate row-order differences but reject set differences
if (!setequal(vb$meta$Geneid, r$meta$Geneid)) {
  stop("Geneid sets differ between the two featureCounts files. ",
       "Re-quantify both cohorts with the same GTF before merging.")
}
if (!identical(vb$meta$Geneid, r$meta$Geneid)) {
  ord <- match(vb$meta$Geneid, r$meta$Geneid)
  r$counts <- r$counts[ord, , drop = FALSE]
  r$meta   <- r$meta[ord, , drop = FALSE]
  message("Re-ordered R-cohort rows to match VB-cohort Geneid order.")
}
if (!identical(vb$meta$Length, r$meta$Length)) {
  warning("Gene Lengths differ between the two featureCounts files. ",
          "Using VB-side lengths for TPM; inspect if differences are large.")
}
gene_meta   <- vb$meta
n_genes_raw <- nrow(gene_meta)
add("[GENES] featureCounts raw: %d", n_genes_raw)

vb_lib_raw <- colSums(vb$counts); r_lib_raw <- colSums(r$counts)
add("[LIB-RAW] mean assigned reads/sample (M): VB %.1f , R %.1f",
    mean(vb_lib_raw)/1e6, mean(r_lib_raw)/1e6)

# ---- 2. sample sheet --------------------------------------------------------
info_all <- as.data.frame(read_excel(info_path, sheet = "total", skip = 1),
                          stringsAsFactors = FALSE)
info_all <- info_all[info_all$seq == "RNA", ]
info_all$ID         <- as.character(info_all$ID)
info_all$individual <- as.character(info_all$individual)
info_all$gender     <- as.character(info_all$gender)

unrel <- as.data.frame(read_excel(info_path, sheet = "unrelated_individuals", skip = 1),
                       stringsAsFactors = FALSE)
stopifnot(nrow(unrel) == 167L)
unrel$ID         <- as.character(unrel$ID)
unrel$individual <- as.character(unrel$individual)
gender_lu        <- setNames(info_all$gender, info_all$ID)
unrel$gender     <- gender_lu[unrel$ID]
unrel$cohort     <- ifelse(grepl("^VB", unrel$individual), "VB", "R")
if (any(is.na(unrel$gender))) stop("Missing gender for some atlas IDs")

add("[ATLAS] unrelated individuals: %d (R=%d, VB=%d)",
    nrow(unrel), sum(unrel$cohort == "R"), sum(unrel$cohort == "VB"))
add("        gender: female=%d, male=%d",
    sum(unrel$gender == "female"), sum(unrel$gender == "male"))
add("        R  cohort: female=%d, male=%d",
    sum(unrel$cohort == "R"  & unrel$gender == "female"),
    sum(unrel$cohort == "R"  & unrel$gender == "male"))
add("        VB cohort: female=%d, male=%d",
    sum(unrel$cohort == "VB" & unrel$gender == "female"),
    sum(unrel$cohort == "VB" & unrel$gender == "male"))

# ---- 3. fold technical replicates -> one column per atlas ID ----------------
bam_indiv_vb <- strip_rep(colnames(vb$counts))
bam_indiv_r  <- strip_rep(colnames(r$counts))

fold_one <- function(target_id, cohort_label) {
  indiv <- strip_rep(target_id)
  if (cohort_label == "VB") { mat <- vb$counts; idx <- bam_indiv_vb }
  else                      { mat <- r$counts;  idx <- bam_indiv_r  }
  cols <- colnames(mat)[idx == indiv]
  if (length(cols) == 0L) return(list(v = NULL,                          n_bam = 0L))
  if (length(cols) == 1L) return(list(v = as.numeric(mat[, cols]),       n_bam = 1L))
  list(v = as.numeric(rowSums(mat[, cols, drop = FALSE])), n_bam = length(cols))
}

merged    <- matrix(NA_real_, nrow = n_genes_raw, ncol = nrow(unrel),
                    dimnames = list(gene_meta$Geneid, unrel$ID))
n_bam_used <- integer(nrow(unrel))
for (i in seq_len(nrow(unrel))) {
  res <- fold_one(unrel$ID[i], unrel$cohort[i])
  if (is.null(res$v))
    stop(sprintf("No BAM found for atlas ID '%s' (cohort %s, individual %s)",
                 unrel$ID[i], unrel$cohort[i], unrel$individual[i]))
  merged[, i]   <- res$v
  n_bam_used[i] <- res$n_bam
}
add("[FOLD] individuals with >1 BAM (raw counts summed): %d", sum(n_bam_used > 1L))
add("        R  cohort: %d , VB cohort: %d",
    sum(n_bam_used > 1L & unrel$cohort == "R"),
    sum(n_bam_used > 1L & unrel$cohort == "VB"))

# ---- 4. library size after folding -----------------------------------------
lib_size <- colSums(merged)
is_VB    <- unrel$cohort == "VB"
is_R     <- !is_VB
add("[LIB] assigned-count library size after folding (M):")
add("       R  median %.1f  range %.1f - %.1f",
    median(lib_size[is_R])/1e6,  min(lib_size[is_R])/1e6,  max(lib_size[is_R])/1e6)
add("       VB median %.1f  range %.1f - %.1f",
    median(lib_size[is_VB])/1e6, min(lib_size[is_VB])/1e6, max(lib_size[is_VB])/1e6)
if (median(lib_size[is_VB]) > 1.5 * median(lib_size[is_R]) ||
    median(lib_size[is_R])  > 1.5 * median(lib_size[is_VB])) {
  warning("Median library size differs >1.5x between cohorts; TMM + INT will ",
          "handle this, but expect cohort to dominate early PCs.")
}

# ---- 5. TPM -----------------------------------------------------------------
gene_len <- gene_meta$Length
rpk      <- merged / gene_len
tpm      <- t( t(rpk) / (colSums(rpk) / 1e6) )

# ---- 6. gene filter ---------------------------------------------------------
expr_vb_ok <- rowMeans(tpm[, is_VB] > 0.1) >= 0.5
expr_r_ok  <- rowMeans(tpm[, is_R ] > 0.1) >= 0.5
keep       <- expr_vb_ok & expr_r_ok
add("[GENE FILTER] TPM>0.1 in >=50%% of samples in both cohorts: %d / %d",
    sum(keep), n_genes_raw)
add("              expressed in VB only: %d ; R only: %d ; both: %d ; neither: %d",
    sum(expr_vb_ok & !expr_r_ok), sum(expr_r_ok & !expr_vb_ok),
    sum(keep), sum(!expr_vb_ok & !expr_r_ok))

counts_f <- merged[keep, , drop = FALSE]
meta_f   <- gene_meta[keep, , drop = FALSE]

# ---- 7. TMM -> log2(CPM+1) -> per-gene INT ----------------------------------
dge      <- DGEList(counts = counts_f)
dge      <- calcNormFactors(dge, method = "TMM")
logcpm   <- cpm(dge, log = TRUE, prior.count = 1)
expr_int <- t(apply(logcpm, 1, INT))
dimnames(expr_int) <- dimnames(logcpm)
add("[NORM] TMM + log2(CPM+1) + per-gene inverse-normal transform")

# ---- 8. PCA -> covariates ---------------------------------------------------
expr_c <- expr_int - rowMeans(expr_int)
pca    <- prcomp(t(expr_c), center = FALSE, scale. = FALSE)
ve     <- pca$sdev^2 / sum(pca$sdev^2)

n_pc_use <- min(sum(ve > 0.01), 20L)
add("[PCA] PCs with >1%% variance: %d (using top %d)", sum(ve > 0.01), n_pc_use)
add("       PC1 %.1f%%  PC2 %.1f%%  PC3 %.1f%%  PC4 %.1f%%  PC5 %.1f%%",
    ve[1]*100, ve[2]*100, ve[3]*100, ve[4]*100, ve[5]*100)

# how much do early PCs separate the two cohorts? both Wilcoxon P and R^2
cohort01 <- as.integer(is_VB)
for (k in 1:min(10, n_pc_use)) {
  p  <- wilcox.test(pca$x[is_VB, k], pca$x[is_R, k])$p.value
  r2 <- summary(lm(pca$x[, k] ~ cohort01))$r.squared
  add("[BATCH] PC%d: cohort Wilcoxon P=%.2e ; R^2(cohort) = %.3f", k, p, r2)
}

# diagnostic plots
pdf(out_scree, width = 6, height = 4)
plot(seq_len(min(30, length(ve))), ve[seq_len(min(30, length(ve)))] * 100,
     type = "b", pch = 19, xlab = "PC", ylab = "Variance explained (%)",
     main = sprintf("Expression PCA scree (top 30 of %d PCs)", length(ve)))
abline(h = 1, col = "red", lty = 2)
dev.off()

pdf(out_pcplot, width = 6, height = 5.5)
col_cohort <- ifelse(is_VB, "#E76F51", "#264653")
plot(pca$x[, 1], pca$x[, 2], col = col_cohort, pch = 19,
     xlab = sprintf("PC1 (%.1f%%)", ve[1] * 100),
     ylab = sprintf("PC2 (%.1f%%)", ve[2] * 100),
     main = "Expression PC1 vs PC2 (after INT)")
legend("topright", c(sprintf("VB cohort (n=%d)", sum(is_VB)),
                     sprintf("R cohort (n=%d)",  sum(is_R))),
       col = c("#E76F51", "#264653"), pch = 19, bty = "n")
dev.off()

# ---- 9. write outputs -------------------------------------------------------
PCmat <- pca$x[, seq_len(n_pc_use), drop = FALSE]
colnames(PCmat) <- paste0("PC", seq_len(n_pc_use))

cov_df <- data.frame(ID = unrel$ID, individual = unrel$individual,
                     cohort = unrel$cohort, sex = unrel$gender,
                     stringsAsFactors = FALSE)
cov_df <- cbind(cov_df, PCmat)

info_df <- data.frame(ID = unrel$ID, individual = unrel$individual,
                      cohort = unrel$cohort, sex = unrel$gender,
                      n_BAMs_folded = n_bam_used,
                      lib_size_assigned = lib_size,
                      stringsAsFactors = FALSE)

fwrite(data.frame(Geneid = rownames(expr_int), expr_int,
                  check.names = FALSE, stringsAsFactors = FALSE),
       file = out_expr, sep = "\t", quote = FALSE)
fwrite(cov_df,  file = out_cov,  sep = "\t", quote = FALSE)
fwrite(info_df, file = out_info, sep = "\t", quote = FALSE)

# gene info with TSS for Phase 2 cis-window mapping
# featureCounts may join multi-region values with ";" -- collapse safely
split_first <- function(z) vapply(strsplit(z, ";", fixed = TRUE),
                                  function(x) x[1], character(1))
split_last  <- function(z) vapply(strsplit(z, ";", fixed = TRUE),
                                  function(x) x[length(x)], character(1))

gene_chr    <- split_first(as.character(meta_f$Chr))
gene_strand <- split_first(as.character(meta_f$Strand))
gene_start  <- as.integer(split_first(as.character(meta_f$Start)))
gene_end    <- as.integer(split_last (as.character(meta_f$End)))
gene_tss    <- ifelse(gene_strand == "+", gene_start, gene_end)

gene_df <- data.frame(Geneid = meta_f$Geneid,
                      chr    = gene_chr,
                      start  = gene_start,
                      end    = gene_end,
                      strand = gene_strand,
                      length = meta_f$Length,
                      TSS    = gene_tss,
                      stringsAsFactors = FALSE)
fwrite(gene_df, file = out_gene, sep = "\t", quote = FALSE)

# manifest: input MD5s + dims + run metadata (reproducibility)
mani <- data.frame(
  item  = c("script", "VB_featurecounts", "R_featurecounts", "sample_info",
            "n_genes_raw", "n_genes_kept", "n_individuals", "n_BAMs_VB",
            "n_BAMs_R", "n_reps_folded_VB", "n_reps_folded_R",
            "n_PCs_exported", "run_time_min"),
  value = c("eSTR_01_expression_prep.R",
            fc_vb_path, fc_r_path, info_path,
            as.character(n_genes_raw), as.character(sum(keep)),
            as.character(nrow(unrel)), as.character(ncol(vb$counts)),
            as.character(ncol(r$counts)),
            as.character(sum(n_bam_used > 1L & unrel$cohort == "VB")),
            as.character(sum(n_bam_used > 1L & unrel$cohort == "R")),
            as.character(n_pc_use),
            sprintf("%.1f", as.numeric(difftime(Sys.time(), t0, units = "mins")))),
  md5   = c(NA_character_,
            tryCatch(unname(md5sum(fc_vb_path)), error = function(e) NA_character_),
            tryCatch(unname(md5sum(fc_r_path)),  error = function(e) NA_character_),
            tryCatch(unname(md5sum(info_path)),  error = function(e) NA_character_),
            rep(NA_character_, 9)),
  stringsAsFactors = FALSE
)
fwrite(mani, file = out_manifest, sep = "\t", quote = FALSE)

writeLines(c(
  sprintf("# eSTR expression preparation revision values (auto-generated %s)",
          format(Sys.time(), "%Y-%m-%d %H:%M")),
  sprintf("# inputs: %s + %s", fc_vb_path, fc_r_path),
  rv
), out_rv)

cat("\nDone. Outputs in 8_eSTR/:\n")
cat("  ", out_expr,     "  -- gene x 167 ID, INT log2(CPM+1)\n")
cat("  ", out_cov,      "  -- covariates for eSTR (sex, cohort, PCs)\n")
cat("  ", out_info,     "  -- per-ID QC after fold\n")
cat("  ", out_gene,     "  -- gene info with TSS (input to Phase 2)\n")
cat("  ", out_manifest, "  -- input/output manifest (reproducibility)\n")
cat("  ", out_rv,       "  -- revision values (for writing list)\n")
cat("  ", out_scree,    "  -- PCA scree\n")
cat("  ", out_pcplot,   "  -- PC1 vs PC2 by cohort\n")
