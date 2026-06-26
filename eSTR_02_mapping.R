# =============================================================================
# eSTR_02_mapping.R   (intronic pSTR main analysis, n up to 167)
#
# Purpose: For each intronic pSTR with an annotated host gene in the expression
#          set, fit
#
#              expr_INT[gene, sample] ~ dosage[locus, sample]
#                                     + sex + PC1 + ... + PC9
#
#          dosage = STR length in bp (sum of two allele lengths). Report BH
#          FDR-adjusted p; emit catalog, motif enrichment, top examples.
#
#          Why this design avoids circularity (R1's main concern):
#            * intronic STR genotype comes from intronic / pre-mRNA reads;
#            * gene expression (featureCounts -t exon) comes from exonic reads;
#            * the two read pools are essentially disjoint, so a positive
#              association is not a self-confirmation.
#
# Inputs:
#   expression_for_eSTR.tsv         (from eSTR_01)
#   eSTR_covariates.tsv             (from eSTR_01; includes PC1..PCk)
#   gene_info_eSTR.tsv              (from eSTR_01; for reference)
#   ../pstr.rna.recode.vcf.gz       (HipSTR pSTR VCF, 167 samples)
#   ../2_annotation/STR_annotation_group.txt   (from annotation.R)
#
# Outputs (8_eSTR/):
#   eSTRs_catalog.tsv               -- one row per tested STR-gene pair
#   eSTR_motif_enrichment.tsv       -- if motif column present in annotation
#   eSTR_examples/                  -- top-3 example scatter PDFs
#   revision_values_eSTR_mapping.txt
#
# Modelling notes:
#   * Batch (the two data-collection rounds) is absorbed by the expression PCs:
#     PC1 is ~97% explained by that axis (R^2=0.975), so retaining PC1..PC9 and
#     omitting an explicit cohort term controls the same structure without a
#     near-collinear pair. (Sensitivity: adding a cohort term, with or without
#     PC1, leaves the eSTRs essentially unchanged; see eSTR_03.)
#   * dosage is computed from the GB FORMAT field (per-allele bp diff from
#     reference), not parsed from REF/ALT sequence; this gives length-based
#     dosage even for iso-alleles (same length, different sequence).
# =============================================================================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# biomaRt is only a *fallback* for ID mapping now (see section 2b); the primary
# path uses the offline org.Hs.eg.db, installed lazily there if/when needed.
if (!requireNamespace("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt", update = FALSE, ask = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(vcfR)
  library(ggplot2)
})

t0 <- Sys.time()

# ---- 0. paths ---------------------------------------------------------------
expr_path      <- "expression_for_eSTR.tsv"
cov_path       <- "eSTR_covariates.tsv"
gene_info_path <- "gene_info_eSTR.tsv"
vcf_path       <- "../pstr.rna.recode.vcf.gz"
ann_path       <- "../2_annotation/STR_annotation_group.txt"

out_catalog      <- "eSTRs_catalog.tsv"
out_motif_enrich <- "eSTR_motif_enrichment.tsv"
out_rv           <- "revision_values_eSTR_mapping.txt"
out_examples_dir <- "eSTR_examples"
dir.create(out_examples_dir, showWarnings = FALSE)

stopifnot(file.exists(expr_path), file.exists(cov_path),
          file.exists(gene_info_path), file.exists(vcf_path),
          file.exists(ann_path))

rv <- character()
add <- function(...) { rv[length(rv)+1L] <<- sprintf(...) }
strip_ver <- function(x) sub("\\.\\d+$", "", x)

# ---- 1. expression + covariates --------------------------------------------
expr_dt  <- fread(expr_path)
genes    <- expr_dt$Geneid
expr_mat <- as.matrix(expr_dt[, -1])
rownames(expr_mat) <- genes

cov <- fread(cov_path)
stopifnot(all(cov$ID %in% colnames(expr_mat)))
expr_mat <- expr_mat[, cov$ID, drop = FALSE]
add("[INPUT] expression: %d genes x %d samples", nrow(expr_mat), ncol(expr_mat))

# strip ENSG version suffixes (cheap; harmless if absent)
rownames(expr_mat) <- strip_ver(rownames(expr_mat))
genes <- rownames(expr_mat)

# ---- 2. STR annotation -----------------------------------------------------
ann <- fread(ann_path)
add("[ANN] columns: %s", paste(colnames(ann), collapse = ", "))

# locate the locus-ID column
locus_col <- NULL
for (cand in c("locus", "ID", "locus_id", "STR_id")) {
  if (cand %in% colnames(ann)) { locus_col <- cand; break }
}
if (is.null(locus_col))
  stop("Cannot find locus ID column in annotation (tried locus / ID / locus_id / STR_id). ",
       "Got: ", paste(colnames(ann), collapse = ", "))
add("[ANN] using '%s' as locus key", locus_col)

gene_col   <- if ("gene_id"     %in% colnames(ann)) "gene_id"     else stop("No gene_id column in annotation")
region_col <- if ("final_group" %in% colnames(ann)) "final_group" else stop("No final_group column in annotation")

# ---- 2b. gene ID format detection + NM -> ENSG mapping (if needed) ---------
# featureCounts used GENCODE GTF -> gene IDs are ENSG######.
# annotation.R extracted VEP's Feature column (col 5) = RefSeq NM_ accessions.
# If gene_id contains NM_/NR_ values they won't match expression; map via biomaRt.
sample_ids <- na.omit(ann[[gene_col]])
sample_ids <- sample_ids[sample_ids != "" & sample_ids != "-"]
is_refseq <- mean(grepl("^NM_|^NR_|^XM_|^XR_", sample_ids)) > 0.5
add("[ANN] gene_id format: %s (example: %s)",
    if (is_refseq) "RefSeq NM/NR" else "other (ENSG or symbol)",
    head(sample_ids, 1))

if (is_refseq) {
  cat("gene_id contains RefSeq accessions (NM_/NR_); mapping to ENSG ",
      "(offline via org.Hs.eg.db, biomaRt archive as fallback)...\n", sep = "")
  nm_query <- unique(strip_ver(sample_ids))
  
  # ---- mapper: NM_/NR_ -> ENSG -------------------------------------------------
  # Preferred path is fully offline (org.Hs.eg.db ships the RefSeq<->Ensembl
  # table), so a live-Ensembl outage / expired-cert mirror no longer aborts the
  # run. biomaRt is attempted only if the offline package is unavailable.
  map_refseq_to_ensg <- function(nm) {
    
    # (1) offline AnnotationDbi -- no network at query time --------------------
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
      tryCatch(BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE),
               error = function(e) message("  could not install org.Hs.eg.db: ",
                                           conditionMessage(e)))
    if (requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
        requireNamespace("AnnotationDbi", quietly = TRUE)) {
      suppressPackageStartupMessages({ library(org.Hs.eg.db); library(AnnotationDbi) })
      # org.Hs.eg.db REFSEQ keys are unversioned and cover NM_/NR_ (and more)
      valid <- intersect(nm, keys(org.Hs.eg.db, keytype = "REFSEQ"))
      if (length(valid) > 0L) {
        sel <- suppressMessages(
          AnnotationDbi::select(org.Hs.eg.db, keys = valid,
                                keytype = "REFSEQ", columns = "ENSEMBL"))
        sel <- as.data.table(sel)[!is.na(ENSEMBL) & ENSEMBL != ""]
        if (nrow(sel) > 0L) {
          sel <- sel[, .(ENSEMBL = ENSEMBL[1L]), by = REFSEQ]   # first ENSG per accession
          message("  mapped via org.Hs.eg.db (offline).")
          return(setNames(sel$ENSEMBL, sel$REFSEQ))
        }
      }
      message("  org.Hs.eg.db returned no matches; trying biomaRt fallback...")
    }
    
    # (2) biomaRt fallback -- try archive hosts (separate infra, often up) -----
    suppressPackageStartupMessages(library(biomaRt))
    hosts <- c("https://nov2023.archive.ensembl.org",
               "https://www.ensembl.org",
               "https://useast.ensembl.org")
    mart <- NULL
    for (h in hosts) {
      mart <- tryCatch(
        useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", host = h),
        error = function(e) { message("  biomaRt host failed: ", h); NULL })
      if (!is.null(mart)) { message("  using biomaRt host: ", h); break }
    }
    if (is.null(mart))
      stop("NM->ENSG mapping failed: org.Hs.eg.db unavailable AND all biomaRt ",
           "hosts unreachable.\n  Fix once (needs any working internet, not ",
           "Ensembl):  BiocManager::install('org.Hs.eg.db')")
    # query mRNA (NM_) and ncRNA (NR_) accessions, then combine
    bm1 <- as.data.table(getBM(c("refseq_mrna",  "ensembl_gene_id"),
                               filters = "refseq_mrna",  values = nm, mart = mart))
    bm2 <- as.data.table(getBM(c("refseq_ncrna", "ensembl_gene_id"),
                               filters = "refseq_ncrna", values = nm, mart = mart))
    setnames(bm1, "refseq_mrna",  "refseq"); setnames(bm2, "refseq_ncrna", "refseq")
    bm <- rbindlist(list(bm1, bm2))
    bm <- bm[refseq != "" & ensembl_gene_id != ""]
    bm <- bm[, .(ensembl_gene_id = ensembl_gene_id[1L]), by = refseq]
    setNames(bm$ensembl_gene_id, bm$refseq)
  }
  
  lu_nm2ensg <- map_refseq_to_ensg(nm_query)
  ann[, gene_id_mapped := lu_nm2ensg[strip_ver(get(gene_col))]]
  n_mapped <- sum(!is.na(ann$gene_id_mapped))
  add("[MAP] NM->ENSG (org.Hs.eg.db / biomaRt fallback): %d / %d IDs mapped (%.1f%%)",
      n_mapped, nrow(ann), 100 * n_mapped / nrow(ann))
  cat(sprintf("  Mapped %d / %d RefSeq IDs (%.1f%%)\n", n_mapped, nrow(ann),
              100 * n_mapped / nrow(ann)))
  # override gene_id_strip with ENSG values where available
  ann[, gene_id_strip_ensg := strip_ver(gene_id_mapped)]
  # use ENSG-mapped ID as the effective gene identifier
  gene_col_eff <- "gene_id_strip_ensg"
} else {
  # already ENSG (or symbol); strip_ver applied later per the existing logic
  gene_col_eff <- gene_col
  ann[, gene_id_strip_ensg := strip_ver(get(gene_col))]
}

# Filter to intronic STRs with gene_id in expressed set
# use gene_id_strip_ensg (ENSG IDs after mapping, or stripped ENSG if already ENSG)
ann[, gene_id_strip := gene_id_strip_ensg]   # unify column name for downstream code
ann_intron <- ann[get(region_col) == "intron"]
add("[FILTER] all intronic pSTRs: %d", nrow(ann_intron))
ann_intron <- ann_intron[!is.na(gene_id_strip) & gene_id_strip != ""]
add("[FILTER] intronic with annotated gene_id: %d", nrow(ann_intron))
ann_intron <- ann_intron[gene_id_strip %in% genes]
add("[FILTER] intronic with gene in expression set: %d", nrow(ann_intron))

# ---- 3. load VCF and build dosage matrix from GB ---------------------------
cat("Reading VCF (may take ~1-2 min)...\n")
vcf      <- read.vcfR(vcf_path, verbose = FALSE)
vcf_loci <- getID(vcf)
ref_bp   <- nchar(getREF(vcf))
add("[VCF] loci: %d", length(vcf_loci))
if (any(is.na(vcf_loci)) || any(vcf_loci == ""))
  warning(sum(is.na(vcf_loci) | vcf_loci == ""), " VCF records have missing ID; ",
          "annotation-side matching uses the ID field.")

cat("Extracting GB and computing dosage...\n")
gb_mat  <- extract.gt(vcf, element = "GB")
add("[VCF] samples: %d", ncol(gb_mat))

# Parse "a|b" or "a/b" -> two integers (bp difference from REF); NA otherwise
gb_vec <- as.vector(gb_mat)
a <- suppressWarnings(as.integer(sub("^(-?\\d+)[|/].+$", "\\1", gb_vec)))
b <- suppressWarnings(as.integer(sub("^.+[|/](-?\\d+)$", "\\1", gb_vec)))
diff_sum <- a + b
diff_mat <- matrix(diff_sum, nrow = nrow(gb_mat), ncol = ncol(gb_mat),
                   dimnames = dimnames(gb_mat))
dosage <- diff_mat + 2L * ref_bp
add("[VCF] dosage matrix built (bp scale)")

# ---- 4. sample alignment ----------------------------------------------------
common <- intersect(colnames(dosage), colnames(expr_mat))
if (length(common) < 100) {
  cat("First VCF samples:  ", head(colnames(dosage),  5), "\n")
  cat("First expr samples: ", head(colnames(expr_mat), 5), "\n")
  stop("Fewer than 100 common samples between VCF and expression. ",
       "Check sample naming.")
}
add("[ALIGN] common samples (VCF n expr): %d", length(common))
dosage   <- dosage[, common, drop = FALSE]
expr_mat <- expr_mat[, common, drop = FALSE]
cov_sub  <- as.data.frame(cov[match(common, cov$ID), ])

# ---- 5. locus alignment (annotation -> VCF) --------------------------------
ann_intron[, vcf_idx := match(get(locus_col), vcf_loci)]
n_lost <- sum(is.na(ann_intron$vcf_idx))
add("[ALIGN] intronic loci not found in VCF: %d (dropped)", n_lost)
ann_intron <- ann_intron[!is.na(vcf_idx)]
add("[ALIGN] intronic loci after VCF match: %d", nrow(ann_intron))

# ---- 6. STR-level QC filter -------------------------------------------------
d_sub    <- dosage[ann_intron$vcf_idx, , drop = FALSE]
n_geno   <- rowSums(!is.na(d_sub))
n_levels <- vapply(seq_len(nrow(d_sub)),
                   function(i) length(unique(na.omit(d_sub[i, ]))),
                   integer(1))
keep_locus <- n_geno >= 0.8 * ncol(d_sub) & n_levels >= 3L
add("[FILTER] loci passing QC (>=80%% genotyping, >=3 dosage levels): %d / %d",
    sum(keep_locus), length(keep_locus))
ann_intron <- ann_intron[keep_locus]
d_sub      <- d_sub[keep_locus, , drop = FALSE]
n_geno     <- n_geno[keep_locus]
n_levels   <- n_levels[keep_locus]
n_tests    <- nrow(ann_intron)
add("[TESTS] STR-gene pairs to test: %d", n_tests)

# ---- 7. fit linear models (fast: .lm.fit + chol2inv) -----------------------
# Build base X (without dosage; same for every locus)
pcs_in_cov <- grep("^PC", colnames(cov_sub), value = TRUE)
pcs_use    <- pcs_in_cov                        # keep ALL PCs (PC1..PCk): they absorb latent technical variation (incl. library-prep / collection-time differences); no explicit cohort term
add("[MODEL] expr ~ dosage + sex + %s  (all PCs retained; PC1 absorbs the technical axis, no cohort term)",
    paste(pcs_use, collapse = " + "))
cov_sub$sex    <- factor(cov_sub$sex)
cov_sub$cohort <- factor(cov_sub$cohort)        # kept only for plot colouring below, NOT a model covariate
fmla   <- as.formula(paste("~ sex +", paste(pcs_use, collapse = " + ")))
X_base <- model.matrix(fmla, data = cov_sub)
add("[MODEL] base design matrix: %d samples x %d columns", nrow(X_base), ncol(X_base))

beta_v <- rep(NA_real_, n_tests)
se_v   <- rep(NA_real_, n_tests)
t_v    <- rep(NA_real_, n_tests)
p_v    <- rep(NA_real_, n_tests)
n_v    <- integer(n_tests)

cat("Fitting", n_tests, "linear models...\n")
report_every <- max(1L, n_tests %/% 10L)
for (i in seq_len(n_tests)) {
  d  <- d_sub[i, ]
  ok <- !is.na(d)
  if (sum(ok) < 30L) next
  g  <- ann_intron$gene_id_strip[i]
  y  <- expr_mat[g, ok]
  if (any(is.na(y))) {
    ok2 <- !is.na(y)
    y   <- y[ok2]
    d_ok <- d[ok][ok2]
    Xi  <- cbind(d = d_ok, X_base[ok, , drop = FALSE][ok2, , drop = FALSE])
  } else {
    Xi  <- cbind(d = d[ok], X_base[ok, , drop = FALSE])
  }
  if (length(y) < 30L) next
  fit <- tryCatch(.lm.fit(Xi, y), error = function(e) NULL)
  if (is.null(fit)) next
  rss  <- sum(fit$residuals^2)
  df.r <- length(y) - ncol(Xi)
  if (df.r < 1L) next
  XtX_inv <- tryCatch(chol2inv(chol(crossprod(Xi))), error = function(e) NULL)
  if (is.null(XtX_inv)) next
  se_d   <- sqrt(XtX_inv[1, 1] * rss / df.r)
  beta_d <- fit$coefficients[1]
  t_d    <- beta_d / se_d
  p_d    <- 2 * pt(abs(t_d), df = df.r, lower.tail = FALSE)
  beta_v[i] <- beta_d
  se_v[i]   <- se_d
  t_v[i]    <- t_d
  p_v[i]    <- p_d
  n_v[i]    <- length(y)
  if (i %% report_every == 0L)
    cat(sprintf("  %d / %d (%.0f%%)\n", i, n_tests, 100 * i / n_tests))
}

# ---- 8. results table + BH FDR ---------------------------------------------
results <- data.table(
  locus            = ann_intron[[locus_col]],
  gene_id          = ann_intron[[gene_col]],      # original annotation ID (RefSeq NM_/NR_ or ENSG)
  gene_ensg        = ann_intron$gene_id_strip,    # ENSG key used to index expr_mat (post NM->ENSG mapping)
  n                = n_v,
  n_dosage_levels  = n_levels,
  genotyping_rate  = n_geno / ncol(dosage),
  beta_bp          = beta_v,
  se               = se_v,
  t                = t_v,
  p                = p_v
)
# carry useful annotation columns through if present
for (ec in c("biotype", "motif", "period", "Chr", "Start", "End",
             "intergenic_subtype")) {
  if (ec %in% colnames(ann_intron)) results[, (ec) := ann_intron[[ec]]]
}
results <- results[!is.na(p)]
results[, fdr := p.adjust(p, method = "BH")]
results[, is_eSTR_005 := fdr < 0.05]
results[, is_eSTR_010 := fdr < 0.10]
setorder(results, p)

n_05 <- sum(results$is_eSTR_005)
n_10 <- sum(results$is_eSTR_010)
add("[RESULT] tests completed: %d", nrow(results))
add("[RESULT] eSTRs at FDR<0.05: %d (%.2f%%)", n_05, 100 * n_05 / nrow(results))
add("[RESULT] eSTRs at FDR<0.10: %d (%.2f%%)", n_10, 100 * n_10 / nrow(results))
add("[RESULT] direction at FDR<0.05: longer-STR-up=%d, longer-STR-down=%d",
    sum(results$is_eSTR_005 & results$beta_bp > 0),
    sum(results$is_eSTR_005 & results$beta_bp < 0))
if (n_05 > 0) {
  abs_b <- abs(results[is_eSTR_005 == TRUE]$beta_bp)
  add("[RESULT] eSTR effect size |beta| (bp^-1, INT-norm scale): median %.3f, IQR %.3f-%.3f",
      median(abs_b), quantile(abs_b, 0.25), quantile(abs_b, 0.75))
}

# ---- 9. motif enrichment (if motif column carried through) -----------------
if ("motif" %in% colnames(results)) {
  bg  <- table(results$motif)
  hit <- table(results[is_eSTR_005 == TRUE]$motif)
  motifs <- names(bg)
  enr <- data.table(motif       = motifs,
                    n_tested    = as.integer(bg),
                    n_eSTR_005  = as.integer(hit[motifs]))
  enr[is.na(n_eSTR_005), n_eSTR_005 := 0L]
  enr[, pct_eSTR := 100 * n_eSTR_005 / n_tested]
  N_eSTR <- sum(enr$n_eSTR_005)
  N_bg   <- sum(enr$n_tested)
  enr[, fisher_p := vapply(seq_len(.N), function(i) {
    a <- n_eSTR_005[i]; b <- n_tested[i] - a
    c <- N_eSTR - a; d <- N_bg - N_eSTR - b
    fisher.test(matrix(c(a, b, c, d), 2, 2), alternative = "greater")$p.value
  }, numeric(1))]
  enr[, fisher_fdr := p.adjust(fisher_p, method = "BH")]
  setorder(enr, fisher_p)
  fwrite(enr, out_motif_enrich, sep = "\t")
  enriched <- enr[fisher_fdr < 0.05]
  add("[ENRICH] motifs enriched in eSTRs at Fisher FDR<0.05: %d", nrow(enriched))
  if (nrow(enriched) > 0)
    add("[ENRICH] top motifs: %s",
        paste(head(enriched$motif, 5), collapse = ", "))
}

# ---- 10. example plots (top 3 by smallest p) -------------------------------
# NB: results$gene_id is the *original* annotation ID (RefSeq NM_/NR_ when the
# annotation came from VEP's Feature column). expr_mat is keyed by ENSG
# (featureCounts/GENCODE), so the locus must be looked up via the ENSG key that
# was actually used during fitting (results$gene_ensg / ann_intron$gene_id_strip),
# NOT via strip_ver(gene_id) -- on an NM_ accession strip_ver is a no-op and the
# lookup fails with "subscript out of bounds".
safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

top_n <- min(3L, nrow(results))
if (top_n > 0) {
  top <- results[seq_len(top_n)]
  for (i in seq_len(top_n)) {
    locus_i <- top$locus[i]
    gene_i  <- top$gene_id[i]            # original ID (labels only)
    ensg_i  <- top$gene_ensg[i]          # ENSG key for expression lookup
    
    idx <- which(ann_intron[[locus_col]] == locus_i &
                   ann_intron[[gene_col]] == gene_i)[1]
    if (is.na(idx)) next
    
    # resolve the expression-matrix key; fall back to the row's gene_id_strip
    g <- ensg_i
    if (is.na(g) || !nzchar(g)) g <- ann_intron$gene_id_strip[idx]
    if (is.na(g) || !(g %in% rownames(expr_mat))) {
      add("[PLOT] skip %s / %s: gene key '%s' not in expression matrix",
          locus_i, gene_i, g)
      next
    }
    
    d  <- d_sub[idx, ]
    y  <- expr_mat[g, ]
    ok <- !is.na(d) & !is.na(y)          # exactly the samples that entered the model
    if (sum(ok) < 2L) next
    
    df_plot <- data.frame(dosage = d[ok], expr = y[ok],
                          cohort = cov_sub$cohort[ok],
                          sex    = cov_sub$sex[ok])
    gene_lab <- if (!is.na(gene_i) && gene_i != g)
      sprintf("%s (%s)", g, gene_i) else g
    pp <- ggplot(df_plot, aes(x = dosage, y = expr)) +
      geom_point(aes(color = cohort), alpha = 0.65, size = 1.8) +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
      scale_color_manual(values = c("R" = "#264653", "VB" = "#E76F51")) +
      labs(x = sprintf("STR length (bp)  —  locus %s", locus_i),
           y = sprintf("Inv-normal expression  —  %s", gene_lab),
           title = sprintf("eSTR #%d   %s ~ %s", i, gene_lab, locus_i),
           subtitle = sprintf("beta = %.3f / bp,  p = %.2e,  FDR = %.2e,  n = %d",
                              top$beta_bp[i], top$p[i], top$fdr[i], top$n[i])) +
      theme_minimal(base_size = 11)
    ggsave(file.path(out_examples_dir,
                     sprintf("eSTR_example_%02d_%s_%s.pdf",
                             i, safe_name(g), safe_name(locus_i))),
           pp, width = 5.2, height = 4.0)
  }
  add("[OUTPUT] example PDFs in %s/", out_examples_dir)
}

# ---- 11. write outputs ------------------------------------------------------
fwrite(results, out_catalog, sep = "\t")
add("[OUTPUT] full catalog: %s", out_catalog)

writeLines(c(
  sprintf("# eSTR mapping revision values (auto-generated %s)",
          format(Sys.time(), "%Y-%m-%d %H:%M")),
  rv,
  sprintf("# Run time: %.1f min",
          as.numeric(difftime(Sys.time(), t0, units = "mins")))
), out_rv)

cat(sprintf("\nDone in %.1f min.\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat("  ", out_catalog,      "  -- one row per STR-gene pair (FDR adjusted)\n")
cat("  ", out_motif_enrich, "  -- motif enrichment (if motif column present)\n")
cat("  ", out_rv,           "  -- revision values for writing list\n")
cat("  ", out_examples_dir, "/  -- top-3 example scatter PDFs\n")

source("eSTR_03_sensitivity_PC1.R")
library(patchwork)
library(ggplot2)
source("eSTR_S18_patch.R")
