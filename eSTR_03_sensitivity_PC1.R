# =============================================================================
# eSTR_03_sensitivity_PC1.R   --  batch-control sensitivity check
#
# Question: do the eSTRs depend on HOW the data-collection/batch effect is handled?
#
# The main analysis (eSTR_02) controls this structure with the expression PCs
# (PC1..PC9 retained, no explicit cohort term; PC1 is ~97.5% explained by the
# batch axis). A reviewer may ask whether residual batch is creating the hits.
# Here we re-fit the SAME 685 STR-gene pairs under three specifications and check
# that the effect sizes, p-values, and the eSTR set are stable across them:
#
#   B  PC1..PC9, no cohort  (PRIMARY -- identical to the paper)
#   A  cohort + PC2..PC9    (explicit cohort term, PC1 dropped)
#   C  cohort + PC1..PC9    (belt-and-suspenders; accepts mild collinearity)
#
# If B/A/C give near-identical betas (Spearman rho ~1) and the FDR<0.10 hit set
# overlaps almost completely, the result is not an artifact of how batch is
# handled -> that is the sentence you put in the rebuttal.
#
# HOW TO RUN
#   Easiest: run this in the SAME R session right after eSTR_02_mapping.R
#   finishes (it reuses the in-memory objects d_sub, expr_mat, cov_sub,
#   ann_intron, results, locus_col). If you start a fresh session, re-run
#   eSTR_02 first (through the fitting section) so those objects exist.
# =============================================================================
suppressPackageStartupMessages(library(data.table))

stopifnot(exists("d_sub"), exists("expr_mat"), exists("cov_sub"),
          exists("ann_intron"), exists("results"), exists("locus_col"))

# covariates available (cov_sub already has sex, cohort as factors + PC1..PCk)
pcs_all <- grep("^PC", colnames(cov_sub), value = TRUE)        # PC1..PC9
if (!is.factor(cov_sub$sex))    cov_sub$sex    <- factor(cov_sub$sex)
if (!is.factor(cov_sub$cohort)) cov_sub$cohort <- factor(cov_sub$cohort)

# ---- refit the dosage coefficient for every pair, given a covariate formula --
# Mirrors eSTR_02 section 7 exactly; only the base design matrix changes.
fit_dosage <- function(fmla, min_n = 30L) {
  X_base <- model.matrix(fmla, data = cov_sub)
  nL     <- nrow(d_sub)
  beta   <- rep(NA_real_, nL); pval <- rep(NA_real_, nL)
  for (i in seq_len(nL)) {
    d  <- d_sub[i, ]; ok <- !is.na(d)
    if (sum(ok) < min_n) next
    g  <- ann_intron$gene_id_strip[i]
    if (is.na(g) || !(g %in% rownames(expr_mat))) next
    y  <- expr_mat[g, ok]
    ok2 <- !is.na(y); y <- y[ok2]
    if (length(y) < min_n) next
    Xi  <- cbind(d = d[ok][ok2], X_base[ok, , drop = FALSE][ok2, , drop = FALSE])
    fit <- tryCatch(.lm.fit(Xi, y), error = function(e) NULL); if (is.null(fit)) next
    rss  <- sum(fit$residuals^2); df.r <- length(y) - ncol(Xi); if (df.r < 1L) next
    iv   <- tryCatch(chol2inv(chol(crossprod(Xi))), error = function(e) NULL); if (is.null(iv)) next
    se   <- sqrt(iv[1, 1] * rss / df.r)
    beta[i] <- fit$coefficients[1]
    pval[i] <- 2 * pt(abs(beta[i] / se), df = df.r, lower.tail = FALSE)
  }
  data.table(locus     = ann_intron[[locus_col]],
             gene_ensg = ann_intron$gene_id_strip,
             beta = beta, p = pval, fdr = p.adjust(pval, "BH"))[!is.na(p)]
}

specs <- list(
  A_cohort_PC2_9   = as.formula(paste("~ sex + cohort +",
                                       paste(setdiff(pcs_all, "PC1"), collapse = " + "))),
  B_PC1_9_noCohort = as.formula(paste("~ sex +", paste(pcs_all, collapse = " + "))),
  C_cohort_PC1_9   = as.formula(paste("~ sex + cohort +", paste(pcs_all, collapse = " + ")))
)

cat("Refitting 685 pairs under 3 batch-control specs...\n")
fits <- lapply(specs, fit_dosage)

# ---- compare each spec to the primary (B = PC1..PC9, no cohort) --------------
base <- fits[["B_PC1_9_noCohort"]]
hits_base <- base[fdr < 0.10]$locus
cmp <- rbindlist(lapply(names(fits), function(nm) {
  f <- fits[[nm]]
  m <- merge(base[, .(locus, gene_ensg, beta0 = beta, p0 = p)],
             f[,    .(locus, gene_ensg, beta1 = beta, p1 = p)],
             by = c("locus", "gene_ensg"))
  hits <- f[fdr < 0.10]$locus
  data.table(
    spec             = nm,
    n_pairs          = nrow(f),
    beta_rho         = round(cor(m$beta0, m$beta1, method = "spearman"), 4),
    neglogp_rho      = round(cor(-log10(m$p0), -log10(m$p1), method = "spearman"), 4),
    n_FDR_lt_0.10    = length(hits),
    overlap_with_primary = length(intersect(hits_base, hits)))
}))

cat("\n================ batch-handling sensitivity ================\n")
print(cmp)
cat("\nPrimary (B, PC1..PC9 no cohort) FDR<0.10 hits:", length(hits_base), "\n")
cat("Interpretation: beta_rho and neglogp_rho near 1.0, plus near-complete\n",
    "overlap with the primary, mean the eSTRs are robust to how batch is handled.\n", sep = "")

fwrite(cmp, "eSTR_sensitivity_PC1_vs_cohort.tsv", sep = "\t")
cat("\n[OUTPUT] eSTR_sensitivity_PC1_vs_cohort.tsv\n")
