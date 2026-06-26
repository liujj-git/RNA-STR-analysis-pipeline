# =============================================================================
# annotation.R  —  clean, revision build
# Manuscript: Section 2.3 ; Figure 2a (region composition), Figures S4/S5 (biotype)
# Pipeline output consumed by scripts 3-6: STR_annotation_group.txt
#
# DESIGN (restricted-revision):
#   * final_group (CDS/5_UTR/3_UTR/intron/intergenic) uses the annotation_order
#     columns `group` + `order`, which are UNCHANGED -> Figure 2a and every
#     region-stratified figure in scripts 3-6 reproduce, no re-code.
#   * NEW: the intergenic class is decomposed via the two columns just added to
#     annotation_order.xlsx -- `subgroup` (gene_proximal / ncRNA_exonic /
#     ncRNA_associated / intergenic_true) and `suborder` (within-intergenic
#     priority). This is authoritative and consistent with final_group, and
#     answers Reviewer 1, point 5, without changing any existing number.
#     (Backed up by an assumption-free raw consequence-term frequency table.)
#   * CHANGED: biotype is now the biotype of the transcript that produced the
#     highest-priority consequence (consistent with final_group) -> fixes
#     Reviewer 2's Figure S5 concern. This shifts the biotype percentages and
#     Figures S4/S5; genetic_character.R must be RE-RUN (not re-coded). The
#     shift is quantified in biotype_reconciliation.txt.
#   * FIXED: the if(grepl(...)) length>1 hazard (R >= 4.2) and the fragile single
#     col-5/col-14 parsing, both handled per-transcript.
#
# REQUIRES the UPDATED annotation_order.xlsx (with subgroup + suborder columns)
# and the VEP-tab output. Verify region_composition.txt against Figure 2a after
# running, and eyeball intergenic_term_frequency.txt.
# =============================================================================

rm(list = ls())

library(readxl)
library(ggplot2)
library(ggsci)

base_theme <- theme_bw() +
  theme(axis.title  = element_text(face = "bold", colour = "black", size = 13),
        axis.text.x = element_text(face = "bold", colour = "black", size = 9, angle = 30, hjust = 1),
        axis.text.y = element_text(face = "bold", colour = "black", size = 10))

region_levels <- c("CDS", "5_UTR", "3_UTR", "intron", "intergenic")

# =============================================================================
# 0. Inputs
# =============================================================================
# annotation_order.xlsx columns: group | annotation (VEP term) | order |
#   subgroup (fine class, only meaningful for the intergenic group) | suborder
annotation_order <- read_excel("annotation_order.xlsx")

# lookups keyed by VEP consequence term
ord_lu  <- setNames(as.numeric(annotation_order$order),                       annotation_order$annotation)
grp_lu  <- setNames(as.character(annotation_order$group),                     annotation_order$annotation)
sub_lu  <- setNames(as.character(annotation_order$subgroup),                  annotation_order$annotation)
subo_lu <- setNames(suppressWarnings(as.numeric(annotation_order$suborder)),  annotation_order$annotation)

# VEP tab output (one line per transcript-consequence). Columns used (as in the
# original): 1 = locus | 5 = gene/feature id | 7 = Consequence | 14 = Extra.
# read.table skips '#'-prefixed header lines; fread() is a faster alternative.
ann <- read.table("pstr.rna.recode.annotated.vcf.gz", sep = "\t", quote = "",
                  comment.char = "#", fill = TRUE, stringsAsFactors = FALSE)
ann <- ann[, c(1, 5, 7, 14)]
colnames(ann) <- c("locus", "gene", "consequence", "extra")

# =============================================================================
# 1. Helpers
# =============================================================================
term_order <- function(tt) { o <- ord_lu[tt]; o[is.na(o)] <- Inf; as.numeric(o) }  # unknown -> lowest
term_group <- function(tt) as.character(grp_lu[tt])

parse_biotype <- function(extra) {                 # one VEP line's Extra field
  if (length(extra) == 0 || is.na(extra) || !grepl("BIOTYPE=", extra)) return(NA_character_)
  kv <- unlist(strsplit(extra, ";"))
  bt <- kv[grepl("BIOTYPE=", kv)]
  if (length(bt) == 0) return(NA_character_)
  sub(".*BIOTYPE=", "", bt[1])
}

classify_locus <- function(rows) {
  terms_per_row <- strsplit(rows$consequence, ",")
  row_min  <- vapply(terms_per_row, function(tt) min(term_order(tt)), numeric(1))
  win      <- which.min(row_min)                   # transcript with top-priority consequence
  wt       <- terms_per_row[[win]]
  wo       <- term_order(wt)
  fgroup   <- term_group(wt[which.min(wo)])        # SAME final_group as the original
  all_terms <- unique(unlist(terms_per_row))
  # intergenic subtype from the authoritative table: among this locus's terms
  # belonging to the intergenic group, take the one with the smallest suborder.
  ig_sub <- NA_character_
  if (identical(fgroup, "intergenic")) {
    igt <- all_terms[grp_lu[all_terms] %in% "intergenic"]
    if (length(igt) > 0) ig_sub <- sub_lu[igt[which.min(subo_lu[igt])]]
  }
  row_biotypes <- vapply(rows$extra, parse_biotype, character(1))
  data.frame(
    locus              = rows$locus[1],
    final_group        = fgroup,
    final_order        = min(wo),
    feature_type       = paste(all_terms, collapse = ","),
    gene_id            = rows$gene[win],
    biotype            = parse_biotype(rows$extra[win]),   # reconciled to winning transcript
    multi_biotype      = length(unique(na.omit(row_biotypes))) > 1,
    intergenic_subtype = ig_sub,
    stringsAsFactors = FALSE, row.names = NULL)
}

# =============================================================================
# 2. Classify every locus  (single rbind; no growing loop)
# =============================================================================
feature_df <- do.call(rbind, lapply(split(ann, ann$locus), classify_locus))
rownames(feature_df) <- NULL

## pipeline output (column NAMES locus/final_group/biotype preserved for 3-6) --
write.table(feature_df[, c("locus", "final_group", "final_order",
                           "feature_type", "gene_id", "biotype", "intergenic_subtype")],
            "STR_annotation_group.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# =============================================================================
# 3. Revision values  ->  written to files for the writing list
# =============================================================================
## region composition  (Figure 2a -- should reproduce 82.41% intron, etc.) ----
region_tab <- as.data.frame(table(factor(feature_df$final_group, levels = region_levels)))
colnames(region_tab) <- c("region", "count")
region_tab$percentage <- round(100 * region_tab$count / sum(region_tab$count), 2)
write.table(region_tab, "region_composition.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

## intergenic decomposition  (NEW -- Reviewer 1, point 5) ---------------------
ig <- feature_df[feature_df$final_group == "intergenic", ]
ig_tab <- as.data.frame(table(ig$intergenic_subtype))
colnames(ig_tab) <- c("subtype", "count")
ig_tab$pct_of_intergenic <- round(100 * ig_tab$count / sum(ig_tab$count), 2)
ig_tab$pct_of_all_pSTR   <- round(100 * ig_tab$count / nrow(feature_df), 2)
ig_tab <- ig_tab[order(ig_tab$count, decreasing = TRUE), ]
write.table(ig_tab, "intergenic_composition.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# assumption-free backup: raw consequence-term frequency among intergenic loci
ig_terms <- as.data.frame(table(unlist(strsplit(ig$feature_type, ","))))
colnames(ig_terms) <- c("consequence_term", "count")
ig_terms <- ig_terms[order(ig_terms$count, decreasing = TRUE), ]
write.table(ig_terms, "intergenic_term_frequency.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

## biotype distribution after reconciliation  (Reviewer 2 -- Figure S5) -------
bt <- feature_df$biotype
bt[is.na(bt)] <- "NoAnnotation"
bt[bt == "protein_coding"] <- "mRNA (protein_coding)"
bt_tab <- as.data.frame(table(bt))
colnames(bt_tab) <- c("biotype", "count")
bt_tab$percentage <- round(100 * bt_tab$count / sum(bt_tab$count), 2)
bt_tab <- bt_tab[order(bt_tab$count, decreasing = TRUE), ]
write.table(bt_tab, "biotype_reconciliation.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

## NEW figure: composition of intergenic-annotated pSTRs ----------------------
ig_plot <- ig_tab
ig_plot$subtype <- factor(ig_plot$subtype, levels = ig_plot$subtype)
intergenic_composition_barplot <- ggplot(ig_plot, aes(subtype, count, fill = subtype)) +
  geom_col(colour = "black") +
  geom_text(aes(label = sprintf("%d (%.1f%%)", count, pct_of_intergenic)), vjust = -0.3, size = 3) +
  scale_fill_lancet() + base_theme + theme(legend.position = "none") +
  labs(title = "Composition of intergenic-annotated pSTRs",
       x = "Intergenic subtype", y = "pSTR count")
ggsave("intergenic_composition_barplot.pdf", intergenic_composition_barplot, width = 8, height = 5)

## master reference file ------------------------------------------------------
rv <- character(0); add <- function(...) rv[[length(rv) + 1]] <<- sprintf(...)
add("# annotation.R revision values  (auto-generated %s)", as.character(Sys.Date()))
add("# loci annotated: %d", nrow(feature_df))
add("")
add("## [REGION COMPOSITION]  (Section 2.3, Figure 2a) -- should reproduce")
for (i in seq_len(nrow(region_tab)))
  add("%s\t%d\t%.2f%%", region_tab$region[i], region_tab$count[i], region_tab$percentage[i])
add("# full table -> region_composition.txt")
add("")
add("## [INTERGENIC DECOMPOSITION]  (NEW; Reviewer 1 point 5) -> intergenic_composition.txt")
for (i in seq_len(nrow(ig_tab)))
  add("%s\t%d\t%.2f%% of intergenic\t%.2f%% of all pSTR",
      ig_tab$subtype[i], ig_tab$count[i], ig_tab$pct_of_intergenic[i], ig_tab$pct_of_all_pSTR[i])
add("# raw consequence terms -> intergenic_term_frequency.txt")
add("")
add("## [BIOTYPE AFTER RECONCILIATION]  (Reviewer 2; Figures S4/S5) -> biotype_reconciliation.txt")
for (i in seq_len(nrow(bt_tab)))
  add("%s\t%d\t%.2f%%", bt_tab$biotype[i], bt_tab$count[i], bt_tab$percentage[i])
add("loci_with_multiple_transcript_biotypes\t%d / %d\t(reconciliation can affect these)",
    sum(feature_df$multi_biotype), nrow(feature_df))
writeLines(rv, "revision_values_annotation.txt")

cat("annotation.R done. final_group unchanged; intergenic decomposed; biotype reconciled.\n")
cat(sprintf("Re-run genetic_character.R afterwards: biotype changed for up to %d multi-biotype loci.\n",
            sum(feature_df$multi_biotype)))




# ============================================================================
# VERIFICATION BLOCK — paste at end of annotation.R, return me the console output
# ============================================================================
cat("\n========== NUMBERS TO CONFIRM (annotation.R) ==========\n")

# (1) Intergenic decomposition — resolve 7.6% vs 7.7%
cat("\n--- Fig S5: intergenic subtypes (exact, 2 decimals) ---\n")
ig_check <- as.data.frame(table(feature_df$intergenic_subtype[feature_df$final_group=="intergenic"]))
colnames(ig_check) <- c("subtype","count")
ig_check$pct_of_intergenic <- 100*ig_check$count/sum(ig_check$count)
ig_check$pct_of_all_pSTR   <- 100*ig_check$count/nrow(feature_df)
print(ig_check, digits=4)
cat(sprintf("intergenic total = %d ; all pSTR = %d\n",
            sum(ig_check$count), nrow(feature_df)))
# the 146 value: print at 3 decimals so we KNOW if it rounds to 7.6 or 7.7
cat(sprintf("ncRNA_exonic pct = %.3f%% of intergenic\n",
            100*ig_check$count[ig_check$subtype=="ncRNA_exonic"]/sum(ig_check$count)))

# (2) CDS count chain — resolve 685 vs 693
cat("\n--- CDS pSTR composition by biotype (resolves 685 vs 693) ---\n")
cds <- feature_df[feature_df$final_group=="CDS", ]
cds_bt <- table(ifelse(is.na(cds$biotype),"NoAnnotation",cds$biotype))
print(cds_bt)
cat(sprintf("TOTAL CDS pSTR (all biotypes, = Fig 2A) = %d\n", nrow(cds)))
cat(sprintf("  of which protein_coding (= Fig S7 left) = %d\n",
            sum(cds$biotype=="protein_coding", na.rm=TRUE)))
cat(sprintf("  of which lncRNA (= Fig S7 right)         = %d\n",
            sum(cds$biotype=="lncRNA", na.rm=TRUE)))
cat(sprintf("  remainder (other biotypes)               = %d\n",
            nrow(cds) - sum(cds$biotype %in% c("protein_coding","lncRNA"), na.rm=TRUE)))

# (3) The 573 link — is Fig S6 NoAnnotation == Fig S5 intergenic_true?
cat("\n--- Does NoAnnotation (Fig S6) == intergenic_true (Fig S5)? ---\n")
n_noann <- sum(is.na(feature_df$biotype))
n_igtrue <- sum(feature_df$intergenic_subtype=="intergenic_true", na.rm=TRUE)
cat(sprintf("NoAnnotation biotype = %d ; intergenic_true = %d ; identical loci? %s\n",
            n_noann, n_igtrue,
            identical(sort(feature_df$locus[is.na(feature_df$biotype)]),
                      sort(feature_df$locus[which(feature_df$intergenic_subtype=="intergenic_true")]))))

# (4) Region composition headline (confirm 82.41 / 7.95 / 5.98 / 2.89 / 0.78)
cat("\n--- Fig 2A region composition ---\n")
print(region_tab)

## ---- Fig S5: intergenic decomposition (annotation.R) ----
library(ggplot2)
theme_jgg <- theme_classic(base_family = "Arial", base_size = 9) +
  theme(axis.title = element_text(colour="black", size=9),
        axis.text  = element_text(colour="black", size=8),
        axis.line  = element_line(colour="black", linewidth=0.4),
        axis.ticks = element_line(colour="black", linewidth=0.4),
        strip.background = element_rect(fill="grey92", colour=NA),
        strip.text = element_text(colour="black", size=8, margin=margin(2,2,2,2)),
        legend.title = element_text(size=8), legend.text = element_text(size=7),
        legend.key.size = unit(3,"mm"),
        plot.tag = element_text(face="bold", size=11))
ig_s5 <- ig_tab
ig_s5$subtype <- factor(ig_s5$subtype, levels=ig_s5$subtype[order(ig_s5$count, decreasing=TRUE)])
sub_lab <- c(gene_proximal="Gene-proximal", ncRNA_exonic="ncRNA-exonic",
             intergenic_true="True intergenic")
figS5 <- ggplot(ig_s5, aes(subtype, count, fill=subtype)) +
  geom_col(colour="black", linewidth=0.3, width=0.7) +
  geom_text(aes(label=sprintf("%d (%.1f%%)", count, pct_of_intergenic)),
            vjust=-0.4, size=2.6) +
  scale_fill_manual(values=c("#2D6DB1","#2E933C","#E69F00"), guide="none") +
  scale_x_discrete(labels=sub_lab) +
  scale_y_continuous(expand=expansion(mult=c(0,0.12))) +
  labs(x="Intergenic subtype", y="pSTR count") + theme_jgg
ggsave("FigureS5.pdf", figS5, width=5, height=4, device=cairo_pdf)
