# RNA-STR-analysis-pipeline

Analysis code for the study **"Transcriptome-Wide Identification and Reliability
Assessment of Polymorphic Short Tandem Repeats in Human Blood."**

The pipeline identifies polymorphic short tandem repeats (pSTRs) directly from
blood RNA-seq, assesses their genotyping reliability against matched
whole-genome sequencing (WGS), characterizes RNA-specific signals (allelic
dropout and imbalance, RNA editing), and maps expression STRs (eSTRs).

This repository contains the **downstream R analysis** that takes per-locus STR
genotypes (from upstream calling) through to the figures and tables in the
manuscript. Upstream read processing and genotype calling are summarized under
[Upstream processing](#upstream-processing) below.

---

## Repository layout

```
scripts/
  01_expression.R               Section 2.1-2.2 | Figure 1, Figures S1-S3, Tables S1-S2
  01b_figS1_S2_S4.R             -> run after 01; Figures S1, S2, S4 (final styling/numbering)

  02_annotation.R               Section 2.3 | locus annotation; input for 03 and 05/06
  03_genetic_character.R        Section 2.3 | Figure 2, Figures S4-S6

  04_concordance.R              Section 2.4 | genotype-level RNA-DNA / RNA-RNA concordance
  04a_ase_analysis.R            -> run after 04; allele-specific expression (ASE) tables
  04b_fig3_concordance.R        -> run after 04; Figure 3, Figures S9-S12
  04c_fig4_ASE.R                -> run after 04 + 04a; Figure 4, Figure S13

  05_polymorphism.R             Section 2.5 | forensic informativeness of pSTRs
  05b_fig5_allele_diversity.R   -> run after 05; Figure 5, Figure S14

  06_disease.R                  Section 2.7 | disease-associated STRs (ExpansionHunter)

  eSTR/
    eSTR_01_expression_prep.R   build the 167-individual gene-expression matrix
    eSTR_02_mapping.R           Section 2.8 | intronic eSTR mapping; eSTR catalog
    eSTR_02b_figS17.R           -> run after eSTR_02; Figure S17 (top-3 eSTR examples)
    eSTR_03_sensitivity_PC1.R   batch/PC1 sensitivity check for the eSTR set
```

### Script-to-figure map

| Manuscript item | Produced by |
|---|---|
| Figure 1, Tables S1-S2 | `01_expression.R` |
| Figure 2 | `02_annotation.R` + `03_genetic_character.R` |
| Figure 3 | `04_concordance.R` + `04b_fig3_concordance.R` |
| Figure 4 | `04_concordance.R` + `04a_ase_analysis.R` + `04c_fig4_ASE.R` |
| Figure 5 | `05_polymorphism.R` + `05b_fig5_allele_diversity.R` |
| Figure 6 | (schematic; no code) |
| Figures S1, S2, S4 | `01b_figS1_S2_S4.R` |
| Figure S3 | `01_expression.R` |
| Figures S4-S6 | `03_genetic_character.R` |
| Figures S9-S12 | `04b_fig3_concordance.R` |
| Figure S13 | `04c_fig4_ASE.R` |
| Figure S14 | `05b_fig5_allele_diversity.R` |
| Figures S15-S16 | `06_disease.R` |
| Figure S17 | `eSTR/eSTR_02b_figS17.R` |
| Tables S3-S10 | written by the corresponding section scripts above |

> **Note on numbering.** A few scripts were originally written as revision
> "patches" and their *internal* comments may reference an earlier figure
> numbering (e.g. "Fig S18", "Fig 4" for the allele-diversity panel). The
> **filenames and the banner at the top of each script reflect the final
> manuscript numbering**, which is authoritative. `fig4_v2.R` (a richer draft
> of Figure 4) was discarded in favour of `04c_fig4_ASE.R`.

---

## Running the pipeline

Scripts are grouped by manuscript section. Within a section, run the main
script first, then `source()` the lettered scripts **in the same R session**,
because they reuse in-memory objects produced by the main script.

```r
# Section 2.1-2.2
source("scripts/01_expression.R")
source("scripts/01b_figS1_S2_S4.R")

# Section 2.3
source("scripts/02_annotation.R")
source("scripts/03_genetic_character.R")

# Section 2.4   (order matters: 04 -> 04a -> 04b/04c)
source("scripts/04_concordance.R")
source("scripts/04a_ase_analysis.R")
source("scripts/04b_fig3_concordance.R")
source("scripts/04c_fig4_ASE.R")

# Section 2.5
source("scripts/05_polymorphism.R")
source("scripts/05b_fig5_allele_diversity.R")

# Section 2.7
source("scripts/06_disease.R")

# Section 2.8  (eSTR)
source("scripts/eSTR/eSTR_01_expression_prep.R")
source("scripts/eSTR/eSTR_02_mapping.R")
source("scripts/eSTR/eSTR_02b_figS17.R")
source("scripts/eSTR/eSTR_03_sensitivity_PC1.R")
```

`02_annotation.R` must be run before `03_genetic_character.R`, `05_polymorphism.R`,
and `06_disease.R`, since it writes the locus-annotation table those scripts read.

---

## Inputs

The scripts expect, per section, the STR genotype call sets and supporting
tables produced by the upstream steps, including:

- HipSTR VCFs for the transcriptome-wide pSTR/mSTR call set (RNA and WGS);
- the relaxed-allele-bias VCF used for the allele-level / ASE analysis;
- ExpansionHunter outputs for the disease-associated loci;
- a gene-level expression matrix (featureCounts) for the eSTR analysis;
- the locus annotation table and the sample sheet.

Paths are set near the top of each script; adjust them to your local layout.
Input data are available under controlled access (see the manuscript's Data
Availability statement); the VCFs and matrices themselves are not redistributed
in this repository.

---

## Upstream processing

Read processing and genotype calling were performed with the following tools
(versions and non-default parameters as reported in the manuscript Methods):

- **Fastp** v0.23.2 - read filtering
- **STAR** v2.7.10 - RNA-seq alignment
- **BWA-MEM** v0.7.17 - WGS alignment
- **HipSTR** v0.6.2 - transcriptome-wide STR genotyping (`--no-rmdup` for RNA-seq)
- **ExpansionHunter** v5.0.0 - disease-associated STR genotyping
- **DumpSTR** (TRTools) v6.1.0 - genotype filtering

---

## Dependencies

R (>= 4.0) with, among others: `vcfR`, `data.table`, `tidyr`, `readxl`,
`ggplot2`, `ggpubr`, `ggrepel`, `ggsci`, `patchwork`, `EnvStats`, `pheatmap`,
`rjson`, `jsonlite`. Figures use an Arial-capable device (`cairo_pdf`).

Install the core plotting/analysis packages with:

```r
install.packages(c("vcfR","data.table","tidyr","readxl","ggplot2","ggpubr",
                   "ggrepel","ggsci","patchwork","EnvStats","pheatmap",
                   "rjson","jsonlite"))
```

---

## Citation

If you use this code, please cite the associated manuscript (citation to be
added upon publication).

## License

Released under the terms of the [LICENSE](LICENSE) file in this repository.
