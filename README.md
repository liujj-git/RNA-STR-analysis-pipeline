```markdown
# RNA-STR analysis pipeline

This repository contains the R scripts used for the transcriptome‑wide characterization, genotyping reliability assessment, and forensic/disease‑related analysis of polymorphic short tandem repeats (STRs) in human blood, as described in:

> **An Atlas of Polymorphic Short Tandem Repeats in the Human Blood Transcriptome**  
> (Manuscript under review)

The pipeline processes RNA‑seq and DNA‑seq data after initial read alignment and STR genotyping (performed with **HipSTR** and **ExpansionHunter**; see Materials and Methods in the manuscript). All custom R code is provided to reproduce the main figures and supplementary results.

---

## Repository structure

```
.
├── 01_expression_analysis.R
├── 02_annotation.R
├── 03_genetic_character.R
├── 04_concordance_analysis.R
├── 05_polymorphism.R
├── 06_disease-associated_STRs.R
└── README.md
```

---

## Dependencies

The scripts were developed and tested with:

- **R** (≥ 4.2.0)
- Required R packages:
  - `tidyverse` (dplyr, tidyr, readxl, ggplot2, purrr)
  - `vcfR`, `data.table`
  - `ggpubr`, `ggsci`, `ggrepel`, `EnvStats`, `pheatmap`, `ggvenn`
  - `rjson`, `jsonlite`

Install missing packages with:

```r
install.packages(c("tidyverse", "vcfR", "data.table", "ggpubr", "ggsci",
                   "ggrepel", "EnvStats", "pheatmap", "ggvenn", "rjson", "jsonlite"))
```

---

## Input data required (before running scripts)

The scripts expect the following files (paths are relative; adjust as needed):

| File | Description |
|------|-------------|
| `../GRCh38.hipstr_reference.refine.bed` | HipSTR reference BED (filtered for 2‑6 bp motifs, ≤50 repeat units) |
| `../VBsampleinfo.xlsx` | Excel file with sample metadata (sheets: `unrelated_individuals`, `repeatability`, `concordance`, `total`) |
| `hipstr.filtered.rna.vcf.gz` | HipSTR VCF output for RNA‑seq data |
| `pstr.rna.recode.vcf.gz` | Re‑coded VCF for RNA‑seq pSTRs (used by `polymorphism.R`) |
| `pstr.rna.recode.annotated.vcf.gz` | VEP‑annotated VCF |
| `../pstr.recode.vcf.gz` | Combined VCF for RNA‑seq and DNA‑seq (used by `concordance.R`) |
| `../2_annotation/STR_annotation_group.txt` | Output from `02_annotation.R` (pSTR annotation groups) |
| `../1_expression/locus_list.polymorphic.txt` | Output from `01_expression_analysis.R` (list of pSTR loci) |
| `disease.rna.vcf`, `disease.dna.vcf` | ExpansionHunter output for 31 disease‑associated STRs |
| `eh.hg38.variant_catalog.disease.json` | ExpansionHunter disease catalog |
| `disease.locusposition.txt`, `locus_intersection.txt`, `disease_threshold.txt` | Custom files for disease‑STR coordinates and thresholds |

> **Note**: The initial steps – read filtering (Fastp), alignment (STAR/BWA‑MEM), STR genotyping (HipSTR/ExpansionHunter), and VEP annotation – were performed using command‑line tools with default parameters except where noted in the manuscript (e.g., `--no-rmdup` for RNA‑seq HipSTR). Those steps are **not** part of these R scripts.

---

## Script descriptions and order of execution

Run the scripts in the order listed below (each script depends on outputs from previous ones).

### 1. `01_expression_analysis.R`
**Purpose**:  
- Compute genotyping rates, classify loci as polymorphic (pSTR) or monomorphic (mSTR).  
- Analyse expression depth (DP) and generate basic plots (boxplots, correlations).  

**Main outputs**:  
- `Loci_withnoexpression.txt` – loci with genotyping rate = 0.  
- `locus_list.txt` – all filtered loci with mean DP, group, coordinates.  
- `locus_list.polymorphic.txt` – list of pSTR loci for downstream analyses.  
- PDF figures: `genotyping_rate_boxplot.pdf`, `genotyping_rate_corplot.pdf`, `expression_boxplot.pdf`, `expression_corplot.pdf`, `pheatmap_polymorphic.pdf`, `pheatmap_monomorphic.pdf`, `gr_cfl.pdf`, `ho_cfl.pdf`, `heatmap_summary.pdf`.

### 2. `02_annotation.R`
**Purpose**:  
- Parse VEP annotation and assign each STR to a final genomic feature (CDS, 5′UTR, 3′UTR, intron, intergenic) using a prioritisation hierarchy.  
- Extract biotype information.  

**Main output**:  
- `STR_annotation_group.txt` – table with locus, final group, order, feature type, gene ID, biotype.

### 3. `03_genetic_character.R`
**Purpose**:  
- Build genetic characterisation table for pSTRs (period size, copy number, motif).  
- Generate bar plots of biotype distribution and pSTR composition (total, mRNA‑derived, lncRNA‑derived).  
- Plot density distributions of reference allele length and copy number across regions.  
- Analyse motif patterns (collapsing reverse complements and circular permutations).  

**Main outputs**:  
- `genetic_character_df.txt` – per‑pSTR characteristics.  
- `biotype_barplot.pdf`, `pSTR_counts_barplot*.pdf`, `ral_distribution_density.pdf`, `ncopy_distribution_density.pdf`.  
- `motif_df.txt`, `motif_barplot.pdf`, `CDS_tri_motif_barplot.pdf`, `intergenic_*_motif_barplot.pdf`.  
- Supplementary files: `ral_mean_by_region.txt`, `ncopy_mean_by_region.txt`, `biotype_pieplot.pdf`.

### 4. `04_concordance_analysis.R`
**Purpose**:  
- Build allele‑sequence ladder (including iso‑allele resolution).  
- Assess genotyping concordance in RNA technical replicates (RNA‑RNA) and RNA‑DNA matched pairs.  
- Quantify step differences for partial‑identical calls and analyse sequence discordance (iso‑allelic variants).  

**Main outputs**:  
- `allele_sequence.concordance.txt`, `concordance_repeated_df.txt`, `concordance_compared_df.txt`.  
- Multiple PDFs: `repeatedsample_GTcounts_barplot.pdf`, `concordance_repeated_sample_plot*.pdf`, `stepdif_halfidentical_plot*.pdf`, `PHPcounts_histogram*.pdf`, `iav_changed_base_pattern_*.pdf`, `concordance_compared_sample_plot*.pdf`, `concordance_region_*.pdf`, `concordance_allele_length_*.pdf`, `Discordance_freq_line.pdf`.

### 5. `05_polymorphism.R`
**Purpose**:  
- Convert HipSTR genotypes to repeat‑count format (ISFG naming).  
- Compute forensic parameters: Na, Ho, He, MP, DP, PE2, PE3, PIC, Ae.  
- Generate allele frequency tables.  
- Analyse major allele frequency, copy‑number differences between major/reference alleles, and correlations with PIC.  

**Main outputs**:  
- `allele_sequence.txt`, `GT.copynumber.txt`.  
- `forensic_parameters.txt`, `allele_freq.txt`.  
- `reported_locus.txt`, `forensic_parameters.reported_locus.txt`.  
- PDFs: `region_na_boxplot.pdf`, `region_pic_boxplot.pdf`, `major_allele_freq_histogram.pdf`, `major_alt_ncopy_dif_barplot.pdf`, `major_reference_ncopy_dif_barplot.pdf`, `major_ncopy_PIC_corplot.pdf`, `reference_ncopy_PIC_corplot.pdf`.

### 6. `06_disease-associated_STRs.R`
**Purpose**:  
- Process ExpansionHunter output for 31 disease‑associated STRs.  
- Compute genotyping rates in RNA and DNA samples.  
- Evaluate concordance in RNA replicates and RNA‑DNA pairs.  
- Identify alleles exceeding pathogenic thresholds.  

**Main outputs**:  
- `disease_threshold.config.txt`, `disease_allele_summary.txt`.  
- PDFs: `genotyping_rate_barplot*.pdf`, `repeated_*_concordance_barplot.pdf`, `compared_*_concordance_barplot.pdf`, `genotyping_rate_barplot3.pdf`.

---

## Usage example

After preparing all input files (see above) and installing required R packages, run the scripts sequentially:

```bash
Rscript 01_expression_analysis.R
Rscript 02_annotation.R
Rscript 03_genetic_character.R
Rscript 04_concordance_analysis.R
Rscript 05_polymorphism.R
Rscript 06_disease-associated_STRs.R
```

> **Path adjustment**: The scripts use relative paths (e.g., `../GRCh38.hipstr_reference.refine.bed`). If you place the scripts in a subfolder (e.g., `scripts/`), ensure that the input files are in the parent directory, or modify the paths accordingly.

---

## License

This code is released under the **MIT License** (see `LICENSE` file in the repository).  
If you use this code in your work, please cite our manuscript (citation details will be added upon publication).

---

## Contact

For questions or issues, please open an issue on this GitHub repository or contact the corresponding author:  
[Your Name / Email – to be added]
```
