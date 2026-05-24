# analysis/

R scripts that reproduce all figures and supplementary tables in the OpenSplice paper.
Each script is self-contained: source it from RStudio or run it with `Rscript` from the
repository root.

---

## Prerequisites

All scripts source `analysis/config.R` at the top, which defines every input/output path.
Open `OpenSplice.Rproj` in RStudio so that `here()` resolves paths from the project root
automatically — no manual path editing needed.

Before running any analysis script, the master table must exist:

```r
source("analysis/00_master_table_creation.R")
```

---

## Script index

### 00 — Master table

| Script | Description |
|--------|-------------|
| `00_master_table_creation.R` | Joins PSI measurements with variant coordinates, computes per-variant statistics (delta-PSI, delta-LogitPSI, z-scores, FDR), assigns genomic regions, and writes `Supplementary_Table4.tsv` (= the master table used by all downstream scripts) and `results/psi_per_variant/coverage_per_exon.tsv`. |

### 01 — Dataset QC and overview

| Script | Description |
|--------|-------------|
| `01.1_replicates_correlation_plots.R` | Replicate correlation (GGally pairs plots) for all mutagenesis and WT screening libraries. |
| `01.2_other_dataset_comparison.R` | Comparison against external datasets: FAS INDEL vs. Baeza-Centurion et al. 2025, single-clone gel validation, and SpliceVarDB. |
| `01.3_number_variant_coverage.R` | Variant counts by mutation type and region; per-exon sequence coverage across libraries. |
| `01.4_overwiev_variant_effect.R` | Distribution of effect sizes (ΔLogitPSI bins) across mutation types and genomic regions; Fisher's exact and chi-squared tests. |

### 02 — Heatmaps

| Script | Description |
|--------|-------------|
| `02.1_heatmap_585_exons.R` | ΔPSI heatmap across all exons, ordered by WT PSI (built in two phases: matrix construction then plotting). |
| `02.2_all_heatmaps_and_wt_logo.R` | Per-exon heatmaps and WT sequence logos (Supplementary Figure 1). |

### 03 — Mutational sensitivity

| Script | Description |
|--------|-------------|
| `03.1_mutational_sensitivity_WT_PSI.R` | Effect of WT PSI on mutational sensitivity, by mutation type and region. |
| `03.2_mutational_sensitivity_exon_length.R` | Effect of exon length on mutational sensitivity; PSI loss below 30 nt. |
| `03.3_mutational_sensitivity_splice_site.R` | Effect of splice site strength (MAxEntScan) on mutational sensitivity |
| `03.4_mutational_sensitivity_region_mut_type.R` | Effect of position/WT nucleotide identity on mutational sensitivity|

### 04 — Branch points (BP)

| Script | Description |
|--------|-------------|
| `04.1_branch_point_processing_data.R` | Process external BP datasets (Mercer 2015, Taggart 2017, LaBranchor). |
| `04.2_branch_point_manual.R` | Mutational signature of manually annotated branch points. |
| `04.3_branch_point_mercer_taggart_labranchor.R` | Mutational signature of preidcted branch points, external dataset comparisons, and Supplementary Table 7. |

### 06 — Splicing cis-regulatory elements (SRE)

| Script | Description |
|--------|-------------|
| `06.1_cis_regulatory_element_mapping.R` | Map SRE positions (enhancer / silencer / overlap / neutral) from max/min ΔLogitPSI profiles. |
| `06.2_preparing_clustering_files_regulatory_elements.R` | Build per-position SRE annotation matrix and ΔLogitPSI heatmap matrices for clustering. |
| `06.3_clustering_regulatory_state.R` | Compute per-region SRE coverage, hierarchically cluster exons (Ward.D2, k = 6), assign cluster IDs. |
| `06.4_state_description_overall.R` | SRE state composition across all exons and PSI groups (all / PSI = 20–80 / PSI < 20 / PSI > 80). |
| `06.5_triangle_plot_regulatory_elements.R` | Ternary plots of per-exon SRE composition; correlation with WT PSI and exon length. |
| `06.6_transition.R` | SRE state-transition rates vs. WT PSI and exon length, per region and transition type. |
| `06.7_clustering_plot1.R` | Heatmaps of the SRE state matrix and ΔLogitPSI matrices ordered by cluster; per-cluster state-coverage bar charts. |
| `06.8_clustering_plot2.R` | Cluster description plots, state-transition analysis per cluster, pairwise tests, summary tables. |
| `06.9.1_SpliceMaps_main.R` | Splice maps for 30 curated example exons arranged in a 6-column grid. |
| `06.9.2_SpliceMaps_all.R` | Splice maps for all exons (≥ 50 % coverage), two-panel per exon (Supplementary Figure 2). |
| `06_shared.R` | Helper functions shared by `06.6` and `06.8` (`calculate_state_transitions()`). Not run directly. |

### 07 — ClinVar benchmarking

| Script | Description |
|--------|-------------|
| `07_clinvar.R` | Bar charts and ROC curves benchmarking splicing effect calls against ClinVar classifications; run for all variants and for high-confidence annotations (status ≥ 2). |

---

## Outputs

| Type | Location |
|------|----------|
| Figures (PNG) | `figures/<NN_group>/<NN.script_name>/` |
| Result tables | `results/analysis/<NN_group>/<NN.script_name>/` |
| Supplementary tables | `results/supplementary_tables/` |

---

## Shared utilities

- **`config.R`** — central path definitions (inputs, outputs, database files). Sourced by every script. Edit only if adding new database paths.
- **`06_shared.R`** — helper functions for the SRE clustering scripts; sourced explicitly by `06.6` and `06.8`.
