# OpenSplice

## Overview

OpenSplice quantifies the impact of >590,000 variants on the alternative splicing of 608 human
exons by massively parallel site-saturation mutagenesis in minigene constructs.
This repository contains all code to reproduce the analyses and figures in the paper.

**Interactive browser ExonExplorer:** https://results.hgi.sanger.ac.uk/OpenSplice/

---

## Repository structure

```
OpenSplice/
├── OpenSplice.Rproj                    ← open this in RStudio
│
├── libraries_design/                   ← library design scripts (run once, before sequencing)
│   ├── 01_wt_screening_libraries/      ← exon selection for the 6k WT screen
│   └── 02_mutagenesis_libraries/       ← saturation mutagenesis oligo design
│
├── psi_calculation_pipeline/           ← HPC sequencing processing pipeline
│   ├── 01_bc_var_association/          ← barcode–variant association (bash + Python)
│   ├── 02_dimsum/                      ← read count aggregation via DiMSum
│   ├── 03_psi_per_barcode/             ← PSI per barcode (bash + Python + R)
│   └── 04_psi_per_variant/             ← aggregate + error model + normalisation (R)
│
├── analysis/                           ← all R analyses and figures (see analysis/README.md)
│   ├── config.R                        ← central path hub — sourced by every script
│   ├── 00_master_table_creation.R      ← build master table from pipeline outputs
│   ├── 01.x – 07_*.R                   ← analysis scripts (see analysis/README.md)
│   └── README.md
│
├── data/
│   ├── raw/                            ← raw sequencing data 
│   ├── processed/                      ← HPC pipeline outputs (add here the PSI per barcode)
│   └── databases/                      ← external reference files (MaxEntScan, ClinVar, …)
│
├── results/
│   ├── psi_per_variant/                ← output of 'psi_calculation_pipeline/04_psi_per_variant/'
│   ├── analysis/                       ← per-script result tables
│   └── supplementary_tables/           ← Supplementary Tables (TSV)
│
├── figures/                            ← stored figures
│
└── envs/
    ├── requirements.txt                ← Python dependencies (pip / conda)
    └── session_info.R                  ← capture R package versions for reproducibility
```

---

## Quick start

### 1. Clone and open the project

```bash
git clone https://github.com/lehner-lab/OpenSplice.git
cd OpenSplice
```

Open `OpenSplice.Rproj` in RStudio. The `here` package resolves all paths automatically
from the project root — no manual path editing is needed.

### 2. Install R dependencies

```r
install.packages("renv")
renv::restore()
```

### 3. Download data

Download raw and processed data from ENA ([PRJEB111846](https://www.ebi.ac.uk/ena/browser/view/PRJEB111846)) and place them in `data/raw/` 
Download PSI per barcode tables from Figshare (https://doi.org/10.6084/m9.figshare.32337414) and place them in `data/processed/03_psi_per_barcode`

The following external database files must be downloaded separately and placed in
`data/databases/` (paths are defined in `analysis/config.R`):

| File | Source |
|------|--------|
| `SpliceVarDB/20250224_splicevardb.download.tsv` | [splicevardb.org](https://splicevardb.org) |
| `other_dms/GSE307247_Processed_PSIs_All_Cells.csv` | [GEO GSE307247](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307247) — file `GSE307247_Processed_PSIs_All_Cells.csv.gz`, decompress before use |
| `branch_point/` — Mercer et al. 2015 Supplementary Table 1 | [Genome Research](https://doi.org/10.1101/gr.182899.114) |
| `branch_point/` — Taggart et al. 2017 Supplementary Table 5 | [Genome Research](https://doi.org/10.1101/gr.202820.115) |
| `branch_point/lstm.gencode_v19.hg19.top.bed.gz` | LaBranchoR — http://bejerano.stanford.edu/labranchor/downloads/dat/lstm.gencode_v19.hg19.top.bed.gz |
| `branch_point/lstm.gencode_v19.hg19.all.tsv.gz` | LaBranchoR — http://bejerano.stanford.edu/labranchor/downloads/dat/lstm.gencode_v19.hg19.all.tsv.gz |

### 4. HPC pipeline

The `psi_calculation_pipeline/` directory contains the compute-intensive processing steps
(barcode–variant association → DiMSum → PSI per barcode → PSI per variant).
See [`psi_calculation_pipeline/README.md`](psi_calculation_pipeline/README.md) for
cluster-specific instructions.

> **Shortcut 1:** if you download the `psi_per_barcode` files from Figshare
> (https://doi.org/10.6084/m9.figshare.32337414) and place them under
> `data/processed/03_psi_per_barcode/`, you can skip steps 01–03 and start
> directly from `psi_calculation_pipeline/04_psi_per_variant/`.

> **Shortcut 2:** if you download the `Supplementary_Table4.tsv` file from the preprint
>  and place it under `results/supplementary_tables/Supplementary_Table4.tsv`,
> you can skip all the psi_calculation_pipeline steps + analysis/00_master_table_creation.R and
> start directly from `analysis/01.1_replicates_correlation_plots.R`.

### 5. Build the master table

Run once to produce `results/supplementary_tables/Supplementary_Table4.tsv`, which is
read by every downstream analysis script:

```r
source("analysis/00_master_table_creation.R")
```

### 6. Run analyses

Each numbered script in `analysis/` is self-contained. Run in order:

```r
source("analysis/01.1_replicates_correlation_plots.R")
source("analysis/01.2_other_dataset_comparison.R")
# ...
source("analysis/07_clinvar.R")
```

See [`analysis/README.md`](analysis/README.md) for a full description of every script.

---

## Data availability

- **DNA/cDNA sequencing data** — European Nucleotide Archive (ENA) under accession [PRJEB111846](https://www.ebi.ac.uk/ena/browser/view/PRJEB111846).
- **PSI values per variant** — Supplementary Table 4 of the paper.
- **Processed prediction scores** — Supplementary Table 12.
- **Barcode-level read counts, exon/variant sequences, flanking genomic regions, and unprocessed predictor scores** — Figshare: https://doi.org/10.6084/m9.figshare.32337414

---

## Citation

> Quarantani G, Clarke J, Thompson M, Sang F, Valcárcel J, Lehner B.
> OpenSplice: the impact of half a million mutations on the alternative splicing of 600 human exons.
> *bioRxiv* 2026. https://doi.org/10.64898/2026.05.22.727141
