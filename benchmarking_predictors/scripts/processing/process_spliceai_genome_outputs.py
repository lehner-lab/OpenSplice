#@title Process SpliceAI genome-mode per-exon outputs into a clean canonical-score table
"""
Process SpliceAI genome-mode in silico mutagenesis outputs.

This script takes the per-exon parquet files produced by the SpliceAI genome-mode
SNV and deletion inference scripts, extracts the canonical acceptor/donor scores,
merges them with variant metadata, computes mutant-minus-WT deltas,
and saves a clean table for downstream model benchmarking.

The scientific logic is intentionally kept close to the original Colab analysis:
- canonical positions are read from each per-exon parquet file;
- the corresponding score columns are selected by position;
- each variant is compared against the WT row for the same exon;
- SNV and deletion outputs are concatenated after the same processing step.
"""

from __future__ import annotations

import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================
# INPUTS / OUTPUTS
# ============================================================

# Repository-relative placeholders.
# Replace these filenames with the final public/supplementary-data filenames.
DATA_DIR = Path("data/input")
RESULTS_DIR = Path("results/spliceai/genome_mode")

META_FILE = DATA_DIR / "opensplice_predictors_benchmarking_variant_metadata.tsv"

DEL_DIR = RESULTS_DIR / "deletions"

SNV_DIR = RESULTS_DIR / "snvs"

OUT_DIR = RESULTS_DIR / "processed"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_FULL_CSV = OUT_DIR / "df_spliceai_genome_canonical_all_raw.csv"
OUT_CLEAN_CSV = OUT_DIR / "df_spliceai_genome_canonical_clean_curated.csv"
OUT_CORR_SUMMARY = OUT_DIR / "df_spliceai_genome_canonical_per_exon_corrs.csv"

# Small test mode for checking file handling on a subset of exons.
TEST = False
N_EXON_FILES = 5


# ============================================================
# LOAD METADATA
# ============================================================

print("Loading variant metadata...")
meta = pd.read_csv(META_FILE, sep="\t")
print("  meta shape:", meta.shape)

required_meta_cols = {"ensembl_exon_id", "variant_id"}
missing = required_meta_cols - set(meta.columns)
if missing:
    raise ValueError(f"Metadata missing required columns: {missing}")

n_meta_exons = meta["ensembl_exon_id"].nunique()
print(f"  Unique exons in metadata: {n_meta_exons}")


# ============================================================
# PROCESS ONE PER-EXON SPLICEAI PARQUET
# ============================================================

def process_spliceai_genome_file(pq_path: str | os.PathLike[str], var_source: str) -> pd.DataFrame:
    """
    Process one SpliceAI genome-mode parquet file for one exon.

    Expected input columns from the inference scripts:
    - ensembl_exon_id
    - Identifier, renamed here to variant_id
    - acceptor_pos0
    - donor_pos0
    - acceptor_score_<acceptor_pos0>
    - donor_score_<donor_pos0>

    Parameters
    ----------
    pq_path
        Path to one per-exon SpliceAI genome-mode parquet file.
    var_source
        Variant source label, usually "SNV" or "DEL".

    Returns
    -------
    pandas.DataFrame
        Metadata-merged dataframe with WT, mutant, and delta canonical scores.
        Returns an empty dataframe if the file fails required checks.
    """

    pq_path = Path(pq_path)
    df = pd.read_parquet(pq_path)

    # Standardise the variant ID column from the inference output.
    if "Identifier" in df.columns:
        df = df.rename(columns={"Identifier": "variant_id"})

    if "ensembl_exon_id" not in df.columns:
        print(f"    WARNING: no ensembl_exon_id in {pq_path.name}; skipping.")
        return pd.DataFrame()

    if "variant_id" not in df.columns:
        print(f"    WARNING: no variant_id in {pq_path.name}; skipping.")
        return pd.DataFrame()

    exon_id = df["ensembl_exon_id"].iloc[0]
    print(f"  Processing exon {exon_id} from {pq_path.name}")

    meta_exon = meta[meta["ensembl_exon_id"] == exon_id].copy()
    if meta_exon.empty:
        print(f"    WARNING: no metadata rows for exon {exon_id}; skipping.")
        return pd.DataFrame()

    # Canonical genome-mode positions should already be stored in each parquet.
    required_pos_cols = {"acceptor_pos0", "donor_pos0"}
    missing_pos = required_pos_cols - set(df.columns)
    if missing_pos:
        print(f"    WARNING: missing position columns {missing_pos} for exon {exon_id}; skipping.")
        return pd.DataFrame()

    acceptor_pos_vals = df["acceptor_pos0"].dropna().unique()
    donor_pos_vals = df["donor_pos0"].dropna().unique()

    if len(acceptor_pos_vals) != 1 or len(donor_pos_vals) != 1:
        print(
            f"    WARNING: expected one acceptor_pos0 and donor_pos0 per exon; skipping {exon_id}.\n"
            f"      acceptor_pos0 values: {acceptor_pos_vals}\n"
            f"      donor_pos0 values: {donor_pos_vals}"
        )
        return pd.DataFrame()

    acceptor_pos0 = int(acceptor_pos_vals[0])
    donor_pos0 = int(donor_pos_vals[0])

    acceptor_col = f"acceptor_score_{acceptor_pos0}"
    donor_col = f"donor_score_{donor_pos0}"

    if acceptor_col not in df.columns:
        print(f"    WARNING: expected canonical acceptor column {acceptor_col} missing; skipping.")
        return pd.DataFrame()

    if donor_col not in df.columns:
        print(f"    WARNING: expected canonical donor column {donor_col} missing; skipping.")
        return pd.DataFrame()

    # Merge SpliceAI predictions with the experimental/variant metadata.
    merged = df.merge(meta_exon, on=["ensembl_exon_id", "variant_id"], how="inner")
    if merged.empty:
        print(f"    WARNING: no variant overlap for exon {exon_id}; skipping.")
        return pd.DataFrame()

    wt_rows = merged[merged["variant_id"].astype(str).str.endswith("_wt")]
    if wt_rows.empty:
        print(f"    WARNING: no WT row after merge for exon {exon_id}; skipping.")
        return pd.DataFrame()

    if len(wt_rows) > 1:
        print(f"    WARNING: multiple WT rows found for exon {exon_id}; using first one.")

    wt_row = wt_rows.iloc[0]
    acc_wt = float(wt_row[acceptor_col])
    don_wt = float(wt_row[donor_col])

    # Canonical WT, mutant, and mutant-minus-WT score columns.
    merged["spliceai_genome_acceptor_wt"] = acc_wt
    merged["spliceai_genome_donor_wt"] = don_wt
    merged["spliceai_genome_acceptor_mut"] = merged[acceptor_col]
    merged["spliceai_genome_donor_mut"] = merged[donor_col]

    merged["spliceai_genome_delta_acceptor"] = merged["spliceai_genome_acceptor_mut"] - acc_wt
    merged["spliceai_genome_delta_donor"] = merged["spliceai_genome_donor_mut"] - don_wt
    merged["spliceai_genome_delta_mean"] = (
        merged["spliceai_genome_delta_acceptor"] + merged["spliceai_genome_delta_donor"]
    ) / 2.0

    # Aliases retained to match earlier downstream processing conventions.
    merged["spliceai_genome_mean_delta_signed"] = merged["spliceai_genome_delta_mean"]
    merged["spliceai_genome_max_signed_delta_donor"] = merged["spliceai_genome_delta_donor"]

    # Provenance columns for auditing which canonical columns were used.
    merged["spliceai_genome_var_source"] = var_source
    merged["spliceai_genome_acceptor_pos0"] = acceptor_pos0
    merged["spliceai_genome_donor_pos0"] = donor_pos0
    merged["spliceai_genome_acceptor_col"] = acceptor_col
    merged["spliceai_genome_donor_col"] = donor_col

    return merged


# ============================================================
# PROCESS SNV AND DELETION FILES
# ============================================================

print("\nProcessing SpliceAI genome-mode SNV parquet files...")
snv_files = sorted(glob.glob(str(SNV_DIR / "*_snv_spliceai_genome_mode.parquet")))
print("  SNV files found:", len(snv_files))
if TEST:
    snv_files = snv_files[:N_EXON_FILES]
    print(f"  TEST mode: restricting SNV files to first {len(snv_files)}")

snv_chunks = [process_spliceai_genome_file(pq, var_source="SNV") for pq in snv_files]
snv_chunks = [df for df in snv_chunks if not df.empty]
snv_master = pd.concat(snv_chunks, ignore_index=True) if snv_chunks else pd.DataFrame()
print("SNV master shape:", snv_master.shape)


print("\nProcessing SpliceAI genome-mode deletion parquet files...")
del_files = sorted(glob.glob(str(DEL_DIR / "*_spliceai_scores_genome_mode.parquet")))
print("  DEL files found:", len(del_files))
if TEST:
    del_files = del_files[:N_EXON_FILES]
    print(f"  TEST mode: restricting DEL files to first {len(del_files)}")

del_chunks = [process_spliceai_genome_file(pq, var_source="DEL") for pq in del_files]
del_chunks = [df for df in del_chunks if not df.empty]
del_master = pd.concat(del_chunks, ignore_index=True) if del_chunks else pd.DataFrame()
print("DEL master shape:", del_master.shape)


# ============================================================
# CONCATENATE AND SAVE FULL TABLE
# ============================================================

frames = []
if not snv_master.empty:
    frames.append(snv_master)
if not del_master.empty:
    frames.append(del_master)

if not frames:
    raise ValueError("No SpliceAI genome-mode SNV or DEL exons were successfully processed.")

spliceai_all = pd.concat(frames, ignore_index=True)
print("\nCombined SpliceAI genome-mode table:", spliceai_all.shape)
print("  Unique exons:", spliceai_all["ensembl_exon_id"].nunique())


spliceai_all.to_csv(OUT_FULL_CSV, index=False)
print(f"Saved full SpliceAI genome canonical dataframe to: {OUT_FULL_CSV}")


# ============================================================
# SAVE CLEAN CURATED TABLE
# ============================================================

keep_cols = [
    "ensembl_exon_id",
    "variant_id",
    "spliceai_genome_var_source",
    "spliceai_genome_acceptor_pos0",
    "spliceai_genome_donor_pos0",
    "spliceai_genome_acceptor_wt",
    "spliceai_genome_donor_wt",
    "spliceai_genome_acceptor_mut",
    "spliceai_genome_donor_mut",
    "spliceai_genome_delta_acceptor",
    "spliceai_genome_delta_donor",
    "spliceai_genome_delta_mean",
]

clean_cols = [c for c in keep_cols if c in spliceai_all.columns]
spliceai_clean = spliceai_all[clean_cols].copy()

print("\nFinal clean SpliceAI genome canonical dataframe shape:", spliceai_clean.shape)
print("Columns in final clean dataframe:")
for col in spliceai_clean.columns:
    print("  ", col)

spliceai_clean.to_csv(OUT_CLEAN_CSV, index=False)
print(f"Saved clean curated SpliceAI genome canonical dataframe to: {OUT_CLEAN_CSV}")
