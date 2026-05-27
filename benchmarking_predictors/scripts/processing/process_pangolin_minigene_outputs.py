#@title Process Pangolin minigene-mode outputs into clean benchmarking tables
"""
Process Pangolin minigene-mode SNV and deletion outputs into final benchmarking tables.

This script starts from the per-exon Pangolin parquet files produced by:

    pangolin_minigene_deletions_inference.py
    pangolin_minigene_snvs_inference.py

For each exon, it:
    1. loads the per-exon Pangolin output parquet,
    2. merges the model output with variant metadata,
    3. identifies the canonical acceptor and donor score columns,
    4. computes WT, mutant, and delta Pangolin scores per tissue,
    5. computes the signed max-absolute canonical Pangolin metric across tissues,
    6. saves both a full table and a clean curated table.

This version intentionally does not compute correlations or make plots, because
experimental delta_psi is not required at this processing stage.
"""
from __future__ import annotations

import glob
from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================
# INPUTS AND OUTPUTS
# ============================================================

DATA_DIR = Path("data/input")
RESULTS_DIR = Path("results/pangolin/minigene_mode")

# Variant-level metadata.
#
# Required here:
#   ensembl_exon_id, variant_id, exon_length
#
# delta_psi is not required by this processing script.
#
# libraries/exon_paper_1/opensplice_predictors_benchmarking_variant_metadata.tsv
#
# NOTE:
# Replace this placeholder filename with the final supplementary-data filename
# during public repository finalisation.
META_FILE = DATA_DIR / "opensplice_predictors_benchmarking_variant_metadata.tsv"

# Per-exon Pangolin minigene-mode output directories.
#
# benchmarking_pangolin/snvs_all/CORRECT_pre_seq
#
# CORRECT_pre_seqbenchmarking_pangolin/dels_all
PANGO_SNV_DIR = RESULTS_DIR / "snvs"
PANGO_DEL_DIR = RESULTS_DIR / "deletions"

# Directory where this processing script writes final tables.
OUTPUT_DIR = RESULTS_DIR / "processed"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_FULL_CSV = OUTPUT_DIR / "df_pangolin_minigene_with_all_cols.csv"
OUT_CLEAN_CSV = OUTPUT_DIR / "df_pangolin_minigene_clean_curated.csv"


# ============================================================
# SETTINGS
# ============================================================

# Tissues encoded in Pangolin output.
TISSUES = ["brain", "heart", "liver", "testis"]

# Test mode can be useful for checking the script on a few exon files.
TEST = False
N_EXON_FILES = 5

# Canonical acceptor coordinate in the original minigene construct.
ACCEPTOR_POS0 = 216

# File patterns from the per-exon Pangolin inference scripts.
SNV_GLOB = "*_pangolin_scores_minigene_snvs.parquet"
DEL_GLOB = "*_pangolin_scores_minigene.parquet"


# ============================================================
# LOAD METADATA
# ============================================================

print("Loading metadata")
meta = pd.read_csv(META_FILE, sep="\t")
print("  meta shape:", meta.shape)

required_meta_cols = {"ensembl_exon_id", "variant_id", "exon_length"}
missing = required_meta_cols - set(meta.columns)
if missing:
    raise ValueError(f"Metadata missing required columns: {missing}")

meta["variant_id"] = meta["variant_id"].astype(str)


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def signed_max_abs(row: pd.Series) -> float:
    """
    Return the value with the largest absolute magnitude, keeping its sign.

    This is used to combine tissue-specific Pangolin deltas into a single
    signed max-absolute acceptor or donor metric.
    """
    vals = row.values.astype(float)
    if np.all(np.isnan(vals)):
        return np.nan

    idx = np.nanargmax(np.abs(vals))
    return vals[idx]


def process_pangolin_file(pq_path: str | Path, var_source: str) -> pd.DataFrame:
    """
    Process one Pangolin minigene-mode parquet for a single exon.

    For each tissue, this creates:
        pangolin_minigene_{tissue}_acceptor_wt
        pangolin_minigene_{tissue}_donor_wt
        pangolin_minigene_{tissue}_acceptor_mut
        pangolin_minigene_{tissue}_donor_mut
        pangolin_minigene_{tissue}_delta_acceptor
        pangolin_minigene_{tissue}_delta_donor
        pangolin_minigene_{tissue}_delta_mean

    The donor coordinate is derived from WT exon length:
        donor_pos0 = 216 - 1 + exon_length
    """
    pq_path = Path(pq_path)
    df = pd.read_parquet(pq_path).copy()

    if "Identifier" in df.columns:
        df = df.rename(columns={"Identifier": "variant_id"})

    if "ensembl_exon_id" not in df.columns:
        print(f"    No 'ensembl_exon_id' in {pq_path.name}; skipping.")
        return pd.DataFrame()

    if "variant_id" not in df.columns:
        print(f"    No 'variant_id' in {pq_path.name}; skipping.")
        return pd.DataFrame()

    df["variant_id"] = df["variant_id"].astype(str)
    exon_id = df["ensembl_exon_id"].iloc[0]

    print(f"  Processing exon {exon_id} from {pq_path.name}")

    # Subset metadata for this exon.
    meta_exon = meta[meta["ensembl_exon_id"] == exon_id].copy()
    if meta_exon.empty:
        print(f"    No metadata rows for exon {exon_id}; skipping.")
        return pd.DataFrame()

    # Use WT metadata to define the expected donor coordinate.
    wt_meta_rows = meta_exon[meta_exon["variant_id"].str.endswith("_wt", na=False)]
    if wt_meta_rows.empty:
        print(f"    No WT metadata row for exon {exon_id}; skipping.")
        return pd.DataFrame()

    wt_meta = wt_meta_rows.iloc[0]
    exon_len = int(wt_meta["exon_length"])
    donor_pos0 = ACCEPTOR_POS0 - 1 + exon_len

    print(f"    WT exon_length = {exon_len}")
    print(f"    donor_pos0 = {donor_pos0}")

    # Merge Pangolin output with metadata.
    m = df.merge(meta_exon, on=["ensembl_exon_id", "variant_id"], how="inner")
    if m.empty:
        print(f"    No variant_id overlap for exon {exon_id}; skipping.")
        return pd.DataFrame()

    # Ensure the WT row is present after merging.
    wt_merged_rows = m[m["variant_id"] == wt_meta["variant_id"]]
    if wt_merged_rows.empty:
        print(f"    WT variant_id from metadata not found in Pangolin output for {exon_id}; skipping.")
        return pd.DataFrame()

    wt_row = wt_merged_rows.iloc[0]

    # Compute WT, mutant, and delta columns tissue-by-tissue.
    for tissue in TISSUES:
        tissue_ps_cols = [c for c in df.columns if c.startswith(f"{tissue}_ps_score_")]
        if not tissue_ps_cols:
            print(f"    No PS-score columns for tissue {tissue} in exon {exon_id}; skipping tissue.")
            continue

        acc_col = f"{tissue}_ps_score_{ACCEPTOR_POS0}"
        if acc_col not in tissue_ps_cols:
            print(f"    Missing {acc_col} for exon {exon_id}; skipping tissue {tissue}.")
            continue

        expected_donor_col = f"{tissue}_ps_score_{donor_pos0}"
        if expected_donor_col in tissue_ps_cols:
            don_col = expected_donor_col
        else:
            # Retain original fallback behaviour from the notebook:
            # if there is exactly one non-acceptor PS column, use it.
            donor_candidates = [c for c in tissue_ps_cols if c != acc_col]
            if len(donor_candidates) == 1:
                don_col = donor_candidates[0]
                print(
                    f"    Expected {expected_donor_col} not found for exon {exon_id}, "
                    f"tissue {tissue}; using {don_col}."
                )
            else:
                print(
                    f"    Ambiguous or missing donor column for exon {exon_id}, "
                    f"tissue {tissue}; candidates: {donor_candidates}. Skipping tissue."
                )
                continue

        if acc_col not in m.columns or don_col not in m.columns:
            print(
                f"    Columns {acc_col} or {don_col} missing after merge for exon "
                f"{exon_id}, tissue {tissue}; skipping tissue."
            )
            continue

        acc_wt_val = wt_row[acc_col]
        don_wt_val = wt_row[don_col]

        m[f"pangolin_minigene_{tissue}_acceptor_wt"] = acc_wt_val
        m[f"pangolin_minigene_{tissue}_donor_wt"] = don_wt_val

        m[f"pangolin_minigene_{tissue}_acceptor_mut"] = m[acc_col]
        m[f"pangolin_minigene_{tissue}_donor_mut"] = m[don_col]

        m[f"pangolin_minigene_{tissue}_delta_acceptor"] = (
            m[f"pangolin_minigene_{tissue}_acceptor_mut"] - acc_wt_val
        )
        m[f"pangolin_minigene_{tissue}_delta_donor"] = (
            m[f"pangolin_minigene_{tissue}_donor_mut"] - don_wt_val
        )

        m[f"pangolin_minigene_{tissue}_delta_mean"] = (
            m[f"pangolin_minigene_{tissue}_delta_acceptor"]
            + m[f"pangolin_minigene_{tissue}_delta_donor"]
        ) / 2.0

    # Provenance and coordinate columns.
    m["pangolin_minigene_var_source"] = var_source
    m["pangolin_minigene_acceptor_pos0"] = ACCEPTOR_POS0
    m["pangolin_minigene_donor_pos0"] = donor_pos0

    return m


# ============================================================
# PROCESS ALL SNV AND DELETION FILES
# ============================================================

print("\nProcessing Pangolin minigene SNV files")
snv_files = sorted(glob.glob(str(PANGO_SNV_DIR / SNV_GLOB)))
print("  SNV files found:", len(snv_files))

if TEST:
    snv_files = snv_files[:N_EXON_FILES]
    print(f"  TEST mode: restricting SNV files to first {len(snv_files)}")

snv_chunks = []
for pq in snv_files:
    snv_chunks.append(process_pangolin_file(pq, var_source="SNV"))

snv_chunks = [df for df in snv_chunks if not df.empty]
snv_master = pd.concat(snv_chunks, ignore_index=True) if snv_chunks else pd.DataFrame()
print("SNV master shape:", snv_master.shape)


print("\nProcessing Pangolin minigene DEL files")
del_files = sorted(glob.glob(str(PANGO_DEL_DIR / DEL_GLOB)))
print("  DEL files found:", len(del_files))

if TEST:
    del_files = del_files[:N_EXON_FILES]
    print(f"  TEST mode: restricting DEL files to first {len(del_files)}")

del_chunks = []
for pq in del_files:
    del_chunks.append(process_pangolin_file(pq, var_source="DEL"))

del_chunks = [df for df in del_chunks if not df.empty]
del_master = pd.concat(del_chunks, ignore_index=True) if del_chunks else pd.DataFrame()
print("DEL master shape:", del_master.shape)


# ============================================================
# COMBINE SNVS AND DELETIONS
# ============================================================

frames = []
if not snv_master.empty:
    frames.append(snv_master)
if not del_master.empty:
    frames.append(del_master)

if not frames:
    raise ValueError("No SNV or DEL exons were successfully processed.")

pangolin_all = pd.concat(frames, ignore_index=True)
print("\nCombined Pangolin minigene table:", pangolin_all.shape)


# ============================================================
# CANONICAL MAX-ABS SIGNED PANGOLIN METRIC
# ============================================================

acc_delta_cols = [
    c
    for c in pangolin_all.columns
    if c.startswith("pangolin_minigene_") and c.endswith("_delta_acceptor")
]
don_delta_cols = [
    c
    for c in pangolin_all.columns
    if c.startswith("pangolin_minigene_") and c.endswith("_delta_donor")
]

if not acc_delta_cols:
    raise ValueError("No Pangolin acceptor delta columns were created.")
if not don_delta_cols:
    raise ValueError("No Pangolin donor delta columns were created.")

print("\nAcceptor delta columns:", acc_delta_cols)
print("Donor delta columns:", don_delta_cols)

pangolin_all["pangolin_minigene_max_signed_delta_acceptor"] = (
    pangolin_all[acc_delta_cols].apply(signed_max_abs, axis=1)
)

pangolin_all["pangolin_minigene_max_signed_delta_donor"] = (
    pangolin_all[don_delta_cols].apply(signed_max_abs, axis=1)
)

pangolin_all["pangolin_minigene_mean_delta_signed"] = (
    pangolin_all["pangolin_minigene_max_signed_delta_acceptor"]
    + pangolin_all["pangolin_minigene_max_signed_delta_donor"]
) / 2.0

print("\nCanonical Pangolin minigene signed mean-delta summary:")
print(pangolin_all["pangolin_minigene_mean_delta_signed"].describe())


# ============================================================
# CLEAN FINAL TABLE
# ============================================================

keep_cols = [
    "ensembl_exon_id",
    "variant_id",
    "exon_length",
    "Type",
    "pangolin_minigene_var_source",
    "pangolin_minigene_acceptor_pos0",
    "pangolin_minigene_donor_pos0",
    "pangolin_minigene_max_signed_delta_acceptor",
    "pangolin_minigene_max_signed_delta_donor",
    "pangolin_minigene_mean_delta_signed",
]

# Keep deletion metadata if present.
optional_extra_cols = [
    "up_5k_len",
    "deletion_start_insert_pos0",
    "deletion_length",
]
for col in optional_extra_cols:
    if col in pangolin_all.columns:
        keep_cols.append(col)

for tissue in TISSUES:
    keep_cols.extend(
        [
            f"pangolin_minigene_{tissue}_acceptor_wt",
            f"pangolin_minigene_{tissue}_donor_wt",
            f"pangolin_minigene_{tissue}_acceptor_mut",
            f"pangolin_minigene_{tissue}_donor_mut",
            f"pangolin_minigene_{tissue}_delta_acceptor",
            f"pangolin_minigene_{tissue}_delta_donor",
            f"pangolin_minigene_{tissue}_delta_mean",
        ]
    )

clean_cols = [c for c in keep_cols if c in pangolin_all.columns]
pangolin_clean = pangolin_all[clean_cols].copy()

print("\nFinal clean Pangolin minigene table:", pangolin_clean.shape)
print("Clean columns:")
for col in pangolin_clean.columns:
    print("  ", col)


# ============================================================
# SAVE OUTPUT TABLES
# ============================================================

pangolin_all.to_csv(OUT_FULL_CSV, index=False)
pangolin_clean.to_csv(OUT_CLEAN_CSV, index=False)

print(f"\nSaved full Pangolin minigene table to: {OUT_FULL_CSV}")
print(f"Saved clean Pangolin minigene table to: {OUT_CLEAN_CSV}")
