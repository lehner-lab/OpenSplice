#@title Process SpliceAI minigene-mode outputs into clean benchmarking tables
"""
Process SpliceAI minigene-mode SNV and deletion outputs into final benchmarking tables.

This script starts from the per-exon SpliceAI parquet files produced by:

    spliceai_minigene_snvs_inference.py
    spliceai_minigene_deletions_inference.py

For each exon, it:
    1. loads the per-exon SpliceAI output parquet,
    2. merges the model output with variant metadata,
    3. identifies the canonical acceptor and donor score columns,
    4. computes WT, mutant, and delta SpliceAI scores,
    5. saves both a full table and a clean curated table.

Plotting code from the original notebook has been removed. The numerical
processing logic is intentionally kept close to the original manuscript-analysis
notebook.
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
RESULTS_DIR = Path("results/spliceai/minigene_mode")

# Variant-level metadata containing experimental delta_psi.
#
# libraries/exon_paper_1/opensplice_predictors_benchmarking_variant_metadata.tsv
#
# NOTE:
# Replace this placeholder filename with the final supplementary-data filename
# during public repository finalisation.
META_FILE = DATA_DIR / "opensplice_predictors_benchmarking_variant_metadata.tsv"

# Per-exon SpliceAI minigene-mode output directories.
#
# benchmarking_spliceai/snvs/CORRECT_pre_seq
#
# benchmarking_spliceai/dels/CORRECT_pre_seq
SNV_DIR = RESULTS_DIR / "snvs"
DEL_DIR = RESULTS_DIR / "deletions"

# Directory where this processing script writes final tables.
OUTPUT_DIR = RESULTS_DIR / "processed"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_FULL_CSV = OUTPUT_DIR / "df_spliceai_minigene_all_raw.csv"
OUT_CLEAN_CSV = OUTPUT_DIR / "df_spliceai_minigene_clean_curated.csv"


# ============================================================
# SETTINGS
# ============================================================

TEST = False
N_EXON_FILES = 5

# Canonical acceptor coordinate in the original minigene construct.
ACCEPTOR_POS0 = 216

# Minimum variants per exon for optional correlation summary.
MIN_VARIANTS_FOR_CORR = 2


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
n_meta_exons = meta["ensembl_exon_id"].nunique()
print(f"  Unique exons in metadata: {n_meta_exons}")


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def process_spliceai_file(pq_path: str | Path, var_source: str) -> pd.DataFrame:
    """
    Process one SpliceAI minigene-mode parquet for a single exon.

    Assumptions from the original minigene analysis:
        - canonical acceptor score column is acceptor_score_216
        - canonical donor position is 216 - 1 + exon_length

    Returns a merged DataFrame for this exon, or an empty DataFrame if required
    information is missing.
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

    acceptor_col = f"acceptor_score_{ACCEPTOR_POS0}"
    expected_donor_col = f"donor_score_{donor_pos0}"

    print(f"    WT exon_length = {exon_len}")
    print(f"    acceptor_col = {acceptor_col}")
    print(f"    expected donor_col = {expected_donor_col}")

    # Merge model output with metadata.
    m = df.merge(meta_exon, on=["ensembl_exon_id", "variant_id"], how="inner")
    if m.empty:
        print(f"    No variant overlap for exon {exon_id}; skipping.")
        return pd.DataFrame()

    if acceptor_col not in m.columns:
        print(f"    {acceptor_col} missing for exon {exon_id}; skipping.")
        return pd.DataFrame()

    donor_cols = [c for c in m.columns if c.startswith("donor_score_")]
    if not donor_cols:
        print(f"    No donor_score_* columns for exon {exon_id}; skipping.")
        return pd.DataFrame()

    if expected_donor_col in donor_cols:
        donor_col = expected_donor_col
    elif len(donor_cols) == 1:
        # Retain original fallback behaviour from the notebook.
        donor_col = donor_cols[0]
        print(
            f"    Expected {expected_donor_col} not found; using {donor_col} "
            f"for exon {exon_id}."
        )
    else:
        print(
            f"    Ambiguous donor columns for exon {exon_id}. "
            f"Expected {expected_donor_col}, found {donor_cols}. Skipping exon."
        )
        return pd.DataFrame()

    # Identify WT row in the merged table.
    wt_merged_rows = m[m["variant_id"] == wt_meta["variant_id"]]
    if wt_merged_rows.empty:
        print(
            f"    WT variant_id from metadata not present in SpliceAI output "
            f"for exon {exon_id}; skipping."
        )
        return pd.DataFrame()

    wt_row = wt_merged_rows.iloc[0]

    acc_wt = float(wt_row[acceptor_col])
    don_wt = float(wt_row[donor_col])

    # WT scores are constant for all variants in this exon.
    m["spliceai_minigene_acceptor_wt"] = acc_wt
    m["spliceai_minigene_donor_wt"] = don_wt

    # Mutant scores are row-specific.
    m["spliceai_minigene_acceptor_mut"] = m[acceptor_col]
    m["spliceai_minigene_donor_mut"] = m[donor_col]

    # Signed delta scores.
    m["spliceai_minigene_delta_acceptor"] = (
        m["spliceai_minigene_acceptor_mut"] - acc_wt
    )
    m["spliceai_minigene_delta_donor"] = (
        m["spliceai_minigene_donor_mut"] - don_wt
    )

    m["spliceai_minigene_delta_mean"] = (
        m["spliceai_minigene_delta_acceptor"]
        + m["spliceai_minigene_delta_donor"]
    ) / 2.0

    # Aliases retained for consistency with downstream benchmarking tables.
    m["spliceai_minigene_max_signed_delta_donor"] = (
        m["spliceai_minigene_delta_donor"]
    )
    m["spliceai_minigene_mean_delta_signed"] = m["spliceai_minigene_delta_mean"]

    # Provenance and coordinate columns.
    m["spliceai_minigene_var_source"] = var_source
    m["spliceai_minigene_acceptor_pos0"] = ACCEPTOR_POS0
    m["spliceai_minigene_donor_pos0"] = donor_pos0
    m["spliceai_minigene_acceptor_col"] = acceptor_col
    m["spliceai_minigene_donor_col"] = donor_col

    return m



# ============================================================
# PROCESS ALL SNV AND DELETION FILES
# ============================================================

print("\nProcessing SpliceAI minigene SNV parquet files")
snv_files = sorted(glob.glob(str(SNV_DIR / "*_snv_spliceai.parquet")))
print("  SNV files found:", len(snv_files))

if TEST:
    snv_files = snv_files[:N_EXON_FILES]
    print(f"  TEST mode: restricting SNV files to first {len(snv_files)}")

snv_chunks = []
for pq in snv_files:
    snv_chunks.append(process_spliceai_file(pq, var_source="SNV"))

snv_chunks = [df for df in snv_chunks if not df.empty]
snv_master = pd.concat(snv_chunks, ignore_index=True) if snv_chunks else pd.DataFrame()

print("SNV master shape:", snv_master.shape)
print("  SNV exons:", snv_master["ensembl_exon_id"].nunique() if not snv_master.empty else 0)


print("\nProcessing SpliceAI minigene DEL parquet files")
del_files = sorted(glob.glob(str(DEL_DIR / "*spliceai_scores_minigene.parquet")))
print("  DEL files found:", len(del_files))

if TEST:
    del_files = del_files[:N_EXON_FILES]
    print(f"  TEST mode: restricting DEL files to first {len(del_files)}")

del_chunks = []
for pq in del_files:
    del_chunks.append(process_spliceai_file(pq, var_source="DEL"))

del_chunks = [df for df in del_chunks if not df.empty]
del_master = pd.concat(del_chunks, ignore_index=True) if del_chunks else pd.DataFrame()

print("DEL master shape:", del_master.shape)
print("  DEL exons:", del_master["ensembl_exon_id"].nunique() if not del_master.empty else 0)


# ============================================================
# COMBINE SNVS AND DELETIONS
# ============================================================

frames = []
if not snv_master.empty:
    frames.append(snv_master)
if not del_master.empty:
    frames.append(del_master)

if not frames:
    raise ValueError("No SpliceAI SNV or DEL exons were successfully processed.")

spliceai_all = pd.concat(frames, ignore_index=True)
print("\nCombined SpliceAI minigene table:", spliceai_all.shape)

n_exons_spliceai = spliceai_all["ensembl_exon_id"].nunique()
print(f"  Unique exons in SpliceAI minigene table: {n_exons_spliceai}")



# ============================================================
# CLEAN FINAL TABLE
# ============================================================

keep_cols = [
    "ensembl_exon_id",
    "variant_id",
    "exon_length",
    "spliceai_minigene_var_source",
    "spliceai_minigene_acceptor_pos0",
    "spliceai_minigene_donor_pos0",
    "spliceai_minigene_acceptor_wt",
    "spliceai_minigene_donor_wt",
    "spliceai_minigene_acceptor_mut",
    "spliceai_minigene_donor_mut",
    "spliceai_minigene_delta_acceptor",
    "spliceai_minigene_delta_donor",
    "spliceai_minigene_delta_mean",
    "spliceai_minigene_mean_delta_signed",
    "spliceai_minigene_max_signed_delta_donor",
]

optional_extra_cols = [
    "Type",
    "up_5k_len",
    "deletion_start_insert_pos0",
    "deletion_length",
]
for col in optional_extra_cols:
    if col in spliceai_all.columns:
        keep_cols.append(col)

clean_cols = [c for c in keep_cols if c in spliceai_all.columns]
spliceai_clean = spliceai_all[clean_cols].copy()

print("\nFinal clean SpliceAI minigene table:", spliceai_clean.shape)
print("Clean columns:")
for col in spliceai_clean.columns:
    print("  ", col)


# ============================================================
# SAVE OUTPUT TABLES
# ============================================================

spliceai_all.to_csv(OUT_FULL_CSV, index=False)
spliceai_clean.to_csv(OUT_CLEAN_CSV, index=False)

print(f"\nSaved full SpliceAI minigene table to: {OUT_FULL_CSV}")
print(f"Saved clean SpliceAI minigene table to: {OUT_CLEAN_CSV}")
