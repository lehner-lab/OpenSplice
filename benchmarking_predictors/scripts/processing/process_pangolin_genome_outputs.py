#@title Process Pangolin genome-mode outputs into clean benchmarking tables
"""
Process Pangolin genome-mode SNV and deletion outputs into final benchmarking tables.

This script starts from the per-exon Pangolin parquet files produced by:

    pangolin_genome_deletions_inference.py
    pangolin_genome_snvs_inference.py

For each exon, it:
    1. loads the per-exon Pangolin output parquet,
    2. identifies the canonical acceptor and donor score columns,
    3. merges the model output with variant metadata,
    4. computes WT, mutant, and delta Pangolin scores per tissue,
    5. computes the signed max-absolute canonical Pangolin metric across tissues,
    6. saves both a full table and a clean curated table.

Plotting and exploratory debugging code from the original notebook have been
removed. The numerical processing logic is intentionally kept close to the
original manuscript-analysis notebook.
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
RESULTS_DIR = Path("results/pangolin/genome_mode")

# Variant-level metadata containing experimental delta_psi.
#
# libraries/exon_paper_1/opensplice_predictors_benchmarking_variant_metadata.tsv
#
# NOTE:
# Replace this placeholder filename with the final supplementary-data filename
# during public repository finalisation.
META_FILE = DATA_DIR / "opensplice_predictors_benchmarking_variant_metadata.tsv"

# Per-exon Pangolin genome-mode output directories.
#
# benchmarking_pangolin/snvs/genome_mode/
#
# benchmarking_pangolin/dels/genome_mode_reverse_exon_order/
PANGO_SNV_DIR = RESULTS_DIR / "snvs"
PANGO_DEL_DIR = RESULTS_DIR / "deletions"

# Directory where this processing script writes final tables.
OUTPUT_DIR = RESULTS_DIR / "processed"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_FULL_CSV = OUTPUT_DIR / "df_pangolin_genome_with_all_cols.csv"
OUT_CLEAN_CSV = OUTPUT_DIR / "df_pangolin_genome_clean_curated.csv"


# ============================================================
# SETTINGS
# ============================================================

# Tissues encoded in Pangolin output.
TISSUES = ["brain", "heart", "liver", "testis"]

# Test mode can be useful for checking the script on a few exon files.
TEST = False
N_EXON_FILES = 5

# File patterns from the per-exon Pangolin inference scripts.
SNV_GLOB = "*_pangolin_scores_genome_mode.parquet"
DEL_GLOB = "*_pangolin_scores_genome_mode.parquet"

# Minimum variants per exon used for the optional correlation-summary table.
MIN_VARIANTS_FOR_CORR = 4


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

# Make variant_id string once up front for safer matching.
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


def get_ref_row(df: pd.DataFrame) -> pd.Series | None:
    """
    Get the REF/WT row from a Pangolin genome-mode parquet-derived table.

    Preferred logic:
        1. Type == REF
        2. variant_id ends with _wt
    """
    if "Type" in df.columns:
        ref_rows = df[df["Type"].astype(str) == "REF"].copy()
        if not ref_rows.empty:
            return ref_rows.iloc[0]

    if "variant_id" in df.columns:
        wt_rows = df[df["variant_id"].astype(str).str.endswith("_wt", na=False)].copy()
        if not wt_rows.empty:
            return wt_rows.iloc[0]

    return None


def get_acceptor_donor_positions(
    df: pd.DataFrame,
) -> tuple[int, int] | tuple[None, None]:
    """
    Read canonical acceptor/donor coordinates directly from a genome-mode file.

    The inference scripts explicitly write:
        acceptor_pos0
        donor_pos0

    These should be constant within a single exon file.
    """
    if "acceptor_pos0" not in df.columns or "donor_pos0" not in df.columns:
        return None, None

    acc_vals = df["acceptor_pos0"].dropna().unique()
    don_vals = df["donor_pos0"].dropna().unique()

    if len(acc_vals) == 0 or len(don_vals) == 0:
        return None, None

    if len(acc_vals) > 1:
        print(f"    Warning: multiple acceptor_pos0 values found: {acc_vals}. Using first.")
    if len(don_vals) > 1:
        print(f"    Warning: multiple donor_pos0 values found: {don_vals}. Using first.")

    return int(acc_vals[0]), int(don_vals[0])


def process_pangolin_file(pq_path: str | Path, var_source: str) -> pd.DataFrame:
    """
    Process one Pangolin genome-mode parquet for a single exon.

    Expected parquet structure:
        - one REF row,
        - many ALT rows,
        - explicit acceptor_pos0 and donor_pos0 columns,
        - score columns such as heart_ps_score_<acceptor_pos0>,
        - optional deletion_start_insert_pos0 and deletion_length for deletion files.

    Returns an empty DataFrame if required information is missing.
    """
    pq_path = Path(pq_path)
    df = pd.read_parquet(pq_path).copy()

    if "Identifier" in df.columns:
        df = df.rename(columns={"Identifier": "variant_id"})

    if "ensembl_exon_id" not in df.columns:
        print(f"    No 'ensembl_exon_id' in {pq_path.name}; skipping.")
        return pd.DataFrame()

    if "variant_id" not in df.columns:
        print(f"    No 'variant_id' after renaming in {pq_path.name}; skipping.")
        return pd.DataFrame()

    df["variant_id"] = df["variant_id"].astype(str)
    exon_id = str(df["ensembl_exon_id"].iloc[0])

    print(f"  Processing exon {exon_id} from {pq_path.name}")

    # Subset metadata for this exon.
    meta_exon = meta[meta["ensembl_exon_id"].astype(str) == exon_id].copy()
    if meta_exon.empty:
        print(f"    No metadata rows for exon {exon_id}; skipping.")
        return pd.DataFrame()

    # Read canonical genome-mode positions from the parquet file.
    acceptor_pos, donor_pos = get_acceptor_donor_positions(df)
    if acceptor_pos is None or donor_pos is None:
        print(f"    Could not read acceptor_pos0/donor_pos0 for exon {exon_id}; skipping.")
        return pd.DataFrame()

    # Merge Pangolin output with experimental metadata.
    m = df.merge(meta_exon, on=["ensembl_exon_id", "variant_id"], how="inner")
    if m.empty:
        print(f"    No variant_id overlap for exon {exon_id}; skipping.")
        return pd.DataFrame()

    # Find WT/REF row after merging.
    ref_row = get_ref_row(m)
    if ref_row is None:
        print(f"    No REF/WT row found after merge for exon {exon_id}; skipping.")
        return pd.DataFrame()

    # Provenance and canonical coordinate columns.
    m["pangolin_genome_var_source"] = var_source
    m["pangolin_genome_acceptor_pos0"] = acceptor_pos
    m["pangolin_genome_donor_pos0"] = donor_pos

    # Compute WT, mutant, and delta scores tissue-by-tissue.
    for tissue in TISSUES:
        acc_col = f"{tissue}_ps_score_{acceptor_pos}"
        don_col = f"{tissue}_ps_score_{donor_pos}"

        if acc_col not in m.columns:
            print(f"    Missing {acc_col} for exon {exon_id}, tissue {tissue}; skipping tissue.")
            continue

        if don_col not in m.columns:
            print(f"    Missing {don_col} for exon {exon_id}, tissue {tissue}; skipping tissue.")
            continue

        acc_wt_val = ref_row[acc_col]
        don_wt_val = ref_row[don_col]

        m[f"pangolin_genome_{tissue}_acceptor_wt"] = acc_wt_val
        m[f"pangolin_genome_{tissue}_donor_wt"] = don_wt_val

        m[f"pangolin_genome_{tissue}_acceptor_mut"] = m[acc_col]
        m[f"pangolin_genome_{tissue}_donor_mut"] = m[don_col]

        m[f"pangolin_genome_{tissue}_delta_acceptor"] = (
            m[f"pangolin_genome_{tissue}_acceptor_mut"] - acc_wt_val
        )
        m[f"pangolin_genome_{tissue}_delta_donor"] = (
            m[f"pangolin_genome_{tissue}_donor_mut"] - don_wt_val
        )

        m[f"pangolin_genome_{tissue}_delta_mean"] = (
            m[f"pangolin_genome_{tissue}_delta_acceptor"]
            + m[f"pangolin_genome_{tissue}_delta_donor"]
        ) / 2.0

    return m



# ============================================================
# PROCESS ALL SNV AND DELETION FILES
# ============================================================

print("\nProcessing SNV Pangolin genome-mode files")
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


print("\nProcessing DEL Pangolin genome-mode files")
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
print("\nCombined Pangolin genome-mode table:", pangolin_all.shape)



# ============================================================
# CANONICAL MAX-ABS SIGNED PANGOLIN METRIC
# ============================================================

acc_delta_cols = [
    c
    for c in pangolin_all.columns
    if c.startswith("pangolin_genome_") and c.endswith("_delta_acceptor")
]
don_delta_cols = [
    c
    for c in pangolin_all.columns
    if c.startswith("pangolin_genome_") and c.endswith("_delta_donor")
]

if not acc_delta_cols:
    raise ValueError("No Pangolin acceptor delta columns were created.")
if not don_delta_cols:
    raise ValueError("No Pangolin donor delta columns were created.")

print("\nAcceptor delta columns:", acc_delta_cols)
print("Donor delta columns:", don_delta_cols)

pangolin_all["pangolin_genome_max_signed_delta_acceptor"] = (
    pangolin_all[acc_delta_cols].apply(signed_max_abs, axis=1)
)

pangolin_all["pangolin_genome_max_signed_delta_donor"] = (
    pangolin_all[don_delta_cols].apply(signed_max_abs, axis=1)
)

pangolin_all["pangolin_genome_mean_delta_signed"] = (
    pangolin_all["pangolin_genome_max_signed_delta_acceptor"]
    + pangolin_all["pangolin_genome_max_signed_delta_donor"]
) / 2.0

print("\nCanonical Pangolin genome-mode signed mean-delta summary:")
print(pangolin_all["pangolin_genome_mean_delta_signed"].describe())


# ============================================================
# CLEAN FINAL TABLE
# ============================================================

keep_cols = [
    "ensembl_exon_id",
    "variant_id",
    "exon_length",
    "Type",
    "pangolin_genome_var_source",
    "pangolin_genome_acceptor_pos0",
    "pangolin_genome_donor_pos0",
    "pangolin_genome_max_signed_delta_acceptor",
    "pangolin_genome_max_signed_delta_donor",
    "pangolin_genome_mean_delta_signed",
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
            f"pangolin_genome_{tissue}_acceptor_wt",
            f"pangolin_genome_{tissue}_donor_wt",
            f"pangolin_genome_{tissue}_acceptor_mut",
            f"pangolin_genome_{tissue}_donor_mut",
            f"pangolin_genome_{tissue}_delta_acceptor",
            f"pangolin_genome_{tissue}_delta_donor",
            f"pangolin_genome_{tissue}_delta_mean",
        ]
    )

clean_cols = [c for c in keep_cols if c in pangolin_all.columns]
pangolin_clean = pangolin_all[clean_cols].copy()

print("\nFinal clean Pangolin genome-mode table:", pangolin_clean.shape)
print("Clean columns:")
for col in pangolin_clean.columns:
    print("  ", col)


# ============================================================
# SAVE OUTPUT TABLES
# ============================================================

pangolin_all.to_csv(OUT_FULL_CSV, index=False)
pangolin_clean.to_csv(OUT_CLEAN_CSV, index=False)

print(f"\nSaved full Pangolin genome-mode table to: {OUT_FULL_CSV}")
print(f"Saved clean Pangolin genome-mode table to: {OUT_CLEAN_CSV}")
