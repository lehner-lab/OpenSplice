#!/usr/bin/env python
#@title Combine AlphaGenome minigene SNV and deletion outputs into one processed table
from __future__ import annotations

"""
Purpose
-------
Process raw AlphaGenome minigene SNV and deletion score files together into
one final variant-level dataframe.

Fix vs previous version
-----------------------
This version only reads RAW AlphaGenome score files and excludes any already
processed outputs such as:
- *_processed.tsv
- *_alphagenome_minigene_processed.tsv
- combined summary files

What it adds
------------
For every ALT row:
- alphagenome_minigene_var_source
- alphagenome_minigene_acceptor_wt
- alphagenome_minigene_donor_wt
- alphagenome_minigene_acceptor_mut
- alphagenome_minigene_donor_mut
- alphagenome_minigene_delta_acceptor
- alphagenome_minigene_delta_donor
- alphagenome_minigene_delta_mean
- alphagenome_minigene_variant_class   ("snv" or "deletion")
- alphagenome_minigene_input_file

It also keeps all original ALT metadata columns from the raw files.
"""

import os
import glob
import gc
import pandas as pd

# ============================================================
# Paths
# ============================================================
SNV_DIR = "results/alphagenome/minigene_mode/snvs/per_exon"
DEL_DIR = "results/alphagenome/minigene_mode/deletions/per_exon"

OUTFILE = (
    "results/alphagenome/"
    "alphagenome_minigene_clean_curated.tsv"
)

# ============================================================
# Helper: keep only RAW AlphaGenome score files
# ============================================================
def get_raw_score_files(input_dir: str) -> list[str]:
    """
    Return only raw AlphaGenome score TSVs, excluding any downstream processed files.
    """
    all_tsvs = sorted(glob.glob(os.path.join(input_dir, "*.tsv")))

    raw_files = []
    for f in all_tsvs:
        base = os.path.basename(f)

        # Must be a raw score file
        if "alphagenome_scores" not in base:
            continue

        # Exclude processed / derived outputs
        exclude_tokens = [
            "_processed.tsv",
            "_alphagenome_minigene_processed.tsv",
            "_with_deltas.tsv",
            "_combined.tsv",
        ]
        if any(tok in base for tok in exclude_tokens):
            continue

        raw_files.append(f)

    return sorted(raw_files)

# ============================================================
# Discover files
# ============================================================
snv_files = get_raw_score_files(SNV_DIR)
del_files = get_raw_score_files(DEL_DIR)

print(f"Found raw SNV files: {len(snv_files)}")
print(f"Found raw deletion files: {len(del_files)}")

if len(snv_files) == 0 and len(del_files) == 0:
    raise FileNotFoundError("No raw AlphaGenome input files found.")

# ============================================================
# Per-file processor
# ============================================================
def process_one_alphagenome_file(file_path: str, variant_class: str) -> pd.DataFrame:
    """
    Process one raw AlphaGenome score file.

    Keeps all original ALT metadata columns and adds standardized
    AlphaGenome minigene WT/MUT/delta columns.
    """
    df = pd.read_csv(file_path, sep="\t")

    required_cols = [
        "Type",
        "Identifier",
        "alphagenome_acceptor_score",
        "alphagenome_donor_score",
    ]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required columns {missing} in file:\n{file_path}"
        )

    ref_df = df[df["Type"] == "REF"].copy()
    alt_df = df[df["Type"] == "ALT"].copy()

    if len(ref_df) != 1:
        raise ValueError(
            f"Expected exactly 1 REF row, found {len(ref_df)} in file:\n{file_path}"
        )

    if len(alt_df) == 0:
        raise ValueError(f"No ALT rows found in file:\n{file_path}")

    wt_row = ref_df.iloc[0]
    wt_acceptor = wt_row["alphagenome_acceptor_score"]
    wt_donor = wt_row["alphagenome_donor_score"]

    alt_df = alt_df.copy()

    # Standardized AlphaGenome minigene columns
    alt_df["alphagenome_minigene_var_source"] = alt_df["Identifier"]

    alt_df["alphagenome_minigene_acceptor_wt"] = wt_acceptor
    alt_df["alphagenome_minigene_donor_wt"] = wt_donor

    alt_df["alphagenome_minigene_acceptor_mut"] = alt_df["alphagenome_acceptor_score"]
    alt_df["alphagenome_minigene_donor_mut"] = alt_df["alphagenome_donor_score"]

    alt_df["alphagenome_minigene_delta_acceptor"] = (
        alt_df["alphagenome_minigene_acceptor_mut"]
        - alt_df["alphagenome_minigene_acceptor_wt"]
    )

    alt_df["alphagenome_minigene_delta_donor"] = (
        alt_df["alphagenome_minigene_donor_mut"]
        - alt_df["alphagenome_minigene_donor_wt"]
    )

    alt_df["alphagenome_minigene_delta_mean"] = (
        alt_df["alphagenome_minigene_delta_acceptor"]
        + alt_df["alphagenome_minigene_delta_donor"]
    ) / 2

    alt_df["alphagenome_minigene_variant_class"] = variant_class
    alt_df["alphagenome_minigene_input_file"] = os.path.basename(file_path)

    return alt_df

# ============================================================
# Process all files
# ============================================================
processed_tables = []

for i, file_path in enumerate(snv_files, start=1):
    print(f"[SNV {i}/{len(snv_files)}] {os.path.basename(file_path)}")
    tmp = process_one_alphagenome_file(file_path, variant_class="snv")
    processed_tables.append(tmp)
    del tmp
    gc.collect()

for i, file_path in enumerate(del_files, start=1):
    print(f"[DEL {i}/{len(del_files)}] {os.path.basename(file_path)}")
    tmp = process_one_alphagenome_file(file_path, variant_class="deletion")
    processed_tables.append(tmp)
    del tmp
    gc.collect()

# ============================================================
# Combine
# ============================================================
final_df = pd.concat(processed_tables, axis=0, ignore_index=True, sort=False)

# ============================================================
# Reorder important columns to front
# ============================================================
front_cols = [
    "ensembl_exon_id",
    "alphagenome_minigene_variant_class",
    "alphagenome_minigene_input_file",
    "Identifier",
    "alphagenome_minigene_var_source",

    "alphagenome_minigene_acceptor_wt",
    "alphagenome_minigene_donor_wt",
    "alphagenome_minigene_acceptor_mut",
    "alphagenome_minigene_donor_mut",
    "alphagenome_minigene_delta_acceptor",
    "alphagenome_minigene_delta_donor",
    "alphagenome_minigene_delta_mean",

    # useful shared construct metadata
    "exon_len",
    "wt_core_len",
    "target_len",
    "wt_left_pad_len",
    "wt_right_pad_len",
    "pre_manual_len",
    "post_manual_len",
    "acceptor_construct_pos0",
    "donor_construct_pos0",

    # SNV-specific metadata
    "snv_start_in_exon_1based",
    "snv_wt_base",
    "snv_alt_base",
    "snv_construct_pos0",
    "snv_padded_pos0",

    # deletion-specific metadata
    "deletion_start_in_nt_seq_1based",
    "deletion_len_bp",
    "deletion_start_construct_pos0",
]

front_cols_existing = [c for c in front_cols if c in final_df.columns]
other_cols = [c for c in final_df.columns if c not in front_cols_existing]
final_df = final_df[front_cols_existing + other_cols]

# ============================================================
# Save
# ============================================================
final_df.to_csv(OUTFILE, sep="\t", index=False)

print("\nDone.")
print("Final shape:", final_df.shape)
print("Saved combined file:")
print(OUTFILE)

print("\nPreview:")
with pd.option_context("display.max_columns", None, "display.width", 100000):
    print(final_df.head())
