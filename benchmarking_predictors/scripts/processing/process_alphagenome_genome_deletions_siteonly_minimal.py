#!/usr/bin/env python
#@title Batch process AlphaGenome DELETION SITEONLY outputs into one merged variant table
"""
Process AlphaGenome deletion SITEONLY_realigned output files.

This reproduces the original deletion processing logic, but removes:
    - experimental delta_psi loading
    - Pearson/Spearman correlations
    - plotting
    - per-exon correlation stats output

It still:
    - zero-fills ALT scores when a canonical splice site is deleted
    - collapses the two canonical-site rows per variant
      into ONE row per variant
    - retains WT and ALT canonical acceptor/donor scores
    - computes:
          delta_acceptor
          delta_donor
          mean_delta_splice
    - writes:
          alphagenome_genome_deletions_clean_curated.csv
"""

from __future__ import annotations

import os
import glob
import pandas as pd


# ============================================================
# SETTINGS
# ============================================================

ALPHAGENOME_DEL_DIR = (
    "results/alphagenome/"
    "exon_608_splice_sites_base_resolution_deletions_strandaware_SITEONLY_realigned_FIXEDWINDOW"
)

TEST_N_EXONS = None  # e.g. 25 for a quick test, or None for all

OUT_MERGED_VARIANTS_CSV = os.path.join(
    ALPHAGENOME_DEL_DIR,
    "alphagenome_genome_deletions_clean_curated.csv",
)

print("ALPHAGENOME_DEL_DIR:", ALPHAGENOME_DEL_DIR)


# ============================================================
# DISCOVER FILES
# ============================================================

all_csv = sorted(
    glob.glob(os.path.join(ALPHAGENOME_DEL_DIR, "*.csv"))
)

print("\nCSV files found in directory:", len(all_csv))
print("First 20 CSVs:")
for fpath in all_csv[:20]:
    print(" -", os.path.basename(fpath))

files = [
    fpath
    for fpath in all_csv
    if "SITEONLY" in os.path.basename(fpath).upper()
]

print("\nCandidate SITEONLY CSVs found:", len(files))

if len(files) == 0:
    raise FileNotFoundError(
        "No SITEONLY CSVs found in ALPHAGENOME_DEL_DIR. "
        "Check the printed filenames above and adjust the filter."
    )

if isinstance(TEST_N_EXONS, int):
    files = files[:TEST_N_EXONS]
    print("TEST MODE: using first", len(files), "files")


# ============================================================
# ZERO-FILL DELETED CANONICAL SITES
# ============================================================

def zero_fill_deleted_sites(df: pd.DataFrame) -> pd.DataFrame:
    """
    If a canonical base is deleted in ALT, define the ALT score at that
    canonical coordinate as 0.0.

    Uses these optional flags if present:
        - acceptor_site_deleted_in_alt
        - donor_site_deleted_in_alt
    """
    out = df.copy()

    if {
        "acceptor_site_deleted_in_alt",
        "acceptor_alt_score_realigned",
    }.issubset(out.columns):
        mask = out["acceptor_site_deleted_in_alt"].astype(bool)
        out.loc[mask, "acceptor_alt_score_realigned"] = 0.0

    if {
        "donor_site_deleted_in_alt",
        "donor_alt_score_realigned",
    }.issubset(out.columns):
        mask = out["donor_site_deleted_in_alt"].astype(bool)
        out.loc[mask, "donor_alt_score_realigned"] = 0.0

    return out


# ============================================================
# COLLAPSE HELPER
# ============================================================

def collapse_alphagenome_del_siteonly_realigned(
    df_ag_raw: pd.DataFrame,
) -> pd.DataFrame:
    """
    Collapse deletion SITEONLY_realigned output into one row per variant_id.

    The input has separate canonical acceptor and canonical donor rows.
    This function preserves canonical WT/ALT scores and computes deltas.
    """
    needed = {
        "variant_id",
        "canonical_role_by_strand",
        "acceptor_ref_score_at_site",
        "acceptor_alt_score_realigned",
        "donor_ref_score_at_site",
        "donor_alt_score_realigned",
    }

    missing = needed - set(df_ag_raw.columns)

    if missing:
        raise ValueError(
            f"Missing AlphaGenome DEL columns: {sorted(missing)}"
        )

    df = df_ag_raw.copy()

    df["canonical_role_by_strand"] = (
        df["canonical_role_by_strand"]
        .astype(str)
        .str.strip()
        .str.lower()
    )

    df = df[
        df["canonical_role_by_strand"].isin(["acceptor", "donor"])
    ].copy()

    score_cols = [
        "acceptor_ref_score_at_site",
        "acceptor_alt_score_realigned",
        "donor_ref_score_at_site",
        "donor_alt_score_realigned",
        "canonical_role_by_strand",
    ]

    base_cols = [
        c for c in df.columns
        if c not in score_cols
    ]

    base = (
        df.sort_values(["variant_id"])
        .drop_duplicates("variant_id", keep="first")[base_cols]
        .set_index("variant_id")
    )

    acc = (
        df[df["canonical_role_by_strand"] == "acceptor"]
        .set_index("variant_id")[[
            "acceptor_ref_score_at_site",
            "acceptor_alt_score_realigned",
        ]]
        .rename(columns={
            "acceptor_ref_score_at_site": "acceptor_ref_at_canonical",
            "acceptor_alt_score_realigned": "acceptor_alt_at_canonical",
        })
    )

    don = (
        df[df["canonical_role_by_strand"] == "donor"]
        .set_index("variant_id")[[
            "donor_ref_score_at_site",
            "donor_alt_score_realigned",
        ]]
        .rename(columns={
            "donor_ref_score_at_site": "donor_ref_at_canonical",
            "donor_alt_score_realigned": "donor_alt_at_canonical",
        })
    )

    out = (
        base
        .join(acc, how="left")
        .join(don, how="left")
        .reset_index()
    )

    out["delta_acceptor"] = (
        out["acceptor_alt_at_canonical"]
        - out["acceptor_ref_at_canonical"]
    )

    out["delta_donor"] = (
        out["donor_alt_at_canonical"]
        - out["donor_ref_at_canonical"]
    )

    out["mean_delta_splice"] = (
        out[["delta_acceptor", "delta_donor"]]
        .mean(axis=1, skipna=True)
    )

    return out


# ============================================================
# PROCESS FILES
# ============================================================

all_merged_rows = []

for i, fpath in enumerate(files, start=1):

    try:
        df_ag_raw = pd.read_csv(fpath)

        # Original processing step: set ALT score to 0 when the canonical
        # splice-site coordinate itself is deleted.
        df_ag_raw = zero_fill_deleted_sites(df_ag_raw)

        # Original processing step: collapse acceptor/donor canonical rows
        # into one row per deletion variant.
        df_ag_one = collapse_alphagenome_del_siteonly_realigned(df_ag_raw)

        keep_cols = [
            "exon_id",
            "strand",
            "variant_id",
            "variant_pos_1based",
            "ref_allele",
            "alt_allele",
            "deletion_span_start_1based",
            "deletion_span_end_1based",
            "deletion_len_bp",
            "acceptor_ref_at_canonical",
            "acceptor_alt_at_canonical",
            "donor_ref_at_canonical",
            "donor_alt_at_canonical",
            "delta_acceptor",
            "delta_donor",
            "mean_delta_splice",
        ]

        keep_cols = [
            c for c in keep_cols
            if c in df_ag_one.columns
        ]

        all_merged_rows.append(df_ag_one[keep_cols])

    except Exception as e:
        print(f"[{i}/{len(files)}] ERROR {os.path.basename(fpath)}: {e}")
        continue

    if i % 25 == 0 or i == 1 or i == len(files):
        print(f"[{i}/{len(files)}] processed")


# ============================================================
# CONCATENATE + SAVE
# ============================================================

df_all = (
    pd.concat(all_merged_rows, ignore_index=True)
    if all_merged_rows
    else pd.DataFrame()
)

print("\nCombined merged variant table shape:", df_all.shape)

df_all.to_csv(OUT_MERGED_VARIANTS_CSV, index=False)

print("\nWrote:", OUT_MERGED_VARIANTS_CSV)

print("\nColumns in merged output:")
for c in df_all.columns:
    print(" -", c)

print("\nDONE.")
