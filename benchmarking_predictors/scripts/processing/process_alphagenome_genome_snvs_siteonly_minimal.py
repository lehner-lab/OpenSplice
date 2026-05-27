#!/usr/bin/env python
#@title Batch process AlphaGenome SITEONLY SNV outputs into one merged variant table
"""
This script reproduces the original AlphaGenome SNV processing logic,
but removes:
    - experimental delta_psi loading
    - Pearson/Spearman correlations
    - plotting

It still:
    - collapses the two canonical-site rows per variant
      into ONE row per variant
    - retains WT and ALT canonical acceptor/donor scores
    - computes:
          delta_acceptor
          delta_donor
          mean_delta_splice
    - writes:
          alphagenome_genome_snvs_clean_curated.csv
"""

from __future__ import annotations

import os
import glob
import pandas as pd


# ============================================================
# SETTINGS
# ============================================================

ALPHAGENOME_DIR = "results/alphagenome/genome_mode/snvs/siteonly_fixedwindow"

TEST_N_EXONS = 650

OUT_MERGED_VARIANTS_CSV = os.path.join(
    ALPHAGENOME_DIR,
    "alphagenome_genome_snvs_clean_curated.csv"
)


# ============================================================
# DISCOVER FILES
# ============================================================

all_csv = sorted(
    glob.glob(os.path.join(ALPHAGENOME_DIR, "*.csv"))
)

print("CSV files found in directory:", len(all_csv))

files = [
    f for f in all_csv
    if "SITEONLY" in os.path.basename(f)
]

print("SITEONLY CSVs found:", len(files))

if len(files) == 0:
    raise FileNotFoundError(
        "No SITEONLY CSVs found in the directory."
    )

if isinstance(TEST_N_EXONS, int):
    files = files[:TEST_N_EXONS]
    print("TEST MODE: using first", len(files), "files")


# ============================================================
# COLLAPSE HELPER
# ============================================================

def collapse_alphagenome_siteonly(df_ag_raw: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse SITEONLY AlphaGenome output from two rows per variant
    into one row per variant.
    """

    needed_ag = {
        "variant_id",
        "canonical_role_by_strand",
        "donor_ref_score_at_site",
        "donor_alt_score_at_site",
        "acceptor_ref_score_at_site",
        "acceptor_alt_score_at_site",
    }

    missing_ag = needed_ag - set(df_ag_raw.columns)

    if missing_ag:
        raise ValueError(
            f"Missing AlphaGenome columns: {sorted(missing_ag)}"
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
        "donor_ref_score_at_site",
        "donor_alt_score_at_site",
        "acceptor_ref_score_at_site",
        "acceptor_alt_score_at_site",
    ]

    base_cols = [
        c for c in df.columns
        if c not in (["canonical_role_by_strand"] + score_cols)
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
            "acceptor_alt_score_at_site",
        ]]
        .rename(columns={
            "acceptor_ref_score_at_site": "acceptor_ref_at_canonical",
            "acceptor_alt_score_at_site": "acceptor_alt_at_canonical",
        })
    )

    don = (
        df[df["canonical_role_by_strand"] == "donor"]
        .set_index("variant_id")[[
            "donor_ref_score_at_site",
            "donor_alt_score_at_site",
        ]]
        .rename(columns={
            "donor_ref_score_at_site": "donor_ref_at_canonical",
            "donor_alt_score_at_site": "donor_alt_at_canonical",
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
        .mean(axis=1)
    )

    return out


# ============================================================
# PROCESS FILES
# ============================================================

all_merged_rows = []

for i, fpath in enumerate(files, start=1):

    try:
        df_ag_raw = pd.read_csv(fpath)

        df_ag_one = collapse_alphagenome_siteonly(df_ag_raw)

        keep_cols = [
            "exon_id",
            "strand",
            "variant_id",
            "variant_pos_1based",
            "ref_allele",
            "alt_allele",
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
