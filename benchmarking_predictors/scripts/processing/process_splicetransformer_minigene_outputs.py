"""
#@title Process SpliceTransformer minigene outputs (SNVs + deletions)

Repository-ready processing script for SpliceTransformer minigene-mode
mutagenesis outputs.

This script:
    - Loads per-exon SpliceTransformer parquet outputs
    - Loads combined WT predictions
    - Applies splice-site deletion coordinate fixes
    - Computes WT vs mutant deltas
    - Preserves full audit columns
    - Exports cleaned SNV and deletion tables

Key correction
--------------
For deletions:
    If a canonical splice site is deleted,
    the mutant splice-site score is forced to 0.

Delta is always:
    mutant_score - WT_score

Expected inputs
---------------
1. Experimental metadata table
2. Exon metadata table
3. Per-exon SpliceTransformer parquet outputs
4. Combined WT parquet file

Outputs
-------
- splicetransformer_minigene_deletions_clean_curated.tsv
- splicetransformer_minigene_snvs_clean_curated.tsv
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

# ============================================================
# SETTINGS
# ============================================================

TEST_MODE = False

# Canonical acceptor position
ACC_POS = 71

# ============================================================
# INPUT FILES
# ============================================================

DATA_DIR = Path("data/input")

EXON_META_FILE = (
    DATA_DIR /
    "opensplice_predictors_benchmarking_exon_metadata.tsv"
)

EXP_VAR_FILE = (
    DATA_DIR /
    "opensplice_predictors_benchmarking_variant_metadata.tsv"
)

# Example original paths:
#
# EXON_META_FILE = (
#     "gioia_exon_sat_mut/libraries/exon_paper_1/"
#     "opensplice_predictors_benchmarking_exon_metadata.tsv"
# )
#
# EXP_VAR_FILE = (
#     "gioia_exon_sat_mut/libraries/exon_paper_1/"
#     "opensplice_predictors_benchmarking_variant_metadata.tsv"
# )

# ============================================================
# SPLICE TRANSFORMER OUTPUTS
# ============================================================

ST_DIR = Path(
    "results/splice_transformer/minigene_mode/inference"
)
# Set ST_DIR to wherever you saved per-exon SpliceTransformer minigene parquet outputs.

WT_COMBINED_FILE = (
    ST_DIR /
    "WT" /
    "all_WT_minigene_splice_transformer_minigene.parquet"
)

# ============================================================
# OUTPUT DIRECTORY
# ============================================================

OUTPUT_DIR = Path(
    "results/splice_transformer/minigene_mode/processed"
)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_DEL = (
    OUTPUT_DIR /
    "splicetransformer_minigene_deletions_clean_curated.tsv"
)

OUT_SNV = (
    OUTPUT_DIR /
    "splicetransformer_minigene_snvs_clean_curated.tsv"
)

# ============================================================
# LOAD INPUT TABLES
# ============================================================

print("=" * 80)
print("Loading input tables")
print("=" * 80)

df_exons = pd.read_csv(
    EXON_META_FILE,
    sep="\t",
)

df_var = pd.read_csv(
    EXP_VAR_FILE,
    sep="\t",
)

wt_all = pd.read_parquet(
    WT_COMBINED_FILE
)

print("Exon metadata:", df_exons.shape)
print("Variant metadata:", df_var.shape)
print("Combined WT predictions:", wt_all.shape)

# ============================================================
# EXON LENGTH LOOKUP
# ============================================================

exon_lengths = (
    df_exons
    .set_index("ensembl_exon_id")["exon_length"]
    .to_dict()
)

# ============================================================
# FIND PARQUET FILES
# ============================================================

st_files = sorted([
    f for f in os.listdir(ST_DIR)
    if (
        f.endswith(".parquet")
        and "all_WT_minigene" not in f
    )
])

if TEST_MODE:
    st_files = st_files[:5]

print("Per-exon parquet files:", len(st_files))

# ============================================================
# OUTPUT CONTAINERS
# ============================================================

all_del_results = []
all_snv_results = []

# ============================================================
# MAIN LOOP
# ============================================================

for file in tqdm(
    st_files,
    desc="Processing exons",
):

    exon_id = file.replace(
        "_splice_transformer_minigene.parquet",
        "",
    )

    st = pd.read_parquet(
        os.path.join(ST_DIR, file)
    )

    # --------------------------------------------------------
    # Exon length
    # --------------------------------------------------------

    exon_len = exon_lengths.get(exon_id)

    if exon_len is None or pd.isna(exon_len):

        print(
            f"Missing exon length for {exon_id}"
        )

        continue

    # --------------------------------------------------------
    # Canonical donor position
    # --------------------------------------------------------

    DON_POS = ACC_POS + int(exon_len) - 1

    # ========================================================
    # WT SCORES
    # ========================================================

    wt_exon = wt_all[
        (
            wt_all["ensembl_exon_id"]
            .astype(str)
            .eq(str(exon_id))
        )
        &
        (
            wt_all["variant"]
            .astype(str)
            .str.contains(
                "WT",
                case=False,
                na=False,
            )
        )
    ].copy()

    wt_acc = wt_exon.loc[
        wt_exon["nt_position"] == ACC_POS,
        "Acceptor",
    ].max()

    wt_don = wt_exon.loc[
        wt_exon["nt_position"] == DON_POS,
        "Donor",
    ].max()

    if pd.isna(wt_acc) or pd.isna(wt_don):

        print(
            f"Missing WT splice-site scores for "
            f"{exon_id}"
        )

        continue

    # ========================================================
    # REMOVE WT ROWS FROM PARQUET
    # ========================================================

    st = st[
        ~st["variant"]
        .astype(str)
        .str.contains(
            "WT",
            case=False,
            na=False,
        )
    ].copy()

    # ========================================================
    # VARIANT METADATA
    # ========================================================

    var_exon = df_var[
        df_var["ensembl_exon_id"] == exon_id
    ].copy()

    if var_exon.empty:
        continue

    # ========================================================
    # DELETIONS
    # ========================================================

    exon_del = var_exon[
        var_exon["variant_id"]
        .astype(str)
        .str.contains(
            "del",
            case=False,
            na=False,
        )
    ].copy()

    if not exon_del.empty:

        # ----------------------------------------------------
        # Deletion coordinates
        # ----------------------------------------------------

        exon_del["del_start"] = (
            exon_del["start"].astype(int)
        )

        exon_del["del_end"] = (
            exon_del["end"].astype(int)
        )

        exon_del["del_len"] = (
            exon_del["length"].astype(int)
        )

        # ----------------------------------------------------
        # Preserve old audit columns
        # ----------------------------------------------------

        exon_del["del_start_global"] = (
            exon_del["del_start"]
        )

        exon_del["del_end_global"] = (
            exon_del["del_end"]
        )

        # ----------------------------------------------------
        # Determine whether canonical sites deleted
        # ----------------------------------------------------

        exon_del["acc_deleted"] = (
            (exon_del["del_start_global"] <= ACC_POS)
            &
            (exon_del["del_end_global"] >= ACC_POS)
        )

        exon_del["don_deleted"] = (
            (exon_del["del_start_global"] <= DON_POS)
            &
            (exon_del["del_end_global"] >= DON_POS)
        )

        # ----------------------------------------------------
        # Coordinate shifting
        # ----------------------------------------------------

        exon_del["shifted_acc"] = exon_del.apply(
            lambda r:
                ACC_POS - r["del_len"]
                if (
                    r["del_end_global"] < ACC_POS
                    and not r["acc_deleted"]
                )
                else ACC_POS,
            axis=1,
        )

        exon_del["shifted_don"] = exon_del.apply(
            lambda r:
                DON_POS - r["del_len"]
                if (
                    r["del_end_global"] < DON_POS
                    and not r["don_deleted"]
                )
                else DON_POS,
            axis=1,
        )

        # ----------------------------------------------------
        # Extract acceptor/donor scores
        # ----------------------------------------------------

        acc_df = (
            st[
                ["variant", "nt_position", "Acceptor"]
            ]
            .rename(
                columns={"nt_position": "pos"}
            )
        )

        don_df = (
            st[
                ["variant", "nt_position", "Donor"]
            ]
            .rename(
                columns={"nt_position": "pos"}
            )
        )

        exon_del = exon_del.merge(
            acc_df,
            how="left",
            left_on=["id", "shifted_acc"],
            right_on=["variant", "pos"],
        ).rename(
            columns={"Acceptor": "acc_score"}
        )

        exon_del = exon_del.merge(
            don_df,
            how="left",
            left_on=["id", "shifted_don"],
            right_on=["variant", "pos"],
        ).rename(
            columns={"Donor": "don_score"}
        )

        # ----------------------------------------------------
        # Corrected deletion logic
        # ----------------------------------------------------

        exon_del[
            "splice_transformer_minigene_acceptor_wt"
        ] = wt_acc

        exon_del[
            "splice_transformer_minigene_donor_wt"
        ] = wt_don

        exon_del[
            "splice_transformer_minigene_acceptor_mut"
        ] = np.where(
            exon_del["acc_deleted"],
            0.0,
            exon_del["acc_score"],
        )

        exon_del[
            "splice_transformer_minigene_donor_mut"
        ] = np.where(
            exon_del["don_deleted"],
            0.0,
            exon_del["don_score"],
        )

        # ----------------------------------------------------
        # Delta calculations
        # ----------------------------------------------------

        exon_del["delta_acc_ST"] = (
            exon_del[
                "splice_transformer_minigene_acceptor_mut"
            ] - wt_acc
        )

        exon_del["delta_don_ST"] = (
            exon_del[
                "splice_transformer_minigene_donor_mut"
            ] - wt_don
        )

        exon_del["mean_delta_ST"] = (
            exon_del[
                ["delta_acc_ST", "delta_don_ST"]
            ]
            .mean(axis=1)
        )

        # ----------------------------------------------------
        # Standardized output columns
        # ----------------------------------------------------

        exon_del[
            "splice_transformer_minigene_delta_acceptor"
        ] = exon_del["delta_acc_ST"]

        exon_del[
            "splice_transformer_minigene_delta_donor"
        ] = exon_del["delta_don_ST"]

        exon_del[
            "splice_transformer_minigene_delta_mean"
        ] = exon_del["mean_delta_ST"]

        exon_del[
            "splice_transformer_minigene_var_source"
        ] = "DEL"

        exon_del["ensembl_exon_id"] = exon_id
        exon_del["exon_key"] = exon_id
        exon_del["variant_source"] = "DEL"

        all_del_results.append(exon_del)

    # ========================================================
    # SNVS
    # ========================================================

    exon_snv = var_exon[
        ~var_exon["variant_id"]
        .astype(str)
        .str.contains(
            "del",
            case=False,
            na=False,
        )
    ].copy()

    if not exon_snv.empty:

        snv = st[
            (
                ~st["variant"]
                .astype(str)
                .str.contains(
                    "del",
                    case=False,
                    na=False,
                )
            )
            &
            (
                st["nt_position"]
                .isin([ACC_POS, DON_POS])
            )
        ].copy()

        # ----------------------------------------------------
        # Extract canonical-site scores
        # ----------------------------------------------------

        acc = snv[
            snv["nt_position"] == ACC_POS
        ][
            ["variant", "Acceptor"]
        ]

        don = snv[
            snv["nt_position"] == DON_POS
        ][
            ["variant", "Donor"]
        ]

        merged = acc.merge(
            don,
            on="variant",
            how="outer",
        )

        # ----------------------------------------------------
        # Delta calculations
        # ----------------------------------------------------

        merged["delta_acc_ST"] = (
            merged["Acceptor"] - wt_acc
        )

        merged["delta_don_ST"] = (
            merged["Donor"] - wt_don
        )

        merged["mean_delta_ST"] = (
            merged[
                ["delta_acc_ST", "delta_don_ST"]
            ]
            .mean(axis=1)
        )

        merged["ensembl_exon_id"] = exon_id

        # ----------------------------------------------------
        # WT + mutant scores
        # ----------------------------------------------------

        merged[
            "splice_transformer_minigene_acceptor_wt"
        ] = wt_acc

        merged[
            "splice_transformer_minigene_donor_wt"
        ] = wt_don

        merged[
            "splice_transformer_minigene_acceptor_mut"
        ] = merged["Acceptor"]

        merged[
            "splice_transformer_minigene_donor_mut"
        ] = merged["Donor"]

        # ----------------------------------------------------
        # Standardized output columns
        # ----------------------------------------------------

        merged[
            "splice_transformer_minigene_delta_acceptor"
        ] = merged["delta_acc_ST"]

        merged[
            "splice_transformer_minigene_delta_donor"
        ] = merged["delta_don_ST"]

        merged[
            "splice_transformer_minigene_delta_mean"
        ] = merged["mean_delta_ST"]

        # ----------------------------------------------------
        # Merge back onto metadata
        # ----------------------------------------------------

        exon_snv = exon_snv.merge(
            merged,
            left_on="id",
            right_on="variant",
            how="left",
        )

        exon_snv[
            "splice_transformer_minigene_var_source"
        ] = "SNV"

        exon_snv["exon_key"] = exon_id
        exon_snv["variant_source"] = "SNV"

        all_snv_results.append(exon_snv)

# ============================================================
# CONCATENATE OUTPUTS
# ============================================================

df_del = (
    pd.concat(
        all_del_results,
        ignore_index=True,
    )
    if all_del_results
    else pd.DataFrame()
)

df_snv = (
    pd.concat(
        all_snv_results,
        ignore_index=True,
    )
    if all_snv_results
    else pd.DataFrame()
)

# ============================================================
# SAVE OUTPUTS
# ============================================================

df_del.to_csv(
    OUT_DEL,
    sep="\t",
    index=False,
)

df_snv.to_csv(
    OUT_SNV,
    sep="\t",
    index=False,
)

print("\nSaved deletion table:")
print(OUT_DEL)
print(df_del.shape)

print("\nSaved SNV table:")
print(OUT_SNV)
print(df_snv.shape)

# ============================================================
# FINISHED
# ============================================================

print("\n" + "=" * 80)
print("Completed SpliceTransformer minigene processing")
print("=" * 80)
