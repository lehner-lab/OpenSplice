#@title Process SpliceTransformer genomic outputs into corrected SNV/deletion delta tables
"""
Process per-exon SpliceTransformer genomic inference parquet files.

This script takes the raw per-position SpliceTransformer outputs generated for
SNVs and deletions, extracts canonical acceptor/donor scores, applies the
splice-site deletion delta fix, and saves full audit tables for downstream
benchmarking.

Key behaviour preserved from the analysis notebook:
    - WT acceptor position is fixed at ACC_POS = 71.
    - WT donor position is ACC_POS + exon_length - 1.
    - For deletions, canonical splice-site coordinates are shifted back onto
      the mutant sequence when sequence before the splice site is deleted.
    - If a deletion removes the canonical acceptor or donor itself, the mutant
      score for that site is set to 0 before calculating delta.
    - Delta is always calculated as mutant_score - WT_score.
    - The output keeps the original experimental/variant metadata columns plus
      audit columns used in earlier analyses.

Repository notes:
    - The file paths below are placeholders. Edit them to match your repo/data
      layout before running.
    - This script intentionally avoids changing the core calculation logic.
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

# Set True for a one-exon smoke test before running the full dataset.
TEST_MODE = False

# In these SpliceTransformer genomic/minigene-aligned outputs, the canonical
# acceptor is at nt_position 71. The donor is inferred from exon_length.
ACC_POS = 71

# Raw SpliceTransformer output filename suffix from the inference script.
ST_FILE_SUFFIX = "_splice_transformer_genomic.parquet"


# ============================================================
# INPUT FILES — PLACEHOLDERS TO EDIT FOR YOUR REPO
# ============================================================

# NOTE: Replace these paths with repository-relative paths when finalising.
# Example repo-style paths are shown first; the original Colab paths are kept
# below as comments for traceability.

EXON_META_FILE = Path("data/input/opensplice_predictors_benchmarking_exon_metadata.tsv")
EXP_VAR_FILE = Path("data/input/opensplice_predictors_benchmarking_variant_metadata.tsv")
ST_DIR = Path("results/splice_transformer/genomic/raw_parquet")
# Set ST_DIR to wherever you saved per-exon SpliceTransformer genomic parquet outputs.

# Original analysis paths:
# EXON_META_FILE = Path(
#     "libraries/exon_paper_1/opensplice_predictors_benchmarking_exon_metadata.tsv"
# )
# EXP_VAR_FILE = Path(
#     "libraries/exon_paper_1/opensplice_predictors_benchmarking_variant_metadata.tsv"
# )
# ST_DIR = Path(
#     "splice_transformer/608_exons_parquet/"
# )


# ============================================================
# OUTPUT FILES — PLACEHOLDERS TO EDIT FOR YOUR REPO
# ============================================================

OUTDIR = Path("results/splice_transformer/genomic/processed")
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_DEL = OUTDIR / "splicetransformer_genome_deletions_clean_curated.tsv"
OUT_SNV = OUTDIR / "splicetransformer_genome_snvs_clean_curated.tsv"


# ============================================================
# HELPER FUNCTIONS
# ============================================================


def require_columns(df: pd.DataFrame, required_cols: list[str], table_name: str) -> None:
    """Fail early with a clear error if an expected input column is missing."""
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(
            f"{table_name} is missing required columns: {missing}\n"
            f"Available columns are: {list(df.columns)}"
        )


# ============================================================
# LOAD INPUT TABLES
# ============================================================

print("Loading exon metadata:", EXON_META_FILE)
df_exons = pd.read_csv(EXON_META_FILE, sep="\t")

print("Loading experimental/variant metadata:", EXP_VAR_FILE)
df_var = pd.read_csv(EXP_VAR_FILE, sep="\t")

require_columns(
    df_exons,
    required_cols=["ensembl_exon_id", "exon_length"],
    table_name="EXON_META_FILE",
)
require_columns(
    df_var,
    required_cols=["ensembl_exon_id", "variant_id", "id", "start", "end", "length"],
    table_name="EXP_VAR_FILE",
)

exon_lengths = df_exons.set_index("ensembl_exon_id")["exon_length"].to_dict()

st_files = sorted([f for f in os.listdir(ST_DIR) if f.endswith(".parquet")])

if TEST_MODE:
    st_files = st_files[:1]
    print("TEST_MODE=True: processing only the first SpliceTransformer parquet file.")

print("Found SpliceTransformer genomic parquet files:", len(st_files))

all_del_results: list[pd.DataFrame] = []
all_snv_results: list[pd.DataFrame] = []


# ============================================================
# MAIN LOOP OVER EXONS
# ============================================================

for file_name in tqdm(st_files, desc="Processing exons"):
    exon_id = file_name.replace(ST_FILE_SUFFIX, "")
    st_path = ST_DIR / file_name

    st = pd.read_parquet(st_path)
    require_columns(
        st,
        required_cols=["variant", "nt_position", "Acceptor", "Donor"],
        table_name=f"SpliceTransformer parquet for {exon_id}",
    )

    exon_len = exon_lengths.get(exon_id)
    if exon_len is None or pd.isna(exon_len):
        print("Missing exon_length for", exon_id)
        continue

    don_pos = ACC_POS + int(exon_len) - 1

    # WT rows are used as the reference score for this exon.
    wt = st[st["variant"].astype(str).str.contains("WT", case=False, na=False)]
    wt_acc = wt.loc[wt["nt_position"] == ACC_POS, "Acceptor"].max()
    wt_don = wt.loc[wt["nt_position"] == don_pos, "Donor"].max()

    if pd.isna(wt_acc) or pd.isna(wt_don):
        print(f"WT splice site scores missing for {exon_id}, skipping.")
        continue

    var_exon = df_var[df_var["ensembl_exon_id"] == exon_id].copy()
    if var_exon.empty:
        continue

    # ========================================================
    # DELETIONS
    # ========================================================

    exon_del = var_exon[
        var_exon["variant_id"].astype(str).str.contains("del", case=False, na=False)
    ].copy()

    if not exon_del.empty:
        exon_del["del_start"] = exon_del["start"].astype(int)
        exon_del["del_end"] = exon_del["end"].astype(int)
        exon_del["del_len"] = exon_del["length"].astype(int)

        # Retain old audit columns exactly.
        exon_del["del_start_global"] = exon_del["del_start"]
        exon_del["del_end_global"] = exon_del["del_end"]

        exon_del["acc_deleted"] = (
            (exon_del["del_start_global"] <= ACC_POS)
            & (exon_del["del_end_global"] >= ACC_POS)
        )

        exon_del["don_deleted"] = (
            (exon_del["del_start_global"] <= don_pos)
            & (exon_del["del_end_global"] >= don_pos)
        )

        # If a deletion is before the splice site, the corresponding site moves
        # left in the mutant sequence by del_len. If the site itself is deleted,
        # keep the original coordinate for auditing but set mutant score to 0.
        exon_del["shifted_acc"] = exon_del.apply(
            lambda row: ACC_POS - row["del_len"]
            if (row["del_end_global"] < ACC_POS and not row["acc_deleted"])
            else ACC_POS,
            axis=1,
        )

        exon_del["shifted_don"] = exon_del.apply(
            lambda row: don_pos - row["del_len"]
            if (row["del_end_global"] < don_pos and not row["don_deleted"])
            else don_pos,
            axis=1,
        )

        acc_df = st[["variant", "nt_position", "Acceptor"]].rename(
            columns={"nt_position": "pos"}
        )
        don_df = st[["variant", "nt_position", "Donor"]].rename(
            columns={"nt_position": "pos"}
        )

        exon_del = exon_del.merge(
            acc_df,
            how="left",
            left_on=["id", "shifted_acc"],
            right_on=["variant", "pos"],
        ).rename(columns={"Acceptor": "acc_score"})

        exon_del = exon_del.merge(
            don_df,
            how="left",
            left_on=["id", "shifted_don"],
            right_on=["variant", "pos"],
        ).rename(columns={"Donor": "don_score"})

        # Corrected logic:
        #   deleted canonical site -> mutant raw score = 0
        #   delta = mutant - WT
        exon_del["splice_transformer_genomic_acceptor_wt"] = wt_acc
        exon_del["splice_transformer_genomic_donor_wt"] = wt_don

        exon_del["splice_transformer_genomic_acceptor_mut"] = np.where(
            exon_del["acc_deleted"],
            0.0,
            exon_del["acc_score"],
        )

        exon_del["splice_transformer_genomic_donor_mut"] = np.where(
            exon_del["don_deleted"],
            0.0,
            exon_del["don_score"],
        )

        exon_del["delta_acc_ST"] = (
            exon_del["splice_transformer_genomic_acceptor_mut"] - wt_acc
        )
        exon_del["delta_don_ST"] = (
            exon_del["splice_transformer_genomic_donor_mut"] - wt_don
        )
        exon_del["mean_delta_ST"] = exon_del[["delta_acc_ST", "delta_don_ST"]].mean(
            axis=1
        )

        exon_del["splice_transformer_genomic_delta_acceptor"] = exon_del[
            "delta_acc_ST"
        ]
        exon_del["splice_transformer_genomic_delta_donor"] = exon_del["delta_don_ST"]
        exon_del["splice_transformer_genomic_delta_mean"] = exon_del["mean_delta_ST"]
        exon_del["splice_transformer_genomic_var_source"] = "DEL"

        exon_del["ensembl_exon_id"] = exon_id
        exon_del["exon_key"] = exon_id
        exon_del["variant_source"] = "DEL"

        all_del_results.append(exon_del)

    # ========================================================
    # SNVs
    # ========================================================

    exon_snv = var_exon[
        ~var_exon["variant_id"].astype(str).str.contains("del", case=False, na=False)
    ].copy()

    if not exon_snv.empty:
        snv = st[
            (~st["variant"].astype(str).str.contains("WT", case=False, na=False))
            & (~st["variant"].astype(str).str.contains("del", case=False, na=False))
            & (st["nt_position"].isin([ACC_POS, don_pos]))
        ].copy()

        acc = snv[snv["nt_position"] == ACC_POS][["variant", "Acceptor"]]
        don = snv[snv["nt_position"] == don_pos][["variant", "Donor"]]

        merged = acc.merge(don, on="variant", how="outer")

        merged["delta_acc_ST"] = merged["Acceptor"] - wt_acc
        merged["delta_don_ST"] = merged["Donor"] - wt_don
        merged["mean_delta_ST"] = merged[["delta_acc_ST", "delta_don_ST"]].mean(axis=1)
        merged["ensembl_exon_id"] = exon_id

        merged["splice_transformer_genomic_acceptor_wt"] = wt_acc
        merged["splice_transformer_genomic_donor_wt"] = wt_don
        merged["splice_transformer_genomic_acceptor_mut"] = merged["Acceptor"]
        merged["splice_transformer_genomic_donor_mut"] = merged["Donor"]
        merged["splice_transformer_genomic_delta_acceptor"] = merged["delta_acc_ST"]
        merged["splice_transformer_genomic_delta_donor"] = merged["delta_don_ST"]
        merged["splice_transformer_genomic_delta_mean"] = merged["mean_delta_ST"]

        exon_snv = exon_snv.merge(
            merged,
            left_on="id",
            right_on="variant",
            how="left",
        )

        exon_snv["splice_transformer_genomic_var_source"] = "SNV"
        exon_snv["exon_key"] = exon_id
        exon_snv["variant_source"] = "SNV"

        all_snv_results.append(exon_snv)


# ============================================================
# SAVE FULL TABLES
# ============================================================

df_del = pd.concat(all_del_results, ignore_index=True) if all_del_results else pd.DataFrame()
df_snv = pd.concat(all_snv_results, ignore_index=True) if all_snv_results else pd.DataFrame()

df_del.to_csv(OUT_DEL, sep="\t", index=False)
df_snv.to_csv(OUT_SNV, sep="\t", index=False)

print("Saved DEL:", OUT_DEL, df_del.shape)
print("Saved SNV:", OUT_SNV, df_snv.shape)
