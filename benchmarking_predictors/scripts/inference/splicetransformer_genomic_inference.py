#@title Run SpliceTransformer genomic-context SNV/deletion inference with auto-resume
"""
Run SpliceTransformer inference for genomic-context SNV/deletion mutagenesis.

This script reads a per-variant input table containing WT and mutant genomic
sequence windows, runs SpliceTransformer for each variant, extracts the variant
sequence window from the model output, and saves one Parquet file per exon.

The original Colab logic is intentionally preserved:
    input_seq = "N" * 4000 + spliceai_seq_genome + "N" * 4000
    start     = 5001
    end       = 5001 + len(nt_seq)

Auto-resume is enabled by default: if an exon-level output file already exists,
that exon is skipped.

Expected input columns
----------------------
The input Parquet file should contain at least:
    - ensembl_exon_id      : exon/group identifier used to split outputs
    - ID                   : variant identifier
    - spliceai_seq_genome  : genomic-context sequence to feed into SpliceTransformer
    - nt_seq               : variant sequence whose length defines output slice

Notes
-----
This script is intended for repository use. Replace the placeholder paths below
with the final public data paths used in your repo.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from tqdm import tqdm
from sptransformer import Annotator


# ============================================================
# USER SETTINGS
# ============================================================

# Set True for a smoke test on the first exon only.
TEST_MODE = False

# Padding used on each side of the genomic-context input sequence.
N_PADDING = 4000

# Original notebook slice offset. Keep unchanged unless the upstream input
# construction changes.
MODEL_OUTPUT_START = 5001

# -------------------------------------------------------------------------
# Placeholder repository paths.
# Replace these with final repo-relative paths before release.
# -------------------------------------------------------------------------

INPUT_FILE = Path("data/input/splice_transformer_genomic_variant_input.parquet")
OUTDIR = Path("results/splice_transformer/genomic/raw_parquet")

# Original analysis paths retained here for reference.
# INPUT_FILE = Path(
#     "libraries/exon_608/splice_transformer_genomic_variant_input.parquet"
# )
# OUTDIR = Path(
#     "608_exons_parquet/"
# )


# ============================================================
# SPLICE TRANSFORMER CHANNEL NAMES
# ============================================================

CHANNEL_NAMES = (
    ["Neither", "Acceptor", "Donor"]
    + [
        "Splicing_" + name
        for name in [
            "Adipose Tissue",
            "Blood",
            "Blood Vessel",
            "Brain",
            "Colon",
            "Heart",
            "Kidney",
            "Liver",
            "Lung",
            "Muscle",
            "Nerve",
            "Small Intestine",
            "Skin",
            "Spleen",
            "Stomach",
        ]
    ]
)


# ============================================================
# HELPER FUNCTIONS
# ============================================================


def validate_input_columns(df: pd.DataFrame) -> None:
    """Check that the input table contains the columns required by the script."""
    required_cols = ["ensembl_exon_id", "ID", "spliceai_seq_genome", "nt_seq"]
    missing_cols = [col for col in required_cols if col not in df.columns]

    if missing_cols:
        raise ValueError(
            "Input file is missing required columns: "
            + ", ".join(missing_cols)
            + f"\nAvailable columns: {list(df.columns)}"
        )



def make_splicetransformer_input(seq_genome: str) -> str:
    """Pad a genomic-context sequence with Ns on both sides."""
    return "N" * N_PADDING + str(seq_genome) + "N" * N_PADDING



def run_variant_inference(
    model,
    variant_id: str,
    seq_genome: str,
    mut_len: int,
) -> pd.DataFrame:
    """Run SpliceTransformer on one variant and return the sliced output table."""
    input_seq = make_splicetransformer_input(seq_genome)

    encoded = model.one_hot_encode(input_seq).T
    tensor = torch.tensor(encoded).unsqueeze(0).float().to(model.device)

    with torch.no_grad():
        out = model.step(tensor)

    out_np = out.squeeze(0).cpu().numpy().T
    df_out = pd.DataFrame(out_np, columns=CHANNEL_NAMES)

    start = MODEL_OUTPUT_START
    end = MODEL_OUTPUT_START + mut_len

    df_out = df_out.iloc[start:end].copy()
    df_out["variant"] = variant_id
    df_out["nt_position"] = np.arange(1, mut_len + 1)

    return df_out



def make_parquet_safe(df: pd.DataFrame) -> pd.DataFrame:
    """Convert object columns to strings so mixed dtypes do not break Parquet saving."""
    df = df.copy()
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].astype(str)
    return df


# ============================================================
# MAIN SCRIPT
# ============================================================


def main() -> None:
    """Run exon-wise SpliceTransformer inference with auto-resume."""
    OUTDIR.mkdir(parents=True, exist_ok=True)

    print(f"Loading input file: {INPUT_FILE}")
    df = pd.read_parquet(INPUT_FILE)
    validate_input_columns(df)

    print("Loaded input:", df.shape)
    print("Unique ensembl_exon_ids:", df["ensembl_exon_id"].nunique())

    print("Initialising SpliceTransformer...")
    annotator = Annotator()
    model = annotator.model
    model.eval()

    all_exons = df["ensembl_exon_id"].dropna().unique()

    if TEST_MODE:
        exon_ids_to_process = all_exons[:1]
        print(f"TEST MODE: processing only: {exon_ids_to_process[0]}")
    else:
        exon_ids_to_process = all_exons
        print(f"FULL RUN: processing {len(all_exons)} exons")

    for exon_id in tqdm(exon_ids_to_process, desc="Processing exons"):
        out_path = OUTDIR / f"{exon_id}_splice_transformer_genomic.parquet"

        # Auto-resume: skip exon if output already exists.
        if out_path.exists():
            print(f"Skipping {exon_id} because output already exists: {out_path}")
            continue

        df_exon = df[df["ensembl_exon_id"] == exon_id]
        all_outputs = []

        for _, row in df_exon.iterrows():
            variant_id = row["ID"]
            seq_genome = row["spliceai_seq_genome"]
            mut_len = len(row["nt_seq"])

            try:
                df_out = run_variant_inference(
                    model=model,
                    variant_id=variant_id,
                    seq_genome=seq_genome,
                    mut_len=mut_len,
                )
                all_outputs.append(df_out)

            except Exception as exc:
                print(f"Error processing variant {variant_id} in exon {exon_id}: {exc}")
                continue

        if not all_outputs:
            print(f"No outputs generated for exon {exon_id}; nothing saved.")
            continue

        exon_df = pd.concat(all_outputs, ignore_index=True)
        exon_df = make_parquet_safe(exon_df)
        exon_df.to_parquet(out_path, index=False)

        print(f"Saved: {out_path} ({len(exon_df)} rows)")

    print("Completed processing. Auto-resume was enabled.")


if __name__ == "__main__":
    main()
