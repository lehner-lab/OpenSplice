"""
#@title Run SpliceTransformer minigene inference (SNVs + deletions + WT together)

Repository-ready script for running SpliceTransformer inference on
minigene-context mutagenesis data.

Key features
------------
- Runs ALL variants together (WT + SNVs + deletions)
- Avoids needing a second WT-only recovery script
- Auto-resume support (skip already processed exons)
- Repository-style structure with placeholder paths
- Preserves original inference logic
- Saves one parquet per exon

Expected input columns
----------------------
Required:
    ID
    ensembl_exon_id
    spliceai_seq_minigene
    nt_seq

Typical inputs include WT, SNV, and deletion rows together.
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
# SETTINGS
# ============================================================

TEST_MODE = False

# ============================================================
# INPUT FILES
# ============================================================

# NOTE:
# Replace these placeholder paths with repository-relative paths
# during final repository cleanup.

DATA_DIR = Path("data/input")

INPUT_FILE = (
    DATA_DIR /
    "splice_transformer_genomic_variant_input.parquet"
)

# Example original path:
# INPUT_FILE = (
#     "gioia_exon_sat_mut/libraries/exon_608/"
#     "splice_transformer_genomic_variant_input.parquet"
# )

# ============================================================
# OUTPUT DIRECTORY
# ============================================================

OUTPUT_DIR = Path(
    "results/splice_transformer/minigene_mode/inference"
)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# LOAD INPUT
# ============================================================

print("=" * 80)
print("Loading input file")
print("=" * 80)

df = pd.read_parquet(INPUT_FILE)

print("Loaded input shape:", df.shape)
print("Unique exons:", df["ensembl_exon_id"].nunique())

# ============================================================
# INITIALIZE SPLICE TRANSFORMER
# ============================================================

print("\nInitializing SpliceTransformer model...")

annotator = Annotator()
model = annotator.model

# ============================================================
# OUTPUT CHANNELS
# ============================================================

channel_names = (
    ["Neither", "Acceptor", "Donor"] +
    [
        "Splicing_" + name for name in [
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

print("Number of output channels:", len(channel_names))

# ============================================================
# DETERMINE EXONS TO PROCESS
# ============================================================

all_exons = df["ensembl_exon_id"].unique()

if TEST_MODE:
    exon_ids_to_process = all_exons[:1]
    print("\nTEST MODE ENABLED")
    print("Processing exon:", exon_ids_to_process[0])

else:
    exon_ids_to_process = all_exons
    print("\nFULL RUN")
    print("Total exons:", len(exon_ids_to_process))

# ============================================================
# MAIN LOOP
# ============================================================

print("\nStarting inference...")

for exon_id in tqdm(exon_ids_to_process, desc="Processing exons"):

    # --------------------------------------------------------
    # Output path
    # --------------------------------------------------------

    out_path = (
        OUTPUT_DIR /
        f"{exon_id}_splice_transformer_minigene.parquet"
    )

    # --------------------------------------------------------
    # AUTO-RESUME
    # --------------------------------------------------------

    if out_path.exists():
        print(f"Skipping {exon_id} (already exists)")
        continue

    # --------------------------------------------------------
    # Extract exon rows
    # --------------------------------------------------------

    df_exon = (
        df[df["ensembl_exon_id"] == exon_id]
        .copy()
    )

    if len(df_exon) == 0:
        print(f"No rows found for {exon_id}")
        continue

    all_outputs = []

    # ========================================================
    # VARIANT LOOP
    # ========================================================

    for _, row in df_exon.iterrows():

        try:

            # ------------------------------------------------
            # Required fields
            # ------------------------------------------------

            variant_id = str(row["ID"])

            seq_minigene = str(
                row["spliceai_seq_minigene"]
            )

            mut_len = len(str(row["nt_seq"]))

            # ------------------------------------------------
            # Construct model input
            # ------------------------------------------------

            # Keep original logic:
            # 4000 N padding on both sides

            input_seq = (
                ("N" * 4000) +
                seq_minigene +
                ("N" * 4000)
            )

            # ------------------------------------------------
            # One-hot encode sequence
            # ------------------------------------------------

            encoded = model.one_hot_encode(input_seq).T

            tensor = (
                torch.tensor(encoded)
                .unsqueeze(0)
                .float()
                .to(model.device)
            )

            # ------------------------------------------------
            # Run inference
            # ------------------------------------------------

            with torch.no_grad():
                out = model.step(tensor)

            # ------------------------------------------------
            # Convert output to DataFrame
            # ------------------------------------------------

            out_np = (
                out.squeeze(0)
                .cpu()
                .numpy()
                .T
            )

            df_out = pd.DataFrame(
                out_np,
                columns=channel_names,
            )

            # ------------------------------------------------
            # Extract mutational window
            # ------------------------------------------------

            # IMPORTANT:
            # This keeps the exact same extraction logic used
            # in the original notebook workflow.

            start = 146
            end = 146 + mut_len

            df_out = df_out.iloc[start:end].copy()

            # ------------------------------------------------
            # Add metadata columns
            # ------------------------------------------------

            df_out["variant"] = variant_id

            df_out["nt_position"] = np.arange(
                1,
                mut_len + 1,
            )

            df_out["ensembl_exon_id"] = exon_id

            # Helpful for downstream processing
            df_out["variant_length"] = mut_len

            # ------------------------------------------------
            # Store output
            # ------------------------------------------------

            all_outputs.append(df_out)

        except Exception as e:

            print(
                f"Error processing "
                f"{variant_id} ({exon_id}): {e}"
            )

            continue

    # ========================================================
    # SAVE OUTPUT
    # ========================================================

    if len(all_outputs) == 0:

        print(
            f"No successful outputs generated for "
            f"{exon_id}"
        )

        continue

    exon_df = pd.concat(
        all_outputs,
        ignore_index=True,
    )

    # --------------------------------------------------------
    # Parquet-safe string conversion
    # --------------------------------------------------------

    for col in exon_df.columns:

        if exon_df[col].dtype == object:
            exon_df[col] = exon_df[col].astype(str)

    # --------------------------------------------------------
    # Save parquet
    # --------------------------------------------------------

    exon_df.to_parquet(
        out_path,
        index=False,
    )

    print(
        f"Saved: {out_path} "
        f"({len(exon_df)} rows)"
    )

# ============================================================
# FINISHED
# ============================================================

print("\n" + "=" * 80)
print("Completed SpliceTransformer minigene inference")
print("=" * 80)
