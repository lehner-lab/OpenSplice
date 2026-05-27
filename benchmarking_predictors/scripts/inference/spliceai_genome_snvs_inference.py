#@title Run SpliceAI SNV inference in genomic-context mode.
"""
Run SpliceAI SNV-only inference in genomic-context mode.

This script scores SNV-containing exon constructs using genomic flanks:

    up_5k + variant_nt_seq + down_5k

For each exon, the script:
  1. loads the five SpliceAI ensemble models,
  2. scores the WT sequence at the canonical acceptor and donor positions,
  3. scores each SNV sequence at the same canonical positions,
  4. saves one parquet file per exon.

The computational logic is intentionally kept close to the original Colab notebook
used for the manuscript analyses. In particular, the canonical acceptor coordinate
is defined as:

    len(up_5k) + 70

because nt_seq is assumed to begin with 70 nt of upstream intronic sequence before
the assayed exon.
"""

from __future__ import annotations

import gc
import os
from pathlib import Path

import numpy as np
import pandas as pd
from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode


# ============================================================
# INPUT FILES
# ============================================================

# Repository-relative placeholder paths.
# These filenames can be updated once the public supplementary data filenames
# are finalised.
DATA_DIR = Path("data/input")
OUTPUT_DIR = Path("results/spliceai/genome_mode/snvs")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

META_FILE = DATA_DIR / "opensplice_predictors_benchmarking_variant_metadata.tsv"

CONTEXT_FILE = DATA_DIR / "opensplice_predictors_benchmarking_exon_metadata.tsv"


# ============================================================
# SETTINGS
# ============================================================

# Same SpliceAI padding logic as the original notebook.
CONTEXT = 10000

# Coordinate assumption:
# nt_seq begins with the 70 nt upstream intronic region before the exon.
# Therefore the canonical acceptor is at:
#   len(up_5k) + 70
UPSTREAM_INTRON_IN_NTSEQ = 70

# Batch size for SNV inference.
BATCH_SIZE = 1000


# ============================================================
# LOAD INPUT TABLES
# ============================================================

meta_df = pd.read_csv(META_FILE, sep="\t")
ctx_df = pd.read_csv(CONTEXT_FILE, sep="\t")

print("Loaded variant metadata:", meta_df.shape)
print("Loaded context metadata:", ctx_df.shape)

required_meta_cols = [
    "ensembl_exon_id",
    "variant_id",
    "nt_seq",
    "exon_length",
]
missing_meta = [c for c in required_meta_cols if c not in meta_df.columns]
if missing_meta:
    raise ValueError(f"META_FILE is missing required columns: {missing_meta}")

required_ctx_cols = [
    "ensembl_exon_id",
    "up_5k",
    "wt_seq",
    "down_5k",
]
missing_ctx = [c for c in required_ctx_cols if c not in ctx_df.columns]
if missing_ctx:
    raise ValueError(f"CONTEXT_FILE is missing required columns: {missing_ctx}")

# Keep one context row per exon.
ctx_df = ctx_df.drop_duplicates(subset=["ensembl_exon_id"]).copy()
ctx_map = ctx_df.set_index("ensembl_exon_id")


# ============================================================
# FILTER TO SNV-CONTAINING EXONS
# ============================================================

# Keep all non-deletion rows. This includes WT and SNVs.
snv_df = meta_df[
    ~meta_df["variant_id"].astype(str).str.contains("_del", na=False)
].copy()

# Exons that have at least one non-WT, non-deletion variant.
snv_only_df = snv_df[
    ~snv_df["variant_id"].astype(str).str.endswith("_wt", na=False)
].copy()

exon_ids = snv_only_df["ensembl_exon_id"].dropna().unique().tolist()
print(f"Will process {len(exon_ids)} exons (SNV-only, genomic context)")


# ============================================================
# MAIN LOOP
# ============================================================

for exon_id in exon_ids:
    out_file = OUTPUT_DIR / f"{exon_id}_snv_spliceai_genome_mode.parquet"

    if os.path.exists(out_file):
        print(f"Skipping {exon_id}, output already exists: {out_file}")
        continue

    print(f"\nProcessing exon: {exon_id}")

    # --------------------------------------------------------
    # 1. Reload models per exon to keep GPU/VRAM usage lower.
    # --------------------------------------------------------
    models = [
        load_model(
            resource_filename("spliceai", f"models/spliceai{x}.h5"),
            compile=False,
        )
        for x in range(1, 6)
    ]

    # --------------------------------------------------------
    # 2. Subset this exon.
    # --------------------------------------------------------
    exon_meta = snv_df[snv_df["ensembl_exon_id"] == exon_id].copy()

    if exon_id not in ctx_map.index:
        print(f"No genomic context row found for {exon_id}; skipping")
        del models
        gc.collect()
        continue

    ctx_row = ctx_map.loc[exon_id]

    up_5k = str(ctx_row["up_5k"])
    ctx_wt_seq = str(ctx_row["wt_seq"])
    down_5k = str(ctx_row["down_5k"])

    # --------------------------------------------------------
    # 3. Find WT row for this exon.
    # --------------------------------------------------------
    wt_candidates = exon_meta[
        exon_meta["variant_id"].astype(str).str.endswith("_wt", na=False)
    ].copy()

    if wt_candidates.empty:
        print(f"No WT row found for {exon_id}; skipping")
        del models
        gc.collect()
        continue

    wt_row = wt_candidates.iloc[0]
    ref_id = wt_row["variant_id"]
    wt_seq = str(wt_row["nt_seq"])
    exon_len = int(wt_row["exon_length"])

    # Optional sanity check: the WT nt_seq in the context metadata should match
    # the WT nt_seq from the variant metadata. We keep using the variant metadata
    # sequence for consistency with the SNV rows.
    if wt_seq != ctx_wt_seq:
        print(
            f"WT nt_seq mismatch for {exon_id}: "
            f"META_FILE nt_seq length={len(wt_seq)}, "
            f"CONTEXT_FILE wt_seq length={len(ctx_wt_seq)}. "
            f"Using META_FILE nt_seq for inference."
        )

    # --------------------------------------------------------
    # 4. Define genomic-context splice-site coordinates.
    # --------------------------------------------------------
    # Old minigene mode used:
    #   ACC_POS = 216
    #   DON_POS = 216 + exon_len - 1
    # because:
    #   146 nt manual prefix + 70 nt upstream intron = 216
    #
    # Genomic-context mode uses:
    #   ACC_POS = len(up_5k) + 70
    #   DON_POS = ACC_POS + exon_len - 1
    upstream_context_len = len(up_5k)
    ACC_POS = upstream_context_len + UPSTREAM_INTRON_IN_NTSEQ
    DON_POS = ACC_POS + exon_len - 1
    pois = [ACC_POS, DON_POS]

    print(f"  up_5k length: {upstream_context_len}")
    print(f"  exon_length: {exon_len}")
    print(f"  acceptor position: {ACC_POS}")
    print(f"  donor position:    {DON_POS}")

    # --------------------------------------------------------
    # 5. REF / WT prediction.
    # --------------------------------------------------------
    ref_core_seq = up_5k + wt_seq + down_5k
    ref_seq = "N" * (CONTEXT // 2) + ref_core_seq + "N" * (CONTEXT // 2)

    ref_oh = one_hot_encode(ref_seq)[None, :]
    ref_pred = np.mean([m.predict(ref_oh, verbose=0) for m in models], axis=0)[0]

    ref_acc = ref_pred[:, 1]
    ref_don = ref_pred[:, 2]

    max_needed_pos = max(pois)
    if max_needed_pos >= len(ref_acc):
        raise ValueError(
            f"{exon_id}: requested position {max_needed_pos} "
            f"is outside REF prediction length {len(ref_acc)}"
        )

    results = [
        {
            "ensembl_exon_id": exon_id,
            "Type": "REF",
            "Identifier": ref_id,
            "up_5k_len": upstream_context_len,
            "acceptor_pos0": ACC_POS,
            "donor_pos0": DON_POS,
            **{f"acceptor_score_{p}": ref_acc[p] for p in pois},
            **{f"donor_score_{p}": ref_don[p] for p in pois},
        }
    ]

    # --------------------------------------------------------
    # 6. ALT SNV predictions.
    # --------------------------------------------------------
    alt_meta = exon_meta[
        ~exon_meta["variant_id"].astype(str).str.endswith("_wt", na=False)
    ].copy()

    print(f"  SNVs to process: {len(alt_meta)}")

    for i in range(0, len(alt_meta), BATCH_SIZE):
        batch = alt_meta.iloc[i : i + BATCH_SIZE]

        encs = []
        batch_rows = []

        for _, row in batch.iterrows():
            alt_nt_seq = str(row["nt_seq"])

            # In genomic mode the variant-specific construct is:
            #   up_5k + variant nt_seq + down_5k
            alt_core_seq = up_5k + alt_nt_seq + down_5k
            alt_seq = "N" * (CONTEXT // 2) + alt_core_seq + "N" * (CONTEXT // 2)

            encs.append(one_hot_encode(alt_seq))
            batch_rows.append(row)

        x_batch = np.stack(encs)
        y_batch = np.mean([m.predict(x_batch, verbose=0) for m in models], axis=0)

        for idx, row in enumerate(batch_rows):
            acc = y_batch[idx, :, 1]
            don = y_batch[idx, :, 2]

            if max_needed_pos >= len(acc):
                raise ValueError(
                    f"{exon_id} / {row['variant_id']}: requested position "
                    f"{max_needed_pos} is outside ALT prediction length {len(acc)}"
                )

            results.append(
                {
                    "ensembl_exon_id": exon_id,
                    "Type": "SNV",
                    "Identifier": row["variant_id"],
                    "up_5k_len": upstream_context_len,
                    "acceptor_pos0": ACC_POS,
                    "donor_pos0": DON_POS,
                    **{f"acceptor_score_{p}": acc[p] for p in pois},
                    **{f"donor_score_{p}": don[p] for p in pois},
                }
            )

        del x_batch, y_batch, encs, batch_rows
        gc.collect()

    # --------------------------------------------------------
    # 7. Save per-exon parquet.
    # --------------------------------------------------------
    res_df = pd.DataFrame(results)
    res_df.to_parquet(out_file, engine="pyarrow", compression="snappy")
    print(f"Saved {len(res_df)} rows: {out_file}")

    # --------------------------------------------------------
    # 8. Cleanup.
    # --------------------------------------------------------
    del models, ref_oh, ref_pred, ref_acc, ref_don, results, res_df
    gc.collect()
