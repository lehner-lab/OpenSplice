#@title Run SpliceAI deletion mutagenesis in genomic-context mode
"""
Run SpliceAI deletion-mutagenesis inference in genomic-context mode.

This script scores deletion variants for each exon using:
    up_5k + variant nt_seq + down_5k
with N-padding on both sides, then realigns deletion predictions back to
WT coordinates by inserting zeros at the deleted interval.

The computational logic is intentionally kept close to the original Colab
notebook used for the manuscript analyses.
"""
from __future__ import annotations

import gc
from pathlib import Path

import numpy as np
import pandas as pd
from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode

# ============================================================
# INPUT FILES
# ============================================================

# NOTE:
# These repository-relative filenames are placeholders for the public paper repo.
# During finalisation, align these with the supplementary-data filenames.
DATA_DIR = Path("data/input")
OUTPUT_DIR = Path("results/spliceai/genome_mode/deletions")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

META_FILE = DATA_DIR / "opensplice_predictors_benchmarking_variant_metadata.tsv"

CONTEXT_FILE = DATA_DIR / "opensplice_predictors_benchmarking_exon_metadata.tsv"

# ============================================================
# SETTINGS
# ============================================================

# SpliceAI N-padding length used in the original analysis.
CONTEXT = 10_000

# Important coordinate assumption:
# nt_seq in the variant metadata starts with the 70 nt upstream intronic
# segment before the assayed exon. Therefore the canonical acceptor site is:
#     len(up_5k) + 70
UPSTREAM_INTRON_IN_NTSEQ = 70

# Batch size used for deletion prediction.
BATCH_SIZE = 250

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
    "start",
    "length",
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

# Keep one row per exon in the genomic-context table.
ctx_df = ctx_df.drop_duplicates(subset=["ensembl_exon_id"]).copy()
ctx_map = ctx_df.set_index("ensembl_exon_id")

# Keep deletion rows for the main loop. WT rows are still retrieved later from
# meta_df within each exon subset.
del_df = meta_df[
    meta_df["variant_id"].astype(str).str.contains("_del", na=False)
].copy()

# Exons to process are those with at least one deletion row.
exon_ids = del_df["ensembl_exon_id"].dropna().unique().tolist()
print(f"Will process {len(exon_ids)} exons with deletion rows")

# ============================================================
# MAIN LOOP
# ============================================================

for exon_id in exon_ids:
    out_file = OUTPUT_DIR / f"{exon_id}_spliceai_scores_genome_mode.parquet"

    if out_file.exists():
        print(f"Skipping {exon_id}, output already exists")
        continue

    print(f"\nProcessing exon: {exon_id}")

    # --------------------------------------------------------
    # 1. Load the five SpliceAI ensemble models fresh per exon.
    #    This mirrors the original memory-saving strategy.
    # --------------------------------------------------------
    models = [
        load_model(
            resource_filename("spliceai", f"models/spliceai{x}.h5"),
            compile=False,
        )
        for x in range(1, 6)
    ]

    # --------------------------------------------------------
    # 2. Subset this exon in the variant table and retrieve context.
    # --------------------------------------------------------
    exon_meta = meta_df[meta_df["ensembl_exon_id"] == exon_id].copy()

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
    # 3. Find the WT row for this exon.
    # --------------------------------------------------------
    wt_candidates = exon_meta[
        exon_meta["variant_id"].astype(str).str.endswith("_wt", na=False)
    ].copy()

    if wt_candidates.empty:
        print(f"No WT row found in META_FILE for {exon_id}; skipping")
        del models
        gc.collect()
        continue

    wt_row = wt_candidates.iloc[0]

    ref_id = wt_row["variant_id"]
    wt_seq = str(wt_row["nt_seq"])
    exon_len = int(wt_row["exon_length"])

    # Optional safety check:
    # The context table also has a wt_seq column. The deletion coordinates are
    # defined on the variant-metadata nt_seq, so that is kept as the source of
    # truth if the two WT strings differ.
    if wt_seq != ctx_wt_seq:
        print(
            f"WT nt_seq mismatch for {exon_id}: "
            f"META_FILE nt_seq length={len(wt_seq)}, "
            f"CONTEXT_FILE wt_seq length={len(ctx_wt_seq)}. "
            f"Using META_FILE nt_seq because deletion coordinates are defined on that sequence."
        )

    # --------------------------------------------------------
    # 4. Define genomic-context splice-site coordinates.
    # --------------------------------------------------------
    # Original minigene mode used:
    #     acceptor = 146 + 70 = 216
    #
    # Genomic-context mode uses:
    #     acceptor = len(up_5k) + 70
    #     donor    = acceptor + exon_len - 1
    upstream_context_len = len(up_5k)
    acceptor_pos = upstream_context_len + UPSTREAM_INTRON_IN_NTSEQ
    donor_pos = acceptor_pos + exon_len - 1
    pois = [acceptor_pos, donor_pos]

    print(f"  up_5k length: {upstream_context_len}")
    print(f"  exon_length: {exon_len}")
    print(f"  acceptor position: {acceptor_pos}")
    print(f"  donor position:    {donor_pos}")

    # --------------------------------------------------------
    # 5. REF / WT prediction.
    # --------------------------------------------------------
    # REF construct:
    #     N-padding + up_5k + WT nt_seq + down_5k + N-padding
    ref_core_seq = up_5k + wt_seq + down_5k
    ref_seq = "N" * (CONTEXT // 2) + ref_core_seq + "N" * (CONTEXT // 2)

    ref_oh = one_hot_encode(ref_seq)[None, :]
    ref_pred = np.mean([m.predict(ref_oh, verbose=0) for m in models], axis=0)[0]

    ref_acc = ref_pred[:, 1]
    ref_don = ref_pred[:, 2]

    # Safety check: requested coordinates must exist inside the prediction track.
    max_needed_pos = max(pois)
    if max_needed_pos >= len(ref_acc):
        raise ValueError(
            f"{exon_id}: requested splice-site position {max_needed_pos} "
            f"is outside REF prediction length {len(ref_acc)}"
        )

    results = [
        {
            "ensembl_exon_id": exon_id,
            "Type": "REF",
            "Identifier": ref_id,
            "up_5k_len": upstream_context_len,
            "acceptor_pos0": acceptor_pos,
            "donor_pos0": donor_pos,
            **{f"acceptor_score_{p}": ref_acc[p] for p in pois},
            **{f"donor_score_{p}": ref_don[p] for p in pois},
        }
    ]

    # --------------------------------------------------------
    # 6. ALT deletion predictions.
    # --------------------------------------------------------
    alt_meta = exon_meta[
        exon_meta["variant_id"].astype(str).str.contains("_del", na=False)
    ].copy()

    if alt_meta.empty:
        print("  No deletion rows for this exon; saving REF only")
    else:
        for del_length, grp in alt_meta.groupby("length"):
            print(f"  - deletions of length {del_length}")

            for i in range(0, len(grp), BATCH_SIZE):
                batch = grp.iloc[i : i + BATCH_SIZE]

                # Build ALT sequences in genomic-context mode:
                #     up_5k + deletion-mutated nt_seq + down_5k
                encs = []
                batch_rows = []

                for _, row in batch.iterrows():
                    alt_nt_seq = str(row["nt_seq"])
                    alt_core_seq = up_5k + alt_nt_seq + down_5k
                    alt_seq = "N" * (CONTEXT // 2) + alt_core_seq + "N" * (CONTEXT // 2)

                    encs.append(one_hot_encode(alt_seq))
                    batch_rows.append(row)

                x_batch = np.stack(encs)
                y_batch = np.mean(
                    [m.predict(x_batch, verbose=0) for m in models],
                    axis=0,
                )

                for idx, row in enumerate(batch_rows):
                    acc = y_batch[idx, :, 1]
                    don = y_batch[idx, :, 2]

                    # Realign ALT predictions back to WT coordinates by inserting
                    # zeros at the deletion start site.
                    #
                    # Original minigene code used:
                    #     ds = 146 + start - 1
                    # because the variable region began after the manual prefix.
                    #
                    # Genomic-context mode uses:
                    #     ds = len(up_5k) + start - 1
                    # because the variable region begins after up_5k.
                    # start is assumed to be 1-based within nt_seq.
                    ds = upstream_context_len + int(row["start"]) - 1
                    deletion_len = int(row["length"])

                    acc = np.insert(acc, ds, [0] * deletion_len)
                    don = np.insert(don, ds, [0] * deletion_len)

                    if max_needed_pos >= len(acc):
                        raise ValueError(
                            f"{exon_id} / {row['variant_id']}: requested position "
                            f"{max_needed_pos} is outside realigned ALT length {len(acc)}"
                        )

                    results.append(
                        {
                            "ensembl_exon_id": exon_id,
                            "Type": "ALT",
                            "Identifier": row["variant_id"],
                            "up_5k_len": upstream_context_len,
                            "deletion_start_insert_pos0": ds,
                            "deletion_length": deletion_len,
                            "acceptor_pos0": acceptor_pos,
                            "donor_pos0": donor_pos,
                            **{f"acceptor_score_{p}": acc[p] for p in pois},
                            **{f"donor_score_{p}": don[p] for p in pois},
                        }
                    )

                del x_batch, y_batch, encs, batch_rows
                gc.collect()

    # --------------------------------------------------------
    # 7. Save one parquet output file per exon.
    # --------------------------------------------------------
    res_df = pd.DataFrame(results)
    res_df.to_parquet(out_file, engine="pyarrow", compression="snappy")

    print(f"Saved {len(res_df)} rows -> {out_file}")

    # --------------------------------------------------------
    # 8. Cleanup model memory before the next exon.
    # --------------------------------------------------------
    del models, ref_oh, ref_pred, ref_acc, ref_don, results, res_df
    gc.collect()
