#@title Run SpliceAI deletion mutagenesis in minigene-context mode
"""
Run SpliceAI deletion-mutagenesis inference in minigene-context mode.

This script scores deletion variants for each exon using the original minigene
construct context:

    N-padding + pre_manual + variant nt_seq + post_manual + N-padding

where:
    pre_manual  = FAS exon 5 + FAS intron 5
    post_manual = FAS intron 6 + FAS exon 7

It extracts SpliceAI acceptor and donor scores at the canonical middle-exon
splice sites using the original minigene coordinate convention:

    acceptor_pos0 = 216
    donor_pos0    = 216 - 1 + exon_length

For deletions, ALT predictions are realigned back to WT coordinates by inserting
zeros at the deletion start position using the original minigene-coordinate rule:

    deletion_start_insert_pos0 = 146 + start - 1

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

DATA_DIR = Path("data/input")
OUTPUT_DIR = Path("results/spliceai/minigene_mode/deletions")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Variant-level metadata used for deletion rows, WT rows, start/length,
# exon_length, and nt_seq.
#
# libraries/exon_paper_1/opensplice_predictors_benchmarking_variant_metadata.tsv
#
# NOTE:
# Replace this placeholder filename with the final supplementary-data filename
# during public repository finalisation.
META_FILE = DATA_DIR / "opensplice_predictors_benchmarking_variant_metadata.tsv"


# ============================================================
# SETTINGS
# ============================================================

# SpliceAI context padding: 5000 N on each side.
CONTEXT = 10000

# Batch size for deletion inference.
BATCH_SIZE = 250

# Optional: skip outputs that already exist.
SKIP_EXISTING = True

# Canonical acceptor coordinate in the minigene construct.
# Original manuscript-analysis convention:
#   146 nt manual prefix + 70 nt upstream intron = 216
ACC_POS = 216

# Deletion realignment offset in minigene mode.
# Original manuscript-analysis convention:
#   deletion_start_insert_pos0 = 146 + start - 1
#
# This reflects the point where the variable nt_seq begins after pre_manual.
PRE_MANUAL_LEN_FOR_DELETION_REALIGNMENT = 146


# ============================================================
# MANUAL MINIGENE FLANKING SEQUENCES
# ============================================================

# These FAS minigene backbone sequences are kept exactly as in the original
# notebook. They define the minigene context used for SpliceAI inference.
FAS_E5 = "ATGTGAACATGGAATCATCAAGGAATGCACACTCACCAGCAACACCAAGTGCAAAGAGGAAG"
FAS_I5 = "GTAATTATTTTTTTACGGTTATATTCTCCTTTCCCCCAACCCCATGGAAAGATGTGAAGAAAAACCAATCACTCTTGATTACTA"
FAS_I6 = "CAGATTGAAATAACTTGGGAAGTAGTTTCTCTTAGTGTGAAAGTATGTTCTCACATGCATTCTACAAGGCTGAGACCTGAGTTGATAAAATTTCTTTGTTCTTTCAG"
FAS_E7 = "TGAAGAGAAAGGAAGTACAGAAAACATGCAGAAAGCACAGAAAGGAA"

PRE_MANUAL = FAS_E5 + FAS_I5
POST_MANUAL = FAS_I6 + FAS_E7


# ============================================================
# LOAD INPUT TABLE
# ============================================================

meta_df = pd.read_csv(META_FILE, sep="\t")

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

# Keep deletion rows only for deciding which exons to process.
del_df = meta_df[
    meta_df["variant_id"].astype(str).str.contains("_del", na=False)
].copy()

exon_ids = del_df["ensembl_exon_id"].dropna().unique().tolist()

print(f"Will process {len(exon_ids)} exons in SpliceAI minigene deletion mode")


# ============================================================
# MAIN LOOP
# ============================================================

for exon_id in exon_ids:
    out_file = OUTPUT_DIR / f"{exon_id}_spliceai_scores_minigene.parquet"

    if SKIP_EXISTING and out_file.exists():
        print(f"Skipping {exon_id}, output already exists: {out_file}")
        continue

    print(f"\nProcessing exon: {exon_id}")

    # --------------------------------------------------------
    # 1. Load SpliceAI ensemble fresh per exon.
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
    exon_meta = meta_df[meta_df["ensembl_exon_id"] == exon_id].copy()

    # --------------------------------------------------------
    # 3. Find the WT row for this exon.
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

    donor_pos = ACC_POS - 1 + exon_len
    pois = [ACC_POS, donor_pos]

    print(f"  exon_length: {exon_len}")
    print(f"  acceptor position: {ACC_POS}")
    print(f"  donor position:    {donor_pos}")

    # --------------------------------------------------------
    # 4. REF / WT prediction.
    # --------------------------------------------------------
    ref_seq = (
        "N" * (CONTEXT // 2)
        + PRE_MANUAL
        + wt_seq
        + POST_MANUAL
        + "N" * (CONTEXT // 2)
    )

    ref_oh = one_hot_encode(ref_seq)[None, :]
    ref_pred = np.mean([m.predict(ref_oh, verbose=0) for m in models], axis=0)[0]

    ref_acc = ref_pred[:, 1]
    ref_don = ref_pred[:, 2]

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
            "acceptor_pos0": ACC_POS,
            "donor_pos0": donor_pos,
            **{f"acceptor_score_{p}": ref_acc[p] for p in pois},
            **{f"donor_score_{p}": ref_don[p] for p in pois},
        }
    ]

    # --------------------------------------------------------
    # 5. ALT deletion predictions.
    # --------------------------------------------------------
    alt_meta = exon_meta[
        exon_meta["variant_id"].astype(str).str.contains("_del", na=False)
    ].copy()

    for del_length, grp in alt_meta.groupby("length"):
        print(f"  deletions of length {del_length}")

        for i in range(0, len(grp), BATCH_SIZE):
            batch = grp.iloc[i:i + BATCH_SIZE]

            encs = []
            batch_rows = []

            for _, row in batch.iterrows():
                alt_seq = (
                    "N" * (CONTEXT // 2)
                    + PRE_MANUAL
                    + str(row["nt_seq"])
                    + POST_MANUAL
                    + "N" * (CONTEXT // 2)
                )

                encs.append(one_hot_encode(alt_seq))
                batch_rows.append(row)

            x_batch = np.stack(encs)
            y_batch = np.mean([m.predict(x_batch, verbose=0) for m in models], axis=0)

            for idx, row in enumerate(batch_rows):
                acc = y_batch[idx, :, 1]
                don = y_batch[idx, :, 2]

                # Realign ALT predictions back to WT coordinates.
                # This keeps the original minigene-mode coordinate rule:
                #   ds = 146 + start - 1
                ds = PRE_MANUAL_LEN_FOR_DELETION_REALIGNMENT + int(row["start"]) - 1
                deletion_length = int(row["length"])

                acc = np.insert(acc, ds, [0] * deletion_length)
                don = np.insert(don, ds, [0] * deletion_length)

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
                        "deletion_start_insert_pos0": ds,
                        "deletion_length": deletion_length,
                        "acceptor_pos0": ACC_POS,
                        "donor_pos0": donor_pos,
                        **{f"acceptor_score_{p}": acc[p] for p in pois},
                        **{f"donor_score_{p}": don[p] for p in pois},
                    }
                )

            del x_batch, y_batch, encs, batch_rows
            gc.collect()

    # --------------------------------------------------------
    # 6. Save per-exon parquet.
    # --------------------------------------------------------
    res_df = pd.DataFrame(results)
    res_df.to_parquet(out_file, engine="pyarrow", compression="snappy")

    print(f"Saved {len(res_df)} rows to {out_file}")

    # --------------------------------------------------------
    # 7. Cleanup model memory.
    # --------------------------------------------------------
    del models, ref_oh, ref_pred, ref_acc, ref_don, results, res_df
    gc.collect()
