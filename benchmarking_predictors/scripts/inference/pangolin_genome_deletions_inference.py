#@title Run Pangolin deletion mutagenesis in genomic-context mode.
"""
Run Pangolin predictions for deletion variants in genomic-context mode.

This script was converted from the original Colab analysis code for the
benchmarking repository. The computational logic is intentionally preserved:
- use per-exon up_5k / down_5k genomic flanks;
- keep the 70 nt upstream-intron offset within nt_seq;
- run all eight Pangolin output heads, averaged across five seeds;
- realign deletion outputs back to WT coordinates by zero insertion;
- save one parquet file per exon.

The paths below use simple repository-style placeholders. The original internal
Google Drive paths used for the manuscript analyses are kept as comments for
provenance and will be mapped to final supplementary-data filenames later.
"""

from __future__ import annotations

import gc
import os
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from pkg_resources import resource_filename

from pangolin.model import AR, L, W, Pangolin


# ============================================================
# INPUT FILES
# ============================================================

DATA_DIR = Path("data/input")
OUTPUT_DIR = Path("results/pangolin/genome_mode/deletions")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

#
# Variant-level metadata used for deletion rows, WT row, start/length,
# exon_length, and nt_seq.
META_FILE = DATA_DIR / "opensplice_predictors_benchmarking_variant_metadata.tsv"

#
# Per-exon metadata providing genomic upstream/downstream context.
CONTEXT_FILE = DATA_DIR / "opensplice_predictors_benchmarking_exon_metadata.tsv"


# ============================================================
# SETTINGS
# ============================================================

# Match the SpliceAI genome-mode script exactly: 5000 N on each side.
CONTEXT = 10000

# Important coordinate assumption:
# nt_seq in the variant metadata is assumed to begin with the 70 nt upstream
# intronic segment before the middle exon starts.
#
# Therefore:
#   acceptor = len(up_5k) + 70
#   donor    = len(up_5k) + 70 + exon_length - 1
UPSTREAM_INTRON_IN_NTSEQ = 70

# Batch size for Pangolin prediction.
BATCH_SIZE = 250

# Use GPU if available.
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using device: {DEVICE}")


# ============================================================
# PANGOLIN MODEL CONFIG
# ============================================================

# Eight Pangolin model heads:
# 0/1 = heart   (P(splice), usage)
# 2/3 = liver   (P(splice), usage)
# 4/5 = brain   (P(splice), usage)
# 6/7 = testis  (P(splice), usage)
MODEL_NUMS = list(range(8))

# Channel indices in the Pangolin output tensor.
INDEX_MAP = {
    0: 1,   # heart_ps
    1: 2,   # heart_usage
    2: 4,   # liver_ps
    3: 5,   # liver_usage
    4: 7,   # brain_ps
    5: 8,   # brain_usage
    6: 10,  # testis_ps
    7: 11,  # testis_usage
}

# Human-readable names for output columns.
MODEL_LABELS = {
    0: "heart_ps",
    1: "heart_usage",
    2: "liver_ps",
    3: "liver_usage",
    4: "brain_ps",
    5: "brain_usage",
    6: "testis_ps",
    7: "testis_usage",
}

# Pangolin one-hot map.
# Output shape should be (4, L), matching Pangolin expectation.
IN_MAP = np.asarray(
    [
        [0, 0, 0, 0],  # N
        [1, 0, 0, 0],  # A
        [0, 1, 0, 0],  # C
        [0, 0, 1, 0],  # G
        [0, 0, 0, 1],  # T
    ],
    dtype=np.float32,
)


def one_hot_pangolin(seq: str) -> np.ndarray:
    """
    One-hot encode a DNA sequence for Pangolin without reverse complementing.

    The input sequences are assumed to already be oriented consistently for this
    benchmark, so this function does not reverse-complement. The returned array
    has shape (4, L), as expected by Pangolin.
    """
    seq = str(seq).upper().translate(str.maketrans("ACGTN", "12340"))
    arr = np.fromiter((int(x) for x in seq), dtype=np.int8, count=len(seq))
    return IN_MAP[arr].T


def load_pangolin_models(model_nums: list[int], device: str):
    """
    Load Pangolin ensembles.

    Models are loaded fresh per exon, used, and then deleted. This mirrors the
    original memory-saving strategy used during the benchmarking runs.
    """
    all_models = {mn: [] for mn in model_nums}

    for mn in model_nums:
        for seed in range(1, 6):
            model = Pangolin(L, W, AR)
            path = resource_filename("pangolin", f"models/final.{seed}.{mn}.3")
            weights = torch.load(path, map_location=device)
            model.load_state_dict(weights)
            model.to(device)
            model.eval()
            all_models[mn].append(model)

    return all_models


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

# Keep one context row per exon and create a fast lookup by exon ID.
ctx_df = ctx_df.drop_duplicates(subset=["ensembl_exon_id"]).copy()
ctx_map = ctx_df.set_index("ensembl_exon_id")

# Keep deletion rows only for deciding which exons to process.
del_df = meta_df[meta_df["variant_id"].astype(str).str.contains("_del", na=False)].copy()

# Exons to process = exons with deletions.
# Reverse the existing order exactly as in the original run.
exon_ids = del_df["ensembl_exon_id"].dropna().unique().tolist()[::-1]

print(f"Will process {len(exon_ids)} exons with deletion rows")
print("Running in reverse exon order")
if exon_ids:
    print("First 5 exon IDs in this run order:", exon_ids[:5])


# ============================================================
# MAIN LOOP
# ============================================================

for exon_i, exon_id in enumerate(exon_ids, start=1):
    out_file = OUTPUT_DIR / f"{exon_id}_pangolin_scores_genome_mode.parquet"

    if out_file.exists():
        print(f"Skipping [{exon_i}/{len(exon_ids)}] {exon_id}; output already exists")
        continue

    print(f"\nProcessing exon [{exon_i}/{len(exon_ids)}]: {exon_id}")

    # --------------------------------------------------------
    # 1. Load Pangolin ensemble fresh per exon.
    # --------------------------------------------------------
    print("Loading Pangolin models ...")
    models = load_pangolin_models(MODEL_NUMS, DEVICE)
    print("Pangolin models loaded")

    # --------------------------------------------------------
    # 2. Subset this exon in both tables.
    # --------------------------------------------------------
    exon_meta = meta_df[meta_df["ensembl_exon_id"] == exon_id].copy()

    if exon_id not in ctx_map.index:
        print(f"No genomic context row found for {exon_id}; skipping")
        del models
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
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
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        continue

    wt_row = wt_candidates.iloc[0]

    ref_id = wt_row["variant_id"]
    wt_seq = str(wt_row["nt_seq"])
    exon_len = int(wt_row["exon_length"])

    # Safety check against the context table WT sequence. We still use the
    # variant metadata nt_seq because deletion coordinates are defined on it.
    if wt_seq != ctx_wt_seq:
        print(
            f"WT nt_seq mismatch for {exon_id}: "
            f"META_FILE nt_seq length={len(wt_seq)}, "
            f"CONTEXT_FILE wt_seq length={len(ctx_wt_seq)}. "
            f"Using META_FILE nt_seq because deletion coordinates are defined on that sequence."
        )

    # --------------------------------------------------------
    # 4. Define genomic-context construct and canonical site coordinates.
    # --------------------------------------------------------
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
    ref_core_seq = up_5k + wt_seq + down_5k
    ref_seq = "N" * (CONTEXT // 2) + ref_core_seq + "N" * (CONTEXT // 2)

    ref_enc = one_hot_pangolin(ref_seq)[None, :, :]
    ref_tensor = torch.tensor(ref_enc, dtype=torch.float32, device=DEVICE)

    ref_scores = {}
    with torch.no_grad():
        for mn in MODEL_NUMS:
            idx = INDEX_MAP[mn]
            preds = []
            for model in models[mn]:
                out = model(ref_tensor)[0, idx, :].detach().cpu().numpy()
                preds.append(out)
            ref_scores[mn] = np.mean(preds, axis=0)

    max_needed_pos = max(pois)
    for mn in MODEL_NUMS:
        if max_needed_pos >= len(ref_scores[mn]):
            raise ValueError(
                f"{exon_id}: requested splice-site position {max_needed_pos} "
                f"is outside REF Pangolin prediction length {len(ref_scores[mn])} "
                f"for model {MODEL_LABELS[mn]}"
            )

    results = []

    row_ref = {
        "ensembl_exon_id": exon_id,
        "Type": "REF",
        "Identifier": ref_id,
        "up_5k_len": upstream_context_len,
        "acceptor_pos0": acceptor_pos,
        "donor_pos0": donor_pos,
    }

    for mn in MODEL_NUMS:
        label = MODEL_LABELS[mn]
        scores = ref_scores[mn]
        for p in pois:
            row_ref[f"{label}_score_{p}"] = float(scores[p])

    results.append(row_ref)

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
            print(f"  • deletions of length {del_length}")

            for i in range(0, len(grp), BATCH_SIZE):
                batch = grp.iloc[i : i + BATCH_SIZE]

                encs = []
                batch_rows = []

                for _, row in batch.iterrows():
                    alt_nt_seq = str(row["nt_seq"])
                    alt_core_seq = up_5k + alt_nt_seq + down_5k
                    alt_seq = "N" * (CONTEXT // 2) + alt_core_seq + "N" * (CONTEXT // 2)

                    encs.append(one_hot_pangolin(alt_seq))
                    batch_rows.append(row)

                batch_arr = np.stack(encs)
                batch_tensor = torch.tensor(batch_arr, dtype=torch.float32, device=DEVICE)

                batch_scores = {}
                with torch.no_grad():
                    for mn in MODEL_NUMS:
                        idx = INDEX_MAP[mn]
                        preds = []
                        for model in models[mn]:
                            out = model(batch_tensor)[:, idx, :].detach().cpu().numpy()
                            preds.append(out)
                        batch_scores[mn] = np.mean(preds, axis=0)

                for bi, row in enumerate(batch_rows):
                    # start is 1-based within nt_seq.
                    ds = upstream_context_len + int(row["start"]) - 1
                    Ldel = int(row["length"])

                    result_row = {
                        "ensembl_exon_id": exon_id,
                        "Type": "ALT",
                        "Identifier": row["variant_id"],
                        "up_5k_len": upstream_context_len,
                        "deletion_start_insert_pos0": ds,
                        "deletion_length": Ldel,
                        "acceptor_pos0": acceptor_pos,
                        "donor_pos0": donor_pos,
                    }

                    for mn in MODEL_NUMS:
                        label = MODEL_LABELS[mn]
                        scr = batch_scores[mn][bi, :]

                        # Realign ALT prediction back to WT coordinates by
                        # inserting zeros at the deletion start.
                        scr_adj = np.insert(scr, ds, [0] * Ldel)

                        if max_needed_pos >= len(scr_adj):
                            raise ValueError(
                                f"{exon_id} / {row['variant_id']}: requested position "
                                f"{max_needed_pos} is outside realigned ALT length {len(scr_adj)} "
                                f"for model {label}"
                            )

                        for p in pois:
                            result_row[f"{label}_score_{p}"] = float(scr_adj[p])

                    results.append(result_row)

                del batch_arr, batch_tensor, batch_scores, encs, batch_rows
                gc.collect()
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()

    # --------------------------------------------------------
    # 7. Save per-exon parquet.
    # --------------------------------------------------------
    res_df = pd.DataFrame(results)
    res_df.to_parquet(out_file, engine="pyarrow", compression="snappy")

    print(f"Saved {len(res_df)} rows -> {out_file}")

    # --------------------------------------------------------
    # 8. Cleanup model memory.
    # --------------------------------------------------------
    del models, ref_enc, ref_tensor, ref_scores, results, res_df
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
