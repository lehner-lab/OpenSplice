#@title Pangolin SNVs in genomic-context mode — use up_5k/down_5k flanks per exon, with optional reverse exon order
from __future__ import annotations

import os
import gc
from pathlib import Path
import numpy as np
import pandas as pd
import torch

from pkg_resources import resource_filename
from pangolin.model import Pangolin, L, W, AR

# ============================================================
# INPUT FILES
# ============================================================

# Variant-level metadata used for SNV rows, WT row, exon_length, nt_seq
DATA_DIR = Path("data/input")
meta_file = DATA_DIR / "opensplice_predictors_benchmarking_variant_metadata.tsv"

# Per-exon metadata providing genomic upstream/downstream context
context_file = DATA_DIR / "opensplice_predictors_benchmarking_exon_metadata.tsv"

# Output directory for Pangolin genomic-context mode SNVs
output_dir = Path("results/pangolin/genome_mode/snvs")
output_dir.mkdir(parents=True, exist_ok=True)

# ============================================================
# SETTINGS
# ============================================================

# Match your genomic-mode deletion script:
# add 5000 N on each side of the full genomic-context construct
context = 10000

# Important coordinate assumption:
# nt_seq in the variant metadata is assumed to begin with the 70 nt upstream
# intronic segment before the middle exon starts.
#
# Therefore, within:
#   up_5k + nt_seq + down_5k
#
# the canonical middle-exon splice sites are at:
#   acceptor = len(up_5k) + 70
#   donor    = len(up_5k) + 70 + exon_length - 1
UPSTREAM_INTRON_IN_NTSEQ = 70

# Batch size for Pangolin prediction
BATCH_SIZE = 1000

# Optional: reverse the order of exon processing
# False = normal order
# True  = reverse order
REVERSE_EXON_ORDER = False

# Optional: skip outputs that already exist
SKIP_EXISTING = True

# Use GPU if available
device = "cuda" if torch.cuda.is_available() else "cpu"
print(f"🖥 Using device: {device}")

# ============================================================
# PANGOLIN MODEL CONFIG
# ============================================================

# 8 Pangolin model heads:
# 0/1 = heart   (P(splice), usage)
# 2/3 = liver   (P(splice), usage)
# 4/5 = brain   (P(splice), usage)
# 6/7 = testis  (P(splice), usage)
model_nums = list(range(8))

# Channel indices in Pangolin output tensor
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

# Human-readable names for output columns
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

# Pangolin one-hot map
# Output shape should be (4, L), matching Pangolin expectation
IN_MAP = np.asarray([
    [0, 0, 0, 0],  # N
    [1, 0, 0, 0],  # A
    [0, 1, 0, 0],  # C
    [0, 0, 1, 0],  # G
    [0, 0, 0, 1],  # T
], dtype=np.float32)

def one_hot_pangolin(seq: str) -> np.ndarray:
    """
    One-hot encode a DNA sequence for Pangolin without reverse complementing.

    Important:
    - Your sequences are already oriented consistently for the model input you use.
    - So, as in your other Pangolin scripts, we do NOT reverse complement.
    - Output shape is (4, L).
    """
    seq = str(seq).upper().translate(str.maketrans("ACGTN", "12340"))
    arr = np.fromiter((int(x) for x in seq), dtype=np.int8, count=len(seq))
    return IN_MAP[arr].T  # shape: (4, L)

def load_pangolin_models(model_nums: list[int], device: str):
    """
    Load all Pangolin ensemble members fresh.

    This mirrors your memory-saving pattern:
    load per exon, run predictions, then free memory.
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

meta_df = pd.read_csv(meta_file, sep="\t")
ctx_df = pd.read_csv(context_file, sep="\t")

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
    raise ValueError(f"meta_file is missing required columns: {missing_meta}")

required_ctx_cols = [
    "ensembl_exon_id",
    "up_5k",
    "wt_seq",
    "down_5k",
]
missing_ctx = [c for c in required_ctx_cols if c not in ctx_df.columns]
if missing_ctx:
    raise ValueError(f"context_file is missing required columns: {missing_ctx}")

# Keep one context row per exon
ctx_df = ctx_df.drop_duplicates(subset=["ensembl_exon_id"]).copy()

# Fast lookup by exon ID
ctx_map = ctx_df.set_index("ensembl_exon_id")

# ============================================================
# IDENTIFY SNV ROWS / EXONS TO PROCESS
# ============================================================

# Keep WT rows
wt_mask = meta_df["variant_id"].astype(str).str.endswith("_wt", na=False)

# Keep SNV-like rows by excluding WT and deletion rows.
# This is the safest generic rule if your SNV naming is not perfectly uniform.
del_mask = meta_df["variant_id"].astype(str).str.contains("_del", na=False)

snv_df = meta_df.loc[~wt_mask & ~del_mask].copy()

# Exons to process = exons with at least one SNV row
exon_ids = snv_df["ensembl_exon_id"].dropna().unique().tolist()

if REVERSE_EXON_ORDER:
    exon_ids = exon_ids[::-1]

print(f"🔬 Will process {len(exon_ids)} exons with SNV rows")
print(f"↕️ Reverse exon order: {REVERSE_EXON_ORDER}")

# ============================================================
# MAIN LOOP
# ============================================================

for exon_id in exon_ids:
    out_file = output_dir / f"{exon_id}_pangolin_scores_genome_mode.parquet"

    if SKIP_EXISTING and out_file.exists():
        print(f"⏭️ Skipping {exon_id}, output already exists")
        continue

    print(f"\nProcessing exon: {exon_id}")

    # --------------------------------------------------------
    # 1. Load Pangolin ensemble fresh per exon
    # --------------------------------------------------------
    print("📦 Loading Pangolin models ...")
    models = load_pangolin_models(model_nums, device)
    print("✅ Pangolin models loaded")

    # --------------------------------------------------------
    # 2. Subset this exon in both tables
    # --------------------------------------------------------
    exon_meta = meta_df[meta_df["ensembl_exon_id"] == exon_id].copy()

    if exon_id not in ctx_map.index:
        print(f"⚠️ No genomic context row found for {exon_id}; skipping")
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
    # 3. Find WT row for this exon
    # --------------------------------------------------------
    wt_candidates = exon_meta[
        exon_meta["variant_id"].astype(str).str.endswith("_wt", na=False)
    ].copy()

    if wt_candidates.empty:
        print(f"⚠️ No WT row found in meta_file for {exon_id}; skipping")
        del models
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        continue

    wt_row = wt_candidates.iloc[0]

    ref_id = wt_row["variant_id"]
    wt_seq = str(wt_row["nt_seq"])
    exon_len = int(wt_row["exon_length"])

    # Optional safety check against context table WT sequence
    if wt_seq != ctx_wt_seq:
        print(
            f"⚠️ WT nt_seq mismatch for {exon_id}: "
            f"meta_file nt_seq length={len(wt_seq)}, "
            f"context_file wt_seq length={len(ctx_wt_seq)}. "
            f"Using meta_file nt_seq because the variant nt_seq rows are defined relative to that sequence."
        )

    # --------------------------------------------------------
    # 4. Define genomic-context construct and splice-site coords
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
    # 5. REF (WT) prediction
    # --------------------------------------------------------
    #
    # Build:
    #   N-padding + up_5k + wt_seq + down_5k + N-padding
    #
    ref_core_seq = up_5k + wt_seq + down_5k
    ref_seq = "N" * (context // 2) + ref_core_seq + "N" * (context // 2)

    ref_enc = one_hot_pangolin(ref_seq)[None, :, :]  # shape: (1, 4, L)
    ref_tensor = torch.tensor(ref_enc, dtype=torch.float32, device=device)

    ref_scores = {}
    with torch.no_grad():
        for mn in model_nums:
            idx = INDEX_MAP[mn]
            preds = []
            for model in models[mn]:
                out = model(ref_tensor)[0, idx, :].detach().cpu().numpy()
                preds.append(out)
            ref_scores[mn] = np.mean(preds, axis=0)

    max_needed_pos = max(pois)
    for mn in model_nums:
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

    for mn in model_nums:
        label = MODEL_LABELS[mn]
        s = ref_scores[mn]
        for p in pois:
            row_ref[f"{label}_score_{p}"] = float(s[p])

    results.append(row_ref)

    # --------------------------------------------------------
    # 6. ALT SNVs
    # --------------------------------------------------------
    #
    # For SNVs, no realignment is needed.
    # ALT construct length matches WT construct length.
    # So we simply read the same canonical positions directly.
    #
    alt_meta = exon_meta.loc[
        ~exon_meta["variant_id"].astype(str).str.endswith("_wt", na=False)
        & ~exon_meta["variant_id"].astype(str).str.contains("_del", na=False)
    ].copy()

    if alt_meta.empty:
        print("  No SNV rows for this exon; saving REF only")
    else:
        print(f"  Number of SNV rows: {len(alt_meta)}")

        for i in range(0, len(alt_meta), BATCH_SIZE):
            batch = alt_meta.iloc[i:i + BATCH_SIZE]

            encs = []
            batch_rows = []

            # Build ALT sequences in genomic-context mode
            for _, row in batch.iterrows():
                alt_nt_seq = str(row["nt_seq"])
                alt_core_seq = up_5k + alt_nt_seq + down_5k
                alt_seq = "N" * (context // 2) + alt_core_seq + "N" * (context // 2)

                encs.append(one_hot_pangolin(alt_seq))
                batch_rows.append(row)

            batch_arr = np.stack(encs)  # shape: (B, 4, L)
            batch_tensor = torch.tensor(batch_arr, dtype=torch.float32, device=device)

            # Predict for all Pangolin models
            batch_scores = {}
            with torch.no_grad():
                for mn in model_nums:
                    idx = INDEX_MAP[mn]
                    preds = []
                    for model in models[mn]:
                        out = model(batch_tensor)[:, idx, :].detach().cpu().numpy()
                        preds.append(out)
                    batch_scores[mn] = np.mean(preds, axis=0)  # shape: (B, Lout)

            # Unpack batch
            for bi, row in enumerate(batch_rows):
                result_row = {
                    "ensembl_exon_id": exon_id,
                    "Type": "ALT",
                    "Identifier": row["variant_id"],
                    "up_5k_len": upstream_context_len,
                    "acceptor_pos0": acceptor_pos,
                    "donor_pos0": donor_pos,
                }

                for mn in model_nums:
                    label = MODEL_LABELS[mn]
                    scr = batch_scores[mn][bi, :]  # shape: (Lout,)

                    if max_needed_pos >= len(scr):
                        raise ValueError(
                            f"{exon_id} / {row['variant_id']}: requested position "
                            f"{max_needed_pos} is outside ALT prediction length {len(scr)} "
                            f"for model {label}"
                        )

                    for p in pois:
                        result_row[f"{label}_score_{p}"] = float(scr[p])

                results.append(result_row)

            del batch_arr, batch_tensor, batch_scores, encs, batch_rows
            gc.collect()
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

    # --------------------------------------------------------
    # 7. Save per-exon parquet
    # --------------------------------------------------------
    res_df = pd.DataFrame(results)
    res_df.to_parquet(out_file, engine="pyarrow", compression="snappy")

    print(f"✅ Saved {len(res_df)} rows → {out_file}")

    # --------------------------------------------------------
    # 8. Cleanup model memory
    # --------------------------------------------------------
    del models, ref_enc, ref_tensor, ref_scores, results, res_df
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
