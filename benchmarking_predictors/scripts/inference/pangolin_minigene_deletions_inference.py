#@title Run Pangolin deletion mutagenesis in minigene-context mode
"""
Run Pangolin deletion-mutagenesis inference in minigene-context mode.

This script scores deletion variants for each exon using the original minigene
construct context:

    N-padding + pre_manual + variant nt_seq + post_manual + N-padding

where:
    pre_manual  = FAS exon 5 + FAS intron 5
    post_manual = FAS intron 6 + FAS exon 7

It extracts Pangolin scores at the canonical middle-exon splice sites using the
original minigene coordinate convention:

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
import torch
from pangolin.model import AR, L, W, Pangolin
from pkg_resources import resource_filename


# ============================================================
# INPUT FILES
# ============================================================

DATA_DIR = Path("data/input")
OUTPUT_DIR = Path("results/pangolin/minigene_mode/deletions")
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

# Pangolin input context: 5000 N on each side.
CONTEXT = 10000

# Batch size used in the original Pangolin minigene deletion run.
BATCH_SIZE = 1000

# Optional: run only the first exon for debugging.
TEST_MODE = False

# Optional: skip outputs that already exist.
SKIP_EXISTING = True

# Canonical acceptor coordinate in the minigene construct.
# Original manuscript-analysis convention:
#   146 nt manual prefix + 70 nt upstream intron = 216
ACC_POS = 216

# Deletion realignment offset in minigene mode.
# Original manuscript-analysis convention:
#   deletion_start_insert_pos0 = 146 + start - 1
PRE_MANUAL_LEN_FOR_DELETION_REALIGNMENT = 146

# Use GPU if available.
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
print(f"Using device: {DEVICE}")


# ============================================================
# MANUAL MINIGENE FLANKING SEQUENCES
# ============================================================

# These FAS minigene backbone sequences are kept exactly as in the original
# notebook. They define the minigene context used for Pangolin inference.
FAS_E5 = "ATGTGAACATGGAATCATCAAGGAATGCACACTCACCAGCAACACCAAGTGCAAAGAGGAAG"
FAS_I5 = "GTAATTATTTTTTTACGGTTATATTCTCCTTTCCCCCAACCCCATGGAAAGATGTGAAGAAAAACCAATCACTCTTGATTACTA"
FAS_I6 = "CAGATTGAAATAACTTGGGAAGTAGTTTCTCTTAGTGTGAAAGTATGTTCTCACATGCATTCTACAAGGCTGAGACCTGAGTTGATAAAATTTCTTTGTTCTTTCAG"
FAS_E7 = "TGAAGAGAAAGGAAGTACAGAAAACATGCAGAAAGCACAGAAAGGAA"

PRE_MANUAL = FAS_E5 + FAS_I5
POST_MANUAL = FAS_I6 + FAS_E7


# ============================================================
# PANGOLIN MODEL CONFIGURATION
# ============================================================

# Eight Pangolin model heads:
#   0/1 = heart  (P(splice), usage)
#   2/3 = liver  (P(splice), usage)
#   4/5 = brain  (P(splice), usage)
#   6/7 = testis (P(splice), usage)
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

# Human-readable names used in output columns.
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
# Output shape from one_hot_pangolin is (4, sequence_length).
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
    One-hot encode a DNA sequence for Pangolin.

    Important:
    - Sequences are already oriented consistently for this analysis.
    - No reverse complementation is applied.
    - Returned array has shape (4, sequence_length).
    """
    seq = str(seq).upper().translate(str.maketrans("ACGTN", "12340"))
    arr = np.fromiter((int(x) for x in seq), dtype=np.int8, count=len(seq))
    return IN_MAP[arr].T


def load_pangolin_models(model_nums: list[int], device: str):
    """
    Load the Pangolin ensemble members.

    This minigene deletion script keeps the original model-loading behaviour:
    all Pangolin models are loaded once before the exon loop.
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
# LOAD MODELS AND INPUT TABLE
# ============================================================

print("Loading Pangolin models")
models = load_pangolin_models(MODEL_NUMS, DEVICE)
print("Pangolin models loaded")

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
if TEST_MODE:
    exon_ids = exon_ids[:1]

print(f"Will process {len(exon_ids)} exons in Pangolin minigene deletion mode")


# ============================================================
# MAIN LOOP
# ============================================================

for exon_id in exon_ids:
    out_file = OUTPUT_DIR / f"{exon_id}_pangolin_scores_minigene.parquet"

    if SKIP_EXISTING and out_file.exists():
        print(f"Skipping {exon_id}, output already exists: {out_file}")
        continue

    print(f"\nProcessing exon: {exon_id}")

    # --------------------------------------------------------
    # 1. Subset this exon.
    # --------------------------------------------------------
    exon_df = meta_df[meta_df["ensembl_exon_id"] == exon_id].copy()

    # --------------------------------------------------------
    # 2. Find WT row and canonical minigene coordinates.
    # --------------------------------------------------------
    wt_candidates = exon_df[
        exon_df["variant_id"].astype(str).str.endswith("_wt", na=False)
    ].copy()

    if wt_candidates.empty:
        print(f"No WT row found for {exon_id}; skipping")
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
    # 3. REF / WT prediction.
    # --------------------------------------------------------
    ref_seq = (
        "N" * (CONTEXT // 2)
        + PRE_MANUAL
        + wt_seq
        + POST_MANUAL
        + "N" * (CONTEXT // 2)
    )

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
        "acceptor_pos0": ACC_POS,
        "donor_pos0": donor_pos,
    }

    for mn in MODEL_NUMS:
        label = MODEL_LABELS[mn]
        scores = ref_scores[mn]

        for p in pois:
            row_ref[f"{label}_score_{p}"] = float(scores[p])

    results.append(row_ref)

    # --------------------------------------------------------
    # 4. ALT deletion predictions.
    # --------------------------------------------------------
    alt_df = exon_df[
        exon_df["variant_id"].astype(str).str.contains("_del", na=False)
    ].copy()

    for del_len, grp in alt_df.groupby("length"):
        print(f"  deletions of length {del_len}, n={len(grp)}")

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
                # Realign ALT predictions back to WT coordinates.
                # This keeps the original minigene-mode coordinate rule:
                #   ds = 146 + start - 1
                ds = PRE_MANUAL_LEN_FOR_DELETION_REALIGNMENT + int(row["start"]) - 1
                deletion_length = int(row["length"])

                result_row = {
                    "ensembl_exon_id": exon_id,
                    "Type": "ALT",
                    "Identifier": row["variant_id"],
                    "deletion_start_insert_pos0": ds,
                    "deletion_length": deletion_length,
                    "acceptor_pos0": ACC_POS,
                    "donor_pos0": donor_pos,
                }

                for mn in MODEL_NUMS:
                    label = MODEL_LABELS[mn]
                    scores = batch_scores[mn][bi, :]

                    scores_adjusted = np.insert(scores, ds, [0] * deletion_length)

                    if max_needed_pos >= len(scores_adjusted):
                        raise ValueError(
                            f"{exon_id} / {row['variant_id']}: requested position "
                            f"{max_needed_pos} is outside realigned ALT length "
                            f"{len(scores_adjusted)} for model {label}"
                        )

                    for p in pois:
                        result_row[f"{label}_score_{p}"] = float(scores_adjusted[p])

                results.append(result_row)

            del batch_arr, batch_tensor, batch_scores, encs, batch_rows
            gc.collect()
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

    # --------------------------------------------------------
    # 5. Save per-exon parquet.
    # --------------------------------------------------------
    out_df = pd.DataFrame(results)
    out_df.to_parquet(out_file, engine="pyarrow", compression="snappy")

    print(f"Saved {len(out_df)} rows to {out_file}")

    del ref_enc, ref_tensor, ref_scores, results, out_df
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
