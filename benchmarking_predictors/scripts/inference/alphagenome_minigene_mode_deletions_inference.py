#!/usr/bin/env python
#@title Run AlphaGenome minigene-mode deletion inference
"""
Purpose
-------
Run AlphaGenome Research locally on custom minigene WT + deletion constructs and extract
the splice-site scores at the intended middle-exon splice sites.

Scoring rule
------------
For every construct, score:
- acceptor at construct_pos0 = 216
- donor at construct_pos0 = 216 + exon_length - 1

Deletion alignment rule
-----------------------
For deletions, ALT construct length is shorter than WT.

So after predicting ALT on the shorter construct, we:
1) extract the ALT donor/acceptor track across the construct only
2) insert zeros into the ALT construct track at the deletion start position
   with length = deletion length
3) this makes the ALT track aligned back to WT construct coordinates

This mirrors your Pangolin / SpliceAI deletion realignment idea.

Important coordinate note
-------------------------
The deletion start used for zero-insertion is:

    deletion_start_construct_pos0 = 146 + start - 1

This matches your previous Pangolin deletion workflow, where `start` is the
1-based position inside the variable `nt_seq` region that begins at construct_pos0 = 146.

Output structure
----------------
One output file per exon:
- one REF row for WT
- one ALT row per deletion

Columns include:
- alphagenome_acceptor_score
- alphagenome_donor_score
- aligned ALT scores after zero-insertion
- deletion metadata

Assumptions
-----------
1) `nt_seq` is already in the assay/transcript orientation.
2) Deletion rows are identified by '_del' in variant_id.
3) AlphaGenome model object `model` has already been loaded in a previous cell.
"""

import os
import gc
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from alphagenome.models import dna_client
from alphagenome.models.dna_output import OutputType

# ──────────────────────────────────────────────────────────────
# CONFIGURATION
# ──────────────────────────────────────────────────────────────

meta_file = "data/input/opensplice_predictors_benchmarking_variant_metadata.tsv"

output_dir = (
    "results/alphagenome/"
    "custom_minigene_deletions_middle_exon_sites_per_exon"
)
os.makedirs(output_dir, exist_ok=True)

TEST_MODE = False
TEST_N_EXONS = 2

# Same manual flanks as your other minigene workflows
fas_e5 = "ATGTGAACATGGAATCATCAAGGAATGCACACTCACCAGCAACACCAAGTGCAAAGAGGAAG"
fas_i5 = "GTAATTATTTTTTTACGGTTATATTCTCCTTTCCCCCAACCCCATGGAAAGATGTGAAGAAAAACCAATCACTCTTGATTACTA"
fas_i6 = "CAGATTGAAATAACTTGGGAAGTAGTTTCTCTTAGTGTGAAAGTATGTTCTCACATGCATTCTACAAGGCTGAGACCTGAGTTGATAAAATTTCTTTGTTCTTTCAG"
fas_e7 = "TGAAGAGAAAGGAAGTACAGAAAACATGCAGAAAGCACAGAAAGGAA"

pre_manual = fas_e5 + fas_i5
post_manual = fas_i6 + fas_e7

# Fixed intended middle-exon acceptor coordinate in construct space
MIDDLE_EXON_ACCEPTOR_CONSTRUCT_POS0 = 216

# Smallest supported AlphaGenome input length
TARGET_LEN = 2**14  # 16384

ontology_terms = ["CL:0002518"]
# kidney epithelial cell (closest available to HEK); note ontology terms mostly
# affect splice-site usage/junction outputs, while this workflow uses splice-site probability.

print("meta_file:", meta_file)
print("output_dir:", output_dir)
print("TARGET_LEN:", TARGET_LEN)
print("middle exon acceptor construct_pos0:", MIDDLE_EXON_ACCEPTOR_CONSTRUCT_POS0)

# AlphaGenome Research offline model.
# Choose the backend that matches your local setup.
# For Kaggle-based weights, replace create_from_huggingface with create_from_kaggle.
from alphagenome_research.model import dna_model

model = dna_model.create_from_huggingface("all_folds")

# ──────────────────────────────────────────────────────────────
# HELPERS
# ──────────────────────────────────────────────────────────────

def clean_seq(seq: str) -> str:
    """
    Upper-case and replace any non-ACGTN character with N.
    """
    seq = str(seq).upper()
    return "".join(base if base in {"A", "C", "G", "T", "N"} else "N" for base in seq)

def build_padded_construct(core_seq: str, target_len: int = TARGET_LEN, pad_base: str = "N"):
    """
    Center a core construct inside a fixed AlphaGenome-supported input length.

    Returns
    -------
    padded_seq : str
        Final padded sequence.
    left_pad_len : int
        Number of N bases added on the left.
    right_pad_len : int
        Number of N bases added on the right.
    """
    core_seq = clean_seq(core_seq)

    if len(core_seq) > target_len:
        raise ValueError(
            f"Core construct length {len(core_seq)} exceeds target length {target_len}."
        )

    left_pad_len = (target_len - len(core_seq)) // 2
    right_pad_len = target_len - len(core_seq) - left_pad_len
    padded_seq = (pad_base * left_pad_len) + core_seq + (pad_base * right_pad_len)

    return padded_seq, left_pad_len, right_pad_len

def get_splice_sites_obj(pred_out):
    """
    Defensive accessor for AlphaGenome splice-site output.
    """
    if hasattr(pred_out, "splice_sites"):
        return pred_out.splice_sites

    for attr in dir(pred_out):
        if attr.startswith("_"):
            continue
        try:
            obj = getattr(pred_out, attr)
        except Exception:
            continue
        if hasattr(obj, "splice_sites"):
            return obj.splice_sites

    raise AttributeError("Could not find splice_sites in AlphaGenome output.")

def donor_acceptor_positive_strand(splice_sites):
    """
    Extract donor and acceptor arrays from positive-strand splice-site tracks.

    Assumption:
    the minigene sequences are already oriented in the assay/transcript direction.
    """
    ss = splice_sites.filter_to_positive_strand()
    arr = np.asarray(ss.values)
    names = [str(n).lower() for n in ss.names]

    donor_idx = names.index("donor")
    acceptor_idx = names.index("acceptor")

    donor = arr[:, donor_idx]
    acceptor = arr[:, acceptor_idx]
    return donor, acceptor

def safe_score(track: np.ndarray, pos0):
    """
    Return track[pos0] if valid, else NaN.
    """
    if pos0 is None:
        return np.nan
    pos0 = int(pos0)
    if 0 <= pos0 < len(track):
        return float(track[pos0])
    return np.nan

def is_del_variant_id(variant_id: str) -> bool:
    """
    Match deletion rows using the same spirit as your earlier workflows.
    """
    return "_del" in str(variant_id)

def extract_construct_track(full_track: np.ndarray, left_pad_len: int, core_len: int):
    """
    Extract the model output track across the construct only (exclude outer N padding).
    """
    start = int(left_pad_len)
    end = start + int(core_len)
    return np.asarray(full_track[start:end], dtype=float)

def realign_deletion_track_by_zero_insertion(alt_construct_track: np.ndarray, del_start_construct_pos0: int, del_len: int):
    """
    Realign a deletion ALT construct track back onto WT construct coordinates by inserting zeros.

    Parameters
    ----------
    alt_construct_track : np.ndarray
        ALT track across the shorter construct only.
    del_start_construct_pos0 : int
        Deletion start in WT construct coordinates.
    del_len : int
        Number of deleted bases.

    Returns
    -------
    aligned_track : np.ndarray
        ALT track with zeros inserted at the deleted span, so it matches WT construct length.
    """
    return np.insert(
        alt_construct_track,
        int(del_start_construct_pos0),
        np.zeros(int(del_len), dtype=float)
    )

# ──────────────────────────────────────────────────────────────
# LOAD METADATA
# ──────────────────────────────────────────────────────────────

meta_df = pd.read_csv(meta_file, sep="\t")

wt_mask = meta_df["variant_id"].astype(str).str.endswith("_wt", na=False)
del_mask = meta_df["variant_id"].astype(str).apply(is_del_variant_id)

work_df = meta_df[wt_mask | del_mask].copy()

exon_ids = (
    work_df.loc[del_mask, "ensembl_exon_id"]
    .dropna()
    .astype(str)
    .unique()
    .tolist()
)

if TEST_MODE:
    exon_ids = exon_ids[:TEST_N_EXONS]

print(f"Will process {len(exon_ids)} exon(s) (TEST_MODE={TEST_MODE})")

# ──────────────────────────────────────────────────────────────
# MAIN EXON LOOP
# ──────────────────────────────────────────────────────────────

for exon_id in tqdm(exon_ids, desc="Exons"):
    out_file = os.path.join(
        output_dir,
        f"{exon_id}_alphagenome_scores_minigene_dels_middle_exon_sites.tsv"
    )

    # Resume / skip completed exons
    if os.path.exists(out_file) and os.path.getsize(out_file) > 0:
        print(f"⏭️ Skipping {exon_id}, output exists")
        continue

    print(f"\n🧬 Processing exon: {exon_id}")

    exon_df = work_df[work_df["ensembl_exon_id"].astype(str) == str(exon_id)].copy()
    if exon_df.empty:
        print(f"[{exon_id}] no rows -> skipping")
        continue

    # ─────────────────────────
    # WT ROW
    # ─────────────────────────
    wt_rows = exon_df[exon_df["variant_id"].astype(str).str.endswith("_wt", na=False)]
    if len(wt_rows) == 0:
        print(f"[{exon_id}] no WT row found -> skipping")
        continue
    if len(wt_rows) > 1:
        print(f"[{exon_id}] warning: multiple WT rows found, taking first")

    wt_row = wt_rows.iloc[0]
    wt_identifier = str(wt_row["variant_id"])
    wt_seq = clean_seq(wt_row["nt_seq"])
    exon_len = int(wt_row["exon_length"])

    # Intended middle-exon site coordinates in WT construct space
    acceptor_construct_pos0 = MIDDLE_EXON_ACCEPTOR_CONSTRUCT_POS0
    donor_construct_pos0 = MIDDLE_EXON_ACCEPTOR_CONSTRUCT_POS0 + exon_len - 1

    # WT construct and padded sequence
    wt_core_construct = clean_seq(pre_manual + wt_seq + post_manual)
    wt_core_len = len(wt_core_construct)

    wt_padded_seq, wt_left_pad_len, wt_right_pad_len = build_padded_construct(
        wt_core_construct,
        target_len=TARGET_LEN
    )

    # WT prediction once per exon
    wt_out = model.predict_sequence(
        sequence=wt_padded_seq,
        organism=dna_client.Organism.HOMO_SAPIENS,
        requested_outputs=[OutputType.SPLICE_SITES],
        ontology_terms=ontology_terms,
        interval=None,
    )

    wt_ss = get_splice_sites_obj(wt_out)
    wt_donor_full, wt_acceptor_full = donor_acceptor_positive_strand(wt_ss)

    # WT construct-only tracks
    wt_donor_construct = extract_construct_track(wt_donor_full, wt_left_pad_len, wt_core_len)
    wt_acceptor_construct = extract_construct_track(wt_acceptor_full, wt_left_pad_len, wt_core_len)

    results = []

    ref_row = {
        "ensembl_exon_id": exon_id,
        "Type": "REF",
        "Identifier": wt_identifier,

        # Geometry
        "exon_len": exon_len,
        "wt_core_len": wt_core_len,
        "target_len": TARGET_LEN,
        "wt_left_pad_len": wt_left_pad_len,
        "wt_right_pad_len": wt_right_pad_len,
        "pre_manual_len": len(pre_manual),
        "post_manual_len": len(post_manual),

        # Fixed intended site coordinates in WT construct space
        "acceptor_construct_pos0": acceptor_construct_pos0,
        "donor_construct_pos0": donor_construct_pos0,

        # Scores at intended middle-exon sites
        "alphagenome_acceptor_score": safe_score(wt_acceptor_construct, acceptor_construct_pos0),
        "alphagenome_donor_score": safe_score(wt_donor_construct, donor_construct_pos0),

        # ALT / deletion-specific metadata blank for REF
        "alt_core_len": np.nan,
        "alt_left_pad_len": np.nan,
        "alt_right_pad_len": np.nan,
        "deletion_start_in_nt_seq_1based": np.nan,
        "deletion_len_bp": np.nan,
        "deletion_start_construct_pos0": np.nan,
    }
    results.append(ref_row)

    # ─────────────────────────
    # ALT DELETION ROWS
    # ─────────────────────────
    del_rows = exon_df[exon_df["variant_id"].astype(str).apply(is_del_variant_id)].copy()

    if len(del_rows) > 0:
        print(f" • Deletions: n={len(del_rows)}")
    else:
        print(f" • No deletions found -> saving REF-only file")

    for _, row in tqdm(del_rows.iterrows(), total=len(del_rows), desc=f"{exon_id} DELs", leave=False):
        alt_identifier = str(row["variant_id"])
        alt_seq = clean_seq(row["nt_seq"])

        # Required deletion metadata from your input table
        if "start" not in row.index:
            raise KeyError(
                f"'start' column required for deletions but missing for exon {exon_id}"
            )
        if "length" not in row.index:
            raise KeyError(
                f"'length' column required for deletions but missing for exon {exon_id}"
            )

        deletion_start_in_nt_seq_1based = int(row["start"])
        deletion_len_bp = int(row["length"])

        # Same deletion-start logic as your Pangolin workflow
        # nt_seq begins at construct_pos0 = 146
        deletion_start_construct_pos0 = 146 + deletion_start_in_nt_seq_1based - 1

        # ALT construct and padded sequence
        alt_core_construct = clean_seq(pre_manual + alt_seq + post_manual)
        alt_core_len = len(alt_core_construct)

        alt_padded_seq, alt_left_pad_len, alt_right_pad_len = build_padded_construct(
            alt_core_construct,
            target_len=TARGET_LEN
        )

        # Predict ALT sequence
        alt_out = model.predict_sequence(
            sequence=alt_padded_seq,
            organism=dna_client.Organism.HOMO_SAPIENS,
            requested_outputs=[OutputType.SPLICE_SITES],
            ontology_terms=ontology_terms,
            interval=None,
        )

        alt_ss = get_splice_sites_obj(alt_out)
        alt_donor_full, alt_acceptor_full = donor_acceptor_positive_strand(alt_ss)

        # ALT construct-only tracks on shorter construct coordinates
        alt_donor_construct = extract_construct_track(alt_donor_full, alt_left_pad_len, alt_core_len)
        alt_acceptor_construct = extract_construct_track(alt_acceptor_full, alt_left_pad_len, alt_core_len)

        # Realign ALT back to WT construct coordinates by inserting zeros at deleted span
        alt_donor_aligned = realign_deletion_track_by_zero_insertion(
            alt_construct_track=alt_donor_construct,
            del_start_construct_pos0=deletion_start_construct_pos0,
            del_len=deletion_len_bp,
        )
        alt_acceptor_aligned = realign_deletion_track_by_zero_insertion(
            alt_construct_track=alt_acceptor_construct,
            del_start_construct_pos0=deletion_start_construct_pos0,
            del_len=deletion_len_bp,
        )

        # Sanity check: after alignment ALT should match WT construct length
        if len(alt_donor_aligned) != wt_core_len or len(alt_acceptor_aligned) != wt_core_len:
            raise ValueError(
                f"Aligned ALT track length does not match WT construct length for {alt_identifier}. "
                f"WT={wt_core_len}, ALT_donor_aligned={len(alt_donor_aligned)}, "
                f"ALT_acceptor_aligned={len(alt_acceptor_aligned)}"
            )

        alt_row = {
            "ensembl_exon_id": exon_id,
            "Type": "ALT",
            "Identifier": alt_identifier,

            # Geometry
            "exon_len": exon_len,
            "wt_core_len": wt_core_len,
            "alt_core_len": alt_core_len,
            "target_len": TARGET_LEN,
            "wt_left_pad_len": wt_left_pad_len,
            "wt_right_pad_len": wt_right_pad_len,
            "alt_left_pad_len": alt_left_pad_len,
            "alt_right_pad_len": alt_right_pad_len,
            "pre_manual_len": len(pre_manual),
            "post_manual_len": len(post_manual),

            # Fixed intended WT construct coordinates
            "acceptor_construct_pos0": acceptor_construct_pos0,
            "donor_construct_pos0": donor_construct_pos0,

            # Deletion metadata
            "deletion_start_in_nt_seq_1based": deletion_start_in_nt_seq_1based,
            "deletion_len_bp": deletion_len_bp,
            "deletion_start_construct_pos0": deletion_start_construct_pos0,

            # Scores at intended middle-exon sites after ALT realignment
            "alphagenome_acceptor_score": safe_score(alt_acceptor_aligned, acceptor_construct_pos0),
            "alphagenome_donor_score": safe_score(alt_donor_aligned, donor_construct_pos0),
        }

        results.append(alt_row)
        gc.collect()

    # ─────────────────────────
    # SAVE ONE FILE FOR THIS EXON
    # ─────────────────────────
    out_df = pd.DataFrame(results)
    out_df.to_csv(out_file, sep="\t", index=False)

    print(f"✅ Saved {len(out_df)} rows → {out_file}")

print("\nDONE.")
