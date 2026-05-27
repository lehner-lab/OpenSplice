#!/usr/bin/env python
#@title Run AlphaGenome minigene-mode SNV inference
"""
Purpose
-------
Run AlphaGenome Research locally on custom minigene WT + SNV constructs and extract
the splice-site scores only at the intended middle-exon splice sites.

Scoring rule
------------
For every construct:
- acceptor score is read at construct_pos0 = 216
- donor score is read at construct_pos0 = 216 + exon_length - 1

Why this works
--------------
In your minigene layout, the intended middle exon starts at construct_pos0 = 216.
So the first nt of the middle exon is always the intended acceptor position, and
the last nt of the middle exon is the intended donor position.

Output structure
----------------
One output file per exon, like your Pangolin workflow:
- one REF row for the WT construct
- one ALT row per SNV construct

Columns include:
- alphagenome_acceptor_score
- alphagenome_donor_score
- fixed construct coordinates used for scoring
- matching padded-sequence coordinates

Assumptions
-----------
1) `nt_seq` is already in the assay/transcript orientation.
2) Non-deletion, non-WT rows are SNVs for this run.
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
    "custom_minigene_snvs_middle_exon_sites_per_exon"
)
os.makedirs(output_dir, exist_ok=True)

TEST_MODE = False
TEST_N_EXONS = 5

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

def is_snv_variant_id(variant_id: str) -> bool:
    """
    Match the same spirit as your Pangolin workflow:
    keep non-deletion, non-WT rows for SNV analysis.
    """
    variant_id = str(variant_id)
    return ("_del" not in variant_id) and (not variant_id.endswith("_wt"))

# ──────────────────────────────────────────────────────────────
# LOAD METADATA
# ──────────────────────────────────────────────────────────────

meta_df = pd.read_csv(meta_file, sep="\t")

wt_mask = meta_df["variant_id"].astype(str).str.endswith("_wt", na=False)
snv_mask = meta_df["variant_id"].astype(str).apply(is_snv_variant_id)

work_df = meta_df[wt_mask | snv_mask].copy()

exon_ids = (
    work_df.loc[snv_mask, "ensembl_exon_id"]
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
        f"{exon_id}_alphagenome_scores_minigene_snvs_middle_exon_sites.tsv"
    )

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

    # Intended middle-exon site coordinates in construct space
    acceptor_construct_pos0 = MIDDLE_EXON_ACCEPTOR_CONSTRUCT_POS0
    donor_construct_pos0 = MIDDLE_EXON_ACCEPTOR_CONSTRUCT_POS0 + exon_len - 1

    # Build WT construct and padded sequence
    wt_core_construct = clean_seq(pre_manual + wt_seq + post_manual)
    wt_padded_seq, wt_left_pad_len, wt_right_pad_len = build_padded_construct(
        wt_core_construct,
        target_len=TARGET_LEN
    )

    # Convert intended construct positions to padded positions
    wt_acceptor_padded_pos0 = wt_left_pad_len + acceptor_construct_pos0
    wt_donor_padded_pos0 = wt_left_pad_len + donor_construct_pos0

    # WT prediction once per exon
    wt_out = model.predict_sequence(
        sequence=wt_padded_seq,
        organism=dna_client.Organism.HOMO_SAPIENS,
        requested_outputs=[OutputType.SPLICE_SITES],
        ontology_terms=ontology_terms,
        interval=None,
    )

    wt_ss = get_splice_sites_obj(wt_out)
    wt_donor_track, wt_acceptor_track = donor_acceptor_positive_strand(wt_ss)

    results = []

    ref_row = {
        "ensembl_exon_id": exon_id,
        "Type": "REF",
        "Identifier": wt_identifier,

        # Geometry
        "exon_len": exon_len,
        "wt_core_len": len(wt_core_construct),
        "target_len": TARGET_LEN,
        "wt_left_pad_len": wt_left_pad_len,
        "wt_right_pad_len": wt_right_pad_len,
        "pre_manual_len": len(pre_manual),
        "post_manual_len": len(post_manual),

        # Fixed intended site coordinates
        "acceptor_construct_pos0": acceptor_construct_pos0,
        "donor_construct_pos0": donor_construct_pos0,
        "acceptor_padded_pos0": wt_acceptor_padded_pos0,
        "donor_padded_pos0": wt_donor_padded_pos0,

        # Scores at intended middle-exon sites
        "alphagenome_acceptor_score": safe_score(wt_acceptor_track, wt_acceptor_padded_pos0),
        "alphagenome_donor_score": safe_score(wt_donor_track, wt_donor_padded_pos0),

        # ALT/SNV-specific metadata blank for REF
        "alt_core_len": np.nan,
        "alt_left_pad_len": np.nan,
        "alt_right_pad_len": np.nan,
        "snv_start_in_exon_1based": np.nan,
        "snv_wt_base": np.nan,
        "snv_alt_base": np.nan,
        "snv_construct_pos0": np.nan,
        "snv_padded_pos0": np.nan,
    }
    results.append(ref_row)

    # ─────────────────────────
    # ALT SNV ROWS
    # ─────────────────────────
    snv_rows = exon_df[exon_df["variant_id"].astype(str).apply(is_snv_variant_id)].copy()

    if len(snv_rows) > 0:
        print(f" • SNVs: n={len(snv_rows)}")
    else:
        print(f" • No SNVs found -> saving REF-only file")

    for _, row in tqdm(snv_rows.iterrows(), total=len(snv_rows), desc=f"{exon_id} SNVs", leave=False):
        alt_identifier = str(row["variant_id"])
        alt_seq = clean_seq(row["nt_seq"])

        # Build ALT construct and padded sequence
        alt_core_construct = clean_seq(pre_manual + alt_seq + post_manual)
        alt_padded_seq, alt_left_pad_len, alt_right_pad_len = build_padded_construct(
            alt_core_construct,
            target_len=TARGET_LEN
        )

        # For SNVs, construct length should be unchanged, so intended construct coordinates stay fixed
        alt_acceptor_padded_pos0 = alt_left_pad_len + acceptor_construct_pos0
        alt_donor_padded_pos0 = alt_left_pad_len + donor_construct_pos0

        # Optional SNV metadata if present in input table
        snv_start_in_exon_1based = row["start"] if "start" in row.index else np.nan
        snv_wt_base = row["wt"] if "wt" in row.index else np.nan
        snv_alt_base = row["mut"] if "mut" in row.index else np.nan

        if pd.notna(snv_start_in_exon_1based):
            snv_start_in_exon_1based = int(snv_start_in_exon_1based)
            snv_construct_pos0 = len(pre_manual) + (snv_start_in_exon_1based - 1)
            snv_padded_pos0 = alt_left_pad_len + snv_construct_pos0
        else:
            snv_construct_pos0 = np.nan
            snv_padded_pos0 = np.nan

        # Predict ALT sequence
        alt_out = model.predict_sequence(
            sequence=alt_padded_seq,
            organism=dna_client.Organism.HOMO_SAPIENS,
            requested_outputs=[OutputType.SPLICE_SITES],
            ontology_terms=ontology_terms,
            interval=None,
        )

        alt_ss = get_splice_sites_obj(alt_out)
        alt_donor_track, alt_acceptor_track = donor_acceptor_positive_strand(alt_ss)

        alt_row = {
            "ensembl_exon_id": exon_id,
            "Type": "ALT",
            "Identifier": alt_identifier,

            # Geometry
            "exon_len": exon_len,
            "wt_core_len": len(wt_core_construct),
            "alt_core_len": len(alt_core_construct),
            "target_len": TARGET_LEN,
            "wt_left_pad_len": wt_left_pad_len,
            "wt_right_pad_len": wt_right_pad_len,
            "alt_left_pad_len": alt_left_pad_len,
            "alt_right_pad_len": alt_right_pad_len,
            "pre_manual_len": len(pre_manual),
            "post_manual_len": len(post_manual),

            # Fixed intended site coordinates
            "acceptor_construct_pos0": acceptor_construct_pos0,
            "donor_construct_pos0": donor_construct_pos0,
            "acceptor_padded_pos0": alt_acceptor_padded_pos0,
            "donor_padded_pos0": alt_donor_padded_pos0,

            # Scores at intended middle-exon sites
            "alphagenome_acceptor_score": safe_score(alt_acceptor_track, alt_acceptor_padded_pos0),
            "alphagenome_donor_score": safe_score(alt_donor_track, alt_donor_padded_pos0),

            # SNV metadata
            "snv_start_in_exon_1based": snv_start_in_exon_1based,
            "snv_wt_base": snv_wt_base,
            "snv_alt_base": snv_alt_base,
            "snv_construct_pos0": snv_construct_pos0,
            "snv_padded_pos0": snv_padded_pos0,
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
