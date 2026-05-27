#!/usr/bin/env python
#@title Run AlphaGenome deletion scoring in genomic-context mode
"""
Fixed-window (per exon) version of the deletion realignment script:

- One fixed interval per exon (length interval_len) centered on exon midpoint.
- WT splice-site tracks computed ONCE via predict_interval(fixed_interval).
- ALT tracks computed per DEL via predict_variant(fixed_interval, variant).
- Report REF donor/acceptor at exon start/end (from WT tracks).
- Report ALT donor/acceptor at exon start/end, realigned to REF coordinate:
    * upstream of deletion: alt_pos = p
    * inside deleted span:  alt_score = NaN and deleted_in_alt=True
    * downstream:          alt_pos = p - deletion_len_bp

Result:
- WT scores at canonical boundaries are constant across variants within an exon.
- ALT values still properly account for coordinate shift from the deletion.
"""

import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from alphagenome import colab_utils
from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.models.dna_output import OutputType

# =========================
# USER PATHS
# =========================
VCF_DIR = "data/input/alphagenome/genome_vcf_per_exon"

EXON_FIRST_PAPER_TXT = "data/input/opensplice_predictors_benchmarking_exon_metadata.tsv"

OUT_DIR = (
    "results/alphagenome/"
    "exon_608_splice_sites_base_resolution_deletions_strandaware_SITEONLY_realigned_FIXEDWINDOW"
)
os.makedirs(OUT_DIR, exist_ok=True)

# =========================
# RUN SETTINGS
# =========================
TEST_MODE = False
TEST_N_EXONS = 2

interval_len = 16384
# Also tested 1_048_576 context; 16,384 gave better correlation with experimental data
# and was used for subsequent analyses.
half_len = interval_len // 2
ontology_terms = ["CL:0002518"]
# kidney epithelial cell (closest available to HEK); note ontology terms mostly
# affect splice-site usage/junction outputs, while this workflow uses splice-site probability.

print("VCF_DIR:", VCF_DIR)
print("EXON_FIRST_PAPER_TXT:", EXON_FIRST_PAPER_TXT)
print("OUT_DIR:", OUT_DIR)
print("interval_len:", interval_len)
print("ontology_terms:", ontology_terms)

# =========================
# AlphaGenome client
# =========================
api_key = colab_utils.get_api_key()
dna_model = dna_client.create(api_key)
organism = dna_client.Organism.HOMO_SAPIENS

# =========================
# Helpers
# =========================
def read_vcf_as_df(path: str) -> pd.DataFrame:
    """Lightweight VCF reader for CHROM/POS/ID/REF/ALT."""
    header, rows = None, []
    with open(path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.lstrip("#").strip().split("\t")
                continue
            if not line.startswith("#") and line.strip():
                rows.append(line.rstrip().split("\t"))

    df = pd.DataFrame(rows, columns=header)
    df["POS"] = df["POS"].astype(int)

    # Standardize CHROM to "chrN"
    df["CHROM"] = df["CHROM"].astype(str)
    df["CHROM"] = df["CHROM"].apply(lambda x: x if x.startswith("chr") else f"chr{x}")

    # Ensure ID exists
    df["ID"] = df["ID"].replace(".", np.nan)
    df["ID"] = df["ID"].fillna(
        df["CHROM"] + ":" + df["POS"].astype(str) + ":" + df["REF"].astype(str) + ">" + df["ALT"].astype(str)
    )
    return df

def is_deletion_by_id(row) -> bool:
    """Your deletion filter: variant ID contains 'del' (case-insensitive)."""
    return "del" in str(row["ID"]).lower()

def interval_centered_on_pos1(chrom: str, center_pos1: int) -> genome.Interval:
    """Centered 0-based half-open interval around a 1-based coordinate."""
    center0 = int(center_pos1) - 1
    start0 = center0 - half_len
    end0 = start0 + interval_len
    return genome.Interval(chromosome=chrom, start=int(start0), end=int(end0))

def exon_midpoint_pos1(start_exon: int, end_exon: int) -> int:
    """Pick a fixed center coordinate for the exon window (1-based)."""
    return int((int(start_exon) + int(end_exon)) // 2)

def get_splice_sites_obj(pred_out):
    """Find object with .splice_sites for either interval or variant outputs."""
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
    raise AttributeError("splice_sites not found in prediction output")

def donor_acceptor_strandaware(splice_sites, strand_symbol: str):
    """Return donor, acceptor arrays after strand filtering."""
    ss = (
        splice_sites.filter_to_positive_strand()
        if strand_symbol == "+"
        else splice_sites.filter_to_negative_strand()
    )
    arr = np.asarray(ss.values)
    names = [str(n).lower() for n in ss.names]
    donor = arr[:, names.index("donor")]
    acceptor = arr[:, names.index("acceptor")]
    return donor, acceptor

def parse_exon_id_from_vcf_path(vcf_path: str) -> str:
    """Extract exon_id from filenames like ENSE00000369924_*.vcf"""
    base = Path(vcf_path).name
    m = re.match(r"^(ENSE\d+)_", base)
    return m.group(1) if m else ""

def norm_strand(x):
    if pd.isna(x):
        return np.nan
    s = str(x).strip()
    if s in {"1", "1.0", "+", "plus"}:
        return "+"
    if s in {"-1", "-1.0", "-", "minus"}:
        return "-"
    return s

def vcf_deleted_span_1based(pos1: int, ref: str, alt: str):
    """
    Canonical VCF deletion semantics (anchor base at POS kept):
      deleted bases are POS+1 .. POS+len(REF)-1 (1-based)
    deletion_len_bp = len(REF) - len(ALT)
    """
    ref = str(ref)
    alt = str(alt)
    deletion_len = len(ref) - len(alt)
    if deletion_len <= 0:
        return None, None, 0
    if len(ref) <= 1:
        return None, None, int(deletion_len)
    del_start1 = pos1 + 1
    del_end1 = pos1 + (len(ref) - 1)
    return int(del_start1), int(del_end1), int(deletion_len)

def refpos_to_altpos_1based(p1: int, del_start1: int | None, del_end1: int | None, del_len: int):
    """Map REF coordinate p1 -> ALT coordinate, or None if deleted."""
    if del_len <= 0 or del_start1 is None or del_end1 is None:
        return int(p1)
    if p1 < del_start1:
        return int(p1)
    if del_start1 <= p1 <= del_end1:
        return None
    return int(p1 - del_len)

def score_at_refpos_with_alt_liftover(ref_track, alt_track, interval_start0: int, p1: int,
                                     del_start1, del_end1, del_len):
    """
    Return (ref_score, alt_score_realigned, alt_pos_used_1based, deleted_in_alt, idx_ref, idx_alt)
    """
    ref_idx = (int(p1) - 1) - int(interval_start0)
    ref_score = float(ref_track[ref_idx]) if 0 <= ref_idx < len(ref_track) else np.nan

    alt_p1 = refpos_to_altpos_1based(int(p1), del_start1, del_end1, del_len)
    if alt_p1 is None:
        return ref_score, np.nan, np.nan, True, ref_idx, np.nan

    alt_idx = (int(alt_p1) - 1) - int(interval_start0)
    alt_score = float(alt_track[alt_idx]) if 0 <= alt_idx < len(alt_track) else np.nan
    return ref_score, alt_score, int(alt_p1), False, ref_idx, alt_idx

# =========================
# 1) Load exon targets
# =========================
df_join = pd.read_csv(EXON_FIRST_PAPER_TXT, sep=None, engine="python")
df_join = df_join[["ensembl_exon_id", "strand", "start_exon", "end_exon"]].copy()
df_join = df_join.rename(columns={"ensembl_exon_id": "exon_id"})

df_join["exon_id"] = df_join["exon_id"].astype(str)
df_join["strand"] = df_join["strand"].apply(norm_strand)
df_join["start_exon"] = pd.to_numeric(df_join["start_exon"], errors="coerce")
df_join["end_exon"] = pd.to_numeric(df_join["end_exon"], errors="coerce")
df_join = df_join.dropna(subset=["strand", "start_exon", "end_exon"]).copy()
df_join["start_exon"] = df_join["start_exon"].astype(int)
df_join["end_exon"] = df_join["end_exon"].astype(int)

df_join = (
    df_join
    .groupby("exon_id", as_index=False)
    .agg({"strand": "first", "start_exon": "min", "end_exon": "max"})
)

print("Final exons with strand+start+end:", df_join["exon_id"].nunique())
print(df_join.head(5).to_string(index=False))

target_exons = set(df_join["exon_id"].unique())
join_lookup = df_join.set_index("exon_id").to_dict(orient="index")

# =========================
# 2) Find VCF files and filter to those exons
# =========================
vcf_paths_all = sorted(str(p) for p in Path(VCF_DIR).glob("*.vcf"))
print("Found VCF files:", len(vcf_paths_all))

exon_to_vcf = {}
for vp in vcf_paths_all:
    exid = parse_exon_id_from_vcf_path(vp)
    if exid and (exid in target_exons) and (exid not in exon_to_vcf):
        exon_to_vcf[exid] = vp

print("VCFs matching target exons:", len(exon_to_vcf))

work_exons = [exid for exid in df_join["exon_id"].tolist() if exid in exon_to_vcf]
if TEST_MODE:
    work_exons = work_exons[:TEST_N_EXONS]
    print(f"TEST_MODE=True -> running {len(work_exons)} exon(s):", work_exons)
else:
    print(f"Running {len(work_exons)} exon(s)")

# =========================
# 3) Main loop (DELETIONS only): FIXED interval per exon + ALT realignment
# =========================
for exon_id in tqdm(work_exons, desc="Exons"):
    out_csv = os.path.join(OUT_DIR, f"{exon_id}_DEL_splice_sites_SITEONLY_realigned_FIXEDWINDOW.csv")
    if os.path.exists(out_csv) and os.path.getsize(out_csv) > 0:
        print(f"\n[{exon_id}] output exists -> skipping: {out_csv}")
        continue

    info = join_lookup[exon_id]
    strand_symbol = info["strand"]
    start_exon = int(info["start_exon"])
    end_exon = int(info["end_exon"])
    vcf_path = exon_to_vcf[exon_id]

    df_vcf = read_vcf_as_df(vcf_path)
    df_del = df_vcf[df_vcf.apply(is_deletion_by_id, axis=1)].copy()

    print(f"\n[{exon_id}] deletions (ID contains 'del'):", len(df_del))
    if len(df_del) == 0:
        pd.DataFrame().to_csv(out_csv, index=False)
        continue

    # ---- FIXED interval for this exon: centered on exon midpoint ----
    chrom_exon = str(df_del["CHROM"].iloc[0])
    center_pos1 = exon_midpoint_pos1(start_exon, end_exon)
    fixed_interval = interval_centered_on_pos1(chrom_exon, center_pos1)

    # ---- WT once per exon ----
    wt_out = dna_model.predict_interval(
        interval=fixed_interval,
        organism=organism,
        requested_outputs=[OutputType.SPLICE_SITES],
        ontology_terms=ontology_terms,
    )
    wt_ss = get_splice_sites_obj(wt_out)
    wt_d, wt_a = donor_acceptor_strandaware(wt_ss, strand_symbol)

    rows = []
    canonical_sites = [
        {"boundary": "start_exon", "site_pos1": start_exon, "role": "acceptor" if strand_symbol == "+" else "donor"},
        {"boundary": "end_exon",   "site_pos1": end_exon,   "role": "donor"    if strand_symbol == "+" else "acceptor"},
    ]

    for _, r in tqdm(df_del.iterrows(), total=len(df_del), desc=f"{exon_id} DELs", leave=False):
        chrom = str(r["CHROM"])
        pos1 = int(r["POS"])
        ref = str(r["REF"])
        alt = str(r["ALT"])
        vid = str(r["ID"])

        if chrom != chrom_exon:
            continue

        del_start1, del_end1, del_len = vcf_deleted_span_1based(pos1, ref, alt)

        var = genome.Variant(
            chromosome=chrom,
            position=pos1,
            reference_bases=ref,
            alternate_bases=alt,
            name=vid
        )

        alt_out = dna_model.predict_variant(
            interval=fixed_interval,
            variant=var,
            organism=organism,
            requested_outputs=[OutputType.SPLICE_SITES],
            ontology_terms=ontology_terms,
        )
        alt_ss = get_splice_sites_obj(alt_out)
        alt_d, alt_a = donor_acceptor_strandaware(alt_ss, strand_symbol)

        for s in canonical_sites:
            p1 = int(s["site_pos1"])

            # donor: WT from fixed window; ALT realigned by liftover
            wt_d_sc, alt_d_sc, alt_d_p1, del_d, idx_ref_d, idx_alt_d = score_at_refpos_with_alt_liftover(
                wt_d, alt_d, fixed_interval.start, p1, del_start1, del_end1, del_len
            )

            # acceptor: WT from fixed window; ALT realigned by liftover
            wt_a_sc, alt_a_sc, alt_a_p1, del_a, idx_ref_a, idx_alt_a = score_at_refpos_with_alt_liftover(
                wt_a, alt_a, fixed_interval.start, p1, del_start1, del_end1, del_len
            )

            rows.append({
                "exon_id": exon_id,
                "strand": strand_symbol,
                "vcf_path": vcf_path,
                "start_exon": start_exon,
                "end_exon": end_exon,

                "fixed_interval_chrom": fixed_interval.chromosome,
                "fixed_interval_start0": int(fixed_interval.start),
                "fixed_interval_end0": int(fixed_interval.end),
                "fixed_interval_center_pos1": int(center_pos1),

                "canonical_boundary": s["boundary"],
                "canonical_role_by_strand": s["role"],
                "canonical_site_1based": p1,

                "variant_id": vid,
                "variant_pos_1based": pos1,
                "ref_allele": ref,
                "alt_allele": alt,

                "deletion_span_start_1based": del_start1,
                "deletion_span_end_1based": del_end1,
                "deletion_len_bp": del_len,

                # donor site-only
                "donor_ref_score_at_site": wt_d_sc,
                "donor_alt_score_realigned": alt_d_sc,
                "donor_alt_pos_used_1based": alt_d_p1,
                "donor_site_deleted_in_alt": bool(del_d),
                "donor_ref_idx": idx_ref_d,
                "donor_alt_idx": idx_alt_d,

                # acceptor site-only
                "acceptor_ref_score_at_site": wt_a_sc,
                "acceptor_alt_score_realigned": alt_a_sc,
                "acceptor_alt_pos_used_1based": alt_a_p1,
                "acceptor_site_deleted_in_alt": bool(del_a),
                "acceptor_ref_idx": idx_ref_a,
                "acceptor_alt_idx": idx_alt_a,
            })

    df_out = pd.DataFrame(rows)
    df_out.to_csv(out_csv, index=False)

    print(f"[{exon_id}] saved {out_csv} (rows={len(df_out)})")
    print(df_out.head(12).to_string(index=False))

print("\nDONE.")

