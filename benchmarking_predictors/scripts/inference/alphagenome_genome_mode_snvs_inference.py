#@title Run AlphaGenome SNV scoring in genomic-context mode
"""
Fixed-window (per exon) version:
- One fixed interval per exon (length interval_len) centered on exon midpoint.
- WT splice-site tracks computed ONCE via predict_interval(fixed_interval).
- ALT tracks computed per SNV via predict_variant(fixed_interval, variant).
- Extract REF+ALT donor/acceptor scores at exon start/end from the SAME fixed interval.

Result:
- WT scores at canonical boundaries are constant across variants within an exon.
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

OUT_DIR = "results/alphagenome/genome_mode/snvs/siteonly_fixedwindow_long_context"
os.makedirs(OUT_DIR, exist_ok=True)

# =========================
# RUN SETTINGS
# =========================
TEST_MODE = False
TEST_N_EXONS = 1

interval_len = 16_384
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
    """Lightweight VCF reader for CHROM/POS/ID/REF/ALT (ignores genotypes)."""
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

    # Fill missing IDs deterministically
    df["ID"] = df["ID"].replace(".", np.nan)
    df["ID"] = df["ID"].fillna(
        df["CHROM"] + ":" + df["POS"].astype(str) + ":" + df["REF"].astype(str) + ">" + df["ALT"].astype(str)
    )
    return df

def is_snv(row) -> bool:
    """True SNV: single-base REF and ALT; exclude multi-allelic ALT."""
    ref = str(row["REF"])
    alt = str(row["ALT"])
    return ("," not in alt) and (len(ref) == 1) and (len(alt) == 1)

def interval_centered_on_pos1(chrom: str, center_pos1: int) -> genome.Interval:
    """Centered 0-based half-open interval around a 1-based coordinate."""
    center0 = int(center_pos1) - 1
    start0 = center0 - half_len
    end0 = start0 + interval_len
    return genome.Interval(chromosome=chrom, start=int(start0), end=int(end0))

def get_splice_sites_obj(pred_out):
    """Find an object with .splice_sites for either interval or variant outputs."""
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

def score_at_site_pos1(track: np.ndarray, interval_start0: int, site_pos1: int) -> float:
    """idx = (site_pos1 - 1) - interval_start0"""
    idx = (int(site_pos1) - 1) - int(interval_start0)
    if 0 <= idx < len(track):
        return float(track[idx])
    return np.nan

def parse_exon_id_from_vcf_path(vcf_path: str) -> str:
    """Extract exon_id from filenames like ENSE00000369924_*.vcf"""
    base = Path(vcf_path).name
    # ✅ FIX: use \d not \\d
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

def exon_midpoint_pos1(start_exon: int, end_exon: int) -> int:
    """Pick a fixed center coordinate for the exon window (1-based)."""
    return int((int(start_exon) + int(end_exon)) // 2)

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
# 3) Main loop (SNVs only): FIXED interval per exon
# =========================
for exon_id in tqdm(work_exons, desc="Exons"):
    out_csv = os.path.join(OUT_DIR, f"{exon_id}_SNV_splice_sites_SITEONLY_FIXEDWINDOW.csv")

    # ---- SKIP if already done ----
    if os.path.exists(out_csv) and os.path.getsize(out_csv) > 0:
        print(f"[SKIP] {exon_id}: output exists -> {out_csv}")
        continue

    info = join_lookup[exon_id]
    strand_symbol = info["strand"]
    start_exon = int(info["start_exon"])
    end_exon = int(info["end_exon"])
    vcf_path = exon_to_vcf[exon_id]

    df_vcf = read_vcf_as_df(vcf_path)

    # SNV filter + exclude deletion IDs by naming convention
    df_snvs = df_vcf[df_vcf.apply(is_snv, axis=1)].copy()
    df_snvs = df_snvs[~df_snvs["ID"].astype(str).str.contains("del", case=False, na=False)].copy()

    if len(df_snvs) == 0:
        print(f"[{exon_id}] No SNVs found -> writing empty file (so we skip next time)")
        pd.DataFrame().to_csv(out_csv, index=False)
        continue

    # ---- FIXED interval for this exon: centered on exon midpoint ----
    chrom_exon = str(df_snvs["CHROM"].iloc[0])
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

    canonical_sites = [
        {"boundary": "start_exon", "site_pos1": start_exon, "role": "acceptor" if strand_symbol == "+" else "donor"},
        {"boundary": "end_exon",   "site_pos1": end_exon,   "role": "donor"    if strand_symbol == "+" else "acceptor"},
    ]

    rows = []

    for _, r in tqdm(df_snvs.iterrows(), total=len(df_snvs), desc=f"{exon_id} SNVs", leave=False):
        chrom = str(r["CHROM"])
        pos1 = int(r["POS"])
        ref = str(r["REF"])
        alt = str(r["ALT"])
        vid = str(r["ID"])

        # should always match, but guard anyway
        if chrom != chrom_exon:
            continue

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
            site_pos1 = int(s["site_pos1"])

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
                "canonical_site_1based": site_pos1,

                "variant_id": vid,
                "variant_pos_1based": pos1,
                "ref_allele": ref,
                "alt_allele": alt,

                # WT (fixed window; constant per exon)
                "donor_ref_score_at_site": score_at_site_pos1(wt_d, fixed_interval.start, site_pos1),
                "acceptor_ref_score_at_site": score_at_site_pos1(wt_a, fixed_interval.start, site_pos1),

                # ALT (same fixed window)
                "donor_alt_score_at_site": score_at_site_pos1(alt_d, fixed_interval.start, site_pos1),
                "acceptor_alt_score_at_site": score_at_site_pos1(alt_a, fixed_interval.start, site_pos1),
            })

    df_out = pd.DataFrame(rows)
    df_out.to_csv(out_csv, index=False)

    print(f"[{exon_id}] SNVs={len(df_snvs)} -> saved {out_csv} (rows={len(df_out)})")
    print(df_out.head(10).to_string(index=False))

print("\nDONE.")

