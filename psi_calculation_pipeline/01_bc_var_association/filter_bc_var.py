"""
filter_bc_var.py — Filter and annotate barcode–variant pairs.

Keeps only barcodes whose variant sequence matches a designed oligo in
variant_mapping_all.tsv AND has >= min_reads supporting reads.
Adds variant annotations (unique identifier, exonic sequence, library ID)
and writes the final barcode dictionary + a per-step summary table.

Usage: python3 filter_bc_var.py
         --lib-id       MUT1
         --bc-var-file  .../MUT1_bc_var_min2reads.tsv.gz
         --mapping-file .../variant_mapping_all.tsv
         --metadata-file .../exon_list_with_metadata.tsv
         --min-reads    5
         --out-dir      .../MUT1/
"""
import argparse
import json
import os
import sys
import pandas as pd

# ── Arguments ──────────────────────────────────────────────────────────────────
p = argparse.ArgumentParser()
p.add_argument("--lib-id",        required=True)
p.add_argument("--bc-var-file",   required=True)
p.add_argument("--mapping-file",  required=True)
p.add_argument("--metadata-file", required=True)
p.add_argument("--min-reads",     type=int, default=5)
p.add_argument("--out-dir",       required=True)
args = p.parse_args()

lib_id    = args.lib_id
min_reads = args.min_reads
print(f"Filtering {lib_id} | min_reads={min_reads}")

# ── Load variant mapping ───────────────────────────────────────────────────────
mapping = pd.read_csv(args.mapping_file, sep="\t", dtype=str)
meta    = pd.read_csv(args.metadata_file, sep="\t", dtype=str,
                      usecols=["ensembl_exon_id", "exon_id", "sat_mutagenesis_library_id"])

mapping = mapping.merge(meta, on="ensembl_exon_id", how="left")
mapping["unique_identifier"] = mapping["exon_id"] + "_" + mapping["identifier"]

# Deduplicate: same nt_seq can appear for multiple positions (synonymous variants)
lookup = mapping.drop_duplicates(subset="nt_seq").set_index("nt_seq")
n_designed = len(lookup)
print(f"  Designed variants loaded: {n_designed:,}")

# ── Load bc_var pairs ──────────────────────────────────────────────────────────
bc = pd.read_csv(args.bc_var_file, sep="\t", dtype={"barcode": str, "variant": str, "count": int})
n_raw = len(bc)
print(f"  bc_var pairs (min 2 reads): {n_raw:,}")

# ── Filter 1: designed variants only ──────────────────────────────────────────
bc_designed = bc[bc["variant"].isin(lookup.index)].copy()
n_designed_match = len(bc_designed)
n_variants_found = bc_designed["variant"].nunique()

# ── Filter 2: min_reads ────────────────────────────────────────────────────────
bc_filt = bc_designed[bc_designed["count"] >= min_reads].copy()
n_filt          = len(bc_filt)
n_variants_kept = bc_filt["variant"].nunique()

print(f"  After designed-variant filter : {n_designed_match:,} barcodes ({n_variants_found:,} unique variants)")
print(f"  After min_reads={min_reads} filter: {n_filt:,} barcodes ({n_variants_kept:,} unique variants)")

# ── Annotate ──────────────────────────────────────────────────────────────────
bc_filt["variant_id"] = bc_filt["variant"].map(lookup["unique_identifier"])
bc_filt["exon"]       = bc_filt["variant"].map(lookup["exonic_seq"])
bc_filt["lib_id"]     = bc_filt["variant"].map(lookup["sat_mutagenesis_library_id"])

bc_out = bc_filt.reset_index(drop=True)
bc_out.insert(0, "barcode_id", range(1, len(bc_out) + 1))
bc_out = bc_out.rename(columns={"count": "barcode_read_count", "variant": "variant_sequence"})
bc_out = bc_out[["barcode_id", "barcode", "barcode_read_count",
                  "variant_id", "exon", "variant_sequence", "lib_id"]]

# ── Write dictionary ──────────────────────────────────────────────────────────
dict_file = os.path.join(args.out_dir, f"{lib_id}_bc_var_dictionary.tsv.gz")
bc_out.to_csv(dict_file, sep="\t", index=False, compression="gzip")
print(f"  Dictionary: {dict_file}")

# ── Collect all stats into final summary ──────────────────────────────────────
rows = []

for step_tag, stats_file in [
    ("step1", os.path.join(args.out_dir, f"{lib_id}_step1_stats.tsv")),
    ("step2", os.path.join(args.out_dir, f"{lib_id}_step2_stats.tsv")),
]:
    if os.path.exists(stats_file):
        df = pd.read_csv(stats_file, sep="\t", dtype=str)
        rows.append(df)

step4_json = os.path.join(args.out_dir, f"{lib_id}_step4_stats.json")
if os.path.exists(step4_json):
    with open(step4_json) as fh:
        s4 = json.load(fh)
    rows.append(pd.DataFrame({
        "step":   "step4_combine",
        "metric": list(s4.keys()),
        "value":  [str(v) for v in s4.values()]
    }))

pct_covered = round(100 * n_variants_kept / n_designed, 1) if n_designed else 0
rows.append(pd.DataFrame({
    "step":   "step5_filter",
    "metric": ["barcodes_before_filter", "barcodes_designed_match",
               "barcodes_min_reads_kept", "unique_variants_kept",
               "pct_designed_variants_covered"],
    "value":  [str(x) for x in [n_raw, n_designed_match, n_filt,
                                  n_variants_kept, pct_covered]]
}))

summary = pd.concat(rows, ignore_index=True)
summary["lib_id"] = lib_id
summary_file = os.path.join(args.out_dir, f"{lib_id}_summary.tsv")
summary.to_csv(summary_file, sep="\t", index=False)
print(f"  Summary: {summary_file}")

print(f"\n── Summary ──────────────────────────────────────────────────────────────")
print(f"  bc_var pairs (≥2 reads) : {n_raw:,}")
print(f"  Barcodes kept (≥{min_reads} reads): {n_filt:,}")
print(f"  Variants covered        : {n_variants_kept:,} / {n_designed:,} ({pct_covered}%)")
