"""
combine_barcode.py — Combine per-chunk barcode pickles for one library.

Usage: python3 combine_barcode.py <lib_id> <chunks_dir> <out_dir> <min_reads>

For each barcode, keeps only the most common variant if it has >= min_reads
supporting reads. Aggregates per-chunk extraction stats into a summary JSON.
"""
import sys
import os
import glob
import json
import pickle
import gzip
import pandas as pd
from collections import Counter

lib_id     = sys.argv[1]
chunks_dir = sys.argv[2]
out_dir    = sys.argv[3]
min_reads  = int(sys.argv[4])

# ── Combine pickles ────────────────────────────────────────────────────────────
pattern   = os.path.join(chunks_dir, f"{lib_id}_chunk*.pkl")
pkl_files = sorted(glob.glob(pattern))

if not pkl_files:
    raise FileNotFoundError(f"No pickle files found: {pattern}")
print(f"Combining {len(pkl_files)} chunks for {lib_id}...")

combined = {}
for fpath in pkl_files:
    with open(fpath, "rb") as fh:
        data = pickle.load(fh)
    for barcode, (var_counter, _) in data.items():
        if barcode not in combined:
            combined[barcode] = Counter()
        combined[barcode].update(var_counter)

# ── Apply min_reads threshold ─────────────────────────────────────────────────
barcodes, variants, counts = [], [], []
for barcode, counter in combined.items():
    top_variant, top_count = counter.most_common(1)[0]
    if top_count >= min_reads:
        barcodes.append(barcode)
        variants.append(top_variant)
        counts.append(top_count)

print(f"Total barcodes: {len(combined):,} | passing min_reads={min_reads}: {len(barcodes):,}")

df = pd.DataFrame({"barcode": barcodes, "variant": variants, "count": counts})

out_file = os.path.join(out_dir, f"{lib_id}_bc_var_min{min_reads}reads.tsv.gz")
df.to_csv(out_file, index=False, sep="\t", compression="gzip")
print(f"Saved: {out_file}")

# ── Aggregate per-chunk extraction stats ──────────────────────────────────────
stats_files = sorted(glob.glob(os.path.join(chunks_dir, f"{lib_id}_chunk*_stats.json")))
agg = {"reads_total": 0, "reads_used": 0, "reads_tossed": 0,
       "tossed_flank": 0, "tossed_barcode": 0, "tossed_length": 0}
for sf in stats_files:
    with open(sf) as fh:
        s = json.load(fh)
    for k in agg:
        agg[k] += s.get(k, 0)

summary = {
    "lib_id":                    lib_id,
    "n_chunks":                  len(pkl_files),
    **agg,
    "barcodes_total":            len(combined),
    f"barcodes_min{min_reads}":  len(barcodes),
    "pct_reads_used":            round(100 * agg["reads_used"] / agg["reads_total"], 2) if agg["reads_total"] else 0,
}
summary_path = os.path.join(out_dir, f"{lib_id}_step4_stats.json")
with open(summary_path, "w") as fh:
    json.dump(summary, fh, indent=2)
print(f"Stats: {summary_path}")
