"""
combine_psi.py — Concatenate per-chunk PSI TSVs into a single file.

Usage: python3 combine_psi.py --lib-id MUT3 --out-dir .../MUT3/
"""
import argparse
import glob
import os
import sys
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--lib-id",  required=True)
p.add_argument("--out-dir", required=True)
args = p.parse_args()

pattern = os.path.join(args.out_dir, "chunks",
                       f"psi_per_barcode_{args.lib_id}_bc*_*.tsv")
files = sorted(glob.glob(pattern),
               key=lambda f: int(f.rsplit("_bc", 1)[1].split("_")[0]))

if not files:
    sys.exit(f"Error: no chunk files found matching {pattern}")

print(f"Combining {len(files)} chunks...")
dfs = [pd.read_csv(f, sep="\t") for f in files]
result = pd.concat(dfs, ignore_index=True)
result = result.sort_values("barcode_id").reset_index(drop=True)

out_file = os.path.join(args.out_dir, f"psi_per_barcode_{args.lib_id}.tsv")
result.to_csv(out_file, sep="\t", index=False)
print(f"  {len(result):,} barcodes → {out_file}")
