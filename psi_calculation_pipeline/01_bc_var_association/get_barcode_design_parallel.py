"""
get_barcode_design_parallel.py — Extract barcode–variant pairs from merged reads.

Usage: python3 get_barcode_design_parallel.py <chunk_file> <upstream_flank>
           <downstream_flank> <design_id> <out_dir>

For each read, the variant sequence is delimited by upstream/downstream flanks
and the barcode is delimited by EXON7_SEQ and BARCODE_FLANK.
Results are saved as a pickle (combined by combine_barcode.py) and a stats JSON.
"""
import sys
import os
import gzip
import json
import pickle
import regex
from collections import Counter
from tqdm.auto import tqdm

# ── Arguments ─────────────────────────────────────────────────────────────────
chunk_file     = sys.argv[1]
upstream_flank = sys.argv[2]
downstream_flank = sys.argv[3]
design_id      = sys.argv[4]   # e.g. MUT1_chunk_03
out_dir        = sys.argv[5]

# ── Constants (minigene construct) ────────────────────────────────────────────
BARCODE_FLANK = "CTACTGATTCGATGCAAGCTT"
EXON7_SEQ     = "GCAGAAAGCACAGAAAGGAA"
BARCODE_LEN   = 38
MAX_ERRS      = 2

# ── Precompile fuzzy patterns ─────────────────────────────────────────────────
upstream_re      = regex.compile(f"({upstream_flank}){{e<={MAX_ERRS}}}")
downstream_re    = regex.compile(f"({downstream_flank}){{e<={MAX_ERRS}}}")
exon7_re         = regex.compile(f"({EXON7_SEQ}){{e<={MAX_ERRS}}}")
barcodeflank_re  = regex.compile(f"({BARCODE_FLANK}){{e<={MAX_ERRS}}}")

# ── Counters ──────────────────────────────────────────────────────────────────
barcodeexondict = {}
total = 0
used  = 0
tossed_flank   = 0   # upstream/downstream flank not found
tossed_barcode = 0   # exon7/barcodeflank not found
tossed_length  = 0   # barcode wrong length

def gen_lines(path):
    with gzip.open(path, "rb") as f:
        for line in f:
            yield line.decode("utf-8").strip()

for line in tqdm(gen_lines(chunk_file), desc=design_id):
    total += 1

    if total % 10_000_000 == 0:
        print(f"  {total:,} reads | used {used:,} | tossed {total - used:,}")

    # ── Extract variant sequence ───────────────────────────────────────────────
    s = line.find(upstream_flank)
    e = line.find(downstream_flank)
    if s != -1 and e != -1 and s < e:
        variant = line[s + len(upstream_flank):e]
    else:
        um = upstream_re.search(line)
        dm = downstream_re.search(line)
        if not um or not dm or um.end() >= dm.start():
            tossed_flank += 1
            continue
        variant = line[um.end():dm.start() + 2]

    # ── Extract barcode ────────────────────────────────────────────────────────
    ei = line.find(EXON7_SEQ)
    bi = line.find(BARCODE_FLANK)
    if ei != -1 and bi != -1 and ei < bi:
        bc_start = ei + len(EXON7_SEQ)
    else:
        em = exon7_re.search(line)
        bm = barcodeflank_re.search(line)
        if not em or not bm or em.end() >= bm.start():
            tossed_barcode += 1
            continue
        bc_start = em.end()

    bc_end = bc_start + BARCODE_LEN
    if bc_end > len(line):
        tossed_length += 1
        continue

    barcode = line[bc_start:bc_end]
    if len(barcode) != BARCODE_LEN:
        tossed_length += 1
        continue

    if barcode not in barcodeexondict:
        barcodeexondict[barcode] = [Counter(), 0]
    barcodeexondict[barcode][0][variant] += 1
    used += 1

# ── Save pickle ───────────────────────────────────────────────────────────────
pkl_path = os.path.join(out_dir, f"{design_id}.pkl")
with open(pkl_path, "wb") as f:
    pickle.dump(barcodeexondict, f)

# ── Save per-chunk stats ───────────────────────────────────────────────────────
stats = {
    "design_id":      design_id,
    "reads_total":    total,
    "reads_used":     used,
    "reads_tossed":   total - used,
    "tossed_flank":   tossed_flank,
    "tossed_barcode": tossed_barcode,
    "tossed_length":  tossed_length,
    "barcodes_found": len(barcodeexondict),
}
stats_path = os.path.join(out_dir, f"{design_id}_stats.json")
with open(stats_path, "w") as f:
    json.dump(stats, f, indent=2)

print(f"[{design_id}] total={total:,} used={used:,} tossed={total-used:,} "
      f"barcodes={len(barcodeexondict):,}")
