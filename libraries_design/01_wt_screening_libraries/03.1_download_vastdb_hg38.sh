#!/bin/bash
# Download and merge VastDB hg38 annotation files
# Usage: bash download_vastdb_hg38.sh /path/to/base_dir

set -e

if [ -z "$1" ]; then
    echo "Usage: $0 /path/to/output_dir"
    exit 1
fi

BASE_DIR="$1"
OUTDIR="$BASE_DIR/libraries_design/01_wt_screening_libraries/necessary_file/VastDB"

mkdir -p "$OUTDIR"
cd "$OUTDIR"


# ------------------------- #
# 1. Download VastDB files
# ------------------------- #
echo "Downloading VastDB hg38 files..."

wget -c https://vastdb.crg.eu/downloads/hg38/SPLICE_SITE_SCORES-hg38.tab.gz
wget -c https://vastdb.crg.eu/downloads/hg38/EVENT_METRICS-hg38.tab.gz
wget -c https://vastdb.crg.eu/downloads/hg38/PSI_TABLE-hg38.tab.gz
wget -c https://vastdb.crg.eu/downloads/hg38/EVENTID_to_GENEID-hg38.tab.gz
wget -c https://vastdb.crg.eu/downloads/hg38/EVENT_INFO-hg38.tab.gz

echo "Uncompressing..."
gunzip -f *.gz

# ------------------------- #
# 2. Reformat SPLICE_SITE_SCORES to 5-column format
# ------------------------- #
echo "Reformatting SPLICE_SITE_SCORES to 5 columns..."
echo -e "EVENT\tss3_seq\tss3_strength\tss5_seq\tss5_strength" > SPLICE_SITE_SCORES_5col-hg38.tab

awk -F'\t' '
NR==1 { next }   # skip header
{
    event=$1
    seq=$2
    score=$3
    len=length(seq)

    if (len==23) {
        ss3_seq[event]=seq
        ss3_strength[event]=score
    }
    else if (len==9) {
        ss5_seq[event]=seq
        ss5_strength[event]=score
    }
}
END {
    PROCINFO["sorted_in"]="@ind_str_asc"  # optional: sort by event
    for (e in ss3_seq) {
        print e, ss3_seq[e], ss3_strength[e], ss5_seq[e], ss5_strength[e]
    }
}
' OFS="\t" SPLICE_SITE_SCORES-hg38.tab >> SPLICE_SITE_SCORES_5col-hg38.tab

# ------------------------- #
# 3. Extract selected EVENT INFO columns
# ------------------------- #

echo "Extracting selected columns from EVENT_INFO..."

# --- Helper: get column index by name ---
get_col_index() {
    local colname="$1"
    local file="$2"
    local idx

    idx=$(head -1 "$file" | tr '\t' '\n' | nl -w1 -s$'\t' \
          | awk -v name="$colname" '$2 == name {print $1}')

    if [[ -z "$idx" ]]; then
        echo "ERROR: Column '$colname' not found in $file" >&2
        exit 1
    fi

    echo "$idx"
}

# --- Detect column indexes dynamically ---
GENE_COL=$(get_col_index "GENE" EVENT_INFO-hg38.tab)
EVENT_COL=$(get_col_index "EVENT" EVENT_INFO-hg38.tab)
COORD_COL=$(get_col_index "COORD_o" EVENT_INFO-hg38.tab)

echo "  GENE column:       $GENE_COL"
echo "  EVENT column:      $EVENT_COL"
echo "  COORD column:      $COORD_COL"


# --- Extract only selected columns ---
awk -v g="$GENE_COL" \
    -v e="$EVENT_COL" \
    -v c="$COORD_COL"  '
    
BEGIN {
    OFS="\t"
    print "GENE","EVENT","COORD"
}
NR > 1 {
    print $g, $e, $c
}
' EVENT_INFO-hg38.tab > EVENT_INFO_selected-hg38.tab

echo "EVENT_INFO_selected-hg38.tab created."

# ------------------------- #
# 4. Extract selected PSI columns
# ------------------------- #

echo "Extracting selected columns from PSI_TABLE..."

# --- Detect column indexes dynamically ---
EVENT_COL=$(get_col_index "EVENT" PSI_TABLE-hg38.tab)
CL293T_COL=$(get_col_index "CL_293T" PSI_TABLE-hg38.tab)
CL293TQ_COL=$(get_col_index "CL_293T-Q" PSI_TABLE-hg38.tab)

echo "  EVENT column:      $EVENT_COL"
echo "  CL_293T column:    $CL293T_COL"
echo "  CL_293T-Q column:  $CL293TQ_COL"

# --- Extract only selected columns ---
awk -v ev="$EVENT_COL" \
    -v t1="$CL293T_COL" \
    -v t2="$CL293TQ_COL" '
BEGIN {
    OFS="\t"
    print "EVENT","CL_293T","CL_293T_Q"
}
NR > 1 {
    print $ev, $t1, $t2
}
' PSI_TABLE-hg38.tab > PSI_TABLE_selected-hg38.tab

echo "PSI_TABLE_selected-hg38.tab created."

echo "removing intermediate files"

gzip PSI_TABLE-hg38.tab EVENT_INFO-hg38.tab 
rm SPLICE_SITE_SCORES-hg38.tab

echo "Done!"


