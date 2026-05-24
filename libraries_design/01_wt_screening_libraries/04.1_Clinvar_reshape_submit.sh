#!/bin/bash
# Download ClinVar and proecess the table
# need Python3 with gzip 
# Usage: bash 04.1_Clinvar_reshape_submit.sh /path/to/base_dir

set -e

if [ -z "$1" ]; then
    echo "Usage: $0 /path/to/output_dir"
    exit 1
fi

BASE_DIR="$1"
OUTDIR="$BASE_DIR/libraries_design/01_wt_screening_libraries/necessary_file/ClinVar"

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# ------------------------- #
# 1. Download VastDB files
# ------------------------- #
echo "Downloading ClinVar version of 20230702"

wget -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2023/clinvar_20230702.vcf.gz

echo "Extracting relevant columns"

cd "$BASE_DIR/libraries_design/01_wt_screening_libraries/"

python3 04.1_Clinvar_reshape.py3 "$OUTDIR"

    