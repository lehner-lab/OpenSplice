#!/bin/bash
# Script to download and prepare hg38 genome fasta + index
# Usage: bash 02.1_download_hg38_genome_fasta_and_index.sh /path/to/base_dir
# need samtools and bgzip

# Exit immediately if a command fails
set -e

# Check if base directory is provided
if [ -z "$1" ]; then
    echo "Usage: $0 /path/to/base_dir"
    exit 1
fi

BASE_DIR="$1"
GENOME_DIR="$BASE_DIR/libraries_design/01_wt_screening_libraries/necessary_file/hg38"

# Create the directory if it doesn't exist
mkdir -p "$GENOME_DIR"

# Move to genome directory
cd "$GENOME_DIR"

# Download the FASTA
echo "Downloading hg38 fasta..."
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# Unzip
echo "Unzipping..."
gunzip -f hg38.fa.gz

# Index fasta with samtools
echo "Creating fasta index..."
samtools faidx hg38.fa

# Compress fasta index with bgzip
echo "Compressing fasta index..."
bgzip -c hg38.fa.fai > hg38.fa.fai.bgz

# Compress fasta with bgzip
echo "Compressing fasta..."
bgzip -c hg38.fa > hg38.fa.bgz

rm hg38.fa hg38.fa.fai

echo "Done. Genome files prepared in $GENOME_DIR"
echo "Proceed with 02_get_exon_intron_sequences.R"
