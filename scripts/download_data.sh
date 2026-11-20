#!/bin/bash
set -euo pipefail

SAMPLESHEET="samplesheet.tsv"
BASE_DIR="raw_fastqs"

mkdir -p "$BASE_DIR"

# parse TSV lines
tail -n +2 "$SAMPLESHEET" | while IFS=$'\t' read -r gsm sample_title diagnosis runs r1_urls r2_urls; do
    # skip empty lines
    [ -z "$gsm" ] && continue

    sample_dir="${BASE_DIR}/${gsm}"
    mkdir -p "$sample_dir"
    echo "=== Processing $gsm ($sample_title | Diagnosis: $diagnosis) ==="

    # Download R1 files 
    for url in $r1_urls; do
        # clean url
        # Strip any carriage returns and ensure ftp:// prefix
        clean_url=$(echo "$url" | tr -d '\r')
        [[ "$clean_url" != ftp://* && "$clean_url" != http* ]] && clean_url="ftp://${clean_url}"
        
        echo "  [R1] Downloading $(basename "$clean_url")..."
        wget -c -q -P "$sample_dir" "$clean_url"
    done

    # R2 files download
    for url in $r2_urls; do
        clean_url=$(echo "$url" | tr -d '\r')
        [[ "$clean_url" != ftp://* && "$clean_url" != http* ]] && clean_url="ftp://${clean_url}"
        
        echo "  [R2] Downloading $(basename "$clean_url")..."
        wget -c -q -P "$sample_dir" "$clean_url"
    done
done

echo "downloads complete!"
