#!/bin/bash
set -euo pipefail

mkdir -p data/merged_fastqs

for dir in data/raw_fastqs/GSM*; do
    [[ -d "$dir" ]] || continue #empty check
    gsm=$(basename "$dir")

    # R1 and R2 files into arrays
    mapfile -t r1_files < <(printf '%s\n' "$dir"/*_1.fastq.gz | sort)
    mapfile -t r2_files < <(printf '%s\n' "$dir"/*_2.fastq.gz | sort)

    # Skip if no matching FASTQs exist
    if [[ ! -e "${r1_files[0]}" || ! -e "${r2_files[0]}" ]]; then
        echo "NO R1/R2 files in ${dir}; skipping." >&2
        continue
    fi

    out_r1="data/merged_fastqs/${gsm}_R1.fastq.gz"
    out_r2="data/merged_fastqs/${gsm}_R2.fastq.gz"

    if [[ ${#r1_files[@]} -eq 1 ]]; then
        echo "Moving single run for ${gsm}..."
        mv "${r1_files[0]}" "$out_r1"
        mv "${r2_files[0]}" "$out_r2"
    else
        echo "Merging ${#r1_files[@]} runs for ${gsm}..."
        cat "${r1_files[@]}" > "$out_r1"
        cat "${r2_files[@]}" > "$out_r2"

        # verify merged exist before removing raw 
        if [[ -s "$out_r1" && -s "$out_r2" ]]; then
            rm -f "${r1_files[@]}" "${r2_files[@]}"
        else
            echo "Error: Merged output empty for ${gsm}; keeping raw files." >&2
            exit 1
        fi
    fi

    # Remove empty GSM directory if nothing remains inside
    rmdir "$dir" 2>/dev/null || true
done
