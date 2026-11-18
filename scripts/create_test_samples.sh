#!/bin/bash
set -eu

# to project root
cd "$(dirname "$0")/.."

RAW_DIR="test_sample/samples"
SUB_DIR="test_sample/subsampled"

# LINES taken by bash :( using READS instead
READS=1000000

# go over all SRR folders and extract the first 250k reads from each fastq.gz file
for SAMPLE_PATH in "${RAW_DIR}"/SRR*; do

    SRR=$(basename "${SAMPLE_PATH}")
    OUT_FOLDER="${SUB_DIR}/${SRR}"
    mkdir -p "${OUT_FOLDER}"

    echo "subsampling begins for ${SRR}"

    # inside SRR folder, loop over _1 and _2 fast.gz files
    for FQ_FILE in "${SAMPLE_PATH}"/*.fastq.gz; do

        BASE=$(basename "${FQ_FILE}" .fastq.gz)
        # save subsampled files as .sub.fastq.gz in subsampled/
        OUT_FILE="${OUT_FOLDER}/${BASE}.sub.fastq.gz"
        zcat "${FQ_FILE}" | head -n ${READS} | gzip > "${OUT_FILE}"
        echo "  Wrote $(zcat "${OUT_FILE}" | wc -l) lines to ${OUT_FILE}"
    done
done