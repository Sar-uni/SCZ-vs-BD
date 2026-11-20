#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/.."

# Paths
GENCODE_DIR="ref/gencode"
GENOME_DIR="ref/star_index"
BASE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50"

FASTA_GZ="${GENCODE_DIR}/GRCh38.primary_assembly.genome.fa.gz"
GTF_GZ="${GENCODE_DIR}/gencode.v50.primary_assembly.annotation.gtf.gz"
FASTA="${GENCODE_DIR}/GRCh38.primary_assembly.genome.fa"
GTF="${GENCODE_DIR}/gencode.v50.primary_assembly.annotation.gtf"
THREADS=16

mkdir -p "${GENCODE_DIR}" "${GENOME_DIR}" 

# Download FASTA and GTF (if not exists)
if [[ ! -f "${FASTA}" && ! -f "${FASTA_GZ}" ]]; then
    echo "--- GRCh38 primary assembly FASTA download"
    wget -c -P "${GENCODE_DIR}" "${BASE_URL}/GRCh38.primary_assembly.genome.fa.gz"
fi

if [[ ! -f "${GTF}" && ! -f "${GTF_GZ}" ]]; then
    echo "--- GENCODE v50 primary assembly GTF download"
    wget -c -P "${GENCODE_DIR}" "${BASE_URL}/gencode.v50.primary_assembly.annotation.gtf.gz"
fi

# Decompress FASTA
if [[ ! -f "${FASTA}" ]]; then
    echo "--- Extracting ${FASTA_GZ}"
    gunzip -k "${FASTA_GZ}"
else
    echo "${FASTA} already extracted skip"
fi

# Decompress GTF 
if [[ ! -f "${GTF}" ]]; then
    echo "-- Extracting ${GTF_GZ}"
    gunzip -k "${GTF_GZ}"
else
    echo "--- ${GTF} already extracted, skip "
fi

#  STAR Genome Index
echo "Building STAR index in ${GENOME_DIR} ==="
STAR \
    --runThreadN "${THREADS}" \
    --runMode genomeGenerate \
    --genomeDir "${GENOME_DIR}" \
    --genomeFastaFiles "${FASTA}" \
    --sjdbGTFfile "${GTF}" \
    --sjdbOverhang 49

echo "index build complete!"
