#!/bin/bash
#SBATCH --job-name=lc2_post
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/lc2_post_%j.log

set -euo pipefail
cd "$(dirname "$0")/.."   

# *requires leafcutter2 environment

# Paths
JUNC_DIR="results/junc"
JUNC_MANIFEST="results/junc_files.txt"
OUT_DIR="results/leafcutter2_test"
GROUPS_FILE="results/leafcutter2/groups_file.txt"
SAMPLESHEET="data/samplesheet.tsv"
SAMPLES_TXT="data/samples.txt"
GROUP_SCRIPT="scripts/build_groupfile.py"

GTF="ref/gencode/chr21.gtf"
FASTA="ref/gencode/chr21.fa"
EXON_FILE="ref/gencode/gencode_v50_exons.txt.gz"

THREADS=16
LC2_PREFIX="scz_vs_bip"
DS_PREFIX="${OUT_DIR}/scz_vs_bip"

mkdir -p "${OUT_DIR}" logs

# feasta index required by lc2
if [[ ! -f "${FASTA}.fai" ]]; then
  echo "samtool indexing FASTA"
  samtools faidx "${FASTA}"
fi

# exon table from GTF (need for gene labels in DS)
if [[ ! -f "${EXON_FILE}" ]]; then
  echo "-Creating exon table"
  leafcutter-gtf-to-exons "${GTF}" "${EXON_FILE%.gz}"
  gzip -f "${EXON_FILE%.gz}"
fi

 
echo "==Run leafcutter2 (full GENCODE v50)"
leafcutter2 \
  -j "${JUNC_MANIFEST}" \
  -r "${OUT_DIR}" \
  -o "${LC2_PREFIX}" \
  -A "${GTF}" \
  -G "${FASTA}" \
  --minclureads 2 \
  --mincluratio 0.0001

COUNTS_FILE="${OUT_DIR}/${LC2_PREFIX}.junction_counts.gz"
[[ -s "${COUNTS_FILE}" ]] || {
  echo "ERROR: missing ${COUNTS_FILE}" >&2
  exit 1
}

# build group file with covariates
python3 ${GROUP_SCRIPT}

echo "Running leafcutter-ds (baseline = Schizophrenia)"
leafcutter-ds \
  --num_threads "${THREADS}" \
  --baseline_group Schizophrenia \
  --min_samples_per_intron 1 \
  --min_samples_per_group 1 \
  --min_coverage 2 \
  --min_unique_vals 2 \
  --exon_file "${EXON_FILE}" \
  --output_prefix "${DS_PREFIX}" \
  "${COUNTS_FILE}" \
  "${GROUPS_FILE}"

 
echo "--leafcutter complete"