#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/.."

# paths
LEAFCUTTER_DIR="${HOME}/tools/leafcutter"
JUNC_MANIFEST="results/junc_files.txt"
OUT_DIR="results/leafcutter_test"
GROUPS_FILE="${OUT_DIR}/groups_file.txt"

mkdir -p "${OUT_DIR}"


# Intron Clustering (lowere m and p for now)
echo "1. clustering introns using leafcutter.regtools"
python3 "${LEAFCUTTER_DIR}/clustering/leafcutter_cluster_regtools.py" \
    -j "${JUNC_MANIFEST}" \
    -r "${OUT_DIR}" \
    -o test_chr21 \
    -l 500000 \
    -m 5 \
    -p 0.0001

COUNTS_FILE="${OUT_DIR}/test_chr21_perind_numers.counts.gz"

# differential splicing 
# test group file
# format (tab seperte): <sample_name> (need to match .junc filenames) <condition>
cat <<EOF > "${GROUPS_FILE}"
SRR3438667	Bipolar
SRR3438668	Schizophrenia
SRR3438633	Schizophrenia
SRR3438631	Bipolar
EOF
# TODO: fix package dependancy issue with leafcutter_ds.R. (leafcutter is incomplete package)
# -i -g set to 2 for 
echo "2. DES leafcutter_ds.R"
Rscript "${LEAFCUTTER_DIR}/scripts/leafcutter_ds.R" \
    --num_threads 4 \
    --min_samples_per_intron 2 \
    --min_samples_per_group 2 \
    --output_prefix "${OUT_DIR}/bip_vs_scz" \
    "${COUNTS_FILE}" \
    "${GROUPS_FILE}"

echo "Cluster significance : ${OUT_DIR}/bip_vs_scz_cluster_significance.txt"
echo "Effect sizes (dPSI)  : ${OUT_DIR}/bip_vs_scz_effect_sizes.txt"
