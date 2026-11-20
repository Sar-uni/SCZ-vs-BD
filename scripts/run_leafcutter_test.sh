#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/.."

# Paths (GONNA CHANGE TO full .GTF and UNZIP .FA *)
JUNC_MANIFEST="results/junc_files.txt"
OUT_DIR="results/leafcutter2_test"
GROUPS_FILE="${OUT_DIR}/groups_file.txt"
GTF_FILE="ref/gencode/chr21.gtf"   
FASTA_FILE="ref/gencode/chr21.fa"   
EXON_FILE="ref/gencode/gencode_v50_exons.txt.gz"    

mkdir -p "${OUT_DIR}"

# Intron Clustering + functional classification: 
#   junctions annotated as productive, unproductive, ambiguous, or intergenic based on GTF and FASTA
echo "Clustering and classifying introns using leafcutter2"
leafcutter2 \
    -j "${JUNC_MANIFEST}" \
    -r "${OUT_DIR}" \
    -o test_chr21 \
    -A "${GTF_FILE}" \
    -G "${FASTA_FILE}" \
    --maxintronlen 500000 \
    --minclureads 5 \
    --mincluratio 0.0001

# Leafcuter2 gives .junction_counts.gz (unlike perind)
COUNTS_FILE="${OUT_DIR}/test_chr21.junction_counts.gz"

# Differential Splicing
# tab-separate <sample_name> (GON CHANGE TO GSM FOR FINAL *) <condition> <covariates*>
cat <<EOF > "${GROUPS_FILE}"
SRR3438667	Bipolar
SRR3438668	Schizophrenia
SRR3438633	Schizophrenia
SRR3438631	Bipolar
EOF

# SMALL VALUES TEST *
# requires exon file to map clusters to gene names
echo "--leafcutter-ds"
leafcutter-ds \
    --num_threads 4 \
    --min_samples_per_intron 2 \
    --min_samples_per_group 2 \
    --min_coverage 5 \
    --min_unique_vals 2 \
    --exon_file "${EXON_FILE}" \
    --output_prefix "${OUT_DIR}/bip_vs_scz" \
    "${COUNTS_FILE}" \
    "${GROUPS_FILE}"

echo "Cluster significance : ${OUT_DIR}/bip_vs_scz_cluster_significance.txt"
echo "Effect sizes (dPSI)  : ${OUT_DIR}/bip_vs_scz_effect_sizes.txt"
