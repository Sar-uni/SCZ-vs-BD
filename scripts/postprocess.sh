#!/bin/bash
#SBATCH --job-name=star_post
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00

set -euo pipefail
cd "$(dirname "$0")/.."

ls results/junc/*.junc | sort > results/junc_files.txt

multiqc results/qc/fastqc results/star_out -o results/qc/multiqc --force
