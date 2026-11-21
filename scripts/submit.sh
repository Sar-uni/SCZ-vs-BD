#!/bin/bash
set -euo pipefail

# start from project root
cd "$(dirname "$0")/.."
mkdir -p logs

# merge into per GSM FASTQs
jid0=$(sbatch --parsable scripts/merge_ena.sh)
echo "Submitted merge_ena.      ID: ${jid0}"

# index build job (--parsable makes it so ID returned)
jid1=$(sbatch --parsable scripts/build_index.sh)
echo "Submitted build_index.    ID: ${jid1}"

# submit STAR array job
jid2=$(sbatch --parsable --dependency=afterok:${jid1} scripts/run_star.sh)
echo "Submitted run_star.       ID: ${jid2}"

# postprcocess job
jid3=$(sbatch --parsable --dependency=afterok:${jid2} scripts/postprocess.sh)
echo "Submitted postprocess.    ID: ${jid3}"