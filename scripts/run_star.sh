#!/bin/bash
#SBATCH --job-name=star_align
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=1-47%6
#SBATCH --output=logs/star_%A_%a.log

set -euo pipefail
cd "$(dirname "$0")/.."

GENOME_DIR="ref/star_index"
MERGED_DIR="data/merged_fastqs"
SAMPLES_FILE="data/samples.txt"
OUT_STAR="results/star_out"
JUNC_DIR="results/junc"
COUNTS_DIR="results/gene_counts"
QC_DIR="results/qc"
THREADS=8

# pick a sample for this array task
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLES_FILE}")
[[ -n "${SAMPLE}" ]] || { echo "No sample at line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }

echo "--- Processing ${SAMPLE} (array task ${SLURM_ARRAY_TASK_ID})"

R1_FILE="${MERGED_DIR}/${SAMPLE}_R1.fastq.gz"
R2_FILE="${MERGED_DIR}/${SAMPLE}_R2.fastq.gz"

# FastQC on merged samples
fastqc \
    -t "${THREADS}" \
    --outdir "${QC_DIR}/fastqc" \
    "${R1_FILE}" \
    "${R2_FILE}"

# STAR alignment
STAR \
    --runThreadN "${THREADS}" \
    --genomeDir "${GENOME_DIR}" \
    --readFilesIn "${R1_FILE}" "${R2_FILE}" \
    --readFilesCommand zcat \
    --twopassMode Basic \
    --quantMode GeneCounts \
    --outSAMtype BAM SortedByCoordinate \
    --limitBAMsortRAM 20000000000 \
    --outSAMstrandField intronMotif \
    --outSAMattributes NH HI AS nM XS \
    --outFileNamePrefix "${OUT_STAR}/${SAMPLE}_"

# copy gene counts to results 
cp "${OUT_STAR}/${SAMPLE}_ReadsPerGene.out.tab" "${COUNTS_DIR}/${SAMPLE}_ReadsPerGene.out.tab"

# index BAM and extract junctions
BAM_FILE="${OUT_STAR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
samtools index -@ "${THREADS}" "${BAM_FILE}"

regtools junctions extract \
    -a 8 -m 50 -M 500000 -s XS \
    "${BAM_FILE}" \
    -o "${JUNC_DIR}/${SAMPLE}.junc"

# add path for LeafCutter manifest
echo "${JUNC_DIR}/${SAMPLE}.junc" >> results/junc_files.txt

# rm large BAM files
if [[ -s "${COUNTS_DIR}/${SAMPLE}_ReadsPerGene.out.tab" && -s "${JUNC_DIR}/${SAMPLE}.junc" ]]; then
    rm -f "${BAM_FILE}" "${BAM_FILE}".bai "${BAM_FILE%.bam}".bai
else
    echo "Oops: no ${SAMPLE}; keeping BAM." >&2
fi

# rm tmp STAR files
rm -rf "${OUT_STAR}/${SAMPLE}__STARpass1" \
       "${OUT_STAR}/${SAMPLE}__STARtmp" \
       "${OUT_STAR}/${SAMPLE}__STARgenome"