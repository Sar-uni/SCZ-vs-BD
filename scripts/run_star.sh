#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/.."

# Paths
GENOME_DIR="test_sample/ref/chr21_star_index"
RAW_DIR="test_sample/samples"
SUB_DIR="test_sample/subsampled"
OUT_STAR="test_sample/star_out"
JUNC_DIR="results/junc"
COUNTS_DIR="results/gene_counts"
QC_DIR="results/qc"

mkdir -p "${SUB_DIR}" "${OUT_STAR}" "${JUNC_DIR}" "${COUNTS_DIR}" "${QC_DIR}/fastqc"

# test samples
SAMPLES=("SRR3438667" "SRR3438668" "SRR3438633" "SRR3438631")

# clear file
> results/junc_files.txt

for SAMPLE in "${SAMPLES[@]}"; do
    echo "=== Processing ${SAMPLE} ==="

    # FastQc on samples
    fastqc \
        -t 4 \
        --outdir "${QC_DIR}/fastqc" \
        "${SUB_DIR}/${SAMPLE}/${SAMPLE}_1.sub.fastq.gz" \
        "${SUB_DIR}/${SAMPLE}/${SAMPLE}_2.sub.fastq.gz"
    
    # star alingment
    STAR \
        --runThreadN 4 \
        --genomeDir "${GENOME_DIR}" \
        --readFilesIn "${SUB_DIR}/${SAMPLE}/${SAMPLE}_1.sub.fastq.gz" "${SUB_DIR}/${SAMPLE}/${SAMPLE}_2.sub.fastq.gz" \
        --readFilesCommand zcat \
        --twopassMode Basic \
        --quantMode GeneCounts \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMstrandField intronMotif \
        --outSAMattributes NH HI AS nM XS \
        --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
        --outFileNamePrefix "${OUT_STAR}/${SAMPLE}_"

    # copy gene counts to results 
    cp "${OUT_STAR}/${SAMPLE}_ReadsPerGene.out.tab" "${COUNTS_DIR}/${SAMPLE}_ReadsPerGene.out.tab"

    # Index BAM and extract junctions
    BAM_FILE="${OUT_STAR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    samtools index -@ 4 "${BAM_FILE}"

    regtools junctions extract \
        -a 8 -m 50 -M 500000 -s XS \
        "${BAM_FILE}" \
        -o "${JUNC_DIR}/${SAMPLE}.junc"

    # add path for LeafCutter manifest
    echo "${JUNC_DIR}/${SAMPLE}.junc" >> results/junc_files.txt
done

# combine FastQC and STAR logs
multiqc \
    "${QC_DIR}/fastqc" \
    "${OUT_STAR}" \
    -o "${QC_DIR}/multiqc" \
    --force
