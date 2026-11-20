#!/bin/bash

set -euo pipefail

cd "$(dirname "$0")/.."

# Paths
GENOME_DIR="ref/star_index"
MERGED_DIR="data/merged_fastqs"
OUT_STAR="results/star_out"
JUNC_DIR="results/junc"
COUNTS_DIR="results/gene_counts"
QC_DIR="results/qc"
THREADS=8

mkdir -p "${OUT_STAR}" "${JUNC_DIR}" "${COUNTS_DIR}" "${QC_DIR}/fastqc" 

# clear LeafCutter manifest file
> results/junc_files.txt

# process merged files
for R1_FILE in "${MERGED_DIR}"/*_R1.fastq.gz; do
    # get GSM (merged_fastqs/GSMXX_R1.fastq.gz -> GSMXX)
    SAMPLE=$(basename "${R1_FILE}" "_R1.fastq.gz")
    R2_FILE="${MERGED_DIR}/${SAMPLE}_R2.fastq.gz"

    echo "---Processing ${SAMPLE}"

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
        --outSAMstrandField intronMotif \
        --outSAMattributes NH HI AS nM XS \
        --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
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
done

# combine FastQC and STAR logs
multiqc \
    "${QC_DIR}/fastqc" \
    "${OUT_STAR}" \
    -o "${QC_DIR}/multiqc" \
    --force
