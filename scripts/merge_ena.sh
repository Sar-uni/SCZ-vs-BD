mkdir -p merged_fastqs
for dir in raw_fastqs/GSM*; do
    gsm=$(basename "$dir")
    cat $(ls "$dir"/*_1.fastq.gz | sort) > "merged_fastqs/${gsm}_R1.fastq.gz"
    cat $(ls "$dir"/*_2.fastq.gz | sort) > "merged_fastqs/${gsm}_R2.fastq.gz"
done
