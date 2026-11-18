conda activate splicing_test

conda create -n splicing_test -c conda-forge -c bioconda \
    star \
    samtools \
    regtools \
    fastqc \
    fastp \
    r-base \
    r-essentials