conda activate splicing_test

conda create -n splicing_test -c conda-forge -c bioconda \
    multiqc \
    star \
    samtools \
    regtools \
    fastqc \
    fastp \
    r-base \
    r-essentials

# leafcutter2 rep
git clone https://github.com/leafcutter2/leafcutter2.git
cd leafcutter2

conda env create -f leafcutter2_env.yml
conda activate leafcutter2

# leafcutter2 package
pip install .
cd ..

git clone https://github.com/leafcutter2/leafcutter-ds.git
cd leafcutter-ds

# leafcutter-ds
TMPDIR=~/pip_tmp python -m pip install --no-cache-dir .
cd ..

mkdir -p ~/pip_tmp
TMPDIR=~/pip_tmp pip install --no-cache-dir leafcutter
rm -rf ~/pip_tmp