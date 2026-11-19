# SCZ-vs-BD

External libraries needed:

- numpy
- pandas
- scikit-learn

## Sequencing Overview & Links
Downloaded raw FASTQ files for DLPFC SCZ (n=29) and BD (n=28) samples
from ENA, then ran STAR tool for alignment

GENCODE fasta and .gtf files (for creating genome index)
- https://www.gencodegenes.org/human/ 

Pre-Paired FastQ files from ENA:
- https://www.ebi.ac.uk/ena/browser/view/PRJNA319583
- paired-end fastq files readily available for each run

GEO SRA Accession #s for DLPFC (B & S) Samples
- https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA319583&f=brain_region_sam_ss%3An%3Adlpfc%3Bclinical_diagnosis_sam_ss%3An%3Aschizophrenia%2Cbipolar%2520disorder&o=bytes_l%3Aa%3Bacc_s%3Aa&s=SRR3438613,SRR3438614,SRR3438629,SRR3438630,SRR3438633,SRR3438634,SRR3438668,SRR3438713,SRR3438714,SRR3438727,SRR3438742,SRR3438764,SRR3438765,SRR3438766,SRR3438776,SRR3438787,SRR3438790,SRR3438824,SRR3438826,SRR3438835,SRR3438855,SRR3438857,SRR3438873,SRR3438878,SRR3438879,SRR3438882,SRR3438883,SRR3438892,SRR3438896,SRR3438571,SRR3438572,SRR3438577,SRR3438578,SRR3438627,SRR3438628,SRR3438631,SRR3438632,SRR3438665,SRR3438667,SRR3438695,SRR3438696,SRR3438718,SRR3438728,SRR3438743,SRR3438753,SRR3438755,SRR3438789,SRR3438800,SRR3438812,SRR3438825,SRR3438845,SRR3438846,SRR3438847,SRR3438866,SRR3438870,SRR3438895,SRR3438906
    - sometimes different runs (SRRXXXX) are done on the same sample (same SRX or GSM #), so those need to be merged. Here, 20 instances (10 samples) where samples are same:

    | SRR | GSM |
    | --- | --- |
    | SRR3438571 | GSM2132213 |
    | SRR3438572 | GSM2132213 |
    | SRR3438577 | GSM2132216 |
    | SRR3438578 | GSM2132216 |
    | SRR3438613 | GSM2132234 |
    | SRR3438614 | GSM2132234 |
    | SRR3438627 | GSM2132241 |
    | SRR3438628 | GSM2132241 |
    | SRR3438629 | GSM2132242 |
    | SRR3438630 | GSM2132242 |
    | SRR3438631 | GSM2132243 |
    | SRR3438632 | GSM2132243 |
    | SRR3438633 | GSM2132244 |
    | SRR3438634 | GSM2132244 |
    | SRR3438695 | GSM2132287 |
    | SRR3438696 | GSM2132287 |
    | SRR3438713 | GSM2132296 |
    | SRR3438714 | GSM2132296 |
    | SRR3438882 | GSM2132460 |
    | SRR3438883 | GSM2132460 |
## Pre-reqs:
- star, samtools, regtools, fastqc, multiqc installed (used conda to create env)
- STAR requires genome index. Can be generated if not exists:
    - requires gencode fasta and gtf (grchh38, primary assembly)
    - via wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_XX/...
    - sjdboverhang = 49 (read length of 50 used in study)

## Test Run Steps
- Genome index creation
    - Downloaded `GRCh38.primary_assembly.genome.fa.gz` and `gencode.v50.primary_assembly.annotation.gtf.gz`
    - Extracted only annotations on chr21 instead of using entire genome (`chr21.gtf`, `chr21.fa`)
    - Generate index with follow command:
    ```
    STAR --runThreadN 4 \
     --runMode genomeGenerate \
     --genomeDir chr21_star_index \
     --genomeFastaFiles chr21.fa \
     --sjdbGTFfile chr21.gtf \
     --sjdbOverhang 49 \
     --genomeSAindexNbases 11
    ```
- Ran `create_test_samples.sh` to converted full .fastq.gz files to smaller .sub.fastq.gz.
Samples stored in test_sample/subsampled
- Aligned samples by running `run_star.sh`
- Matrices stored in `results/`

## cluster notes
- .junc files appending immediately after each sample (may conflict). can  do after
- paths (relative used in test)
- downloads for paired fastq files using wget from ENA needs merging logic
- thread counts can be increased

## Post aligntment
- LeafCutter
    - git leafcutter repo: git clone https://github.com/davidaknowles/leafcutter
