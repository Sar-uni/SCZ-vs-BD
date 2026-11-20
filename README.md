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
    - `data/sample_to_srr.tsv` maps GSM to SRR runs

## Pre-reqs:
- star, samtools, regtools, fastqc, multiqc installed (used conda to create env)
- STAR requires genome index in `ref/star_index`
    - requires **unzipped** gencode fasta and gtf (`ref/gencode`)
    - Unzip command:
    ```
    gunzip GRCh38.primary_assembly.genome.fa.gz gencode.v50.primary_assembly.annotation.gtf.gz
    ```
- merged .fastq.gz files in `data/merged_fastqs` (created by `scripts/download_data.sh`)

## Alignment Steps
1. start with working directory in project root or `scripts/` folder
2. download dataset via `scripts/download_data.sh`
3. Merge .fastq files with saem GSM IDs with `scripts/merge_ena.sh`
4. if index doesn't exist, create index by running `scripts/build_index.sh`
5. Begin alignment by running `scripts/run_star.sh`

## cluster notes
- paths (relative used in test)
- thread counts

## Project Directory Structure

```text
project_root/
├── data/
│   ├── samplesheet.tsv          # GSM mapped to ENA ftp links
│   ├── raw_fastqs/              
│   └── merged_fastqs/           # Merged R1 R2 FASTQs
├── ref/
│   ├── GRCh38.primary_assembly.genome.fa
│   ├── gencode.v50.primary_assembly.annotation.gtf
│   └── star_index/              # Full STAR genome index
├── scripts/
│   ├── download_data.sh         # data/raw_fastqs/GSM*/
│   ├── merge_ena.sh             # data/merged_fastqs/
│   ├── build_index.sh           # ref/star_index/
│   └── run_star.sh              # Main
├── results/
│   ├── star_out/                # Sorted BAMs and STAR logs
│   ├── junc/                    # splice junctions ({GSM}.junc)
│   ├── gene_counts/             # STAR gene-level read counts ({GSM}_ReadsPerGene.out.tab)
│   ├── qc/                      
│   └── junc_files.txt
```

## Post aligntment
- LeafCutter
    - git leafcutter repo: git clone https://github.com/davidaknowles/leafcutter
