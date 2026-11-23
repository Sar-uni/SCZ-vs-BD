# SCZ-vs-BD

External libraries needed:

- numpy
- pandas
- scikit-learn

## Sequencing Overview & Links
Downloaded raw FASTQ files for DLPFC SCZ (n=29) and BD (n=28) samples
from ENA, then ran STAR tool for alignment

GENCODE fasta and .gtf (gene transfer format) files (for creating genome index)
- https://www.gencodegenes.org/human/ 
- GTF headers: chr | source | feature type | start | end | confidence | strand | reading frame | metadata pairs 

Pre-Paired FastQ files from ENA:
- https://www.ebi.ac.uk/ena/browser/view/PRJNA319583
- paired-end fastq files readily available for each run

GEO SRA Accession #s for DLPFC (B & S) Samples
- https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA319583&f=brain_region_sam_ss%3An%3Adlpfc%3Bclinical_diagnosis_sam_ss%3An%3Aschizophrenia%2Cbipolar%2520disorder&o=bytes_l%3Aa%3Bacc_s%3Aa&s=SRR3438613,SRR3438614,SRR3438629,SRR3438630,SRR3438633,SRR3438634,SRR3438668,SRR3438713,SRR3438714,SRR3438727,SRR3438742,SRR3438764,SRR3438765,SRR3438766,SRR3438776,SRR3438787,SRR3438790,SRR3438824,SRR3438826,SRR3438835,SRR3438855,SRR3438857,SRR3438873,SRR3438878,SRR3438879,SRR3438882,SRR3438883,SRR3438892,SRR3438896,SRR3438571,SRR3438572,SRR3438577,SRR3438578,SRR3438627,SRR3438628,SRR3438631,SRR3438632,SRR3438665,SRR3438667,SRR3438695,SRR3438696,SRR3438718,SRR3438728,SRR3438743,SRR3438753,SRR3438755,SRR3438789,SRR3438800,SRR3438812,SRR3438825,SRR3438845,SRR3438846,SRR3438847,SRR3438866,SRR3438870,SRR3438895,SRR3438906
    - `data/sample_to_srr.tsv` maps GSM to SRR runs

## Pre-reqs:
- star, samtools, regtools, fastqc, multiqc installed (used conda to create env)
- STAR requires genome index
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
LeafCutter2
- leafcutter2 repos: `https://github.com/leafcutter2/leafcutter2.git` and `https://github.com/leafcutter2/leafcutter-ds.git`
- .GTF file can stay gz, but fa needs to be unzipped and .fai index needs to be created: `samtools faidx GRCh38.primary_assembly.genome.fa`
- default, leafcutter-ds uses the first group appearing in groups_file.txt

Files Table formats
- `.junc` files: chr | start | end | strand | read_count
- `.junction_counts.gz` (splicing feature matrix): intron id | sample columns
- `cluster_significance.txt`: cluster | status | loglr (lg likelihood) | df (degree freddom) | p | p.adj | genes  | annotations (functional tags present in cluster)
- `effect_sizes.txt`: intron | logef_<group> (lg effect size of non-baseline group) | psi_<baseline> | psi_<group> | deltapsi_<group> 
- `.ratios.gz` - intron id | ratios (reads supporting junction / total cluster reads)
    - 0/0 should be NaN (no supporting cluster coverage), not 0 (intron not used)

DESeq2.R 
- `ReadPerGene.out.tab`
    - Summary (first 4 rows): N_unmapped, N_multimapping, N_noFeature, N_ambiguous
    - filter low counts
- GENE ID | Unstranded | Forward | Backward
- `deseq2_scz_vs_bip.csv`: baseMean | log2FoldChange | lfcSE | stat | pvalue | padj

Merging Results
- `leafcutter2_deseq2_integrated.csv:` cluster | status | loglr | df | p | lc_padj | genes | annotations | intron | logef_Schizophrenia | psi_Bipolar | psi_Schizophrenia | deltapsi | func_type | gene_symbol | gene_id | baseMean | log2FoldChange | lfcSE | stat | pvalue | deseq_padj
- compare w/ leafcutter's UP clusters. Increased unproductive/NMD splciing often means downregulation at total gene expression. For a gene that is transcriptionally downregulated in SCZ (lgFoldChange < 0), is there also a shift toward an unproductive (UP) isoform (func_type == "UP", deltapsi > 0)?

## Classification
- create metadata, expression matrix (R), and splicing matrix standardized for GSM
- produce counts with prepare_counts.py -> then prepare_expression.R -> then Classifier.py
    - issue: preprocessing should be using training samples in CV fold (maybe later) 
- compute AUROC for each repeat instead of fold (only 9 samples per fold)
- permutation test: shuffle diagnosis. Peferable to testing if AUROC is diff from 0.5 bc it accounts for the feature-selection and model fitting process.
- DLPFC (brain region) has high alternative splicing, so expect have 3-4x more features
- The 20 repeats reuse the 47 samples. They provide infoo about sensitivity to data partitioning, but they are not 20 independent validation cohorts. mean and standard deviation are useful descriptive summaries. Not treating them as independent observations for significance testing.
- feature selection: a feature selected across many folds is more interesting