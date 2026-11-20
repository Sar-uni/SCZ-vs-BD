# Merges the cluster_significance.txt and effect_sizes.txt outputs from 
# Leafcutter2 with the DESeq2 results, and filters for NMD downregulation candidates.

import gzip
import re
import pandas as pd

# load leafcutter-ds outputs
sig = pd.read_csv("results/leafcutter2_test/bip_vs_scz_cluster_significance.txt", sep="\t")
eff = pd.read_csv("results/leafcutter2_test/bip_vs_scz_effect_sizes.txt", sep="\t")

# rename 'deltapsi_smth' (ex. deltapsi_Schizophrenia) to 'deltapsi'
dpsi_cols = [c for c in eff.columns if c.startswith("deltapsi")]
if not dpsi_cols:
    raise KeyError(f"No deltapsi column found in: {list(eff.columns)}")
eff = eff.rename(columns={dpsi_cols[0]: "deltapsi"})

# extract matching cluster ID and lc2 tag (PR, UP, NE, IN)
# eg. intron str: "chr21:33031935:33036103:clu_104_+:UP"
parts = eff["intron"].str.split(":", expand=True)
eff["cluster"] = parts[0] + ":" + parts[3]
eff["func_type"] = parts[4] if parts.shape[1] > 4 else None

# join signif + effect sizes, keep only success tests
lc_merged = pd.merge(sig, eff, on="cluster", how="inner").query("status == 'Success'")

if "genes" not in lc_merged.columns or lc_merged["genes"].isna().all():
    raise ValueError(
        "The 'genes' column is empty or missing in cluster_significance.txt. "
    )

lc_merged = (
    lc_merged.dropna(subset=["genes"])
    .assign(genes=lambda df: df["genes"].astype(str).str.split(","))
    .explode("genes")
    .assign(gene_symbol=lambda df: df["genes"].str.strip())
    .rename(columns={"p.adjust": "lc_padj"})
)

# Map Ensembl gene_id to gene_symbol from gencode GTF (col 8 is metadata pairs)
gtf_records = []
with gzip.open("ref/gencode/gencode.v50.primary_assembly.annotation.gtf.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        cols = line.split("\t")
        if cols[2] == "gene":
            gid = re.search(r'gene_id "([^"]+)"', cols[8])
            gname = re.search(r'gene_name "([^"]+)"', cols[8])
            if gid and gname:
                gtf_records.append((gid.group(1), gname.group(1)))

gtf_df = pd.DataFrame(gtf_records, columns=["gene_id", "gene_symbol"]).drop_duplicates()
print(f"Found {len(gtf_df)} unique gene_id to gene_symbol mappings")

# load DESeq2 .CSV and attach gene symbols
deseq = (
    pd.read_csv("results/deseq2_scz_vs_bip.csv", index_col=0)
    .reset_index()
    .rename(columns={"index": "gene_id", "padj": "deseq_padj"})
    .merge(gtf_df, on="gene_id", how="inner")
)

# Merge lc2 and deseq2 by gene_symbol
combined = pd.merge(lc_merged, deseq, on="gene_symbol", how="inner")

# filter nmd downregulation events
nmd_candidates = combined[
    (combined["lc_padj"] < 0.05) &
    (combined["func_type"] == "UP") &
    (combined["deltapsi"] > 0.05) &
    (combined["log2FoldChange"] < 0)
].sort_values("lc_padj")

combined.to_csv("results/leafcutter2_deseq2_integrated.csv", index=False)
nmd_candidates.to_csv("results/nmd_downregulated_candidates.csv", index=False)
print(f"Merged {len(combined)} total rows; found {len(nmd_candidates)} NMD candidates.")
