# produces RAW counts for expression. 
# Run before prepare_expression_matrix.py
from pathlib import Path
import pandas as pd

def load_star_counts(meta,count_dir="results/gene_counts"):
    """Read STAR gene.out.tab files (unstranded col).
       Returns samples x genes.
    """

    count_dir = Path(count_dir)
    samples = {}

    for sample in meta.index:
        path = count_dir / f"{sample}_ReadsPerGene.out.tab"
        if not path.exists():
            continue

        df = pd.read_csv(path,sep="\t",skiprows=4,header=None,index_col=0)

        counts = pd.to_numeric(df.iloc[:, 0],errors="raise")

        samples[sample] = counts

    matrix = pd.DataFrame(samples).T
    matrix.index.name = "sample"
    matrix.columns.name = "gene_id"

    return matrix

if __name__ == "__main__":
    meta = pd.read_csv("results/classification/metadata.csv", index_col=0)
    meta.index = meta.index.astype(str).str.strip()

    counts = load_star_counts(meta)

    counts.to_csv("results/classification/expression_counts.csv")
