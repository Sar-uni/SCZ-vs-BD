import pandas as pd

# from metadata, keep only samples present in junction_counts.gz
counts_cols = pd.read_csv("results/leafcutter2/scz_vs_bip.junction_counts.gz", sep=r"\s+", nrows=0).columns[1:]
meta = pd.read_csv("data/SraRunTable.csv")

# deduplicate
meta = meta.drop_duplicates(subset=["Sample Name"]).copy()

# Map diagnosis and filter to scz vs bd + samples in counts
diag_map = {"Schizophrenia": "Schizophrenia", "Bipolar Disorder": "Bipolar"}
meta["condition"] = meta["clinical_diagnosis"].map(diag_map)
meta = meta[meta["condition"].notna() & meta["Sample Name"].isin(counts_cols)]

# select condition & covariates 
cols = ["Sample Name", "condition", "gender", "age_at_death", "post-mortem_interval", "Brain_pH"]
out = meta[cols].dropna()

out.to_csv("results/leafcutter2/groups_file.txt", sep="\t", header=False, index=False)
