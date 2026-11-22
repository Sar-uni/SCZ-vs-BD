import pandas as pd

meta = pd.read_csv("data/SraRunTable.csv")

# Keep DLPFC and diag
meta = meta[
    (meta["brain_region"] == "DLPFC") &
    (meta["clinical_diagnosis"].isin([
        "Schizophrenia",
        "Bipolar Disorder"
    ]))
].copy()

# Rename cols (simplfy)
meta = meta.rename(columns={
    "GEO_Accession (exp)": "sample",
    "clinical_diagnosis": "diagnosis",
    "age_at_death": "age",
    "gender": "sex",
    "post-mortem_interval": "pmi",
    "Brain_pH": "ph"
})

# Standard labels
meta["diagnosis"] = meta["diagnosis"].replace({
    "Schizophrenia": "SCZ",
    "Bipolar Disorder": "BD"
})

meta["sex"] = meta["sex"].str.upper()

# One row per sample
meta = meta.drop_duplicates("sample").set_index("sample")

# Verify sample counts
print(meta["diagnosis"].value_counts())

assert len(meta) == 47
assert (meta["diagnosis"] == "SCZ").sum() == 24
assert (meta["diagnosis"] == "BD").sum() == 23

meta.to_csv(
    "results/classification/metadata.csv"
)

print("Metadata prepared :)")
