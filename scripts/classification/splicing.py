import pandas as pd
import numpy as np

def parse_leafcutter_ratios(path: str):
  raw = pd.read_csv(path, sep=r"\s+", index_col=0)

  # Stack into a single Series, split once into two columns, then unstack
  parts = raw.stack().str.split("/", expand=True).astype(int)

  numerator = parts[0].unstack()
  denominator = parts[1].unstack()
  psi = numerator.div(denominator.replace(0, np.nan))

  return psi.T, numerator.T, denominator.T


psi, numerator, coverage = parse_leafcutter_ratios(
    "results/leafcutter2/scz_vs_bip.cluster_ratios.gz"
)

print(psi.shape)
print(psi.head())
