"""
Mostly followed tanmay's notes
"""

import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression


# 1.
def load_data(meta_path="metadata.csv", expr_path="expr_logcpm.csv",
spl_path="splicing_psi.csv"):
    '''TO DO
    Idk file format this will use exactly so its csv for now
    '''
    pass

def check_groups(cov, y):
    ''' 
    TO DO
    Compare SCZ vs BD on age, sex, PMI, pH.
    Return a small DataFrame of group means and p-values.
    '''
    pass

# 2.
# idk how to even start
def align_all(sample_table):
    '''
    TO DO
    '''
    pass

# 3.
def filter_features(table, minV=0.0, maxM=0.2):
    '''
    TO DO
    minV = minimum varinace
    maxM = max missing
    drop features with too much missing data or near-zero variance
    and return that table
    '''
    pass


def count_genes(bam_files, gtf):
    pass

# 4.

def model_factory(c, k=None):
    '''
    TO DO
    make pipeline for table
    '''
    steps = []

    return Pipeline(steps)


def make_baseline_model(C=0.1):
    """Covariates-only model (age, sex, PMI, pH): no feature selection."""
    return Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
        ("clf", LogisticRegression(C=C, max_iter=5000)),
    ])


# 5.
def CVpredictions(X, y, k=None, repeats=10, seed=0):
    '''
    TO DO
    Returns one averaged out-of-fold
    probability per sample.
    '''
    pass


# 6.

def top_features(X, y, k=100, n_top=10):
    '''
    TO DO
    Fit make_model(k) on all data (or per CV fold?) and 
    count how often each feature is selected across fold
    '''
    pass


def translate_isoforms(upstream, exon, downstream):
    '''
    TO DO
    Build the CDS with and without the alternative exon

    translate both, and report:
      - frameshift if len(exon) % 3 != 0
      - premature stop codon ('*' before the end of the protein)
    '''
    pass

def find_pfam_domains(protein_seq, pfam_path="?"):
    '''
    TO DO
    run pyhmmer hmmsearch on the protein and return the Pfam domains hit
    '''
    pass

def main():
    print("hello")

if __name__ == "__main__":
    main()