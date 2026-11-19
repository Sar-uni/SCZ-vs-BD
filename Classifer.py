"""
Mostly followed tanmay's notes
"""

import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import roc_auc_score
from scipy import stats
from Bio.Seq import Seq


# 1.
def load_data(meta_path="metadata.csv", expr_path="expr_logcpm.csv",
spl_path="splicing_psi.csv"):
    '''TO DO
    Idk file format this will use exactly so its csv for now
    '''
    meta = pd.read_csv(meta_path, index_col=0)
    expr = pd.read_csv(expr_path, index_col=0)
    spl = pd.read_csv(spl_path, index_col=0)

    ids = meta.index.intersection(expr.index).intersection(spl.index)
    meta, expr, spl = meta.loc[ids], expr.loc[ids], spl.loc[ids]

    y = (meta["diagnosis"] == "SCZ").astype(int).values
    cov = meta[["age", "sex", "pmi", "ph"]].copy()
    cov["sex"] = (cov["sex"] == "M").astype(int)
    return y, expr, spl, cov

def check_groups(cov, y):
    ''' 
    TO DO
    Compare SCZ vs BD on age, sex, PMI, pH.
    Return a small DataFrame of group means and p-values.
    '''
    rows = []
    for col in ["age", "pmi", "ph"]:
        scz, bd = cov[col][y == 1], cov[col][y == 0]
        p = stats.mannwhitneyu(scz, bd).pvalue
        rows.append((col, scz.mean(), bd.mean(), p))

    table = pd.crosstab(cov["sex"], y)
    p_sex = stats.fisher_exact(table)[1]
    rows.append(("sex (fraction male)", cov["sex"][y == 1].mean(),
                 cov["sex"][y == 0].mean(), p_sex))
    return pd.DataFrame(rows, columns=["variable", "SCZ", "BD", "p"])

# 2.
# alignment done seperately

# 3.
def filter_features(table, minV=0.0, maxM=0.2):
    '''
    TO DO
    minV = minimum varinace
    maxM = max missing
    drop features with too much missing data or near-zero variance
    and return that table
    '''
    keep_missing = table.isna().mean() <= maxM
    keep_var = table.var() > minV
    return table.loc[:, keep_missing & keep_var]


def count_genes(bam_files, gtf):
    pass

# 4.

def model_factory(c=0.1, k=None):
    '''
    TO DO
    make pipeline for table
    '''
    steps = [("impute", SimpleImputer(strategy="median")),
             ("scale", StandardScaler())]
    if k is not None:
        steps.append(("select", SelectKBest(f_classif, k=k)))
    steps.append(("clf", LogisticRegression(C=c, max_iter=5000)))

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
    X = np.asarray(X, dtype=float)
    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=repeats, random_state=seed)
    total = np.zeros(len(y))
    count = np.zeros(len(y))
    for train, test in cv.split(X, y):
        kk = None if k is None else min(k, X.shape[1])
        model = model_factory(k=kk).fit(X[train], y[train])
        total[test] += model.predict_proba(X[test])[:, 1]
        count[test] += 1
    return total / count


# 6.

def top_features(X, y, k=100, n_top=10, repeats=5, seed=0):
    '''
    TO DO
    Fit make_model(k) on all data (or per CV fold?) and 
    count how often each feature is selected across fold
    '''
    names = np.array(X.columns)
    Xv = X.values.astype(float)
    k = min(k, Xv.shape[1])
    counts = pd.Series(0, index=names)
    cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=repeats, random_state=seed)
    for train, _ in cv.split(Xv, y):
        model = model_factory(k=k).fit(Xv[train], y[train])
        chosen = model.named_steps["select"].get_support()
        counts[names[chosen]] += 1
    return counts.sort_values(ascending=False).head(n_top)

def filter_features(table, minV=0.0, maxM=0.2):
    '''
    TO DO
    Drop bad features. Does NOT use diagnosis labels.

    table = DataFrame, samples as rows, features as columns
    minV  = minimum variance (drop features at or below this)
    maxM  = max fraction of missing values allowed (drop features above this)

    Steps:
      1. missing_fraction = fraction of NaN in each column  (table.isna().mean())
      2. variance of each column                            (table.var())
      3. keep columns with missing_fraction <= maxM AND variance > minV
      4. return table with only those columns

    Returns the filtered DataFrame
    '''
    pass

def permutation_test(X, y, k=None, n_perm=200, repeats=2, seed=3):
    '''
    TO DO
    Is the model's AUROC better than luck?

    Steps:
      1. real_auc = AUROC of CVpredictions on the real labels
      2. repeat n_perm times:
           - shuffle y (np.random.default_rng(seed).permutation(y))
           - run CVpredictions on the shuffled labels
           - save the AUROC (scored against the shuffled labels)
      3. p_value = (1 + number of shuffled AUROCs >= real_auc) / (1 + n_perm)

    Returns (real_auc, p_value)
    Note: use fewer repeats than the main run, this reruns the whole pipeline
    n_perm times and gets slow.
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

    result = {"exon_len": len(exon), "frameshift_if_skipped": len(exon) % 3 != 0}
    for label, dna in [("included", upstream + exon + downstream),
                       ("skipped", upstream + downstream)]:
        dna = dna[: len(dna) - len(dna) % 3]
        protein = str(Seq(dna).translate())
        result[label + "_protein"] = protein
        result[label + "_premature_stop"] = "*" in protein.rstrip("*")
    return result


def find_pfam_domains(protein_seq, pfam_path="?"):
    '''
    TO DO
    does this protein contain any known functional domains?
    run pyhmmer hmmsearch on the protein and return the Pfam domains hit

    Steps:
      1. remove the '*' stop symbol from protein_seq
      2. turn the sequence into a pyhmmer sequence (TextSequence, digitize with
         the amino alphabet)
      3. open Pfam-A.hmm with HMMFile
      4. run pyhmmer.hmmsearch(hmms, [seq], E=evalue)
      5. keep the domains whose hits are "included"

    Returns a list of Pfam domain names, e.g. ["Pkinase", "SH3_1"]
    Needs: pip install pyhmmer, and the Pfam-A.hmm file.

    Use it on both isoforms and compare the two lists:
      lost = set(domains_included) - set(domains_skipped)
    '''
    pass




def make_example_data(n=40, n_genes=300, n_introns=500, seed=0):
    '''
    Writes three small fake csv files so the pipeline can be tested
    Samples are rows. diagnosis is SCZ or BD. A few features are given a
    small real signal so the models have something to find
    '''
    rng = np.random.default_rng(seed)
    ids = [f"s{i:02d}" for i in range(n)]
    diag = np.array(["SCZ"] * (n // 2) + ["BD"] * (n - n // 2))
    y = (diag == "SCZ").astype(int)

    meta = pd.DataFrame({
        "diagnosis": diag,
        "age": rng.normal(55, 10, n).round(),
        "sex": rng.choice(["M", "F"], n),
        "pmi": rng.normal(20, 6, n).round(1),
        "ph": rng.normal(6.4, 0.3, n).round(2),
    }, index=ids)

    expr = pd.DataFrame(rng.normal(5, 1, (n, n_genes)),
                        index=ids, columns=[f"gene{i}" for i in range(n_genes)])
    expr.iloc[:, :10] += y[:, None] * 0.6

    spl = pd.DataFrame(rng.beta(2, 2, (n, n_introns)),
                       index=ids, columns=[f"chr1:{i}:{i + 100}:clu{i}" for i in range(n_introns)])
    spl.iloc[:, :10] += y[:, None] * 0.15
    spl = spl.clip(0, 1)
    spl = spl.mask(rng.random(spl.shape) < 0.05)   # some missing PSI values

    meta.to_csv("metadata.csv")
    expr.to_csv("expr_logcpm.csv")
    spl.to_csv("splicing_psi.csv")


def bootstrap_delta_auc(y, p_a, p_b, n_boot=2000, seed=1):
    '''
    delta AUROC = AUROC(a) - AUROC(b), resampling the samples.
    Returns (mean delta, ci_low, ci_high).
    CI above 0: a better. CI below 0: b better. Includes 0: no difference.
    '''
    rng = np.random.default_rng(seed)
    deltas = []
    while len(deltas) < n_boot:
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) < 2:
            continue
        deltas.append(roc_auc_score(y[idx], p_a[idx]) - roc_auc_score(y[idx], p_b[idx]))
    lo, hi = np.percentile(deltas, [2.5, 97.5])
    return np.mean(deltas), lo, hi

def main():
    make_example_data()   # remove this when using real data


    y, expr, spl, cov = load_data()
    print("Samples:", len(y), "| SCZ:", y.sum(), "| BD:", (1 - y).sum())

    print("\nGroup comparison:")
    print(check_groups(cov, y).round(3).to_string(index=False))

    expr = filter_features(expr)
    spl = filter_features(spl)
    print("\nFeatures kept:", expr.shape[1], "genes,", spl.shape[1], "introns")

    K = 50   # same number of features for both tables
    p_expr = CVpredictions(expr.values, y, k=K)
    p_spl = CVpredictions(spl.values, y, k=K)
    p_cov = CVpredictions(cov.values, y, k=None)

    print("\nAUROC")
    print(f"  expression: {roc_auc_score(y, p_expr):.3f}")
    print(f"  splicing:   {roc_auc_score(y, p_spl):.3f}")
    print(f"  covariates: {roc_auc_score(y, p_cov):.3f}  (baseline)")

    delta, lo, hi = bootstrap_delta_auc(y, p_spl, p_expr)
    print(f"\ndelta (splicing - expression) = {delta:.3f}, 95% CI [{lo:.3f}, {hi:.3f}]")
    if lo > 0:
        print("-> splicing better (also check it beats the covariate baseline)")
    elif hi < 0:
        print("-> expression better")
    else:
        print("-> no detectable difference")

    print("\nTop splicing events (times selected across folds):")
    print(top_features(spl, y, k=K))

if __name__ == "__main__":
    main()