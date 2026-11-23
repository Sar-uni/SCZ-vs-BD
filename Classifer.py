"""
Mostly followed tanmay's notes
"""

import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.feature_selection import VarianceThreshold
from sklearn.exceptions import ConvergenceWarning
from sklearn.model_selection import cross_validate
import warnings
from scipy import stats
# from Bio.Seq import Seq
from pathlib import Path


# 1.
def load_metadata(path="results/classification/metadata.csv"):
    """Load a row per sample. Returns samples x covarites"""

    meta = pd.read_csv(path, index_col=0)
    meta.index = meta.index.astype(str).str.strip()

    meta["diagnosis"] = (
        meta["diagnosis"].astype(str).str.strip()
        .replace({
            "Schizophrenia": "SCZ",
            "Bipolar Disorder": "BD"
        })
    )
    meta["sex"] = meta["sex"].astype(str).str.upper().str.strip()
    for col in ["age", "pmi", "ph"]:
        meta[col] = pd.to_numeric(meta[col], errors="coerce")

    return meta


def load_expression_matrix(path="results/classification/expression_matrix.csv"):
    """Load VSTnormalized expression matrix.
       Returns samples x genes.
    """
    expr = pd.read_csv(path, index_col=0)
    expr.index = expr.index.astype(str).str.strip()
    expr.index.name = "sample"

    # want numeric expression values
    expr = expr.apply(pd.to_numeric, errors="raise")

    return expr


def load_leafcutter_ratios(path="results/leafcutter2/scz_vs_bip.cluster_ratios.gz", min_coverage=10):
    """Parse Leafcutter2 num/denom values.
       Returns samples x intron events.
    """
    raw = pd.read_csv(path, sep=r"\s+", index_col=0, dtype=str)

    # Split all "num/den" strings
    parts = raw.stack().str.split("/", n=1, expand=True)

    numerator = pd.to_numeric(parts[0], errors="coerce").unstack()
    denominator = pd.to_numeric(parts[1], errors="coerce").unstack()

    # intron usage proportion
    ratios = numerator.div(denominator.where(denominator > 0))

    # Mask low coverage stuff
    ratios = ratios.mask(denominator < min_coverage)
    ratios = ratios.replace([np.inf, -np.inf], np.nan)

    # Transpose: samples x intron events
    return ratios.T, denominator.T


def prepare_matrices():
    """Create metadata, expression, and splicing matrices with matching samples
        Returns:
            meta: samples x covariates
            expr: samples x genes
            spl: samples x intron events
    """
    meta = load_metadata()
    expr = load_expression_matrix()
    # spl, coverage = load_leafcutter_ratios()
    
    # match same samples
    common = (
        meta.index
        .intersection(expr.index)
        # .intersection(spl.index)
    )

    meta = meta.loc[common]
    expr = expr.loc[common]
    # spl = spl.loc[common]

    # Ensure sample ordering
    assert meta.index.equals(expr.index)
    # assert meta.index.equals(spl.index)

    out = Path("results/classification")
    out.mkdir(parents=True, exist_ok=True)

    meta.to_csv(out / "metadata_matched.csv")
    expr.to_csv(out / "expression_matrix_matched.csv")
    # spl.to_csv(out / "splicing_matrix.csv")
    # coverage.loc[common].to_csv(
    #     out / "splicing_coverage.csv"
    # )
    ############# TESTING ONLY
    spl = pd.read_csv("results/classification/splicing_matrix.csv", index_col=0)

    print(meta["diagnosis"].value_counts())
    print("exp features:", expr.shape[1])
    print("spl features:", spl.shape[1])

    return meta, expr, spl

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
    keep_missing = table.isna().mean() <= maxM
    keep_var = table.var() > minV
    return table.loc[:, keep_missing & keep_var]


def count_genes(bam_files, gtf):
    pass

# 4.
def model_factory(c=0.1, k=None):
    '''
    Create a pipeline for table. Independently fitted in each CV fold.
    '''
    steps = [
        ("impute", SimpleImputer(strategy="median")), # missing ratios
        ("variance",VarianceThreshold(threshold=1e-8)), #remove constant features
        ("scale", StandardScaler()) #standardize features
        ]
    if k is not None:
        steps.append(("select", SelectKBest(score_func=f_classif, k=k))) # best features
    steps.append(("clf", LogisticRegression( 
        solver="saga",
        l1_ratio=1.0,   # L1 logistic regression (0 for l2)
        C=c,
        class_weight="balanced",
        max_iter=5000,
        random_state=42
    )))

    return Pipeline(steps)


def make_baseline_model(C=0.1):
    """Covariates-only model (age, sex, PMI, pH): no feature selection."""
    return Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
        ("clf", LogisticRegression(
            l1_ratio=1.0,
            solver="saga",
            C=C, max_iter=5000, 
            class_weight="balanced", 
            random_state=42)),
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

def repeated_cv_auc(feature_sets,y, k=50, c=0.1, n_splits=5, n_repeats=20, seed=42):
    """
    In each repeat, every sample is held out exactly once. Combines those out-of-fold 
    predictions and calculate 1 AUROC for that repeat. All models use the same 
    train/test splits. This is important because we want differences in performance to 
    reflect the models, instead differences in how the data was divided.
    """
    y = np.asarray(y)
    cv = RepeatedStratifiedKFold(n_splits=n_splits,n_repeats=n_repeats,random_state=seed)
    splits = list(cv.split(np.zeros(len(y)), y))
    results = []

    for repeat in range(n_repeats):
        # One prediction per sample for this repeat
        predictions = {
            name: np.full(len(y), np.nan)
            for name in feature_sets
        }
        # Five folds for this repeat
        repeat_splits = splits[repeat * n_splits:(repeat + 1) * n_splits]

        for train_idx, test_idx in repeat_splits:
            for name, X in feature_sets.items():
                # Covar not need feature selection
                model_k = None if name == "covariates" else k
                model = model_factory(c=c, k=model_k)
                model.fit(X.iloc[train_idx],y[train_idx])
                predictions[name][test_idx] = (
                    model.predict_proba(X.iloc[test_idx])[:, 1]
                )

        # Calculate AUROC
        for name, pred in predictions.items():
            auc = roc_auc_score(y, pred)
            results.append({
                "repeat": repeat + 1,
                "model": name,
                "auc": auc
            })
    return pd.DataFrame(results)

def paired_model_comparison(results):
    """
    Compare model AUROCs using the same CV repeats.
    Delta AUROC = Model A - Model B.
    """

    pivot = results.pivot(
        index="repeat",
        columns="model",
        values="auc"
    )

    comparisons = [
        ("splicing", "expression"),
        ("combined", "expression"),
        ("combined", "splicing"),
        ("expression", "covariates"),
        ("splicing", "covariates")
    ]

    rows = []

    for model_a, model_b in comparisons:
        if model_a not in pivot.columns or model_b not in pivot.columns:
            raise ValueError(
                f"Missing model: {model_a} or {model_b}"
            )

        paired = pivot[[model_a, model_b]].dropna()
        delta = paired[model_a] - paired[model_b]
        rows.append({
            "comparison": f"{model_a} - {model_b}",
            "mean_delta_auc": delta.mean(),
            "median_delta_auc": delta.median(),
            "sd_delta_auc": delta.std(),
            "min_delta_auc": delta.min(),
            "max_delta_auc": delta.max(),
            "model_a_wins": (delta > 0).sum(),
            "model_b_wins": (delta < 0).sum(),
            "ties": (delta == 0).sum(),
            "repeats": len(delta)
        })
    return pd.DataFrame(rows)

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

def permutation_test(feature_sets,y,observed_results,k=50,c=0.1,n_splits=5,n_repeats=20,n_perm=100,seed=3):
    '''
    Is the model's AUROC better than luck?

    Steps:
      1. Mean AUROC for each model across all repeats (observed_results)
      2. repeat n_perm times:
           - shuffle y (np.random.default_rng(seed).permutation(y))
           - run CVpredictions on the shuffled labels
           - save the AUROC (scored against the shuffled labels)
      3. p_value = (1 + number of shuffled AUROCs >= real_auc) / (1 + n_perm)
    '''
    y = np.asarray(y)
    rng = np.random.default_rng(seed)
    observed = (
        observed_results.groupby("model")["auc"]
        .mean()
        .to_dict()
    )

    null_scores = {name: [] for name in feature_sets}

    for i in range(n_perm):
        print(f"Permutation {i + 1}/{n_perm}")
        # Shuffle diagnosis labels
        y_perm = rng.permutation(y)
        # Re-run all models using shuffled labels
        perm_results = repeated_cv_auc(
            feature_sets=feature_sets,
            y=y_perm,
            k=k,
            c=c,
            n_splits=n_splits,
            n_repeats=n_repeats,
            seed=seed
        )

        perm_means = perm_results.groupby("model")["auc"].mean()
        for name in feature_sets:
            null_scores[name].append(perm_means[name])
    rows = []

    # Compare observed AUROC to null distribution for each model
    for name in feature_sets:
        null = np.asarray(null_scores[name])
        # One sided permutation p-value
        p_value = (1 + np.sum(null >= observed[name])) / (n_perm + 1)

        rows.append({
            "model": name,
            "observed_auc": observed[name],
            "null_mean_auc": null.mean(),
            "null_sd_auc": null.std(),
            "permutation_p": p_value
        })

    return pd.DataFrame(rows), pd.DataFrame(null_scores)


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
    meta, expr, spl = prepare_matrices()

    y = (meta["diagnosis"] == "SCZ").astype(int).values
    cov = meta[["age", "sex", "pmi", "ph"]].copy()
    cov["sex"] = (cov["sex"] == "MALE").astype(int)
    print("Samples:", len(y), "| SCZ:", y.sum(), "| BD:", (1 - y).sum())

    # Combine expression and splicing features
    combined = pd.concat([expr.add_prefix("expr__"),spl.add_prefix("spl__")],axis=1)

    feature_sets = {
        "expression": expr,
        "splicing": spl,
        "combined": combined,
        "covariates": cov
    }

    print("\nFeature counts:")
    for name, X in feature_sets.items():
        print(f"{name}: {X.shape[1]}")

    results = repeated_cv_auc(
        feature_sets=feature_sets,
        y=y,
        k=50,
        c=0.1,
        n_splits=5,
        n_repeats=20,
        seed=42
    )

    print("\nAUROC results:")
    summary = (
        results.groupby("model")["auc"]
        .agg(["mean", "std", "min", "max"])
        .sort_values("mean", ascending=False)
    )
    print(summary.round(3))
    results.to_csv("results/classification/classification_cv_results.csv",index=False)
    summary.to_csv("results/classification/classification_cv_summary.csv")

    print("\nPermutation test--")
    perm_summary, perm_null = permutation_test(
        feature_sets=feature_sets,
        y=y,
        observed_results=results,
        k=50,
        c=0.1,
        n_splits=5,
        n_repeats=5,
        n_perm=2,
        seed=123
    )

    print("\nPermutation test results:")
    print(perm_summary.round(3))
    perm_summary.to_csv("results/classification/permutation_test_summary.csv",index=False)
    perm_null.to_csv("results/classification/permutation_null_distributions.csv",index=False)

    print("\nPaired model comparison:")
    paired_results = paired_model_comparison(results)
    print(paired_results.round(3).to_string(index=False))

    paired_results.to_csv("results/classification/paired_model_comparisons.csv",index=False)

    print("\nGroup comparison:")
    print(check_groups(cov, y).round(3).to_string(index=False))

    # expr = filter_features(expr) # removed to prevent test data from influencing feature selection
    # print("\nFeatures kept:", expr.shape[1], "genes,", spl.shape[1], "introns")

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