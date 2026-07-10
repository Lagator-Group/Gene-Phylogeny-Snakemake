"""
Continuous / numerical analogue of gene_trait.py.

Where gene_trait.py treats a metadata column as *categorical* (binarising each
unique value and running purity + a binary Mantel test), this script treats the
column as a single *continuous* variable and measures how strongly the numeric
value is structured by the gene tree ("phylogenetic signal").

For each (gene, trait) it reports four complementary statistics:

  * Blomberg's K   - variance of phylogenetically-independent contrasts vs. the
                     Brownian-motion expectation. K ~ 1 => Brownian, K > 1 =>
                     more clustered than Brownian, K ~ 0 => no signal.
                     p-value from tip-label permutation.
  * Pagel's lambda - ML scaling of the off-diagonal covariance. lambda ~ 1 =>
                     strong signal, lambda ~ 0 => none. p-value from a likelihood
                     ratio test against lambda = 0.
  * Moran's I      - phylogenetic autocorrelation using 1/cophenetic-distance
                     weights. p-value from tip-label permutation.
  * Mantel r       - Pearson correlation between the cophenetic distance matrix
                     and the pairwise |trait_i - trait_j| distance matrix. This
                     is the direct continuous analogue of the binary Mantel test
                     in gene_trait.py. p-value from skbio's Mantel permutation.

Wired into Snakemake exactly like gene_trait.py:
    snakemake.input[0]      -> gene .treefile
    snakemake.input[1]      -> metadata csv
    snakemake.params.gene_col -> accession column (must match treefile leaves)
    snakemake.wildcards.gene  -> gene id
    snakemake.wildcards.trait -> numeric metadata column
    snakemake.output[0]     -> per (gene, trait) csv (single row)
"""

import warnings

import numpy as np
import pandas as pd

from itertools import combinations

from ete3 import Tree
from scipy.optimize import minimize_scalar
from scipy.stats import chi2
from skbio import DistanceMatrix
from skbio.stats.distance import mantel

# Minimum number of tree leaves with a numeric value required to run a test.
MIN_ACCESSIONS = 10
# Number of permutations for the permutation-based p-values (K, Moran's I) and
# for the skbio Mantel test.
N_PERMUTATIONS = 999
SEED = 42


def clean_id(x):
    """Remove duplicate suffixes like _dupelabel1, _dupelabel2, etc.

    Coerces to str first: pandas 3.0's `astype(str)` uses the new ``str`` dtype,
    which leaves NaN as a float rather than the string 'nan', so a bare
    ``x.split`` would raise on missing ids.
    """
    return str(x).split("_dupelabel")[0]


def preprocess_tree(tree):
    """Attach the metadata-space (de-duplicated) id to every leaf."""
    for leaf in tree.iter_leaves():
        leaf.add_feature("clean_id", clean_id(leaf.name))
    return tree


def get_numeric_trait(meta, id_col, trait_col):
    """
    Return a Series of {clean_id: numeric value} for a continuous trait.

    Non-numeric cells are coerced to NaN and dropped. When an accession appears
    on multiple metadata rows the values are averaged, so each accession
    contributes exactly one value (avoiding phylogenetic pseudo-replication).
    """
    sub = meta[[id_col, trait_col]].copy()
    sub = sub.dropna(subset=[id_col])
    sub[id_col] = sub[id_col].astype(str).map(clean_id)
    sub[trait_col] = pd.to_numeric(sub[trait_col], errors="coerce")
    sub = sub.dropna(subset=[trait_col])
    return sub.groupby(id_col)[trait_col].mean()


def build_value_map(tree, trait_series):
    """
    Map one representative original leaf name -> numeric value.

    Collapses _dupelabel copies to a single representative leaf per accession so
    that identical duplicated tips do not artificially inflate the signal.
    Returns {original_leaf_name: value}.
    """
    value_map = {}
    seen = set()
    for leaf in tree.iter_leaves():
        cid = clean_id(leaf.name)
        if cid in seen or cid not in trait_series.index:
            continue
        seen.add(cid)
        value_map[leaf.name] = float(trait_series.loc[cid])
    return value_map


def midpoint_root(tree):
    """Midpoint-root in place so the variance-covariance matrix is meaningful."""
    try:
        outgroup = tree.get_midpoint_outgroup()
        if outgroup is not None and outgroup is not tree:
            tree.set_outgroup(outgroup)
    except Exception as e:  # pragma: no cover - defensive
        warnings.warn(f"[gene_trait_numeric] midpoint rooting failed: {e}")
    return tree


def vcv_matrix(tree, leaf_names):
    """
    Brownian-motion variance-covariance matrix C.

    C[i, j] is the shared root-to-MRCA path length of leaves i and j, built by
    adding each edge's branch length to every pair of its descendant leaves.
    The diagonal is the root-to-tip distance.
    """
    idx = {name: i for i, name in enumerate(leaf_names)}
    n = len(leaf_names)
    C = np.zeros((n, n))
    for node in tree.traverse():
        if node.is_root():
            continue
        b = node.dist
        if b == 0:
            continue
        desc = [idx[l.name] for l in node.iter_leaves() if l.name in idx]
        if desc:
            C[np.ix_(desc, desc)] += b
    return C


def cophenetic_matrix(tree, leaf_names):
    """Pairwise cophenetic (patristic) distance matrix, symmetric, zero diag."""
    leaves = {l.name: l for l in tree.iter_leaves()}
    n = len(leaf_names)
    D = np.zeros((n, n))
    for i, j in combinations(range(n), 2):
        d = tree.get_distance(leaves[leaf_names[i]], leaves[leaf_names[j]])
        if np.isnan(d):
            d = 0.0
        D[i, j] = D[j, i] = d
    return D


def _phylo_mean(C_inv, x):
    ones = np.ones(len(x))
    return (ones @ C_inv @ x) / (ones @ C_inv @ ones)


def _k_ratio(C_inv, x):
    """Observed MSE0 / MSE ratio used inside Blomberg's K."""
    a = _phylo_mean(C_inv, x)
    u = x - a
    mse0 = u @ u
    mse = u @ C_inv @ u
    return mse0 / mse if mse != 0 else np.nan


def blombergs_k(C, x, n_perm=N_PERMUTATIONS, rng=None):
    """Blomberg's K plus a tip-permutation p-value (one-sided, K high => signal)."""
    n = len(x)
    C_inv = np.linalg.pinv(C)
    ones = np.ones(n)
    expected = (np.trace(C) - n / (ones @ C_inv @ ones)) / (n - 1)
    obs_ratio = _k_ratio(C_inv, x)
    if not np.isfinite(obs_ratio) or expected == 0:
        return np.nan, np.nan
    k = obs_ratio / expected

    if rng is None:
        rng = np.random.default_rng(SEED)
    ge = 1  # +1 for the observed value
    for _ in range(n_perm):
        perm_ratio = _k_ratio(C_inv, rng.permutation(x))
        if np.isfinite(perm_ratio) and perm_ratio >= obs_ratio:
            ge += 1
    p = ge / (n_perm + 1)
    return k, p


def _lambda_transform(C, lam):
    Cl = C * lam
    np.fill_diagonal(Cl, np.diag(C))
    return Cl


def _neg_log_lik(lam, C, x):
    n = len(x)
    Cl = _lambda_transform(C, lam)
    sign, logdet = np.linalg.slogdet(Cl)
    if sign <= 0:
        return np.inf
    Cl_inv = np.linalg.pinv(Cl)
    a = _phylo_mean(Cl_inv, x)
    u = x - a
    sigma2 = (u @ Cl_inv @ u) / n
    if sigma2 <= 0:
        return np.inf
    ll = -0.5 * (n * np.log(2 * np.pi * sigma2) + logdet + n)
    return -ll


def pagels_lambda(C, x):
    """
    ML estimate of Pagel's lambda plus a likelihood-ratio p-value against
    lambda = 0 (no phylogenetic signal). Returns (lambda_hat, p_value).
    """
    res = minimize_scalar(
        _neg_log_lik, bounds=(0.0, 1.0), args=(C, x), method="bounded"
    )
    if not res.success:
        return np.nan, np.nan
    lam_hat = float(res.x)
    ll_hat = -res.fun
    ll_null = -_neg_log_lik(0.0, C, x)
    stat = 2 * (ll_hat - ll_null)
    if not np.isfinite(stat) or stat < 0:
        stat = 0.0
    p = chi2.sf(stat, df=1)
    return lam_hat, p


def _morans_i(W, z, s0):
    return (len(z) / s0) * (z @ W @ z) / (z @ z)


def morans_i(D, x, n_perm=N_PERMUTATIONS, rng=None):
    """
    Distance-based Moran's I (weights = 1/cophenetic distance) plus a
    tip-permutation p-value (one-sided, positive autocorrelation => signal).
    """
    n = len(x)
    with np.errstate(divide="ignore"):
        W = np.where(D > 0, 1.0 / D, 0.0)
    np.fill_diagonal(W, 0.0)
    s0 = W.sum()
    if s0 == 0:
        return np.nan, np.nan
    z = x - x.mean()
    if z @ z == 0:
        return np.nan, np.nan
    obs = _morans_i(W, z, s0)

    if rng is None:
        rng = np.random.default_rng(SEED)
    ge = 1
    for _ in range(n_perm):
        zp = rng.permutation(x)
        zp = zp - zp.mean()
        if _morans_i(W, zp, s0) >= obs:
            ge += 1
    p = ge / (n_perm + 1)
    return obs, p


def continuous_mantel(D, x):
    """Mantel r/p between cophenetic distance and |trait_i - trait_j|."""
    names = [str(i) for i in range(len(x))]
    tree_dm = DistanceMatrix(D, names)
    trait_dm = DistanceMatrix(np.abs(x[:, None] - x[None, :]), names)
    r, p, _ = mantel(
        tree_dm, trait_dm, method="pearson", permutations=N_PERMUTATIONS
    )
    return r, p


def run_numeric_test(tree, value_map):
    """Prune the tree to annotated leaves and run all continuous-signal tests."""
    if len(value_map) < MIN_ACCESSIONS:
        raise ValueError(
            f"only {len(value_map)} accessions with numeric values "
            f"(need >= {MIN_ACCESSIONS})"
        )

    tree.prune(list(value_map.keys()), preserve_branch_length=True)
    tree = midpoint_root(tree)

    leaf_names = [l.name for l in tree.iter_leaves()]
    x = np.array([value_map[name] for name in leaf_names], dtype=float)

    if np.var(x) == 0:
        raise ValueError("trait has zero variance across annotated leaves")

    C = vcv_matrix(tree, leaf_names)
    D = cophenetic_matrix(tree, leaf_names)
    rng = np.random.default_rng(SEED)

    k, k_p = blombergs_k(C, x, rng=rng)
    lam, lam_p = pagels_lambda(C, x)
    mi, mi_p = morans_i(D, x, rng=rng)
    mantel_r, mantel_p = continuous_mantel(D, x)

    return {
        "n_accessions": len(x),
        "trait_mean": float(np.mean(x)),
        "trait_sd": float(np.std(x, ddof=1)),
        "blomberg_k": k,
        "blomberg_k_p": k_p,
        "pagel_lambda": lam,
        "pagel_lambda_p": lam_p,
        "morans_i": mi,
        "morans_i_p": mi_p,
        "mantel_r": mantel_r,
        "mantel_p": mantel_p,
    }


def _row(gene, trait, status, error="", result=None):
    base = {
        "Gene": gene,
        "Trait": trait,
        "N_Accessions": np.nan,
        "Trait_Mean": np.nan,
        "Trait_SD": np.nan,
        "Status": status,
        "Blomberg K": np.nan,
        "Blomberg K p-value": np.nan,
        "Pagel lambda": np.nan,
        "Pagel lambda p-value": np.nan,
        "Morans I": np.nan,
        "Morans I p-value": np.nan,
        "Mantel r": np.nan,
        "Mantel p-value": np.nan,
        "Error": error,
    }
    if result is not None:
        base.update(
            {
                "N_Accessions": result["n_accessions"],
                "Trait_Mean": result["trait_mean"],
                "Trait_SD": result["trait_sd"],
                "Blomberg K": result["blomberg_k"],
                "Blomberg K p-value": result["blomberg_k_p"],
                "Pagel lambda": result["pagel_lambda"],
                "Pagel lambda p-value": result["pagel_lambda_p"],
                "Morans I": result["morans_i"],
                "Morans I p-value": result["morans_i_p"],
                "Mantel r": result["mantel_r"],
                "Mantel p-value": result["mantel_p"],
            }
        )
    return base


def main():
    t_path = snakemake.input[0]
    metadata = pd.read_csv(snakemake.input[1], low_memory=False)
    id_col = snakemake.params.gene_col
    trait_col = snakemake.wildcards.trait
    gene = snakemake.wildcards.gene
    path_out = snakemake.output[0]

    trait_series = get_numeric_trait(metadata, id_col, trait_col)

    print(
        f"[gene_trait_numeric] gene={gene} trait={trait_col} "
        f"numeric_accessions={len(trait_series)}"
    )

    if len(trait_series) == 0:
        row = _row(
            gene, trait_col, "insufficient_data",
            error="no numeric values in trait column",
        )
        pd.DataFrame([row]).to_csv(path_out, index=False)
        print(f"[gene_trait_numeric] wrote 1 row to {path_out}")
        return

    base_tree = preprocess_tree(Tree(t_path))
    value_map = build_value_map(base_tree, trait_series)

    if len(value_map) < MIN_ACCESSIONS:
        row = _row(
            gene, trait_col, "insufficient_data",
            error=(
                f"only {len(value_map)} tree leaves overlap numeric metadata "
                f"(need >= {MIN_ACCESSIONS})"
            ),
        )
        pd.DataFrame([row]).to_csv(path_out, index=False)
        print(f"[gene_trait_numeric] wrote 1 row to {path_out}")
        return

    try:
        test_tree = base_tree.copy("newick-extended")
        result = run_numeric_test(test_tree, value_map)
        row = _row(gene, trait_col, "ok", result=result)
        print(
            f"  n={result['n_accessions']} "
            f"K={result['blomberg_k']:.3g} (p={result['blomberg_k_p']:.3g}) "
            f"lambda={result['pagel_lambda']:.3g} (p={result['pagel_lambda_p']:.3g}) "
            f"I={result['morans_i']:.3g} (p={result['morans_i_p']:.3g}) "
            f"mantel_r={result['mantel_r']:.3g} (p={result['mantel_p']:.3g})"
        )
    except Exception as e:
        warnings.warn(f"[gene_trait_numeric] {gene}/{trait_col} failed: {e}")
        row = _row(gene, trait_col, "error", error=str(e))

    pd.DataFrame([row]).to_csv(path_out, index=False)
    print(f"[gene_trait_numeric] wrote 1 row to {path_out}")


if __name__ == "__main__":
    main()
