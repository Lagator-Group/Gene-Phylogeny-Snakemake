from ete3 import Tree
import pandas as pd
import numpy as np
from skbio.stats.distance import mantel, DistanceMatrix
import random


### Mantel test
def cophenetic_matrix(tree, labels):

    n = len(labels)
    mat = np.zeros((n, n), dtype=float)

    for i, a in enumerate(labels):

        for j, b in enumerate(labels):

            if i == j:
                mat[i, j] = 0.0

            else:
                d = tree.get_distance(a, b)

                if np.isnan(d):
                    d = 0.0

                mat[i, j] = d

    # force perfect symmetry
    mat = (mat + mat.T) / 2

    # enforce exact zero diagonal
    np.fill_diagonal(mat, 0)

    return mat


def mantel_calc(tree1_path, tree2_path):
    t1 = Tree(tree1_path)
    t2 = Tree(tree2_path)

    shared = sorted(set(t1.get_leaf_names()) & set(t2.get_leaf_names()))

    t1.prune(shared)
    t2.prune(shared)

    t1.unroot()
    t2.unroot()

    m1 = cophenetic_matrix(t1, shared)
    m2 = cophenetic_matrix(t2, shared)

    dm1 = DistanceMatrix(m1, ids=shared)
    dm2 = DistanceMatrix(m2, ids=shared)

    r, p, n = mantel(dm1, dm2, method="pearson", permutations=999)

    print(f"Mantel r = {r:.4f}")
    print(f"p = {p:.6f}")
    print(f"n = {n}")

    df = pd.DataFrame(columns=["Gene 1", "Gene 2", "Mantel r", "p", "n"])

    df.loc[0] = [gene1, gene2, r, p, n]

    return df


### Robinson-Foulds
def get_shared_leaves(t1, t2):
    """Return shared leaf names."""
    leaves1 = set(t1.get_leaf_names())
    leaves2 = set(t2.get_leaf_names())
    return leaves1 & leaves2


def prune_to_shared(t1, t2):
    """
    Return copies of trees pruned to shared taxa.
    """
    shared = get_shared_leaves(t1, t2)

    if len(shared) < 4:
        raise ValueError(
            f"Only {len(shared)} shared taxa found. RF distance requires >=4 taxa."
        )

    t1p = Tree(t1.write())
    t2p = Tree(t2.write())

    t1p.prune(shared, preserve_branch_length=True)
    t2p.prune(shared, preserve_branch_length=True)

    return t1p, t2p, shared


def compare_trees(tree1, tree2, unrooted=True):
    """
    Compare two trees using RF distance.

    Returns
    -------
    dict containing:
        rf               : RF distance
        max_rf           : maximum possible RF
        norm_rf          : normalized RF
        shared_leaves    : set of shared taxa
        unique_to_tree1  : bipartitions unique to tree1
        unique_to_tree2  : bipartitions unique to tree2
    """

    t1, t2, shared = prune_to_shared(tree1, tree2)

    rf_results = t1.robinson_foulds(t2, unrooted_trees=unrooted)

    (
        rf,
        max_rf,
        common_leaves,
        parts_t1,
        parts_t2,
        discarded_edges_t1,
        discarded_edges_t2,
    ) = rf_results

    norm_rf = rf / max_rf if max_rf > 0 else 0

    unique_to_tree1 = parts_t1 - parts_t2
    unique_to_tree2 = parts_t2 - parts_t1

    return {
        "rf": rf,
        "max_rf": max_rf,
        "norm_rf": norm_rf,
        "shared_leaves": shared,
        "unique_to_tree1": unique_to_tree1,
        "unique_to_tree2": unique_to_tree2,
        "tree1_pruned": t1,
        "tree2_pruned": t2,
    }


def randomize_leaf_labels(tree):
    """
    Randomly permute tip labels while preserving topology.
    """
    t = Tree(tree.write())

    leaves = t.get_leaves()
    labels = [leaf.name for leaf in leaves]

    shuffled = labels[:]
    random.shuffle(shuffled)

    for leaf, label in zip(leaves, shuffled):
        leaf.name = label

    return t


def rf_permutation_test(tree1, tree2, n_permutations=1000, unrooted=True, seed=None):
    """
    Test whether observed RF is smaller than expected by chance.

    Randomizes labels of tree2 while preserving topology.
    """

    if seed is not None:
        random.seed(seed)

    observed = compare_trees(tree1, tree2, unrooted=unrooted)

    obs_rf = observed["norm_rf"]

    null_dist = []

    for _ in range(n_permutations):

        permuted_tree2 = randomize_leaf_labels(tree2)

        result = compare_trees(tree1, permuted_tree2, unrooted=unrooted)

        null_dist.append(result["norm_rf"])

    p_value = (sum(x <= obs_rf for x in null_dist) + 1) / (n_permutations + 1)

    observed["null_distribution"] = null_dist
    observed["p_value"] = p_value

    return observed


def compare_many_pairs(tree_pairs, n_permutations=None, unrooted=True):
    """
    Efficient wrapper for many tree comparisons.

    Parameters
    ----------
    tree_pairs : iterable of (name, treefile1, treefile2)
    """

    results = []

    for pair_name, f1, f2 in tree_pairs:

        t1 = Tree(f1)
        t2 = Tree(f2)

        res = compare_trees(t1, t2, unrooted=unrooted)

        if n_permutations:
            perm = rf_permutation_test(
                t1, t2, n_permutations=n_permutations, unrooted=unrooted
            )
            res["p_value"] = perm["p_value"]

        res["pair"] = pair_name
        results.append(res)

    return results


def robinson_fould(gene1_path, gene2_path):

    tree1 = Tree(gene1_path)
    tree2 = Tree(gene2_path)

    result = compare_trees(tree1, tree2)

    print("RF:", result["rf"])
    print("Max RF:", result["max_rf"])
    print("Normalized RF:", result["norm_rf"])
    print("Shared taxa:", len(result["shared_leaves"]))

    print("\nUnique branches in tree1:")
    for b in result["unique_to_tree1"]:
        print(b)

    print("\nUnique branches in tree2:")
    for b in result["unique_to_tree2"]:
        print(b)

    perm_result = rf_permutation_test(tree1, tree2, n_permutations=1000, seed=42)
    print("\nPermutation p-value:", perm_result["p_value"])

    cols = [
        "Gene 1",
        "Gene 2",
        "RF",
        "Max RF",
        "Normalized RF",
        "Shared taxa",
        "Permutation p-value",
    ]
    df = pd.DataFrame(columns=cols)

    df.loc[0] = [
        gene1,
        gene2,
        result["rf"],
        result["max_rf"],
        result["norm_rf"],
        len(result["shared_leaves"]),
        perm_result["p_value"],
    ]

    return df


### Execution

gene1_path = snakemake.input[0]
gene2_path = snakemake.input[1]
path_out = snakemake.output[0]

gene1 = snakemake.wildcards.gene1
gene2 = snakemake.wildcards.gene2

print(f"{gene1} + {gene2}")

#### Mantel Exec and Save
try:
    mantel_df = mantel_calc(gene1_path, gene2_path)
    print("done")
except Exception as e:
    print(e)
    mantel_df = pd.DataFrame(columns=["Gene 1", "Gene 2", "Mantel r", "p", "n"])
    mantel_df.loc[0] = [gene1, gene2, np.nan, np.nan, np.nan]


#### Robinson-Foulds Exec and Save
try:
    rf_df = robinson_fould(gene1_path, gene2_path)
    print("done")
except Exception as e:
    print(e)
    rf_df = pd.DataFrame(
        columns=[
            "Gene 1",
            "Gene 2",
            "RF",
            "Max RF",
            "Normalized RF",
            "Shared taxa",
            "Permutation p-value",
        ]
    )
    rf_df.loc[0] = [gene1, gene2, np.nan, np.nan, np.nan, np.nan, np.nan]

merged_df = pd.DataFrame(
    {
        "Gene1": gene1,
        "Gene2": gene2,
        "Mantel r": mantel_df["Mantel r"],
        "p": mantel_df["p"],
        "n": mantel_df["n"],
        "RF": rf_df["RF"],
        "Max RF": rf_df["Max RF"],
        "Normalized RF": rf_df["Normalized RF"],
        "Shared taxa": rf_df["Shared taxa"],
        "Permutation p-value": rf_df["Permutation p-value"],
    }
)

merged_df.to_csv(path_out, index=False)
