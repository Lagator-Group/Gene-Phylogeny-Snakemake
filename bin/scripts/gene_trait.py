import pandas as pd
import numpy as np
from collections import Counter
import random
from ete3 import Tree
from itertools import combinations
from skbio import DistanceMatrix
from skbio.stats.distance import mantel


def build_metadata_dict(meta, id_col, trait_col):
    meta = meta.copy()

    meta = meta.dropna(subset=[id_col])

    meta[trait_col] = meta[trait_col].fillna("Unknown").astype(str).str.strip()

    return dict(zip(meta[id_col], meta[trait_col]))


def compute_tree_purity(tree, metadata_dict):
    """
    Clustering of traits on tree.
    """

    purities = []

    for node in tree.traverse():

        leaves = node.get_leaves()
        traits = [metadata_dict.get(l.name) for l in leaves]
        traits = [t for t in traits if t is not None]

        if len(traits) == 0:
            continue

        most_common = Counter(traits).most_common(1)[0][1]
        purity = most_common / len(traits)

        purities.append(purity)

    return np.mean(purities)


def permutation_test(tree, metadata_dict, n_iter=1000, seed=42):
    """
    Randomises metadata labels to test significance.
    """

    random.seed(seed)

    keys = list(metadata_dict.keys())
    values = list(metadata_dict.values())

    observed = compute_tree_purity(tree, metadata_dict)

    null_scores = []

    for _ in range(n_iter):

        shuffled_values = random.sample(values, len(values))
        shuffled_dict = dict(zip(keys, shuffled_values))

        score = compute_tree_purity(tree, shuffled_dict)
        null_scores.append(score)

    null_scores = np.array(null_scores)

    p_value = (np.sum(null_scores >= observed) + 1) / (n_iter + 1)

    return {
        "observed_purity": observed,
        "null_mean": null_scores.mean(),
        "null_std": null_scores.std(),
        "p_value": p_value,
        "null_distribution": null_scores,
    }


def tree_to_distance_matrix(tree):
    """
    Returns skbio DistanceMatrix from tree branch lengths.
    """

    leaves = tree.get_leaves()
    names = [l.name for l in leaves]

    n = len(leaves)
    dist = np.zeros((n, n))

    # compute pairwise distances
    for i, j in combinations(range(n), 2):

        d = tree.get_distance(leaves[i], leaves[j])

        dist[i, j] = d
        dist[j, i] = d

    return DistanceMatrix(dist, names)


def trait_to_distance_matrix(metadata_dict):

    names = list(metadata_dict.keys())
    n = len(names)

    dist = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):

            if metadata_dict[names[i]] == metadata_dict[names[j]]:
                d = 0
            else:
                d = 1

            dist[i, j] = d
            dist[j, i] = d

    return DistanceMatrix(dist, names)


def run_mantel(tree_dm, trait_dm, method="pearson", permutations=999):

    corr, p_value, n = mantel(
        tree_dm, trait_dm, method=method, permutations=permutations
    )

    return {"mantel_r": corr, "p_value": p_value, "permutations": permutations}


def run_full_phylo_association(tree_path, meta, id_col, trait_col, n_iter=1000):

    tree = Tree(tree_path)

    # metadata
    metadata_dict = build_metadata_dict(meta, id_col, trait_col)

    purity_result = permutation_test(
        tree=tree, metadata_dict=metadata_dict, n_iter=n_iter
    )

    # align tips (important!)
    tree_names = set(l.name for l in tree.iter_leaves())
    metadata_dict = {k: v for k, v in metadata_dict.items() if k in tree_names}

    # prune tree to match metadata (important)
    tree.prune(list(metadata_dict.keys()), preserve_branch_length=True)

    # distance matrices
    tree_dm = tree_to_distance_matrix(tree)
    trait_dm = trait_to_distance_matrix(metadata_dict)

    # Mantel test
    mantel_result = run_mantel(tree_dm, trait_dm)

    return  purity_result, mantel_result


t_path = snakemake.input[0]
metadata = pd.read_csv(snakemake.input[1])
id_col = snakemake.params.gene_col
trait_col = snakemake.wildcards.trait

path_out = snakemake.output[0]

try:
    purity_result, mantel_result = run_full_phylo_association(
        tree_path=t_path,
        meta=metadata,
        id_col=id_col,
        trait_col=trait_col,
    )

    print("Observed purity:", purity_result["observed_purity"])
    print("Null mean:", purity_result["null_mean"])
    print("Purity P-value:", purity_result["p_value"])
    print("Mantel r:", mantel_result["mantel_r"])
    print("Mantel P-value:", mantel_result["p_value"])

    df = pd.DataFrame(
        {
            "Gene": snakemake.wildcards.gene,
            "Trait": trait_col,
            "Observed Purity": [purity_result["observed_purity"]],
            "Null Mean": [purity_result["null_mean"]],
            "Purity p-value": [purity_result["p_value"]],
            "Mantel r": [mantel_result["mantel_r"]],
            "Mantel p-value": [mantel_result["p_value"]],
        }
    )

except Exception as e:
    print(e)
    df = pd.DataFrame(
        {
            "Gene": snakemake.wildcards.gene,
            "Trait": trait_col,
            "Observed Purity": [np.nan],
            "Null Mean": [np.nan],
            "Purity p-value": [np.nan],
            "Mantel r": [np.nan],
            "Mantel p-value": [np.nan],
        }
    )

df.to_csv(path_out, index=False)