import numpy as np
import pandas as pd
import random
import warnings

from collections import Counter
from itertools import combinations

from ete3 import Tree
from skbio import DistanceMatrix
from skbio.stats.distance import mantel


# Minimum number of accessions annotated with a trait value required to run
# the Mantel + purity test. Values below this threshold get an
# `insufficient_data` row instead of a meaningless p-value.
MIN_TRAIT_COUNT = 10

# Both classes must have at least this many accessions for Mantel to be
# meaningful. Otherwise the binary trait is too imbalanced.
MIN_CLASS_COUNT = 3


def clean_id(x):
    """
    Removes duplicate suffixes like _dupelabel1, _dupelabel2, etc.
    """
    return x.split("_dupelabel")[0]


def get_unique_trait_values(meta, id_col, trait_col):
    """
    Returns the set of unique non-null values in the trait column, after
    deduplicating (id_col, trait_col) pairs. Assumes the metadata is already
    fully exploded (one value per cell in the trait column).
    """
    sub = meta[[id_col, trait_col]].dropna(subset=[trait_col]).copy()
    sub[id_col] = sub[id_col].astype(str).apply(clean_id)
    sub[trait_col] = sub[trait_col].astype(str).str.strip()
    sub = sub.drop_duplicates(subset=[id_col, trait_col])
    values = sorted(
        v for v in sub[trait_col].unique()
        if v and v.lower() != "unknown"
    )
    return values, sub


def preprocess_tree(tree):
    # unify labels with metadata space
    for leaf in tree.iter_leaves():
        leaf.add_feature("clean_id", clean_id(leaf.name))
    return tree


def compute_tree_purity(tree, metadata_dict):
    purities = []
    for node in tree.traverse():
        leaves = node.get_leaves()
        traits = [metadata_dict.get(l.name) for l in leaves]
        traits = [t for t in traits if t is not None]
        if len(traits) == 0:
            continue
        most_common = Counter(traits).most_common(1)[0][1]
        purities.append(most_common / len(traits))
    return np.mean(purities)


def permutation_test(tree, metadata_dict, n_iter=1000, seed=42):
    random.seed(seed)
    keys = list(metadata_dict.keys())
    values = list(metadata_dict.values())
    observed = compute_tree_purity(tree, metadata_dict)
    null_scores = []
    for _ in range(n_iter):
        shuffled = random.sample(values, len(values))
        shuffled_dict = dict(zip(keys, shuffled))
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
    leaves = tree.get_leaves()
    names = [l.name for l in leaves]
    n = len(leaves)
    dist = np.zeros((n, n))
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
            dist[i, j] = 0 if metadata_dict[names[i]] == metadata_dict[names[j]] else 1
            dist[j, i] = dist[i, j]
    return DistanceMatrix(dist, names)


def run_mantel(tree_dm, trait_dm, method="pearson", permutations=999):
    r, p_value, _ = mantel(tree_dm, trait_dm, method=method, permutations=permutations)
    return {"mantel_r": r, "p_value": p_value, "permutations": permutations}


def build_test_dict(tree, sub_meta, id_col, trait_col, value):
    """
    Build a binary {original_leaf_name: "1" or "0"} dict over the leaves of
    the tree that are also present in the metadata. Keys are ORIGINAL leaf
    names (with any _dupelabel suffix intact) so that ete3's tree.prune()
    can find them. Values come from the cleaned-id metadata join.

    Returns (dict, n_positive, n_negative).
    """
    # Map: cleaned leaf id -> original leaf name (for the leaves that have
    # a unique cleaned form; if multiple leaves share the same cleaned id,
    # we keep them all — they're truly duplicate entries in the tree).
    leaf_clean_to_originals = {}
    for leaf in tree.iter_leaves():
        cid = clean_id(leaf.name)
        leaf_clean_to_originals.setdefault(cid, []).append(leaf.name)

    meta_ids = set(sub_meta[id_col].astype(str))
    annotated = set(
        sub_meta.loc[sub_meta[trait_col] == value, id_col].astype(str)
    )

    test_dict = {}
    for cid, originals in leaf_clean_to_originals.items():
        if cid not in meta_ids:
            # leaf has no metadata at all -> drop from test (tree-only scope)
            continue
        label = "1" if cid in annotated else "0"
        for orig in originals:
            test_dict[orig] = label

    n_pos = sum(1 for v in test_dict.values() if v == "1")
    n_neg = len(test_dict) - n_pos
    return test_dict, n_pos, n_neg


def run_value_test(tree, metadata_dict):
    """
    Run purity + permutation + Mantel on a single binary trait dict.
    Caller is responsible for cloning the tree (each value test mutates it
    via .prune()).
    """
    if len(metadata_dict) < 3:
        raise ValueError(
            f"only {len(metadata_dict)} accessions in test (need >= 3)"
        )
    positives = sum(1 for v in metadata_dict.values() if v == "1")
    negatives = len(metadata_dict) - positives
    if positives < MIN_CLASS_COUNT or negatives < MIN_CLASS_COUNT:
        raise ValueError(
            f"class imbalance too severe: positives={positives}, "
            f"negatives={negatives} (need >= {MIN_CLASS_COUNT} each)"
        )

    tree.prune(list(metadata_dict.keys()), preserve_branch_length=True)

    purity_result = permutation_test(tree=tree, metadata_dict=metadata_dict)
    tree_dm = tree_to_distance_matrix(tree)
    trait_dm = trait_to_distance_matrix(metadata_dict)
    mantel_result = run_mantel(tree_dm, trait_dm)

    return {
        "n_accessions": len(metadata_dict),
        "n_positive": positives,
        "n_negative": negatives,
        "observed_purity": purity_result["observed_purity"],
        "null_mean": purity_result["null_mean"],
        "purity_p_value": purity_result["p_value"],
        "mantel_r": mantel_result["mantel_r"],
        "mantel_p_value": mantel_result["p_value"],
    }


# -----------------------
# Snakemake entry point
# -----------------------
def main():
    t_path = snakemake.input[0]
    metadata = pd.read_csv(snakemake.input[1], low_memory=False)
    id_col = snakemake.params.gene_col
    trait_col = snakemake.wildcards.trait

    path_out = snakemake.output[0]

    # discover the universe of values for this trait column
    trait_values, sub_meta = get_unique_trait_values(metadata, id_col, trait_col)

    # Load + preprocess the tree ONCE per gene. We clone it per value test
    # because .prune() mutates the tree in place.
    base_tree = Tree(t_path)
    base_tree = preprocess_tree(base_tree)

    print(f"[gene_trait] gene={snakemake.wildcards.gene} trait={trait_col} "
          f"values={len(trait_values)}")

    rows = []
    for value in trait_values:
        # Count unique accessions in the FULL metadata annotated with this
        # value. Used only to decide whether to short-circuit with
        # insufficient_data. The actual test candidate space is built per
        # value inside build_test_dict (it depends on which leaves of the tree
        # have metadata at all).
        n_unique_annotated = int(
            sub_meta.loc[sub_meta[trait_col] == value, id_col].nunique()
        )

        if n_unique_annotated < MIN_TRAIT_COUNT:
            rows.append({
                "Gene": snakemake.wildcards.gene,
                "Trait": trait_col,
                "Trait_Value": value,
                "N_Accessions": np.nan,
                "N_Positive": np.nan,
                "N_Negative": np.nan,
                "Status": "insufficient_data",
                "Observed Purity": np.nan,
                "Null Mean": np.nan,
                "Purity p-value": np.nan,
                "Mantel r": np.nan,
                "Mantel p-value": np.nan,
                "Error": "",
            })
            continue

        try:
            metadata_dict, n_pos, n_neg = build_test_dict(
                base_tree, sub_meta, id_col, trait_col, value
            )
            if len(metadata_dict) == 0:
                raise ValueError(
                    "no tree leaves overlap with metadata for this value"
                )

            test_tree = base_tree.copy("newick-extended")
            result = run_value_test(test_tree, metadata_dict)
            rows.append({
                "Gene": snakemake.wildcards.gene,
                "Trait": trait_col,
                "Trait_Value": value,
                "N_Accessions": result["n_accessions"],
                "N_Positive": result["n_positive"],
                "N_Negative": result["n_negative"],
                "Status": "ok",
                "Observed Purity": result["observed_purity"],
                "Null Mean": result["null_mean"],
                "Purity p-value": result["purity_p_value"],
                "Mantel r": result["mantel_r"],
                "Mantel p-value": result["mantel_p_value"],
                "Error": "",
            })
            print(f"  - {value}: n={result['n_accessions']} "
                  f"(+{result['n_positive']}/-{result['n_negative']}) "
                  f"purity_p={result['purity_p_value']:.4g} "
                  f"mantel_p={result['mantel_p_value']:.4g}")
        except Exception as e:
            warnings.warn(f"[gene_trait] {trait_col}={value} failed: {e}")
            # Report the actually-matched counts when available, NaN otherwise.
            try:
                n_pos_safe = n_pos
                n_neg_safe = n_neg
                n_acc_safe = n_pos_safe + n_neg_safe
            except NameError:
                n_pos_safe = np.nan
                n_neg_safe = np.nan
                n_acc_safe = np.nan
            rows.append({
                "Gene": snakemake.wildcards.gene,
                "Trait": trait_col,
                "Trait_Value": value,
                "N_Accessions": n_acc_safe,
                "N_Positive": n_pos_safe,
                "N_Negative": n_neg_safe,
                "Status": "error",
                "Observed Purity": np.nan,
                "Null Mean": np.nan,
                "Purity p-value": np.nan,
                "Mantel r": np.nan,
                "Mantel p-value": np.nan,
                "Error": str(e),
            })

    df = pd.DataFrame(rows)
    df.to_csv(path_out, index=False)
    print(f"[gene_trait] wrote {len(df)} rows to {path_out}")


if __name__ == "__main__":
    main()
