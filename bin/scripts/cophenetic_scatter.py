from ete3 import Tree
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression


def build_distance_matrix(tree, taxa):

    n = len(taxa)

    mat = np.zeros((n, n))

    name2node = {
        leaf.name: leaf
        for leaf in tree.iter_leaves()
    }

    for i, t1 in enumerate(taxa):

        for j, t2 in enumerate(taxa):

            if j < i:
                mat[i, j] = mat[j, i]
                continue

            d = tree.get_distance(
                name2node[t1],
                name2node[t2]
            )

            mat[i, j] = d
            mat[j, i] = d

    return pd.DataFrame(
        mat,
        index=taxa,
        columns=taxa
    )

def cophenetic_scatter(tree1_path, tree2_path, save_path):
    tree1 = Tree(tree1_path)
    tree2 = Tree(tree2_path)

    # ============================================================
    # KEEP SHARED TAXA
    # ============================================================

    shared = sorted(
        set(tree1.get_leaf_names()) &
        set(tree2.get_leaf_names())
    )

    print(f"{len(shared)} shared taxa")

    tree1.prune(shared)
    tree2.prune(shared)

    tree1.unroot()
    tree2.unroot()

    # ============================================================
    # BUILD DISTANCE MATRICES
    # ============================================================



    dist1 = build_distance_matrix(tree1, shared)
    dist2 = build_distance_matrix(tree2, shared)

    # ============================================================
    # FLATTEN UPPER TRIANGLES
    # ============================================================

    tri_idx = np.triu_indices_from(dist1, k=1)

    vals1 = dist1.values[tri_idx]
    vals2 = dist2.values[tri_idx]

    # ============================================================
    # CORRELATION
    # ============================================================

    corr, p = pearsonr(vals1, vals2)

    print(f"Cophenetic correlation = {corr:.4f}")
    print(f"P-value = {p:.3e}")

    # ============================================================
    # LINEAR REGRESSION
    # ============================================================

    X = vals1.reshape(-1, 1)

    model = LinearRegression()
    model.fit(X, vals2)

    xfit = np.linspace(vals1.min(), vals1.max(), 200)
    yfit = model.predict(xfit.reshape(-1, 1))

    # ============================================================
    # PLOT
    # ============================================================

    plt.figure(figsize=(8, 8))

    plt.scatter(
        vals1,
        vals2,
        s=5,
        alpha=0.3
    )

    # regression line
    plt.plot(
        xfit,
        yfit,
        linewidth=2
    )

    # identity line
    mn = min(vals1.min(), vals2.min())
    mx = max(vals1.max(), vals2.max())

    plt.plot(
        [mn, mx],
        [mn, mx],
        linestyle="--"
    )

    plt.xlabel(f"{rep} pairwise distances")
    plt.ylabel(f"{gene} pairwise distances")

    plt.title(
        "Cophenetic Distance Comparison\n"
        f"Pearson r = {corr:.3f}"
    )

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()

    # plt.show()

rep_path = snakemake.input[0]
tf_path = snakemake.input[1]
path_out = snakemake.output[0]

rep = snakemake.wildcards.gene1
gene = snakemake.wildcards.gene2

print(f'{rep} + {gene}')

try:
    cophenetic_scatter(rep_path,tf_path,path_out)
    print('done')

except Exception as e:
    print(e)
    fig, ax = plt.subplots()
    ax.text(0.5, 0.5, "No data / plot failed",
            ha="center", va="center")
    ax.set_axis_off()
    fig.savefig(path_out, bbox_inches="tight")
    plt.close(fig)
    print('passing')
