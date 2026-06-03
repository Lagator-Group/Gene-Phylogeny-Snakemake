from ete3 import Tree
import pandas as pd
import numpy as np
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr
from skbio.stats.distance import mantel, DistanceMatrix

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

def mantel_calc(tree1_path,tree2_path):
    t1 = Tree(tree1_path)
    t2 = Tree(tree2_path)

    # ============================================================
    # KEEP SHARED TAXA ONLY
    # ============================================================

    shared = sorted(
        set(t1.get_leaf_names()) &
        set(t2.get_leaf_names())
    )

    t1.prune(shared)
    t2.prune(shared)

    t1.unroot()
    t2.unroot()

    # ============================================================
    # BUILD COPHENETIC MATRIX
    # ============================================================



    m1 = cophenetic_matrix(t1, shared)
    m2 = cophenetic_matrix(t2, shared)

    # ============================================================
    # CONVERT TO DISTANCEMATRIX
    # ============================================================

    dm1 = DistanceMatrix(m1, ids=shared)
    dm2 = DistanceMatrix(m2, ids=shared)

    # ============================================================
    # MANTEL TEST
    # ============================================================

    r, p, n = mantel(
        dm1,
        dm2,
        method='pearson',
        permutations=999
    )

    print(f"Mantel r = {r:.4f}")
    print(f"p = {p:.6f}")
    print(f"n = {n}")

    df = pd.DataFrame(columns=['Gene 1', 'Gene 2', 'Mantel r', 'p', 'n'])

    df.loc[0] = [rep, gene, r, p, n]

    return df

rep_path = snakemake.input[0]
tf_path = snakemake.input[1]
path_out = snakemake.output[0]
values_notes = snakemake.output[0]

rep = snakemake.wildcards.rep
gene = snakemake.wildcards.gene

print(f'{rep} + {gene}')

try:
    df = mantel_calc(rep_path,tf_path)
    print('done')
except Exception as e:
    print(e)
    df = pd.DataFrame(columns=['Gene 1', 'Gene 2', 'Mantel r', 'p', 'n'])
    df.loc[0] = [rep, gene, np.nan, np.nan, np.nan]

df.to_csv(path_out, index=False)