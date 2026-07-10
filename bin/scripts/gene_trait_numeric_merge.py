"""
Merge the per (gene, trait) continuous-signal tables produced by
gene_trait_numeric.py into a single summary and apply Benjamini-Hochberg
correction, within each Trait, to every p-value column.
"""

import os

import numpy as np
import pandas as pd

try:
    from scipy.stats import false_discovery_control as _bh
    _USE_NEW_API = True
except ImportError:
    from statsmodels.stats.multitest import multipletests as _multipletests
    _USE_NEW_API = False


def benjamini_hochberg(p):
    """Return BH-adjusted p-values. NaN inputs pass through as NaN."""
    p = np.asarray(p, dtype=float)
    mask = ~np.isnan(p)
    out = np.full_like(p, np.nan)
    if mask.sum() == 0:
        return out
    if _USE_NEW_API:
        out[mask] = _bh(p[mask], method="bh")
    else:
        _, adj, _, _ = _multipletests(p[mask], alpha=0.05, method="fdr_bh")
        out[mask] = adj
    # Return a plain ndarray so groupby(...).transform assigns positionally.
    # (Returning a Series with its own RangeIndex would misalign here, because
    # every per-(gene,trait) input file contributes a single row with index 0.)
    return out


dir_in = "data/gene_trait_numeric"

cols = [
    "Gene", "Trait", "N_Accessions", "Trait_Mean", "Trait_SD", "Status",
    "Blomberg K", "Blomberg K p-value",
    "Pagel lambda", "Pagel lambda p-value",
    "Morans I", "Morans I p-value",
    "Mantel r", "Mantel p-value", "Error",
]
df = pd.DataFrame(columns=cols)

for file in os.listdir(dir_in):
    if file == "summary.csv":
        continue
    df_temp = pd.read_csv(os.path.join(dir_in, file))
    df = pd.concat([df, df_temp], axis=0)

df = df.reset_index(drop=True)

# (p-value column, BH output column) pairs to correct, grouped within Trait.
pval_cols = [
    ("Blomberg K p-value", "Blomberg K p-value BH"),
    ("Pagel lambda p-value", "Pagel lambda p-value BH"),
    ("Morans I p-value", "Morans I p-value BH"),
    ("Mantel p-value", "Mantel p-value BH"),
]

if len(df) > 0:
    for src, dst in pval_cols:
        df[dst] = df.groupby("Trait")[src].transform(benjamini_hochberg)

    bh_cols = [dst for _, dst in pval_cols]
    df["Significant BH (q<0.05)"] = (
        (df[bh_cols].astype(float) < 0.05).any(axis=1).fillna(False).astype(bool)
    )

df.to_csv(snakemake.output[0], index=False)
