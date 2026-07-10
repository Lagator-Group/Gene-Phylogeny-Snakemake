import os
import pandas as pd
import numpy as np

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
        # scipy 1.11+: returns adjusted p-values directly
        out[mask] = _bh(p[mask], method="bh")
    else:
        # statsmodels: returns (reject, adjusted_p, _, _)
        _, adj, _, _ = _multipletests(p[mask], alpha=0.05, method="fdr_bh")
        out[mask] = adj
    # Wrap in a float Series so the result is a clean numeric dtype and
    # plays nicely with downstream boolean comparisons (e.g. < 0.05).
    return pd.Series(out, dtype=float)


dir_in = "data/gene_trait"

cols = [
    "Gene", "Trait", "Trait_Value", "N_Accessions", "N_Positive", "N_Negative",
    "Status", "Observed Purity", "Null Mean", "Purity p-value",
    "Parsimony Score", "Parsimony Null Mean", "Parsimony p-value",
    "Mantel r", "Mantel p-value", "Error",
]
df = pd.DataFrame(columns=cols)

for file in os.listdir(dir_in):
    if file == "summary.csv":
        continue
    df_temp = pd.read_csv(os.path.join(dir_in, file))
    df = pd.concat([df, df_temp], axis=0)

if len(df) > 0:
    df["Mantel p-value BH"] = (
        df.groupby(["Gene", "Trait"])["Mantel p-value"]
        .transform(benjamini_hochberg)
    )
    df["Purity p-value BH"] = (
        df.groupby(["Gene", "Trait"])["Purity p-value"]
        .transform(benjamini_hochberg)
    )
    df["Parsimony p-value BH"] = (
        df.groupby(["Gene", "Trait"])["Parsimony p-value"]
        .transform(benjamini_hochberg)
    )
    df["Significant BH (q<0.05)"] = (
        (df["Mantel p-value BH"].astype(float) < 0.05)
        | (df["Purity p-value BH"].astype(float) < 0.05)
        | (df["Parsimony p-value BH"].astype(float) < 0.05)
    ).fillna(False).astype(bool)

df.to_csv(snakemake.output[0], index=False)
