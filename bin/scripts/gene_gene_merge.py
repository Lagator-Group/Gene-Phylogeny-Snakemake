import os
import pandas as pd

def merge_df(dir_in, df):
    for file in os.listdir(dir_in):
        if file == 'summary.csv':
            continue
        df_temp = pd.read_csv(os.path.join(dir_in,file))
        df = pd.concat([df,df_temp],axis=0)
    return df

dir_in = 'data/alisim'

cols = ['Gene1', 'Gene2', 'Mantel r', 'p', 'n', 'Mantel null mean', 'Mantel null sd', 'Mantel Sim p','RF', 'Max RF', 'Normalized RF', 'Shared taxa', 'Permutation p-value', 'nRF null mean','nRF null sd','nRF Sim p']
df = pd.DataFrame(columns=cols)
df = merge_df(dir_in, df)

df.to_csv(snakemake.output[0], index=False)