import os
import pandas as pd

def merge_df(dir_in, df):
    for file in os.listdir(dir_in):
        if file == 'summary.csv':
            continue
        df_temp = pd.read_csv(os.path.join(dir_in,file))
        df = pd.concat([df,df_temp],axis=0)
    return df

### Mantel
dir_in = 'data/mantel'

cols = ['Gene 1', 'Gene 2', 'Mantel r', 'p', 'n']
df = pd.DataFrame(columns=cols)
df = merge_df(dir_in, df)

df.to_csv(os.path.join(dir_in,'summary.csv'), index=False)

### Mantel
dir_in = 'data/rf'

cols = ['Gene 1', 'Gene 2', 'RF', 'Max RF', 'Normalized RF', 'Shared taxa', 'Permutation p-value']
df = pd.DataFrame(columns=cols)
df = merge_df(dir_in, df)

df.to_csv(os.path.join(dir_in,'summary.csv'), index=False)