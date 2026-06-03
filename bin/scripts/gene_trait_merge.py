import os
import pandas as pd

dir_in = 'data/gene_trait'

cols = ["Gene","Trait","Observed Purity","Null Mean", "Purity p-value","Mantel r","Mantel p-value"]
df = pd.DataFrame(columns=cols)

for file in os.listdir(dir_in):
    if file == 'summary.csv':
        continue
    df_temp = pd.read_csv(os.path.join(dir_in,file))
    df = pd.concat([df,df_temp],axis=0)

df.to_csv(os.path.join(dir_in,'summary.csv'), index=False)