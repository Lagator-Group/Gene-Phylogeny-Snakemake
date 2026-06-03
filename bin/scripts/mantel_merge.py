import os
import pandas as pd

dir_in = 'data/mantel'

cols = ['Gene 1', 'Gene 2', 'Mantel r', 'p', 'n']
df = pd.DataFrame(columns=cols)

for file in os.listdir(dir_in):
    if file == 'summary.csv':
        continue
    df_temp = pd.read_csv(os.path.join(dir_in,file))
    df = pd.concat([df,df_temp],axis=0)

df.to_csv(os.path.join(dir_in,'summary.csv'), index=False)