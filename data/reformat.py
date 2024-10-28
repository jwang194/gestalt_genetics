import pandas as pd
import numpy as np

from functools import reduce
from collections import Counter

df = pd.read_csv('depression_neuroticism.pheno.gz',sep='\t')
df.columns = [s.split('.')[1] for s in df.columns]
df.insert(0,'FID',df.eid)
df.columns = ['FID','IID']+list(df.columns[2:])

quant = {'20433': [-121.,-818.],
         '20434': [-121.,-818.],
         '20127': [],
         '20442': [-818.,-999.]}
for q,b in quant.items():
    df[q].iloc[np.where([v in b for v in df[q]])[0]] = np.nan

binaries = {'20441': {-818.: np.nan, 1.: 2, 0.: 1}, 
            '20446': {-818.: np.nan, 1.: 2, 0.: 1}, 
            '20447': {-818.: np.nan, 1.: 2, 0.: 1}, 
            '20445': {-818.: np.nan, -313.: np.nan, -121.: np.nan, 1.: 2, 0.: 1}, 
            '20532': {-818.: np.nan, -121.: np.nan, 1.: 2, 0.: 1},
            '20435': {-818.: np.nan, -121.: np.nan, 1.: 2, 0.: 1},
            '20449': {-818.: np.nan, -121.: np.nan, 1.: 2, 0.: 1},
            '20450': {-818.: np.nan, -121.: np.nan, 1.: 2, 0.: 1},
            '20448': {-818.: np.nan, -121.: np.nan, 1.: 2, 0.: 1},
            '20437': {-818.: np.nan, -121.: np.nan, 1.: 2, 0.: 1},
            '20534': {1.: 2, 0.: 1},
            '20533': {1.: 2, 0.: 1}, 
            '20535': {1.: 2, 0.: 1}
}           

for b,dct in binaries.items():
    for o,n in dct.items():
        df[b].loc[df[b] == o] = n
    df[b] = df[b].astype('Int64')

likert = {'20518': {-818.: np.nan, 1.: 1, 2.: 2, 3.: 2, 4.: 2},
          '20510': {-818.: np.nan, 1.: 1, 2.: 2, 3.: 2, 4.: 2},
          '20507': {-818.: np.nan, 1.: 1, 2.: 2, 3.: 2, 4.: 2},
          '20519': {-818.: np.nan, 1.: 1, 2.: 2, 3.: 2, 4.: 2},
          '20514': {-818.: np.nan, 1.: 1, 2.: 2, 3.: 2, 4.: 2},
          '20511': {-818.: np.nan, 1.: 1, 2.: 2, 3.: 2, 4.: 2},
          '20513': {-818.: np.nan, 1.: 1, 2.: 2, 3.: 2, 4.: 2},
          '20508': {-818.: np.nan, 1.: 1, 2.: 2, 3.: 2, 4.: 2},
          '20517': {-818.: np.nan, 1.: 1, 2.: 2, 3.: 2, 4.: 2},
}

for b,dct in likert.items():
    df[b+'_quant'] = df[b]
    for o,n in dct.items():
        if np.isnan(n):
            df[b+'_quant'].loc[df[b] == o] = n
        else:
            df[b+'_quant'].loc[df[b] == o] = o
        df[b].loc[df[b] == o] = n
    df[b] = df[b].astype('Int64')

semiquant = {'20438': {-818.: np.nan,
                        1.: 0.,
                        2.: 1.5,
                        3.: 3.5,
                        4.: 9.,
                        5.: 18.,
                        6.: 24.},
             '20436': {-818.: np.nan,
                       -121.: np.nan,
                        1.: 0.25,
                        2.: 0.5,
                        3.: 0.75,
                        4.: 1.},
             '20439': {-818.: np.nan,
                       -121.: np.nan,
                        1.: 0.25,
                        2.: 0.75,
                        3.: 1},
             '20440': {-818.: np.nan,
                        0.: 0.,
                        1.: 1.,
                        2.: 2.,
                        3.: 3.},
}

for b,dct in semiquant.items():
    for o,n in dct.items():
        df[b].loc[df[b] == o] = n

del df['20536']

neur = {'1920': {-1.: np.nan, -3.: np.nan, 1.: 2, 0.: 1.},
        '1930': {-1.: np.nan, -3.: np.nan, 1.: 2, 0.: 1.},
        '1940': {-1.: np.nan, -3.: np.nan, 1.: 2, 0.: 1.},
        '1950': {-1.: np.nan, -3.: np.nan, 1.: 2, 0.: 1.},
        '1960': {-1.: np.nan, -3.: np.nan, 1.: 2, 0.: 1.},
        '1970': {-1.: np.nan, -3.: np.nan, 1.: 2, 0.: 1.},
        '1980': {-1.: np.nan, -3.: np.nan, 1.: 2, 0.: 1.},
        '1990': {-1.: np.nan, -3.: np.nan, 1.: 2, 0.: 1.},
        '2000': {-1.: np.nan, -3.: np.nan, 1.: 2, 0.: 1.},
        '2010': {-1.: np.nan, -3.: np.nan, 1.: 2, 0.: 1.},
        '2020': {-1.: np.nan, -3.: np.nan, 1.: 2, 0.: 1.},
        '2030': {-1.: np.nan, -3.: np.nan, 1.: 2, 0.: 1.}
}

for b,dct in neur.items():
    for o,n in dct.items():
        df[b].loc[df[b] == o] = n
    df[b] = df[b].astype('Int64')

non_quant = reduce(lambda x,y: x+y, [list(a.keys()) for a in [neur,binaries,likert]], [])

for n in non_quant:
    print(n)
    print(Counter(df[n]))

df.to_csv('../gwas/data/phens/depression_neuroticism.pheno',sep='\t',index=False,na_rep='NA')
np.savetxt('../gwas/data/phens/colnames.txt',df.columns,'%s')

phq = df[['FID','IID'] + [s+'_quant' for s in likert.keys()]]
phq -= 1
phq['PHQ9_Score'] = phq[phq.columns[2:]].sum(1)
phq['Mild_Exact'] = ((phq.PHQ9_Score > 5) & (phq.PHQ9_Score < 10)).astype('Int64')+1
phq['Moderate_Exact'] = ((phq.PHQ9_Score > 10) & (phq.PHQ9_Score < 15)).astype('Int64')+1
phq['Moderately_Severe_Exact'] = ((phq.PHQ9_Score > 15) & (phq.PHQ9_Score < 20)).astype('Int64')+1
phq['Severe'] = (phq.PHQ9_Score > 20).astype('Int64')+1
phq['Mild_Bound'] = (phq.PHQ9_Score > 5).astype('Int64')+1
phq['Moderate_Bound'] = (phq.PHQ9_Score > 10).astype('Int64')+1
phq['Moderately_Severe_Bound'] = (phq.PHQ9_Score > 15).astype('Int64')+1

phq.to_csv('../gwas/data/phens/phq.pheno',sep='\t',index=False,na_rep='NA')
np.savetxt('../gwas/data/phens/phq_colnames.txt',phq.columns[2:],'%s')
