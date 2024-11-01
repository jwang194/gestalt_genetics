import sys
import numpy as np
import pandas as pd

phenotype_files = sys.argv[1] # prefix to P and MG phenotype files 
train_prefix = sys.argv[2] # path to training fold plink files

P = pd.read_csv(phenotype_files + '_P.txt', sep='\t')
MG = pd.read_csv(phenotype_files + '_MG.txt', sep='\t')
fam = pd.read_csv(train_prefix + '.fam', sep='\t', header = None)

P_filtered = P[P['IID'].isin(fam[1])]
MG_filtered = MG[MG['IID'].isin(fam[1])]

scaled_P = P_filtered.copy()
scaled_P.iloc[:, 2:] = (P_filtered.iloc[:, 2:] - P_filtered.iloc[:, 2:].mean()) / P_filtered.iloc[:, 2:].std()
scaled_MG = MG_filtered.copy()
scaled_MG.iloc[:, 2:] = (MG_filtered.iloc[:, 2:] - MG_filtered.iloc[:, 2:].mean()) / MG_filtered.iloc[:, 2:].std()

scaled_P.to_csv(train_prefix + '_P.txt',sep='\t',index=False)
scaled_MG.to_csv(train_prefix + '_MG.txt',sep='\t',index=False)

