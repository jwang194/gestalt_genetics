import sys
import numpy as np
import pandas as pd

phenotype_files = sys.argv[1] # prefix to P and MG phenotype files 
train_prefix = sys.argv[2] # path to training fold plink files
sums = int(sys.argv[3]) # 1 (true) / 0 (false) indicating whether or not to filter the _SUM.txt file too

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

if sums > 0:
    SUM = pd.read_csv(phenotype_files + '_SUM.txt', sep = '\t')
    SUM_filtered = SUM[SUM['IID'].isin(fam[1])]
    scaled_SUM = SUM_filtered.copy()
    scaled_SUM.iloc[:, 2:] = (scaled_SUM.iloc[:, 2:] - scaled_SUM.iloc[:, 2:].mean()) / scaled_SUM.iloc[:, 2:].std()
    scaled_SUM.to_csv(train_prefix + '_SUM.txt',sep='\t',index=False)
