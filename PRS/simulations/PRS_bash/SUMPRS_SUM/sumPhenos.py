import sys
import os
import numpy as np
import pandas as pd
from itertools import chain

# sum traits and save a new file 

phenotype_files = sys.argv[1] # prefix to P and MG phenotype files 

P = pd.read_csv(phenotype_files + '_P.txt', sep='\t')
P_sum = P.iloc[:,range(2)].copy()
P_sum['PHENSUM'] = P.iloc[:,range(2,len(P.columns))].sum(axis=1)
P_sum['PHENSUM'] = (P_sum['PHENSUM'] - P_sum['PHENSUM'].mean()) / P_sum['PHENSUM'].std()
P_sum.to_csv(phenotype_files + '_SUM.txt',sep='\t',index=False)


