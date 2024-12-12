import sys
import os
import numpy as np
import pandas as pd
from itertools import chain

# sum PRS and save a new file 

prs_files = sys.argv[1] # dir to prs files ex. /u/home/k/kaia/GESTALT/data/sim/PRS/N50kM10k_5traits_rg0.1_re-0.1_SUM/
nphen = int(sys.argv[2]) # number of phenotypes
p_threshold = str(sys.argv[3])
rep = str(sys.argv[4]) # repetition number

PRS_predictions = [pd.read_csv(f"{prs_files}prs_predictions_pheno{phenonum}_thresh{p_threshold}rep{rep}_P.txt", sep='\t') for phenonum in range(1, nphen + 1)]
PRS_predictions_df = PRS_predictions[0].rename(columns={'PRS': 'PRS_1' })
for i in range(1, nphen):
    PRS_predictions[i] = PRS_predictions[i].rename(columns={'PRS': ('PRS_' + str(i+1)) })
    PRS_predictions_df = pd.merge(PRS_predictions_df, PRS_predictions[i], on=["FID", "IID"], how="inner")
PRS_sum = PRS_predictions_df.iloc[:,range(2)].copy()
PRS_sum['SUMPRS'] = PRS_predictions_df.iloc[:,range(2,len(PRS_predictions_df.columns))].sum(axis=1)
PRS_sum.to_csv(f"{prs_files}prs_predictions_thresh{p_threshold}rep{rep}_SUM_sumPRS.txt",sep='\t',index=False)


