import sys
import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from itertools import chain
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler 

# fit optimal linear combination of PRS to predict maxh phenotype

phenotype_files = sys.argv[1] # prefix to MG phenotype files 
prs_files = sys.argv[2] # prefix to PRS predictions of each phenotype 
traintest_fam = sys.argv[3] # fam file of traintest individuals 
gwas_prefix = sys.argv[4] 
nphen = int(sys.argv[5])
p_thresholds = ["0.000005", "0.001", "0.05"]

# read in true MG phenotype and filter individuals 
MG = pd.read_csv(phenotype_files + '_MG.txt', sep='\t')
fam = pd.read_csv(traintest_fam, sep='\t', header = None)
MG_filtered = MG[MG['IID'].isin(fam[1])]

# for each p-value threshold
for p in p_thresholds:
    # read in each phenotype PRS prediction and combine to dataframe, reorder to match roworder of true MG phenotype 
    PRS_predictions = [pd.read_csv(f"{prs_files}{phenonum}_thresh{p}_P.txt", sep='\t') for phenonum in range(1, nphen + 1)]
    PRS_predictions_df = PRS_predictions[0].rename(columns={'PRS': 'PHENO_1' })
    for i in range(1, nphen):
        PRS_predictions[i] = PRS_predictions[i].rename(columns={'PRS': ('PHENO_' + str(i+1)) })
        PRS_predictions_df = pd.merge(PRS_predictions_df, PRS_predictions[i], on=["FID", "IID"], how="inner")
    PRS_predictions_df_reordered = PRS_predictions_df.set_index('IID').reindex(MG_filtered['IID']).reset_index()
    # fit linear model and extract coefficients 
    scaler = StandardScaler() 
    model = LinearRegression()
    X, y = scaler.fit_transform(PRS_predictions_df_reordered.iloc[:,2:(nphen+2)]), MG_filtered.iloc[:,(nphen+1)]
    y = (y - y.mean())/y.std()
    model.fit(X, y)
    coeffs = model.coef_
    # read in GWAS for each phenotype and filter at this p-value threshold 
    GWAS_betas = [pd.read_csv(f"{gwas_prefix}.PHEN{phenonum}.glm.linear", sep='\t')[['ID','BETA','P']] for phenonum in range(1, nphen + 1)]
    GWAS_betas[0].loc[GWAS_betas[0]['P'] >= float(p), 'BETA'] = 0
    GWAS_betas_df = GWAS_betas[0].loc[:,['ID', 'BETA']].rename(columns={'BETA': 'BETA_1'})
    for i in range(1, nphen):
        GWAS_betas[i].loc[GWAS_betas[i]['P'] >= float(p), 'BETA'] = 0
        GWAS_betas[i] = GWAS_betas[i].loc[:,['ID', 'BETA']].rename(columns={'BETA': ('BETA_' + str(i+1)) })
        GWAS_betas_df = pd.merge(GWAS_betas_df, GWAS_betas[i], on=["ID"], how="inner")
    # linearly combine Betas using extracted coefficients
    BETAS_COMBINED = GWAS_betas_df.iloc[:,1:] @ coeffs
    GWAS_betas_df = pd.read_csv(f"{gwas_prefix}.PHEN1.glm.linear", sep='\t').iloc[:,[2,3,11]]
    GWAS_betas_df['BETA'] = BETAS_COMBINED
    # save Betas (these the SNP weights for MaxH phenotype, at a specific p-value threshold)
    GWAS_betas_df.to_csv(phenotype_files + '_thresh_' + p + '_SUM_PRS_BETAS.txt',sep='\t',index=False)
