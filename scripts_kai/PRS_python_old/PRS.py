import sys
import os
import numpy as np
import xarray as xr
import pandas as pd
from scipy.stats import pearsonr,spearmanr
from pandas_plink import read_plink1_bin
from PRS_functions import *


plink_genotype_prefix = sys.argv[1] # path to plink genotype files 
gwas_path = sys.argv[2] # path to directory gwas results ex. '/u/home/k/kaia/GESTALT/data/sim/gwas/testing/'
pheno_path = sys.argv[3] # path to directory gwas results ex. '/u/home/k/kaia/GESTALT/data/sim/phenos/testing/'
output_path = sys.argv[4] # path to output directory ex. '/u/home/k/kaia/GESTALT/data/sim/results/testing/'
C = sys.argv[5] # settings name ex. 'testing'
nphen = int(sys.argv[6]) # number of phenotypes ex. '5'
shared = float(sys.argv[7]) # proportion of trait shared variants ex. '0.25'
specific = float(sys.argv[8]) # proportion of trait specific variants ex. '0.15'
N = int(sys.argv[9]) # number of individuals in this simulation ex. '1000'
M = int(sys.argv[10]) # number of total variants ex. '10000'

#R = 25
R = ['1234']

shared_m = int(shared*M)
specific_m = int(specific*M)

# read in genotypes (people x snps)
genotypes = read_plink1_bin(*[plink_genotype_prefix+e for e in ['.bed','.bim','.fam']])
G = np.array(genotypes.sel(sample=['per%i'%i for i in range(N)]))[:N,:M]
G -= G.mean(0)
G /= G.std(0)

#for r in range(1,R+1):
for r in R:
    true_pheno_dir = pheno_path + '%s_%s_%s_rep%s_P.txt'%(N,M,C,r)
    if not os.path.isfile(true_pheno_dir):
        print('replicate %i missing!'%r)
        continue
    truth_pheno = pd.read_csv(true_pheno_dir, sep='\t')

    true_pheno_maxh_dir = pheno_path + '%s_%s_%s_rep%s_MG.txt'%(N,M,C,r)
    if not os.path.isfile(true_pheno_maxh_dir):
        print('replicate %i missing!'%r)
        continue
    truth_pheno_maxh = pd.read_csv(true_pheno_maxh_dir, sep='\t')

    estimates = [pd.read_csv(gwas_path + '%s_%s_%s_rep%s_P.PHEN%i.glm.linear'%(N,M,C,r,i),sep='\t')[['ID','BETA','P']] for i in range(1, (nphen+1))]
    for i in range(nphen):
        estimates[i].index = estimates[i].ID

    maxh_estimates = [pd.read_csv(gwas_path + '%s_%s_%s_rep%s_MG.PHEN%i.glm.linear'%(N,M,C,r,i),sep='\t')[['ID','BETA','P']] for i in range(1,(nphen+1))]
    for i in range(nphen):
        maxh_estimates[i].index = maxh_estimates[i].ID

    # filter snps to compute PRS with 

    beta_matrix = pd.DataFrame({f'BETA_{i}': df['BETA'].values for i, df in enumerate(estimates)}).to_numpy()
    r_list_base = PRS_evaluate_multi(truth_pheno.iloc[:, 2:], G, beta_matrix)
    print(r_list_base)

    beta_matrix_maxh = pd.DataFrame({f'BETA_{i}': df['BETA'].values for i, df in enumerate(maxh_estimates)}).to_numpy()
    r_list_maxh = PRS_evaluate_multi(truth_pheno_maxh.iloc[:, 2:], G, beta_matrix_maxh)
    print(r_list_maxh)

# 5-fold CV