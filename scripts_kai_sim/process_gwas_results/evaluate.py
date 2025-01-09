import sys
import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from itertools import chain

### This script was edited because the GWAS simulation framework was changed ### 
    ### The genotypes have LD and random subsets of SNPs were chosen as causal ###
    ### Therefore, causal SNPs for each trait have to be retrieved from the assignments file ###

true_betas_path = sys.argv[1] # path to directory with true betas ex. '/u/scratch/k/kaia/GESTALT/sim/phenos/traits_10_causal_0.01_shared_0.05_uniform_rg_0.01_random_re_-0.1'
gwas_path =sys.argv[2] # path to directory gwas results ex. '/u/scratch/k/kaia/GESTALT/sim/gwas/traits_10_causal_0.01_shared_0.05_uniform_rg_0.01_random_re_-0.1'
output_path = sys.argv[3] # path to output directory ex. '/u/home/k/kaia/GESTALT/sim/gwas_results/traits_10_causal_0.01_shared_0.05_uniform_rg_0.01_random_re_-0.1'
setting = sys.argv[4] # settings name ex. 'traits_10_causal_0.01_shared_0.05_uniform_rg_0.01_random_re_-0.1'
N = sys.argv[5] # ex. 50000
M = sys.argv[6] # ex. 98163
nphen = int(sys.argv[7]) # number of phenotypes ex. '10'
R = int(sys.argv[8])
configs = sys.argv[9]

gg = np.loadtxt('%s/%s/%s_gg.txt'%(configs, setting, setting), delimiter=' ', dtype=float)

results_power = np.empty((0, 5))
# rep   heritability   model   effect   power
results_corr = np.empty((0, 5)) # corr with true beta for each phenotype
# rep   heritability    model   effect   corr

for r in range(1,R+1):
    true_betas_file = true_betas_path + '/%s_%s_%s_rep%s_B.txt'%(N,M,setting,r)
    if not os.path.isfile(true_betas_file):
        print('replicate %i missing!'%r)
        continue
    true_betas = pd.read_csv(true_betas_file, sep='\t')
    true_betas.index = true_betas['SNP']

    estimates = [pd.read_csv(gwas_path + '/%s_%s_%s_rep%s_P.PHEN%i.glm.linear'%(N,M,setting,r,i),sep='\t')[['ID','BETA','P']] for i in range(1, (nphen+1))]
    for i in range(nphen):
        estimates[i].index = estimates[i].ID
    
    maxh_estimates = [pd.read_csv(gwas_path + '/%s_%s_%s_rep%s_MG.PHEN%i.glm.linear'%(N,M,setting,r,i),sep='\t')[['ID','BETA','P']] for i in range(1,(nphen+1))]
    for i in range(nphen):
        maxh_estimates[i].index = maxh_estimates[i].ID

    sum_estimates = pd.read_csv(gwas_path + '/%s_%s_%s_rep%s_SD.PHEN%i.glm.linear'%(N,M,setting,r,1),sep='\t')[['ID','BETA','P']]
    sum_estimates.index = sum_estimates.ID

    shared_snps = true_betas.loc[true_betas['assignment'] == ','.join(str(i) for i in range(nphen)), :].index

    # shared snps power for maxh and sum 
    maxh_shared_powr = (maxh_estimates[(nphen-1)].loc[shared_snps, 'P'] < 0.05/float(M)).mean()
    sum_shared_powr = (sum_estimates.loc[shared_snps, 'P'] < 0.05/float(M)).mean()
    results_power = np.vstack([results_power, np.array([r, "MAXH", "MAXH", "SHARED", maxh_shared_powr])])
    results_power = np.vstack([results_power, np.array([r, "SUM", "SUM", "SHARED", sum_shared_powr])])

    for i in range(nphen):
        base_shared_corr = pearsonr(estimates[i].loc[shared_snps, 'BETA'],true_betas.loc[shared_snps, 'pheno'+str(i) ])[0]
        maxh_shared_corr = pearsonr(maxh_estimates[(nphen-1)].loc[shared_snps, 'BETA'],true_betas.loc[shared_snps, 'pheno'+str(i) ])[0]
        sum_shared_corr = pearsonr(sum_estimates.loc[shared_snps, 'BETA'],true_betas.loc[shared_snps, 'pheno'+str(i) ])[0]
        
        base_shared_powr = (estimates[i].loc[shared_snps, 'P'] < 0.05/float(M)).mean()

        specific_snps = true_betas.loc[true_betas['assignment'] == str(i), :].index
        base_specific_corr = pearsonr(estimates[i].loc[specific_snps, 'BETA'],true_betas.loc[specific_snps, 'pheno'+str(i) ])[0] 
        maxh_specific_corr = pearsonr(maxh_estimates[(nphen-1)].loc[specific_snps, 'BETA'],true_betas.loc[specific_snps, 'pheno'+str(i) ])[0]
        sum_specific_corr = pearsonr(sum_estimates.loc[specific_snps, 'BETA'],true_betas.loc[specific_snps, 'pheno'+str(i) ])[0]

        base_specific_powr = (estimates[i].loc[specific_snps, 'P'] < 0.05/float(M)).mean() 
        maxh_specific_powr = (maxh_estimates[(nphen-1)].loc[specific_snps, 'P'] < 0.05/float(M)).mean()
        sum_specific_powr = (sum_estimates.loc[specific_snps, 'P'] < 0.05/float(M)).mean()

        results_power = np.vstack([results_power, np.array([r, gg[i,i], "BASE", "SHARED", base_shared_powr]) ])
        results_power = np.vstack([results_power, np.array([r, gg[i,i], "BASE", "SPECIFIC", base_specific_powr])])
        results_power = np.vstack([results_power, np.array([r, gg[i,i], "MAXH", "SPECIFIC", maxh_specific_powr])])
        results_power = np.vstack([results_power, np.array([r, gg[i,i], "SUM", "SPECIFIC", sum_specific_powr])])

        results_corr = np.vstack([results_corr, np.array([r, gg[i,i], "BASE", "SHARED", base_shared_corr])])
        results_corr = np.vstack([results_corr, np.array([r, gg[i,i], "MAXH", "SHARED", maxh_shared_corr])])
        results_corr = np.vstack([results_corr, np.array([r, gg[i,i], "SUM", "SHARED", sum_shared_corr])])
        results_corr = np.vstack([results_corr, np.array([r, gg[i,i], "BASE", "SPECIFIC", base_specific_corr])])
        results_corr = np.vstack([results_corr, np.array([r, gg[i,i], "MAXH", "SPECIFIC", maxh_specific_corr])])
        results_corr = np.vstack([results_corr, np.array([r, gg[i,i], "SUM", "SPECIFIC", sum_specific_corr])])


output_power = pd.DataFrame(results_power)
output_power.columns = ['rep','heritability','model','effect', 'power'] 
# heritability indicates heritability of target trait
    # 1, 0.1, MAXH, SPECIFIC, POWER 
        # how well does GWAS on MAXH identify snps specific to a trait with heritability 0.1 in rep 1
output_corr = pd.DataFrame(results_corr)
output_corr.columns = ['rep','heritability','model','effect', 'corr']
# heritability indicates heritability of target trait
    # 1, 0.1, MAXH, SHARED, CORR 
        # how well does GWAS effect sizes on MAXH correlate with true shared snp effect sizes on a trait with heritability 0.1 in rep 1

output_power.to_csv(output_path+'/%s_%s_%s_POWER.txt'%(N,M,setting),sep='\t',index=None)
output_corr.to_csv(output_path+'/%s_%s_%s_CORR.txt'%(N,M,setting),sep='\t',index=None)