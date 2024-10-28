import pandas as pd
from scipy.stats import pearsonr

nphen = 5
shared = 0.25
specific = 0.15

r = 1
N = 5*10**4
M = 10**4
C = 'all_and_ind_overlaps_uniform_gg_random_ge_0.1'

shared_m = int(shared*M)
specific_m = int(specific*M)

truth = pd.read_csv('../data/sim/phenos/%s_%s_%s_rep%s_B.txt'%(N,M,C,r),sep=' ',header=None)
truth.index = ['SNP_%i'%i for i in range(M)]

estimates = [pd.read_csv('../gwas/sim/%s_%s_%s_rep%s_P.PHEN%i.glm.linear'%(N,M,C,r,i),sep='\t')[['ID','BETA','P']] for i in range(1,6)]
for i in range(nphen):
    estimates[i].index = estimates[i].ID

estimates = [e.merge(truth[i],left_index=True,right_index=True) for e,i in zip(estimates,range(5))]

maxh_estimates = [pd.read_csv('../gwas/sim/%s_%s_%s_rep%s_MG.PHEN%i.glm.linear'%(N,M,C,r,i),sep='\t')[['ID','BETA','P']] for i in range(1,6)]
for i in range(nphen):
    maxh_estimates[i].index = maxh_estimates[i].ID

for i in range(nphen):
    print('SHARED EFFECT CORRELATION')
    pearsonr(estimates[i].iloc[:shared_m].BETA,estimates[i].iloc[:shared_m][i])
    pearsonr(maxh_estimates[4].iloc[:shared_m].BETA,truth.iloc[:shared_m][i])
    print('SHARED EFFECT POWER')
    (estimates[i].iloc[:shared_m].P < 0.05/M).mean()
    (maxh_estimates[4].iloc[:shared_m].P < 0.05/M).mean()

for i in range(nphen):
    print('SPECIFIC EFFECT CORRELATION')
    pearsonr(estimates[i].iloc[(shared_m+i*specific_m):(shared_m+i*specific_m+specific_m)].BETA,estimates[i].iloc[(shared_m+i*specific_m):(shared_m+i*specific_m+specific_m)][i])
    pearsonr(maxh_estimates[4].iloc[(shared_m+i*specific_m):(shared_m+i*specific_m+specific_m)].BETA,truth.iloc[(shared_m+i*specific_m):(shared_m+i*specific_m+specific_m)][i])
    print('SPECIFIC EFFECT POWER')
    (estimates[i].iloc[(shared_m+i*specific_m):(shared_m+i*specific_m+specific_m)].P < 5e-6).mean()
    (maxh_estimates[4].iloc[(shared_m+i*specific_m):(shared_m+i*specific_m+specific_m)].P < 5e-6).mean()
