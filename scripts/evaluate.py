import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

nphen = 5
shared = 0.25
specific = 0.15

R = 25
N = 5*10**4
M = 10**4
C = 'all_and_ind_overlaps_uniform_gg_random_ge_0.1'

shared_m = int(shared*M)
specific_m = int(specific*M)

gg = np.loadtxt('configs/%s_gg.txt'%C)

output = []

for r in range(1,R+1):
    truth_dir = '../data/sim/phenos/%s_%s_%s_rep%s_B.txt'%(N,M,C,r)
    if not os.path.isfile(truth_dir):
        print('replicate %i missing!'%r)
        continue
    truth = pd.read_csv(truth_dir,sep=' ',header=None)
    truth.index = ['SNP_%i'%i for i in range(M)]

    estimates = [pd.read_csv('../gwas/sim/%s_%s_%s_rep%s_P.PHEN%i.glm.linear'%(N,M,C,r,i),sep='\t')[['ID','BETA','P']] for i in range(1,6)]
    for i in range(nphen):
        estimates[i].index = estimates[i].ID

    estimates = [e.merge(truth[i],left_index=True,right_index=True) for e,i in zip(estimates,range(5))]

    maxh_estimates = [pd.read_csv('../gwas/sim/%s_%s_%s_rep%s_MG.PHEN%i.glm.linear'%(N,M,C,r,i),sep='\t')[['ID','BETA','P']] for i in range(1,6)]
    for i in range(nphen):
        maxh_estimates[i].index = maxh_estimates[i].ID

    running_output = []
    for i in range(nphen):
        base_shared_corr = pearsonr(estimates[i].iloc[:shared_m].BETA,estimates[i].iloc[:shared_m][i])[0]
        maxh_shared_corr = pearsonr(maxh_estimates[4].iloc[:shared_m].BETA,truth.iloc[:shared_m][i])[0]
        base_shared_powr = (estimates[i].iloc[:shared_m].P < 0.05/M).mean()
        maxh_shared_powr = (maxh_estimates[4].iloc[:shared_m].P < 0.05/M).mean()
        base_specific_corr = pearsonr(estimates[i].iloc[(shared_m+i*specific_m):(shared_m+i*specific_m+specific_m)].BETA,estimates[i].iloc[(shared_m+i*specific_m):(shared_m+i*specific_m+specific_m)][i])[0] 
        maxh_specific_corr = pearsonr(maxh_estimates[4].iloc[(shared_m+i*specific_m):(shared_m+i*specific_m+specific_m)].BETA,truth.iloc[(shared_m+i*specific_m):(shared_m+i*specific_m+specific_m)][i])[0]
        base_specific_powr = (estimates[i].iloc[(shared_m+i*specific_m):(shared_m+i*specific_m+specific_m)].P < 5e-6).mean() 
        maxh_specific_powr = (maxh_estimates[4].iloc[(shared_m+i*specific_m):(shared_m+i*specific_m+specific_m)].P < 5e-6).mean()
        running_output.append([base_shared_corr,'Base GWAS: Shared Effects','Correlation',gg[i,i]])
        running_output.append([base_specific_corr,'Base GWAS: Specific Effects','Correlation',gg[i,i]])
        running_output.append([maxh_shared_corr,'MaxH GWAS: Shared Effects','Correlation',gg[i,i]])
        running_output.append([maxh_specific_corr,'MaxH GWAS: Specific Effects','Correlation',gg[i,i]])
        running_output.append([base_shared_powr,'Base GWAS: Shared Effects','Power',gg[i,i]])
        running_output.append([base_specific_powr,'Base GWAS: Specific Effects','Power',gg[i,i]])
        running_output.append([maxh_shared_powr,'MaxH GWAS: Shared Effects','Power',gg[i,i]])
        running_output.append([maxh_specific_powr,'MaxH GWAS: Specific Effects','Power',gg[i,i]])

    output.extend(running_output)

output = pd.DataFrame(output)
output.columns = ['Value','Source','Metric','Heritability']
output.to_csv('../out/%s_%s_%s_RESULTS.txt'%(N,M,C),sep='\t',index=None)
