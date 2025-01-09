import sys
import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from itertools import chain

true_betas_path = sys.argv[1] # path to directory with true betas ex. '/u/home/k/kaia/GESTALT/data/sim/phenos/testing/'
gwas_path =sys.argv[2] # path to directory gwas results ex. '/u/home/k/kaia/GESTALT/data/sim/gwas/testing/'
output_path = sys.argv[3] # path to output directory ex. '/u/home/k/kaia/GESTALT/data/sim/results/testing/'
C = sys.argv[4] # settings name ex. 'testing'
nphen = int(sys.argv[5]) # number of phenotypes ex. '5'
shared = float(sys.argv[6]) # proportion of trait shared variants ex. '0.25'
specific = float(sys.argv[7]) # proportion of trait specific variants ex. '0.15'
N= int(sys.argv[8]) # number of individuals in this simulation ex. '1000'
M = int(sys.argv[9]) # number of total variants ex. '10000'

R = 25

shared_m = int(shared*M)
specific_m = int(specific*M)

gg = np.loadtxt('configs/%s_gg.txt'%C)

output = []

for r in range(1,R+1):
#for r in R:
    truth_dir = true_betas_path + '%s_%s_%s_rep%s_B.txt'%(N,M,C,r)
    if not os.path.isfile(truth_dir):
        print('replicate %i missing!'%r)
        continue
    truth = pd.read_csv(truth_dir, sep=' ',header=None)
    truth.index = ['SNP_%i'%i for i in range(M)]

    estimates = [pd.read_csv(gwas_path + '%s_%s_%s_rep%s_P.PHEN%i.glm.linear'%(N,M,C,r,i),sep='\t')[['ID','BETA','P']] for i in range(1, (nphen+1))]
    for i in range(nphen):
        estimates[i].index = estimates[i].ID

    estimates = [e.merge(truth[i],left_index=True,right_index=True) for e,i in zip(estimates,range(nphen))] # merges true betas to dataframe with estimated betas

    maxh_estimates = [pd.read_csv(gwas_path + '%s_%s_%s_rep%s_MG.PHEN%i.glm.linear'%(N,M,C,r,i),sep='\t')[['ID','BETA','P']] for i in range(1,(nphen+1))]
    for i in range(nphen):
        maxh_estimates[i].index = maxh_estimates[i].ID

    running_output = []
    for i in range(nphen):
        base_shared_corr = pearsonr(estimates[i].iloc[:shared_m].BETA,estimates[i].iloc[:shared_m][i])[0]
        maxh_shared_corr = pearsonr(maxh_estimates[(nphen-1)].iloc[:shared_m].BETA,truth.iloc[:shared_m][i])[0]
        base_shared_powr = (estimates[i].iloc[:shared_m].P < 0.05/M).mean()
        maxh_shared_powr = (maxh_estimates[(nphen-1)].iloc[:shared_m].P < 0.05/M).mean()
        base_specific_corr = pearsonr(estimates[i].iloc[(shared_m+i*specific_m):(shared_m+i*specific_m+specific_m)].BETA,estimates[i].iloc[(shared_m+i*specific_m):(shared_m+i*specific_m+specific_m)][i])[0] 
        maxh_specific_corr = pearsonr(maxh_estimates[(nphen-1)].iloc[(shared_m+i*specific_m):(shared_m+i*specific_m+specific_m)].BETA,truth.iloc[(shared_m+i*specific_m):(shared_m+i*specific_m+specific_m)][i])[0]
        base_specific_powr = (estimates[i].iloc[(shared_m+i*specific_m):(shared_m+i*specific_m+specific_m)].P < 0.05/M).mean() 
        maxh_specific_powr = (maxh_estimates[(nphen-1)].iloc[(shared_m+i*specific_m):(shared_m+i*specific_m+specific_m)].P < 0.05/M).mean()

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
output.to_csv(output_path+'%s_%s_%s_RESULTS.txt'%(N,M,C),sep='\t',index=None)
