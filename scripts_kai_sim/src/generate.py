import sys
from sim_kai import *
from pandas_plink import read_plink1_bin
import pandas as pd
import numpy as np

N,M = [int(s) for s in sys.argv[1:3]] # N = number of samples, M = total number of SNPs
alleles = 2
assignments = np.loadtxt(sys.argv[3],dtype=str) # sorted so that it is the same order as .bim 
num_causal = len(assignments[:,0])
gg = np.loadtxt(sys.argv[4])
ge = np.loadtxt(sys.argv[5])
rep = int(sys.argv[7])
out = sys.argv[9]

genotypes = read_plink1_bin(*[sys.argv[6]+e for e in ['.bed','.bim','.fam']]) # already filtered to causal variants to improve memory and runtime efficiency
G = np.array(genotypes.sel(sample=['per%i'%i for i in range(N)]))[:N,:]
G -= G.mean(0)
G /= G.std(0)

# few improvements for memory and runtime efficiency 
# cholesky decomposition of covariance matrix for faster sampling (can be used for genotype sampling with LD as well)
# parallelize using numba 
B, PG, PE = simulate_phenotypes(N, num_causal, G, assignments, gg, ge) 
P = PG + PE
BETAS = np.column_stack((assignments, B))
BETAS_df = pd.DataFrame(BETAS, columns = ["SNP", "assignment"] + [f"pheno{i}" for i in range(1, gg.shape[0]+1)])

BETAS_df.to_csv(out + '%i_%i_%s_rep%i_B.txt'%(N,M,sys.argv[8],rep), sep = '\t', index = False) # save true effect sizes
np.savetxt(out + '%i_%i_%s_rep%i_P.txt'%(N,M,sys.argv[8],rep),P) # save phenotypes
#np.savetxt(out + '%i_%i_%s_rep%i_B.txt'%(N,M,sys.argv[8],rep),BETAS) # save true effect sizes
np.savetxt(out + '%i_%i_%s_rep%i_PE.txt'%(N,M,sys.argv[8],rep),PE) # save environmental component of phenotypes
