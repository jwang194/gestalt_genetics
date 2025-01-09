import sys
import numpy as np
import pandas as pd

N,M = [int(s) for s in sys.argv[1:3]]
O = sys.argv[3]
R = int(sys.argv[4])
pheno_dir = sys.argv[5]

P = np.loadtxt(pheno_dir + '%i_%i_%s_rep%i_P.txt'%(N,M,O,R))
MG = np.loadtxt(pheno_dir + '%i_%i_%s_rep%i_MG.txt'%(N,M,O,R))
SD = np.loadtxt(pheno_dir + '%i_%i_%s_rep%i_SD.txt'%(N,M,O,R))

ids = np.array(['per%i'%i for i in range(N)])[:,None]
zeros = np.zeros((N, 1), dtype=int)
P = pd.DataFrame(np.hstack((zeros,ids,P)))
MG = pd.DataFrame(np.hstack((zeros,ids,MG)))
SD = pd.DataFrame(np.hstack((zeros,ids,SD)))
P.columns = MG.columns = ['FID','IID'] + ['PHEN%i'%i for i in range(1,P.shape[1]-1)]
SD.columns = ['FID','IID'] + ['PHEN%i'%i for i in range(1,3)]
P.to_csv(pheno_dir + '%i_%i_%s_rep%i_P.txt'%(N,M,O,R),sep='\t',index=False)
MG.to_csv(pheno_dir + '%i_%i_%s_rep%i_MG.txt'%(N,M,O,R),sep='\t',index=False)
SD.to_csv(pheno_dir + '%i_%i_%s_rep%i_SD.txt'%(N,M,O,R),sep='\t',index=False)
