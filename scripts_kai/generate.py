import sys
from sim import *
from pandas_plink import read_plink1_bin

N,M = [int(s) for s in sys.argv[1:3]]
alleles = 2
assignments = np.loadtxt(sys.argv[3],dtype=str)
gg = np.loadtxt(sys.argv[4])
ge = np.loadtxt(sys.argv[5])
rep = int(sys.argv[7])
out = sys.argv[9]

B,PE = simulate(N,M,alleles,assignments,gg,ge,sample_genotypes=False)
#genotypes,P,B,gen_comp,PE = simulate(N,M,alleles,assignments,gg,ge,sample_genotypes=True)
#print(np.cov(gen_comp, rowvar=False))
#print(np.cov(PE, rowvar=False))

genotypes = read_plink1_bin(*[sys.argv[6]+e for e in ['.bed','.bim','.fam']])
G = np.array(genotypes.sel(sample=['per%i'%i for i in range(N)]))[:N,:M]
G -= G.mean(0)
G /= G.std(0)

P = G @ B + PE

#np.savetxt('../data/sim/phenos/%i_%i_%s_rep%i_P.txt'%(N,M,sys.argv[8],rep),P)
#np.savetxt('../data/sim/phenos/%i_%i_%s_rep%i_B.txt'%(N,M,sys.argv[8],rep),B)
#np.savetxt('../data/sim/phenos/%i_%i_%s_rep%i_PE.txt'%(N,M,sys.argv[8],rep),PE)

np.savetxt(out + '%i_%i_%s_rep%i_P.txt'%(N,M,sys.argv[8],rep),P) # save phenotypes
np.savetxt(out + '%i_%i_%s_rep%i_B.txt'%(N,M,sys.argv[8],rep),B) # save true effect sizes
np.savetxt(out + '%i_%i_%s_rep%i_PE.txt'%(N,M,sys.argv[8],rep),PE) # save environmental component of phenotypes
