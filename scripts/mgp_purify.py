import sys
import numpy as np
from scipy.linalg import null_space
from maxgcp import fit_heritability

P = np.loadtxt(sys.argv[1])
PE = np.loadtxt(sys.argv[2])
G = P - PE

GC = G.T @ G
PC = P.T @ P
V = GC[:-1,-1:]
Z = null_space(V.T)

W = fit_heritability(GC[:-1,:-1], PC[:-1,:-1])
W = np.hstack((W, Z@fit_heritability(Z.T @ GC[:-1,:-1] @ Z, Z.T @ PC[:-1,:-1] @ Z)))
sums_and_diffs = np.ones((P.shape[1]-1,2))
sums_and_diffs[::2,1] = -1

np.savetxt(sys.argv[1].split('_P.txt')[0]+'_W.txt',W)
np.savetxt(sys.argv[3],P[:,:-1] @ W)
np.savetxt(sys.argv[4],P[:,:-1] @ sums_and_diffs)
