import time
import numpy as np
import numba as nb

from functools import reduce
from itertools import product
from scipy.stats import matrix_normal
from scipy.stats import pearsonr,spearmanr

# @nb.jit(nopython=True)
def simulate(n,m,alleles,assignments,gg,ge,spikes=None,sample_genotypes=True):
    nphen = gg.shape[0]
    if sample_genotypes:
        mafs = np.random.uniform(0.01,0.5,m) # randomly drawing allele frequencies for each snp (between 0.01 and 0.5)
        genotypes = np.zeros((n,m))
        for a in range(alleles):
            genotypes += np.random.rand(n,m) < mafs
        for c in range(m):
            site = genotypes[:,c]
            genotypes[:,c] = (site - site.mean())/(1 if site.std() == 0 else site.std())
    subsets = set(assignments) # get unique elements in assignments
    betas = np.zeros((m,nphen)) # effect sizes of variants (m x number phenotypes) = each variant has an effect size for each phenotype
    counts = np.zeros((nphen,nphen))
    for s in subsets:
        indices = np.where(np.array(assignments) == s)[0] # get indices of snps with the same assignment (ex. shared snps)
        s = s.split(',')
        for p in product(s,s):
            p = [int(pe) for pe in p]
            counts[p[0],p[1]] += indices.shape[0] # counting number of shared/partially shared/specific snps for traits
    # gg /= counts 
    start_time = time.perf_counter()
    for s in subsets:
        indices = np.where(np.array(assignments) == s)[0]
        sub_m = indices.shape[0] # number of snps in that category (ex. specific to phenotype 0)
        s = [int(se) for se in s.split(',')]
        sub_p = len(s)
        sub_gg = gg[np.ix_(s,s)] # corresponding element in the genetic covariance matrix (diag = additive genetic variance of trait, offdiag = genetic correlation between traits)
        tmp_beta = np.zeros((m,sub_p)) # stores the effect sizes of snps in this category 
        tmp_beta[indices] = matrix_normal(np.zeros((sub_m,sub_p)),np.eye(sub_m),sub_gg/sub_m).rvs(1) # note: for snps that effect multiple traits, this line samples effect sizes jointly for these traits with the genetic covariance matrix sub_gg
        betas[:,s] += tmp_beta 
    # betas.std(0) = standard deviation of genetic effects for each phenotype
    betas *= np.sqrt(np.diag(gg)/m)/betas.std(0) # making sure the additive genetic variance matches the gg matrix diags
    print('betas completed in %f seconds'%(time.perf_counter() - start_time))
    start_time = time.perf_counter()
    env_comp = matrix_normal(np.zeros((n,nphen)),np.eye(n),ge).rvs(1)
    print('noise completed in %f seconds'%(time.perf_counter() - start_time))
    if sample_genotypes:
        gen_comp = genotypes @ betas
        return genotypes,gen_comp+env_comp,betas,gen_comp,env_comp
    else:
        return betas,env_comp

# @nb.jit(nopython=True)
def simulate_pairwise(n,m,alleles,sg,se,overlaps,overlaps_gg,overlaps_ge,spikes=None):
    nphen = overlaps.shape[0]
    overlaps = np.triu(overlaps)
    for i in range(nphen):
        overlaps[i,i] -= overlaps[i,(i+1):].sum() + overlaps[:i,i].sum()
    overlaps = np.round(overlaps,2)
    mafs = {}
    genotypes = {}
    betas = {}
    gen_comp = {}
    env_comp = {}
    for i in range(nphen):
        mafs[i] = {}
        genotypes[i] = {}
        betas[i] = {}
        gen_comp[i] = {}
        env_comp[i] = {}
        for j in range(i,nphen):
            print('%i -- %i'%(i,j))
            mp = int(m*overlaps[i,j])
            mafs[i][j] = np.random.uniform(0.01,0.5,mp)
            genotypes[i][j] = np.zeros((n,mp))
            for a in range(alleles):
                genotypes[i][j] += np.random.rand(n,mp) < mafs[i][j]
            for c in range(mp):
                site = genotypes[i][j][:,c]
                genotypes[i][j][:,c] = (site - site.mean())/(1 if site.std() == 0 else site.std())
            gg = overlaps_gg[i,j]
            ge = overlaps_ge[i,j]
            if i == j:
                betas[i][j] = [np.random.randn(mp)*np.sqrt(sg[i]/m)] if mp != 0 else [np.zeros(0)]
            else:
                # beta_cov = np.block([
                #     [sg[i]*np.diag(np.ones(mp)),        gg*np.diag(np.ones(mp))],
                #     [gg*np.diag(np.ones(mp)).T,      sg[j]*np.diag(np.ones(mp))],
                # ])/m
                # E = np.linalg.eigh(beta_cov)[0]
                # print(list(E[np.where(E<1e-8)[0]]))
                # b = np.random.multivariate_normal(np.zeros(2*mp),beta_cov) if mp != 0 else np.zeros(2*mp)
                b1 = np.sqrt(sg[i])*np.random.randn(mp)/np.sqrt(m)
                b2 = (gg*b1/b1.std() + np.sqrt(1-gg)*np.random.randn(mp))
                b2 /= b2.std()*np.sqrt(m/sg[j])
                betas[i][j] = [b1,b2]
            gen_comp[i][j] = [genotypes[i][j]@b[:,None] for b in betas[i][j]]
            # eps_cov = np.block([
            #     [se[i]*np.diag(np.ones(n)),     ge*np.diag(np.ones(n))],
            #     [ge*np.diag(np.ones(n)).T,      se[j]*np.diag(np.ones(n))],
            # ])/nphen
            # e = np.random.multivariate_normal(np.zeros(2*n),eps_cov)
            # env_comp[i][j] = [eps[:n],eps[n:]]
            # e1 = np.random.randn(n)*np.sqrt(se[i]/nphen)
            # e2 = (np.sqrt(ge)*e1 + np.random.randn(n)*np.sqrt(1-ge*se[i]/nphen))*np.sqrt(se[j]/nphen)
            env_comp[i][j] = [np.random.randn(n)*np.sqrt(se[i]/nphen),np.random.randn(n)*np.sqrt(se[j]/nphen)]
    G = np.zeros((n,1))
    B = [np.zeros(0) for i in range(nphen)]
    P = [np.zeros(n) for i in range(nphen)]
    GC = [np.zeros(n) for i in range(nphen)]
    EC = [np.zeros(n) for i in range(nphen)]
    for i in range(nphen):
        for j in range(i,nphen):
            print('%i -- %i'%(i,j))
            G = np.hstack((G,genotypes[i][j]))
            print(genotypes[i][j].shape)
            P[i] += gen_comp[i][j][0][:,0] + env_comp[i][j][0]
            B[i] = np.concatenate((B[i],betas[i][j][0]))
            GC[i] += gen_comp[i][j][0][:,0]
            EC[i] += env_comp[i][j][0]
            if j != i:
                P[j] += gen_comp[i][j][1][:,0] + env_comp[i][j][1]
                B[j] = np.concatenate((B[j],betas[i][j][1]))
                GC[j] += gen_comp[i][j][1][:,0]
                EC[j] += env_comp[i][j][1]
    return(G[:,1:],P,B,GC,EC,genotypes,betas,gen_comp,env_comp)

def heritability(G,kinship=True):
    n,m = G.shape
    if kinship:
        K = G
    else:
        K = (G @ G.T)/m
    quad = K - np.eye(n)
    trace = np.trace(K @ K)
    def loss(P):
        P -= P.mean()
        P /= P.std()
        numerator = P.T @ quad @ P
        return numerator/(trace-n)
    return loss

def genetic_correlation(G):
    n,m = G.shape
    K = (G @ G.T)/m
    quad = K - np.eye(n)
    trace = np.trace(K @ K)
    her = heritability(K)
    def loss(y1,y2):
        y1 -= y1.mean()
        y1 /= y1.std()
        y2 -= y2.mean()
        y2 /= y2.std()
        gamma = y1.T @ quad @ y2
        return gamma / (np.sqrt(her(y2)) * np.sqrt(her(y1))) / (trace - n)
    return her,loss
