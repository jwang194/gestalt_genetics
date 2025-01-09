import time
import numpy as np

from numba import jit, prange, get_num_threads
from functools import reduce
from itertools import product
from scipy.stats import matrix_normal
from scipy.stats import pearsonr,spearmanr

def timeit_decorator(func):
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        end = time.perf_counter()
        print(f"{func.__name__} execution time: {end - start:.6f} seconds")
        return result
    return wrapper

def sample_genotypes(n, m, alleles = 2):
    # n = number of individuals, m = number of snps
    # note that this function does not simulate any LD structure
        # this can be parameterized by reference LD as done in twas_sim tool to approximate LD for a given ancestry
    mafs = np.random.uniform(0.01,0.5,m) # randomly drawing allele frequencies for each snp (between 0.01 and 0.5)
    genotypes = np.zeros((n,m))
    for a in range(alleles):
        genotypes += np.random.rand(n,m) < mafs
    for c in range(m):
        site = genotypes[:,c]
        genotypes[:,c] = (site - site.mean())/(1 if site.std() == 0 else site.std())
    return genotypes

@timeit_decorator
@jit(nopython=True, parallel=True)
def generate_matrix_normal_parallel(n, m, S): # about 19x faster than generate_matrix_normal without parallel sampling for 10k snps
    # n = number of rows, m = number of columns
    # S  = covariance matrix of columns 
    # this function assumes that the rows are independent
        # if rows are independent, rows can be sampled in parallel using numba package 
        # the cholesky decomposition of S is used to enforce the random normal samples to have the desired covariance structure
            # since the S is the same for every row, the cholesky decomposition only has to be performed once, reducing overhead
    
    if m == 1:
        Z = np.zeros((n, m))
        sd = np.sqrt(S.item())
        for i in prange(n): # for each row in parallel
            Z[i, :] = np.random.normal(0, sd, m) # mean, standard deviation, number of samples
        return Z

    L = np.linalg.cholesky(S)
    Z = np.zeros((n, m))
    for i in prange(n): # for each row in parallel
        Z[i, :] = np.random.normal(0, 1, m) # Z is a vector of m samples from standard normal distribution
    samples = Z @ L.T
    return samples

@timeit_decorator
def generate_matrix_normal(n, m, S):
    # n = number of rows, m = number of columns
    # S  = covariance matrix of columns 
    # this function assumes that the rows are independent
        # if rows are independent, rows can be sampled in parallel using numba package 
        # the cholesky decomposition of S is used to enforce the random normal samples to have the desired covariance structure
            # since the S is the same for every row, the cholesky decomposition only has to be performed once, reducing overhead
    
    if m == 1:
        Z = np.zeros((n, m))
        sd = np.sqrt(S.item())
        for i in range(n): # for each row in parallel
            Z[i, :] = np.random.normal(0, sd, m) # mean, standard deviation, number of samples
        return Z
    
    L = np.linalg.cholesky(S)
    Z = np.zeros((n, m))
    for i in range(n): # for each row, not parallel 
        Z[i, :] = np.random.normal(0, 1, m) # Z is a vector of m samples from standard normal distribution
    samples = Z @ L.T
    return samples

def simulate_phenotypes(n,m_causal,genotypes,assignments_matrix,gg,ge):
    # n = number of individuals, m_causal = number of causal snps
    # genotypes = n x m_causal matrix of normalized genotypes 
    # assignments_matrix = numpy matrix with m_causal rows and 2 columns (column 1 = variant index, column 2 = indicates phenotype this snp effects)
    # gg = nphen x nphen numpy matrix of genetic covariance between phenotypes
    # ge = n x nphen numpy matrix of environmental covariance between phenotypes
    print('Number of threads for simulate_phenotypes: %s'%(get_num_threads()))
    nphen = gg.shape[0]
    assignments = assignments_matrix[:,1]
    subsets = set(assignments) # get unique elements in assignments
    print('Simulating betas')
    betas = np.zeros((m_causal,nphen)) # effect sizes of variants (m x number phenotypes) = each variant has an effect size for each phenotype
    for s in subsets:
        indices = np.where(np.array(assignments) == s)[0]
        sub_m = indices.shape[0] # number of snps in that category (ex. specific to phenotype 0)
        s = [int(se) for se in s.split(',')]
        sub_p = len(s)
        sub_gg = gg[np.ix_(s,s)] # corresponding element in the genetic covariance matrix (diag = additive genetic variance of trait, offdiag = genetic correlation between traits)
        tmp_beta = np.zeros((m_causal,sub_p)) # stores the effect sizes of snps in this category 
        #tmp_beta[indices] = matrix_normal(np.zeros((sub_m,sub_p)),np.eye(sub_m),sub_gg/sub_m).rvs(1) # note: for snps that effect multiple traits, this line samples effect sizes jointly for these traits with the genetic covariance matrix sub_gg
        tmp_beta[indices] = generate_matrix_normal_parallel(sub_m,sub_p,sub_gg/sub_m)
        betas[:,s] += tmp_beta 
    betas *= np.sqrt(np.diag(gg)/m_causal)/betas.std(0) # making sure the additive genetic variance per phenotype matches the gg matrix diags
        # note: since the shared and specific components both have additive genetic variance of the diagonals in gg
            # the total additive genetic variance of the phenotype is double the diagonal of gg
            # scaling the betas like this forces the total additive genetic variance of the phenotype to match what we want but halves the covariances between phenotypes
        # come back to this after talking to Jeremy?
    gen_comp = genotypes @ betas
    print('Covariance between genetic component of phenotypes: \n ', np.cov(gen_comp, rowvar=False))
    print('Simulating environmental noise')
    env_comp = generate_matrix_normal_parallel(n, nphen, ge)
    return betas,gen_comp,env_comp