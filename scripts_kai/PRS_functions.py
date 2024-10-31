import numpy as np
import pandas as pd
from scipy.stats import pearsonr,spearmanr


def PRS_evaluate(true_phenotypes, genotypes, gwas_betas):
    # true_phenotypes = vector of true phenotypes (in same order as genotypes)
    # genotypes = (people x snps) matrix of genotypes
    # gwas_betas = vector of GWAS effect sizes for each snp (in same order as genotypes)
    ### goal: compute PRS scores and return R^2, variance of true phenotype explained
    PRS = genotypes @ gwas_betas
    ss_res = np.sum((true_phenotypes - PRS) ** 2)
    ss_tot = np.sum((true_phenotypes - np.mean(true_phenotypes)) ** 2)
    r2 = 1 - (ss_res / ss_tot)
    return(r2)


def PRS_evaluate_multi(true_phenotypes, genotypes, gwas_betas):
    # true_phenotypes = (people x phenotypes) matrix of true phenotypes (in same order as genotypes)
    # genotypes = (people x snps) matrix of genotypes
    # gwas_betas = (snps, phenotypes) matrix of GWAS effect sizes for each snp (in same order as genotypes)
    ### goal: compute PRS scores for each trait and return R^2, variance of true phenotype explained
    PRS = genotypes @ gwas_betas
    r2_values = []
    for i in range(true_phenotypes.shape[1]):
        ss_res = np.sum((true_phenotypes[:, i] - PRS[:, i]) ** 2)
        ss_tot = np.sum((true_phenotypes[:, i] - np.mean(true_phenotypes[:, i])) ** 2)
        r2 = 1 - (ss_res / ss_tot)
        r2_values.append(r2)
    return(r2)

