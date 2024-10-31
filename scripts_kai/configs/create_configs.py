import sys
import numpy as np
from scipy.stats import matrix_normal
from scipy.stats import pearsonr,spearmanr


heritability_string = sys.argv[1] # comma separated list of heritability values for traits. e.g 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
genetic_correlation = sys.argv[2] # uniform genetic correlation between triats e.g. 0.05
environmental_correlation = sys.argv[3] # mean of randomly drawn environmental covariance between traits
total_causal = sys.argv[4] # total number of causal variants (10000)
shared_causal_fraction = sys.argv[5] # fraction of causal variants that are shared (0.2)
specific_causal_fraction = sys.argv[6] # fraction of causal variants that are specific to each phenotype (0.08)

heritability_list = np.array( [float(x) for x in heritability_string.split(',')] ) 
num_traits = len(heritability_list)
print('Number of phenotypes: ', num_traits)

setting_name = 'traits_%i_shared_%s_specific_%s_uniform_rg_%s_random_re_%s'%(num_traits,shared_causal_fraction,specific_causal_fraction,genetic_correlation,environmental_correlation)
print(setting_name)

genetic_covariance = np.full((num_traits, num_traits), float(genetic_correlation) )
np.fill_diagonal(genetic_covariance, heritability_list)
print('Genetic covariance: ')
print(genetic_covariance)
np.savetxt(setting_name + '_gg.txt', genetic_covariance, fmt='%.12f')

environmental_covariance = np.random.uniform( float(environmental_correlation), 0, size=(num_traits, num_traits))
np.fill_diagonal(environmental_covariance, [ (1-x) for x in heritability_list] )
print('Environmental covariance: ')
print(environmental_covariance)
np.savetxt(setting_name + '_ge.txt', environmental_covariance, fmt='%.12f')

num_shared = int(float(total_causal) * float(shared_causal_fraction))
num_specific_per_pheno = int(float(total_causal) * float(specific_causal_fraction))

variant_assignments = np.concatenate([
    np.repeat(','.join(map(str, range(0, num_traits))), num_shared), 
    [pheno for pheno in list(range(0, num_traits)) for _ in range(num_specific_per_pheno)]
])

np.savetxt(setting_name + '_assignments.txt', variant_assignments, fmt='%s')


