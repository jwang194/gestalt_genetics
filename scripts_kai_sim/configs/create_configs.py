import sys
import numpy as np
import pandas as pd
from scipy.stats import matrix_normal
from scipy.stats import pearsonr,spearmanr
import subprocess

heritability_string = sys.argv[1] # comma separated list of heritability values for traits. e.g 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
bim_file = sys.argv[2] # bim file from genotype data used to simulate phenotypes
genetic_covariance = sys.argv[3] # uniform genetic covariance between triats e.g. 0.05
environmental_covariance = sys.argv[4] # mean of randomly drawn environmental covariance between traits e.g. -0.1
total_snps = sys.argv[5] # total number of variants e.g. 1,000,000 HapMap3 variants
causal_fraction = sys.argv[6]
shared_causal_fraction = sys.argv[7] # fraction of causal variants that are shared (0.2)
outdir = sys.argv[8]

# collecting args and printing settings
heritability_list = np.array( [float(x) for x in heritability_string.split(',')] ) 
num_traits = len(heritability_list)
print('Number of phenotypes: ', num_traits)
specific_causal_fraction = (1 - float(shared_causal_fraction))/num_traits
num_causal = int(float(total_snps) * float(causal_fraction))
print('Number of causal variants: ', num_causal)
setting_name = 'traits_%i_causal_%s_shared_%s_uniform_rg_%s_random_re_%s'%(num_traits,causal_fraction,shared_causal_fraction,genetic_covariance,environmental_covariance)
print(setting_name)
cmd = ['mkdir', '-p', '%s/%s'%(outdir,setting_name)]
subprocess.run(cmd, text= True)

# create genetic covariance matrix 
genetic_covariance_matrix = np.full((num_traits, num_traits), float(genetic_covariance) )
np.fill_diagonal(genetic_covariance_matrix, heritability_list)
print('Genetic covariance: ')
print(genetic_covariance_matrix)
np.savetxt( '%s/%s/%s_gg.txt'%(outdir,setting_name, setting_name), genetic_covariance_matrix, fmt='%.12f')

# create environmental covariance matrix 
environmental_covariance_matrix = np.random.uniform( float(environmental_covariance), 0, size=(num_traits, num_traits))
np.fill_diagonal(environmental_covariance_matrix, [ (1-x) for x in heritability_list] )
print('Environmental covariance: ')
print(environmental_covariance_matrix)
np.savetxt( '%s/%s/%s_ge.txt'%(outdir,setting_name, setting_name), environmental_covariance_matrix, fmt='%.12f')

# create variant assignment file
    # maps snp number (randomly chosen causal variants) <-> assignment
num_shared = int(float(num_causal) * float(shared_causal_fraction))
num_specific_per_pheno = int(float(num_causal) * float(specific_causal_fraction))

# randomly shuffle variant ids
bim = pd.read_csv(bim_file, sep='\t', header=None)
variant_ids = bim.iloc[:,1].to_numpy()
np.random.shuffle(variant_ids) # shuffles in place

# assign shared and specific causal variants 
assignments = np.concatenate([
    np.repeat(','.join(map(str, range(0, num_traits))), num_shared), 
    [pheno for pheno in list(range(0, num_traits)) for _ in range(num_specific_per_pheno)]
])

# save causal variant assignments
variant_assignments = np.column_stack((variant_ids[:len(assignments)], assignments))
sorted_variant_assignments = variant_assignments[ np.array([x.replace('snp', '') for x in variant_assignments[:, 0]]).astype(int).argsort()]
np.savetxt('%s/%s/%s_assignments.txt'%(outdir,setting_name, setting_name), sorted_variant_assignments, fmt='%s', delimiter='\t')
