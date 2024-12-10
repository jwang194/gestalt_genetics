import sys
import numpy as np
import pandas as pd
from scipy.stats import matrix_normal
from scipy.stats import pearsonr,spearmanr

heritability_string = sys.argv[1] # comma separated list of heritability values for traits. e.g 0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1 OR file with 2 columns (Phenotype, rg)
genetic_correlation = sys.argv[2] # uniform genetic correlation between triats e.g. 0.05 OR file with 3 columns (Phenotype 1, Phenotype 2, rg) with no repeated values
environmental_correlation = sys.argv[3] # mean of randomly drawn environmental covariance between traits
total_causal = sys.argv[4] # total number of causal variants (10000)
shared_causal_fraction = sys.argv[5] # fraction of causal variants that are shared (0.2)
specific_causal_fraction = sys.argv[6] # fraction of causal variants that are specific to each phenotype (0.08)
custom = sys.argv[7] # True or False, indicating the use of custom file for heritability and genetic correlation

if not custom: 
    heritability_list = np.array( [float(x) for x in heritability_string.split(',')] ) 
else:
    heritability_df = pd.read_csv(heritability_string, sep = '\t')
    heritability_df.columns = ['Trait', 'Heritability']
    heritability_map = dict()
    for _,row in heritability_df.iterrows():
        heritability_map[row['Trait']] = row['Heritability']
    heritability_map = dict(sorted(heritability_map.items()))
    heritability_list = np.array(list(heritability_map.values()))

num_traits = len(heritability_list)
print('Number of phenotypes: ', num_traits)

if not custom:
    setting_name = 'traits_%i_shared_%s_specific_%s_uniform_rg_%s_random_re_%s'%(num_traits,shared_causal_fraction,specific_causal_fraction,genetic_correlation,environmental_correlation)
    print(setting_name)
else:
    setting_name = 'traits_%i_shared_%s_specific_%s_uniform_rg_custom_random_re_%s'%(num_traits,shared_causal_fraction,specific_causal_fraction,environmental_correlation)
    print(setting_name)

if not custom:
    genetic_covariance = np.full((num_traits, num_traits), float(genetic_correlation) )
    np.fill_diagonal(genetic_covariance, heritability_list)
    print('Genetic covariance: ')
    print(genetic_covariance)
    np.savetxt(setting_name + '_gg.txt', genetic_covariance, fmt='%.12f')
else:
    print('Reading genetic correlation file')
    df_rg = pd.read_csv(genetic_correlation, sep='\t')
    phenotypes = sorted(set(df_rg["Phenotype 1"]).union(df_rg["Phenotype 2"]))
    genetic_covariance = pd.DataFrame(np.nan, index=phenotypes, columns=phenotypes)
    for _, row in df_rg.iterrows():
        covariance = row["rg"]*np.sqrt(heritability_map[row["Phenotype 1"]])*np.sqrt(heritability_map[row["Phenotype 2"]])
        genetic_covariance.loc[row["Phenotype 1"], row["Phenotype 2"]] = covariance
        genetic_covariance.loc[row["Phenotype 2"], row["Phenotype 1"]] = covariance
    genetic_covariance = genetic_covariance.to_numpy()
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


