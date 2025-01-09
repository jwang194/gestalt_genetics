# PRS predictions 

> plink conveniently has flags for PRS score calculations 

**PRS_plinkCV.sh** = main pipeline file 

1. split individuals for cross validation 
2. compute GWAS Betas for each training fold, compute PRS for testing fold 
3. output PRS predictions 

**filter_phenotype.py** = intermediate file to filter phenotype file to individuals in training fold 

**collect_results.R** = R script to collect PRS prediction results and summarize accuracy with R-squared

directories added after:

- SUMPRS_MAXH = directory to evaluate sum of PRS fit to maxh 
- SUMPRS_SUM = directory to evaluate sum of PRS to predict sum of phenotype 
