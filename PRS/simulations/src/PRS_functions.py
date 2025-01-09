from PLINK2_class import *
import numpy as np 
import pandas as pd
import sys
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler 

def PRS_heldout(discovery_fold, testing_fold, freq_file, ld_ref, P_df, intermediate_dir, p_thresholds, plink_validation=None):
    ### PRS on held-out testing cohort ###
    ### P+T, default clumping/pruning parameters and p-value thresholding optionally tuned ###

    ### INPUTS ###
        # discovery_fold: Plink2 object for discovery cohort in which to perform GWAS 
        # testing_fold: Plink2 object for testing cohort in which to build PRS and calculate R-squared
        # freq_file: plink <.afreq/.acount/.gcount/.freq/.frq/.frq.count/.frqx> filename
            # uses an agreed set of minor alleles
        # ld_ref: Plink2 object to use as LD reference (can be the same cohort or external like 1KG)
        # P_df: dataframe of phenotype, first two columns should be FID/IID (contains all individuals)
        # intermediate_dir: directory to store intermediate files 
        # p_thresholds: list of p-value thresholds to test
        # plink_validation: Plink2 object for validation cohort, required if p_thresholds contains more than 1 value
            # None if only 1 p-value

    ### OUTPUTS ###
        # p_thresholds_used: dictionary mapping phenotype name (in P_df) <-> p-value threshold used
        # P_predictions_df_filtered: pandas dataframe of PRS predictions of all phenotypes in P_df file, only for testing cohort samples

    # I/O checks
    validation_boolean = len(p_thresholds) > 1
    if validation_boolean: 
        if plink_validation is None:
            print("PRS_CV: plink_validation is required when p_thresholds contain more than 1 value. See plink_validation_split function")
            sys.exit(1)
        else:
            plink_validation.filter_phenotype(P_df) # validation

    # Prepare output files
    P_predictions_df = pd.DataFrame(np.full( P_df.shape, np.nan), columns = P_df.columns).astype(object)
        # P_df.shape includes validation cohort individuals, so predictions for those individuals will be nan
        # can filter after populating the entire dataframe
    P_predictions_df.iloc[:,0:2] = P_df.iloc[:,0:2].astype(object)
    p_thresholds_used = {} # keeping track of p-value thresholds used for each phenotype in P_df
        # phenotype <-> (fold, p-value threshold)

    print("Building PRS")

    # filter phenotype
    discovery_fold.filter_phenotype(P_df) # discovery
    testing_fold.filter_phenotype(P_df) # test

    # run GWAS (.acount file in main_plink2.bfile + '.acount')
    # if real data, compute genotype PC's (~ 6 top eigenvectors, depends on cohort and sample size) and use as covariates for GWAS 
    results_map = discovery_fold.gwas(freq_file)

    for key, value in results_map.items(): # for each phenotype in pheno file

        # clumping/pruning, using in-sample LD and default parameters
        # 'ID', 'P' are snp and p-value column names in gwas results
        valid_snps = ld_ref.clump(value, 'ID', 'P', '%s.%s'%(discovery_fold.bfile, key))

        # if tuning p-values, generate PRS for validation cohort using various p-value thresholds 
        if validation_boolean: 
            prs_map_validation = plink_validation.prs(value, 'ID', 'P', 'A1', 'BETA', valid_snps, p_thresholds, '%s_validation.%s'%(discovery_fold.bfile, key) )
            # identify the best p-value threshold by comparing to plink_validation.pheno
            r2_validation = []
            pheno_validation = pd.read_csv(plink_validation.pheno, sep = '\t')
            y_actual = pheno_validation[key].values
            for p_val in p_thresholds:
                y_predicted = prs_map_validation[p_val].iloc[:,2].values.reshape(-1,1)
                r2_validation.append(compute_r2(y_predicted, y_actual))
            # define p_best as an array with 1 entry (so we can use the same Plink2.prs function )
            p_best = [p_thresholds[np.nanargmax(r2_validation)]]
        else: # if only 1 p-value is provided, use that
            p_best = p_thresholds
            
        # run prs on testing with p_best
        prs_test = testing_fold.prs(value, 'ID', 'P', 'A1', 'BETA', valid_snps, p_best, '%s.%s'%(testing_fold.bfile, key) )
        p_thresholds_used[(key, 1)] = p_best[0]
            # no cross validation here but other functions/code expect p_thresholds_used to be a dictionary with (phenotype, cv_count) as key. here just default cv_count to 1
            # in the cross validation version, the p-value threshold may be different for each fold hence why I designed the dictionary this way
        prs_test_df = list(prs_test.values())[0]
        indices = P_predictions_df.index[P_predictions_df.iloc[:,1].isin(prs_test_df.iloc[:,1]) ].tolist()
        P_predictions_df.loc[indices, key] = prs_test_df.iloc[:, 2].values

    # filter out validation cohort individuals 
    P_predictions_df_filtered = P_predictions_df.dropna()

    return p_thresholds_used, P_predictions_df_filtered

def PRS_CV(plink_folds, freq_file, ld_ref, P_df, intermediate_dir, p_thresholds, plink_validation=None):
    ### PRS cross-validation ###
    ### P+T, default clumping/pruning parameters and p-value thresholding optionally tuned ###

    ### INPUTS ###
        # plink_folds: array of tuples, each containing two Plink2 objects for each fold: (discovery_plink2, test_plink2)
            # see plink_cv_split function below
        # freq_file: plink <.afreq/.acount/.gcount/.freq/.frq/.frq.count/.frqx> filename
            # uses an agreed set of minor alleles
        # ld_ref: Plink2 object to use as LD reference (can be the same cohort or external like 1KG)
        # P_df: dataframe of phenotype, first two columns should be FID/IID (contains all individuals)
        # intermediate_dir: directory to store intermediate files for cross-validation 
        # p_thresholds: list of p-value thresholds to test
        # plink_validation: Plink2 object for validation cohort, required if p_thresholds contains more than 1 value
            # None if only 1 p-value

    ### OUTPUTS ###
        # p_thresholds_used: dictionary mapping phenotype name (in P_df) <-> p-value threshold used
        # P_predictions_df_filtered: pandas dataframe of PRS predictions of all phenotypes in P_df file (via cross-validation), excluding validation samples

    # I/O checks
    validation_boolean = len(p_thresholds) > 1
    if validation_boolean: 
        if plink_validation is None:
            print("PRS_CV: plink_validation is required when p_thresholds contain more than 1 value. See plink_validation_split function")
            sys.exit(1)
        else:
            plink_validation.filter_phenotype(P_df) # validation

    # Prepare output files
    P_predictions_df = pd.DataFrame(np.full( P_df.shape, np.nan), columns = P_df.columns).astype(object)
        # P_df.shape includes validation cohort individuals, so predictions for those individuals will be nan
        # can filter after populating the entire dataframe
    P_predictions_df.iloc[:,0:2] = P_df.iloc[:,0:2].astype(object)
    p_thresholds_used = {} # keeping track of p-value thresholds used for each phenotype in P_df
        # phenotype <-> (fold, p-value threshold)

    # for each validation / testing fold 
    cv_count = 1
    for plink_fold in plink_folds:

        # counter
        print("Beginning cross-validation fold %i"%(cv_count))

        # filter phenotype
        plink_fold[0].filter_phenotype(P_df) # discovery
        plink_fold[1].filter_phenotype(P_df) # test

        # run GWAS (.acount file in main_plink2.bfile + '.acount')
        # if real data, compute genotype PC's (~ 6 top eigenvectors, depends on cohort and sample size) and use as covariates for GWAS 
        results_map = plink_fold[0].gwas(freq_file)

        for key, value in results_map.items(): # for each phenotype in pheno file

            # clumping/pruning, using in-sample LD and default parameters
            # 'ID', 'P' are snp and p-value column names in gwas results
            valid_snps = ld_ref.clump(value, 'ID', 'P', '%s.%s'%(plink_fold[0].bfile, key))

            # if tuning p-values, generate PRS for validation cohort using various p-value thresholds 
            if validation_boolean: 
                prs_map_validation = plink_validation.prs(value, 'ID', 'P', 'A1', 'BETA', valid_snps, p_thresholds, '%s_validation.%s'%(plink_fold[0].bfile, key) )
                # identify the best p-value threshold by comparing to plink_validation.pheno
                r2_validation = []
                pheno_validation = pd.read_csv(plink_validation.pheno, sep = '\t')
                y_actual = pheno_validation[key].values
                for p_val in p_thresholds:
                    y_predicted = prs_map_validation[p_val].iloc[:,2].values.reshape(-1,1)
                    r2_validation.append(compute_r2(y_predicted, y_actual))
                # define p_best as an array with 1 entry (so we can use the same Plink2.prs function )
                p_best = [p_thresholds[np.nanargmax(r2_validation)]]
            else: # if only 1 p-value is provided, use that
                p_best = p_thresholds
            
            # run prs on testing with p_best
            prs_test = plink_fold[1].prs(value, 'ID', 'P', 'A1', 'BETA', valid_snps, p_best, '%s.%s'%(plink_fold[1].bfile, key) )
            p_thresholds_used[(key, cv_count)] = p_best[0]
            prs_test_df = list(prs_test.values())[0]
            indices = P_predictions_df.index[P_predictions_df.iloc[:,1].isin(prs_test_df.iloc[:,1]) ].tolist()
            P_predictions_df.loc[indices, key] = prs_test_df.iloc[:, 2].values

        cv_count+=1

    # filter out validation cohort individuals 
    P_predictions_df_filtered = P_predictions_df.dropna()

    return p_thresholds_used, P_predictions_df_filtered

def plink_validation_test_split(PLINK2_object, intermediate_dir, fraction, out_prefix):
    ### Split plink file into discovery, testing and validation cohorts ###

    # inputs 
        # PLINK2_object: main PLINK2_class object
        # intermediate_dir: directory to store discovery, testing, validation cohort plink files (scratch)
        # fraction: fraction of individuals to set aside as validation and testing (same fraction for each)
        # out_prefix: prefix of filtered plink files
    # outputs 
        # a tuple with three Plink2 objects, discovery, testing, validation
        
    fam = pd.read_csv(PLINK2_object.bfile + '.fam', sep = '\t', header = None)
    validation_indices = np.arange(fam.shape[0])
    np.random.shuffle(validation_indices)
    validation_fam = fam.iloc[validation_indices[0:int(fam.shape[0]*fraction)],:]
    validation_fam.to_csv( '%s/foldvalidation.txt'%(intermediate_dir) , sep = '\t', index = False, header = None)
    validation_plink2 = PLINK2_object.keep_ids('%s/foldvalidation.txt'%(intermediate_dir), '%s/%s_validation'%(intermediate_dir, out_prefix))
    discovery_test_plink2 = PLINK2_object.remove_ids('%s/foldvalidation.txt'%(intermediate_dir), '%s/%s_discovery_test'%(intermediate_dir, out_prefix))
    testing_fam = fam.iloc[validation_indices[int(fam.shape[0]*fraction):int(fam.shape[0]*fraction)*2],:]
    testing_fam.to_csv( '%s/foldtesting.txt'%(intermediate_dir) , sep = '\t', index = False, header = None)
    testing_plink2 = discovery_test_plink2.keep_ids('%s/foldtesting.txt'%(intermediate_dir), '%s/%s_testing'%(intermediate_dir, out_prefix))
    discovery_plink2 = discovery_test_plink2.remove_ids('%s/foldtesting.txt'%(intermediate_dir), '%s/%s_discovery'%(intermediate_dir, out_prefix))
    return (discovery_plink2, validation_plink2, testing_plink2)

def plink_test_split(PLINK2_object, intermediate_dir, fraction, out_prefix):
    ### Split plink file into discovery and testing cohorts ###

    # inputs 
        # PLINK2_object: main PLINK2_class object
        # intermediate_dir: directory to store discovery/testing cohort plink files (scratch)
        # fraction: fraction of individuals to set aside as testing
        # out_prefix: prefix of filtered plink files
    # outputs 
        # a tuple with two Plink2 objects, one for discovery and one for testing
        
    fam = pd.read_csv(PLINK2_object.bfile + '.fam', sep = '\t', header = None)
    testing_indices = np.arange(fam.shape[0])
    np.random.shuffle(testing_indices)
    testing_fam = fam.iloc[testing_indices[0:int(fam.shape[0]*fraction)],:]
    testing_fam.to_csv( '%s/foldtesting.txt'%(intermediate_dir) , sep = '\t', index = False, header = None)
    testing_plink2 = PLINK2_object.keep_ids('%s/foldtesting.txt'%(intermediate_dir), '%s/%s_testing'%(intermediate_dir, out_prefix))
    discovery_plink2 = PLINK2_object.remove_ids('%s/foldtesting.txt'%(intermediate_dir), '%s/%s_discovery'%(intermediate_dir, out_prefix))
    return (discovery_plink2, testing_plink2)

def plink_validation_split(PLINK2_object, intermediate_dir, fraction, out_prefix):
    ### Split plink file into discovery+testing and validation cohorts ###

    # inputs 
        # PLINK2_object: main PLINK2_class object
        # intermediate_dir: directory to store discovery+testing/validation cohort plink files (scratch)
        # fraction: fraction of individuals to set aside as validation
        # out_prefix: prefix of filtered plink files
    # outputs 
        # a tuple with two Plink2 objects, one for discovery+testing and one for validation
        
    fam = pd.read_csv(PLINK2_object.bfile + '.fam', sep = '\t', header = None)
    validation_indices = np.arange(fam.shape[0])
    np.random.shuffle(validation_indices)
    validation_fam = fam.iloc[validation_indices[0:int(fam.shape[0]*fraction)],:]
    validation_fam.to_csv( '%s/foldvalidation.txt'%(intermediate_dir) , sep = '\t', index = False, header = None)
    validation_plink2 = PLINK2_object.keep_ids('%s/foldvalidation.txt'%(intermediate_dir), '%s/%s_validation'%(intermediate_dir, out_prefix))
    discovery_test_plink2 = PLINK2_object.remove_ids('%s/foldvalidation.txt'%(intermediate_dir), '%s/%s_discovery_test'%(intermediate_dir, out_prefix))
    return (discovery_test_plink2, validation_plink2)

def plink_cv_split(PLINK2_object, intermediate_dir, cv_folds):
    ### Split plink file for cross validation ###

    # inputs
        # PLINK2_object: main PLINK2_class object
        # intermediate_dir: directory to store intermediate files for cross-validation 
        # cv_folds: number of cross-validation folds

    # output 
        # array of tuples, each tuple contains two Plink2 objects (training fold, testing fold)

    folds = []
    fam = pd.read_csv(PLINK2_object.bfile + '.fam', sep = '\t', header = None)
    if fam.shape[0] % cv_folds != 0:
        print("Individuals in validation+testing (N*(1-validation_fraction)) must be divisible by cv_folds to ensure equal counts.")
        sys.exit(1)
    assignments = np.repeat(np.arange(1, cv_folds + 1), int(fam.shape[0] / cv_folds), axis = 0 )
    np.random.shuffle(assignments)
    
    for fold in np.arange(1, cv_folds+1):
        fold_indices = np.where(assignments == fold)[0]
        fold_fam = fam.iloc[fold_indices,:]
        fold_fam.to_csv('%s/fold%itest.txt'%(intermediate_dir, fold), sep = '\t', index = False, header = None)
        test_plink2 = PLINK2_object.keep_ids( '%s/fold%itest.txt'%(intermediate_dir, fold), '%s/fold%itest'%(intermediate_dir, fold) )
        discovery_plink2 = PLINK2_object.remove_ids( '%s/fold%itest.txt'%(intermediate_dir, fold), '%s/fold%idiscovery'%(intermediate_dir, fold) )
        folds.append((discovery_plink2, test_plink2))

    return folds


def evaluate_predictions(pheno_predicted, pheno_actual):
    ### computes the r-squared between each phenotype column and returns it as a dataframe with phenotypes as columns and one row of r-squared values
    # pheno_predicted: phenotype file with first two columns FID/IID with predictions
    # pheno_actual: phenotype file with first two columns FID/IID with real values

    nphen = pheno_predicted.shape[1] - 2
    pheno_names =  list(pheno_predicted.columns)

    # shallow copy dataframe so we can change column name
    pheno_actual = pheno_actual.copy()
    pheno_actual.columns = list(pheno_actual.columns)[0:2] + [ col + '_actual' for col in list(pheno_actual.columns)[2:] ]
    pheno_merged = pd.merge(pheno_predicted, pheno_actual, on = ['FID', 'IID'], how = "inner")
    r2_map = {}
    for i in range(2, nphen + 2):
        r2_map[ pheno_names[i] ] = [compute_r2( np.array(pheno_merged.iloc[:,i].values.reshape(-1,1), dtype=float) , pheno_merged.iloc[:,(i+nphen)].values )]
    return(pd.DataFrame(r2_map))


def compute_r2(y_predicted, y_actual):
    if np.isnan(y_predicted).all():
        return np.nan
    else:
        lreg = LinearRegression().fit( y_predicted ,  y_actual )
        return lreg.score( y_predicted ,  y_actual )

def fit_lm(X, y):
    scaler = StandardScaler() 
    model = LinearRegression()
    X = scaler.fit_transform(X)
    y = (y - y.mean())/y.std()
    model.fit(X, y)
    coeffs = model.coef_
    return coeffs


def clean_up(prefix):
    cmd = "rm %s*"%(prefix)
    subprocess.run(cmd, text=True, shell=True)



