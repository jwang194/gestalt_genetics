import argparse
from PLINK2_class import *
from PRS_functions import *
import numpy as np 
import pandas as pd
import os

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Python script for polygenic risk score analysis of GESTALT trait project simulations."
    )

    ### Required Arguments ###
    parser.add_argument(
        "--plink2_path", type=str, required=True, 
        help="Path to the PLINK2 executable."
    )
    parser.add_argument(
        "--genotypes", type=str, required=True,
        help="Path and prefix to the PLINK genotype files."
    )
    parser.add_argument(
        "--phenotype_prefix", type=str, required=True,
        help="Path and prefix for the phenotype files ( '_P.txt' and '_MG.txt' files)."
    )
    parser.add_argument(
        "--intermediate_dir", type=str, required=True,
        help="Directory for intermediate cross-validation files."
    )
    parser.add_argument(
        "--output_prefix", type=str, required=True,
        help="Path and prefix for the output files."
    )

    ### Optional Arguments ###
    parser.add_argument(
        "--model_type", type=str, default="P,MG,SUMPRS_MG,PRSSUM,SUMPRS",
        help=(
            "Comma-separated list of model/pheno types to evaluate. "
            "Choices: ['P', 'MG', 'SUMPRS_MG', 'PRSSUM', 'SUMPRS']."
        )
    )
    parser.add_argument(
        "--threads", type=int, default=1,
        help="Number of threads to use. Default is 1."
    )
    parser.add_argument(
        "--evaluate_r2", action="store_true",
        help="Evaluate R-squared between true and predicted phenotypes."
    )
    parser.add_argument(
        "--p_values", type=str, default="5e-8,5e-6,5e-4,0.05",
        help="Comma-separated list of p-value thresholds. Default is '5e-8,5e-6,5e-4,0.05' and the best will be chosen via prediction in validation fold."
    )
    parser.add_argument(
        "--validation_test_fraction", type=float, default=0.1, 
        help="Fraction of cohort in which to evaluate PRS accuracy with different p-value thresholds"
    )
    parser.add_argument(
        "--ldref", type=str, default=None,
        help="Path to LD reference file for clumping/pruning. If not provided, it defaults to the value of --genotypes."
    )

    return parser.parse_args()

# functions like this with variable names specific to this script will be kept in the same file
def read_phenotype_files(phenotype_prefix, fam_file, models):
    # reads phenotype files, filters for individuals in fam file, and returns a dictionary of dataframes

    fam = pd.read_csv(fam_file, header = None, sep = '\t')

    # keys: 'P', 'MG', 'SUM'
    phenos = {}
    if 'P' in models or 'PRSSUM' in models or 'SUMPRS' in models:
        if not os.path.exists(phenotype_prefix + '_P.txt'):
            print("Cannot perform PRS for [P, PRSSUM, SUMPRS] without _P.txt file ")
            sys.exit(1)
        P_df = pd.read_csv(phenotype_prefix + '_P.txt', sep='\t')
        P_df_filtered = P_df[P_df.iloc[:,1].isin(fam.iloc[:,1])]
        P_df_filtered = P_df_filtered.set_index(P_df_filtered.columns[1], drop = False)
        P_df_filtered = P_df_filtered.reindex(fam.iloc[:, 1]).reset_index(drop = True)
        P_df_filtered.columns = ["FID", "IID"] + list(P_df_filtered.columns[2:])
        phenos['P'] = P_df_filtered
        #phenos['P'].iloc[:,0] = phenos['P'].iloc[:,1] # ONLY FOR TESTING DATASET
    if 'MG' in models or 'SUMPRS_MG' in models:
        if not os.path.exists(phenotype_prefix + '_MG.txt'):
            print("Cannot perform PRS for [MG, SUMPRS_MG] without _MG.txt file ")
            sys.exit(1)
        MG_df = pd.read_csv(phenotype_prefix + '_MG.txt', sep='\t').iloc[:,[0,1,-1]]
        MG_df_filtered = MG_df[MG_df.iloc[:,1].isin(fam.iloc[:,1])]
        MG_df_filtered = MG_df_filtered.set_index(MG_df_filtered.columns[1], drop = False)
        MG_df_filtered = MG_df_filtered.reindex(fam.iloc[:, 1]).reset_index(drop = True)
        phenos['MG'] = MG_df_filtered
        phenos['MG'].columns = ['FID', 'IID', 'MAXH']
        #phenos['MG'].iloc[:,0] = phenos['MG'].iloc[:,1] # ONLY FOR TESTING DATASET
    if 'PRSSUM' in models or 'SUMPRS' in models:
        phenos['SUM'] = pd.concat( [ phenos['P'].iloc[:,0:2], phenos['P'].iloc[:,2:phenos['P'].shape[1]].sum(axis=1) ], axis=1)
        phenos['SUM'].columns = ['FID', 'IID', 'SUM']
        #phenos['SUM'].iloc[:,0] = phenos['SUM'].iloc[:,1] # ONLY FOR TESTING DATASET
    return phenos

def get_memory_info():
    with open("/proc/meminfo", "r") as f:
        meminfo = f.readlines()
    meminfo = {line.split(":")[0]: line.split(":")[1].strip() for line in meminfo}
    total_mem_kb = int(meminfo["MemTotal"].split()[0])
    free_mem_kb = int(meminfo["MemAvailable"].split()[0])
    return total_mem_kb / 1024, free_mem_kb / 1024  # Convert to MB

def main():
    args = parse_arguments()

    ### Required arguments ###
    print("Path to PLINK2:", args.plink2_path)
    print("Genotype files:", args.genotypes)
    print("Phenotype prefix:", args.phenotype_prefix)
    print("Intermediate directory:", args.intermediate_dir)
    print("Output prefix:", args.output_prefix)

    ### Optional arguments ###
    models = args.model_type.split(",")
    print("Model/Pheno types:", models)
    print("Number of threads:", args.threads)
    print("Evaluate R-squared:", args.evaluate_r2)
    p_thresholds = [float(value) for value in args.p_values.split(",")]
    print("p-value thresholds:", p_thresholds)
    print("Validation and testing fraction:", args.validation_test_fraction)
    if args.ldref is None:
        args.ldref = args.genotypes
    print("LD reference file:", args.ldref)

    ### Initialize objects and data ###

    # initialize main Plink2 object
    main_plink2 = Plink2(args.plink2_path, args.genotypes, args.threads)

    # initialize ldref Plink2 object 
    ldref_plink2 = Plink2(args.plink2_path, args.ldref, args.threads)

    # read in phenotype files
    PHENOS = read_phenotype_files(args.phenotype_prefix, main_plink2.bfile + '.fam' , models) # dictionary of dataframes, plink fam format phenotype data
    
    ### Start PRS evaluation ###
    ### Discovery = train GWAS
    ### Validation = tune p-value thresholds or compute coefficients to combine PRS
    ### Testing = PRS prediction 

    # memory check
    total_memory, available_memory = get_memory_info()
    print(f"Total Memory: {total_memory:.2f} MB")
    print(f"Available Memory: {available_memory:.2f} MB")

    # check if validation cohort is necessary to tune across multiple p-values or fit PRS to predict MAXH phenotype
    if len(p_thresholds) > 1 or 'SUMPRS_MG' in models:
        discovery_test_validation_plink2 = plink_validation_test_split(main_plink2, args.intermediate_dir, args.validation_test_fraction, "main") 
            # "main" is prefix of _discovery.* and _test.* and _validation.* plink files
        discovery_plink2 = discovery_test_validation_plink2[0]
        validation_plink2 = discovery_test_validation_plink2[1]
        testing_plink2 = discovery_test_validation_plink2[2]
    else:
        discovery_test_plink2 = plink_test_split(main_plink2, args.intermediate_dir, args.validation_test_fraction, "main") 
            # "main" is prefix of _discovery.* and _test.* plink files
        discovery_plink2 = discovery_test_plink2[0]
        testing_plink2 = discovery_test_plink2[1]
        validation_plink2 = None

    # PRS section for P phenotypes
    if 'P' in models:
        print("PRS for P model")
        p_thresholds_P, predictions_P = PRS_heldout(discovery_plink2, testing_plink2, main_plink2.bfile + '.acount', ldref_plink2, PHENOS['P'], args.intermediate_dir, p_thresholds, validation_plink2)
        predictions_P.to_csv(args.output_prefix + "_P.txt", index = False, sep = '\t')
        if args.evaluate_r2:
            r2_df_P = evaluate_predictions(predictions_P, PHENOS['P'] )
            r2_df_P.to_csv(args.output_prefix + "_P_r2.txt", index = False, sep = '\t')

    # PRS section for MG phenotypes
    if 'MG' in models:
        print("PRS for MG model")
        p_thresholds_MG, predictions_MG = PRS_heldout(discovery_plink2, testing_plink2, main_plink2.bfile + '.acount', ldref_plink2, PHENOS['MG'], args.intermediate_dir, p_thresholds, validation_plink2)
        predictions_MG.to_csv(args.output_prefix + "_MG.txt", index = False, sep = '\t')
        if args.evaluate_r2:
            r2_df_MG = evaluate_predictions(predictions_MG, PHENOS['MG'] )
            r2_df_MG.to_csv(args.output_prefix + "_MG_r2.txt", index = False, sep = '\t')

    # PRS section for PRSSUM phenotypes
    if 'PRSSUM' in models:
        print("PRS for PRSSUM model")
        p_thresholds_PRSSUM, predictions_PRSSUM = PRS_heldout(discovery_plink2, testing_plink2, main_plink2.bfile + '.acount', ldref_plink2, PHENOS['SUM'], args.intermediate_dir, p_thresholds, validation_plink2)
        predictions_PRSSUM.to_csv(args.output_prefix + "_PRSSUM.txt", index = False, sep = '\t')
        if args.evaluate_r2:
            r2_df_PRSSUM = evaluate_predictions(predictions_PRSSUM, PHENOS['SUM'] )
            r2_df_PRSSUM.to_csv(args.output_prefix + "_PRSSUM_r2.txt", index = False, sep = '\t')

    # PRS section for SUMPRS 
    if 'SUMPRS' in models:
        print("PRS for SUMPRS model")
        if 'P' not in models: # if P not in models, we have to run prs predictions on _P.txt so we can sum the results
            p_thresholds_P, predictions_P = PRS_heldout(discovery_plink2, testing_plink2, main_plink2.bfile + '.acount', ldref_plink2, PHENOS['P'], args.intermediate_dir, p_thresholds, validation_plink2)
        # sum the PRS prediction results 
        predictions_SUMPRS = pd.concat( [predictions_P.iloc[:, 0:2], predictions_P.iloc[:, 2:].sum(axis = 1, skipna = False)], axis = 1)
        predictions_SUMPRS.columns = list(predictions_P.columns)[0:2] + ['SUMPRS']
        predictions_SUMPRS.to_csv(args.output_prefix + "_SUMPRS.txt", index = False, sep = '\t')
        if args.evaluate_r2:
            r2_df_SUMPRS = evaluate_predictions(predictions_SUMPRS, PHENOS['SUM'] )
            r2_df_SUMPRS.to_csv(args.output_prefix + "_SUMPRS_r2.txt", index = False, sep = '\t')

    # fit SUMPRS_MG using validation cohort
    if 'SUMPRS_MG' in models:
        print("PRS for SUMPRS_MG model")
        if 'P' not in models:
            p_thresholds_P, predictions_P = PRS_heldout(discovery_plink2, testing_plink2, main_plink2.bfile + '.acount', ldref_plink2, PHENOS['P'], args.intermediate_dir, p_thresholds, validation_plink2)
        
        # fit coefficients to predict MAXH phenotype from PRS of each trait using validation cohort 
        validation_trainfit_plink2 = plink_validation_split(validation_plink2, args.intermediate_dir, 0.5, "validation") # split validation cohort into two
        # ran into an interesting bug here when I used the best p-value threshold found in discovery-validation cohort 
            # often the p-value threshold is pretty small, and gwas using half the validation cohort does not find any snps at that threshold 
            # and therefore the resulting PRS score is NA, and the linear weights cannot be found with fit_lm 
            # for now I use the 0.05 p-value threshold to avoid this error
        p_thresholds_validation, validation_predictions = PRS_heldout(validation_trainfit_plink2[0], validation_trainfit_plink2[1], main_plink2.bfile + '.acount', ldref_plink2, PHENOS['P'], args.intermediate_dir, [0.05], None)
        validation_plink2.filter_phenotype(PHENOS['MG'])
        validation_MG = pd.read_csv(validation_plink2.pheno, sep = '\t')
        validation_MG_filtered = validation_MG[validation_MG['IID'].isin(validation_predictions['IID'])]
        validation_predictions_reordered = validation_predictions.set_index('IID', drop = False).reindex(validation_MG_filtered['IID']).reset_index(drop = True)
        if (validation_predictions_reordered.isna().any().any()):
            print('SUMPRS_MG cannot be calculated because GWAS on half of the validation cohort did not find any snps p < 0.05. Try increasing the validation_test_fraction argument')
            clean_up(args.intermediate_dir + '/')
            sys.exit(1)
        coefficients = fit_lm(validation_predictions_reordered.iloc[:,2:], validation_MG_filtered.iloc[:,2] ) # fit coefficients with all individuals in validation cohort
        print('Printing coefficients for SUMPRS_MG: ',coefficients)
        predictions_SUMPRS_MG = pd.concat( [predictions_P.iloc[:, 0:2], predictions_P.iloc[:, 2:] @ coefficients ], axis = 1)
        predictions_SUMPRS_MG.columns = list(predictions_P.columns)[0:2] + ['SUMPRS_MG']
        predictions_SUMPRS_MG.to_csv(args.output_prefix + "_SUMPRS_MG.txt", index = False, sep = '\t')
        if args.evaluate_r2:
            r2_df_SUMPRS_MG = evaluate_predictions(predictions_SUMPRS_MG, PHENOS['MG'] )
            r2_df_SUMPRS_MG.to_csv(args.output_prefix + "_SUMPRS_MG_r2.txt", index = False, sep = '\t')

    # clean up intermediate directory
    clean_up(args.intermediate_dir + '/')


if __name__ == "__main__":
    main()